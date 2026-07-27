module elasticWave3D

using Base.Threads

export ElasticThreePointConfig, ElasticWaveState3D
export prepare_elastic_wave_3d, add_ricker_source!, step_elastic_wave_3d!
export displacement_magnitude, stable_elastic_timestep

"""
    ElasticThreePointConfig

Second-order explicit isotropic elastic scheme. `pointsInSpace = 3`,
`pointsInTime = 3`, and `supplementaryOrder = 0` are deliberately fixed so
the numerical method recorded by a simulation cannot disagree with the
implemented stencil.
"""
Base.@kwdef struct ElasticThreePointConfig
    pointsInSpace::Int = 3
    pointsInTime::Int = 3
    supplementaryOrder::Int = 0
    cfl::Float64 = 0.42
    cerjan_width::Int = 12
    cerjan_strength::Float64 = 0.015
end

mutable struct ElasticWaveState3D{T<:AbstractFloat,A<:AbstractArray{T,3}}
    ux_previous::A
    uy_previous::A
    uz_previous::A
    ux::A
    uy::A
    uz::A
    ux_next::A
    uy_next::A
    uz_next::A
    vp2::A
    vs2::A
    solid::BitArray{3}
    damping::A
    spacing::NTuple{3,T}
    dt::T
    time::T
    step::Int
    config::ElasticThreePointConfig
end

function _validate_config(config::ElasticThreePointConfig)
    config.pointsInSpace == 3 ||
        throw(ArgumentError("this solver implements exactly 3 spatial points"))
    config.pointsInTime == 3 ||
        throw(ArgumentError("this solver implements exactly 3 time levels"))
    config.supplementaryOrder == 0 ||
        throw(ArgumentError("this solver requires supplementaryOrder = 0"))
    0 < config.cfl < 1 || throw(ArgumentError("cfl must lie in (0, 1)"))
    config.cerjan_width >= 0 ||
        throw(ArgumentError("cerjan_width must be non-negative"))
    config.cerjan_strength >= 0 ||
        throw(ArgumentError("cerjan_strength must be non-negative"))
    return config
end

"""
    stable_elastic_timestep(vp, spacing; cfl=0.42)

Conservative three-dimensional CFL time step. Velocities and spacings must use
the same length unit per second.
"""
function stable_elastic_timestep(
    vp::AbstractArray,
    spacing::NTuple{3,<:Real};
    cfl::Real=0.42,
    mask=trues(size(vp)),
)
    size(mask) == size(vp) || throw(DimensionMismatch("mask and vp differ"))
    vmax = maximum(vp[mask])
    isfinite(vmax) && vmax > 0 ||
        throw(ArgumentError("the material mask contains no positive finite Vp"))
    dx, dy, dz = spacing
    minimum(spacing) > 0 || throw(ArgumentError("grid spacings must be positive"))
    return Float64(cfl) /
           (Float64(vmax) * sqrt(inv(dx^2) + inv(dy^2) + inv(dz^2)))
end

function _cerjan_damping(
    shape::NTuple{3,Int},
    width::Int,
    strength::Real,
    ::Type{T},
) where {T<:AbstractFloat}
    damping = ones(T, shape)
    width == 0 && return damping
    nx, ny, nz = shape
    Threads.@threads for k in 1:nz
        # The upper z side is the physical topography and is not absorbed.
        dk = max(width - (k - 1), 0)
        for j in 1:ny, i in 1:nx
            di = max(width - min(i - 1, nx - i), 0)
            dj = max(width - min(j - 1, ny - j), 0)
            distance = max(di, dj, dk)
            damping[i, j, k] = exp(-T(strength) * T(distance^2))
        end
    end
    return damping
end

"""
    prepare_elastic_wave_3d(seismic_model, spacing; material_mask, dt=nothing)

Create a `Float32` three-level displacement state. `Vp` and `Vs` are expected
in km/s, as returned by `getParamsAndTopo`; grid spacing is in metres.
Air and fluids are excluded by default (`ρ > 0.01` and `Vs > 0`).
"""
function prepare_elastic_wave_3d(
    seismic_model,
    spacing::NTuple{3,<:Real};
    material_mask=nothing,
    dt=nothing,
    config::ElasticThreePointConfig=ElasticThreePointConfig(),
    T::Type{<:AbstractFloat}=Float32,
)
    _validate_config(config)
    vp = seismic_model.Vpv
    vs = seismic_model.Vsv
    ρ = seismic_model.ρ
    size(vp) == size(vs) == size(ρ) ||
        throw(DimensionMismatch("ρ, Vp and Vs must have the same dimensions"))
    ndims(vp) == 3 || throw(DimensionMismatch("the elastic solver needs 3D fields"))

    solid = isnothing(material_mask) ?
            BitArray((ρ .> 0.01) .& (vs .> 0) .& isfinite.(vp) .& isfinite.(vs)) :
            BitArray(material_mask .& (vs .> 0))
    size(solid) == size(vp) ||
        throw(DimensionMismatch("material_mask and seismic model differ"))

    # Convert km/s to m/s once. Squared velocities are sufficient for the
    # locally heterogeneous isotropic displacement operator.
    vp2 = T.((vp .* 1_000).^2)
    vs2 = T.((vs .* 1_000).^2)
    spacingT = T.(spacing)
    timestep = isnothing(dt) ?
               stable_elastic_timestep(
                   sqrt.(vp2),
                   Tuple(spacingT);
                   cfl=config.cfl,
                   mask=solid,
               ) :
               Float64(dt)
    timestep > 0 || throw(ArgumentError("dt must be positive"))

    zero_field = zeros(T, size(vp))
    damping = _cerjan_damping(size(vp), config.cerjan_width,
                              config.cerjan_strength, T)
    return ElasticWaveState3D(
        copy(zero_field), copy(zero_field), copy(zero_field),
        copy(zero_field), copy(zero_field), copy(zero_field),
        copy(zero_field), copy(zero_field), copy(zero_field),
        vp2, vs2, solid, damping, Tuple(spacingT), T(timestep),
        zero(T), 0, config,
    )
end

@inline function _neighbor(field, solid, i, j, k, ic, jc, kc)
    return solid[i, j, k] ? field[i, j, k] : field[ic, jc, kc]
end

"""
    step_elastic_wave_3d!(state)

Advance one explicit three-level time step. At the irregular material/air
interface, missing neighbors use a zero-normal-gradient ghost value. This is a
stable topography-aware preliminary boundary, not yet an exact traction-free
stress condition.
"""
function step_elastic_wave_3d!(state::ElasticWaveState3D{T}) where {T}
    (; ux, uy, uz, ux_previous, uy_previous, uz_previous,
       ux_next, uy_next, uz_next, vp2, vs2, solid, damping) = state
    dx, dy, dz = state.spacing
    idx2, idy2, idz2 = inv(dx^2), inv(dy^2), inv(dz^2)
    i4dxdy, i4dxdz, i4dydz =
        inv(4dx * dy), inv(4dx * dz), inv(4dy * dz)
    dt2 = state.dt^2
    nx, ny, nz = size(ux)

    fill!(ux_next, zero(T))
    fill!(uy_next, zero(T))
    fill!(uz_next, zero(T))
    Threads.@threads for k in 2:(nz - 1)
        for j in 2:(ny - 1), i in 2:(nx - 1)
            solid[i, j, k] || continue

            ux0, uy0, uz0 = ux[i, j, k], uy[i, j, k], uz[i, j, k]
            ux_xx = (_neighbor(ux, solid, i+1,j,k,i,j,k) - 2ux0 +
                     _neighbor(ux, solid, i-1,j,k,i,j,k)) * idx2
            ux_yy = (_neighbor(ux, solid, i,j+1,k,i,j,k) - 2ux0 +
                     _neighbor(ux, solid, i,j-1,k,i,j,k)) * idy2
            ux_zz = (_neighbor(ux, solid, i,j,k+1,i,j,k) - 2ux0 +
                     _neighbor(ux, solid, i,j,k-1,i,j,k)) * idz2
            uy_xx = (_neighbor(uy, solid, i+1,j,k,i,j,k) - 2uy0 +
                     _neighbor(uy, solid, i-1,j,k,i,j,k)) * idx2
            uy_yy = (_neighbor(uy, solid, i,j+1,k,i,j,k) - 2uy0 +
                     _neighbor(uy, solid, i,j-1,k,i,j,k)) * idy2
            uy_zz = (_neighbor(uy, solid, i,j,k+1,i,j,k) - 2uy0 +
                     _neighbor(uy, solid, i,j,k-1,i,j,k)) * idz2
            uz_xx = (_neighbor(uz, solid, i+1,j,k,i,j,k) - 2uz0 +
                     _neighbor(uz, solid, i-1,j,k,i,j,k)) * idx2
            uz_yy = (_neighbor(uz, solid, i,j+1,k,i,j,k) - 2uz0 +
                     _neighbor(uz, solid, i,j-1,k,i,j,k)) * idy2
            uz_zz = (_neighbor(uz, solid, i,j,k+1,i,j,k) - 2uz0 +
                     _neighbor(uz, solid, i,j,k-1,i,j,k)) * idz2

            uy_xy = (uy[i+1,j+1,k] - uy[i+1,j-1,k] -
                     uy[i-1,j+1,k] + uy[i-1,j-1,k]) * i4dxdy
            uz_xz = (uz[i+1,j,k+1] - uz[i+1,j,k-1] -
                     uz[i-1,j,k+1] + uz[i-1,j,k-1]) * i4dxdz
            ux_xy = (ux[i+1,j+1,k] - ux[i+1,j-1,k] -
                     ux[i-1,j+1,k] + ux[i-1,j-1,k]) * i4dxdy
            uz_yz = (uz[i,j+1,k+1] - uz[i,j+1,k-1] -
                     uz[i,j-1,k+1] + uz[i,j-1,k-1]) * i4dydz
            ux_xz = (ux[i+1,j,k+1] - ux[i+1,j,k-1] -
                     ux[i-1,j,k+1] + ux[i-1,j,k-1]) * i4dxdz
            uy_yz = (uy[i,j+1,k+1] - uy[i,j+1,k-1] -
                     uy[i,j-1,k+1] + uy[i,j-1,k-1]) * i4dydz

            cp2, cs2 = vp2[i,j,k], vs2[i,j,k]
            cross = cp2 - cs2
            ax = cp2 * ux_xx + cs2 * (ux_yy + ux_zz) +
                 cross * (uy_xy + uz_xz)
            ay = cp2 * uy_yy + cs2 * (uy_xx + uy_zz) +
                 cross * (ux_xy + uz_yz)
            az = cp2 * uz_zz + cs2 * (uz_xx + uz_yy) +
                 cross * (ux_xz + uy_yz)
            damp = damping[i,j,k]
            ux_next[i,j,k] = damp * (2ux0 - ux_previous[i,j,k] + dt2 * ax)
            uy_next[i,j,k] = damp * (2uy0 - uy_previous[i,j,k] + dt2 * ay)
            uz_next[i,j,k] = damp * (2uz0 - uz_previous[i,j,k] + dt2 * az)
        end
    end

    state.ux_previous, state.ux, state.ux_next =
        state.ux, state.ux_next, state.ux_previous
    state.uy_previous, state.uy, state.uy_next =
        state.uy, state.uy_next, state.uy_previous
    state.uz_previous, state.uz, state.uz_next =
        state.uz, state.uz_next, state.uz_previous
    state.step += 1
    state.time += state.dt
    return state
end

"""
    add_ricker_source!(state, index; f0, t0=1.5/f0, amplitude=1)

Add a vertical point-force acceleration to the newly computed displacement
level. Call this immediately after `step_elastic_wave_3d!`.
"""
function add_ricker_source!(
    state::ElasticWaveState3D,
    index::CartesianIndex{3};
    f0::Real,
    t0::Real=1.5 / f0,
    amplitude::Real=1,
    component::Symbol=:z,
)
    checkbounds(state.uz, index)
    state.solid[index] || throw(ArgumentError("source must be inside elastic material"))
    τ = state.time - t0
    a = π * f0 * τ
    ricker = (1 - 2a^2) * exp(-a^2)
    field = component === :x ? state.ux :
            component === :y ? state.uy :
            component === :z ? state.uz :
            throw(ArgumentError("component must be :x, :y or :z"))
    field[index] += eltype(field)(amplitude * state.dt^2 * ricker)
    return state
end

displacement_magnitude(state::ElasticWaveState3D) =
    @. sqrt(state.ux^2 + state.uy^2 + state.uz^2)

end
