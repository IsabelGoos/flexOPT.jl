module elasticWave2D

using Base.Threads

export ElasticThreePointConfig2D, ElasticWaveState2D
export prepare_elastic_wave_2d, step_elastic_wave_2d!, add_ricker_source!
export displacement_magnitude, elastic_wave_coordinates

Base.@kwdef struct ElasticThreePointConfig2D
    pointsInSpace::Int = 3
    pointsInTime::Int = 3
    supplementaryOrder::Int = 2
    cfl::Float64 = 0.38
end

mutable struct ElasticWaveState2D{T,A<:AbstractMatrix{T}}
    ux_previous::A
    uz_previous::A
    ux::A
    uz::A
    ux_next::A
    uz_next::A
    ρ::A
    λ::A
    μ::A
    σxx::A
    σzz::A
    σxz::A
    material::BitMatrix
    damping::A
    spacing::NTuple{2,T}
    padding::Matrix{Int}
    dt::T
    time::T
    step::Int
    config::ElasticThreePointConfig2D
end

function _validate(config)
    config.pointsInSpace == 3 ||
        throw(ArgumentError("the 2D solver uses exactly three spatial points"))
    config.pointsInTime == 3 ||
        throw(ArgumentError("the 2D solver uses exactly three time levels"))
    config.supplementaryOrder == 2 ||
        throw(ArgumentError("the 2D development case requires supplementaryOrder=2"))
    0 < config.cfl < inv(sqrt(2)) ||
        throw(ArgumentError("cfl must lie in (0, 1/√2)"))
end

function _extend(array, padding)
    lower = padding[1, :]
    nx, nz = size(array)
    out = similar(array, nx + sum(padding[:, 1]), nz + sum(padding[:, 2]))
    @threads for k in axes(out, 2)
        source_k = clamp(k - lower[2], 1, nz)
        for i in axes(out, 1)
            out[i, k] = array[clamp(i - lower[1], 1, nx), source_k]
        end
    end
    out
end

function _extended_material(material, padding)
    out = _extend(material, padding)
    # Padding is elastic continuation of the physical side/bottom boundary.
    # There is normally no upper-z padding because that side is the free surface.
    BitMatrix(out)
end

function _cerjan(shape, padding, strength, ::Type{T}) where T
    result = ones(T, shape)
    lower = padding[1, :]
    physical = (lower[1] + 1):(shape[1] - padding[2, 1]),
               (lower[2] + 1):(shape[2] - padding[2, 2])
    @threads for k in axes(result, 2)
        for i in axes(result, 1)
            dx = i < first(physical[1]) ? first(physical[1]) - i :
                 i > last(physical[1]) ? i - last(physical[1]) : 0
            dz = k < first(physical[2]) ? first(physical[2]) - k :
                 k > last(physical[2]) ? k - last(physical[2]) : 0
            result[i, k] = exp(-T(strength) * T(max(dx, dz)^2))
        end
    end
    result
end

"""
    prepare_elastic_wave_2d(model, spacing; material_mask, boundary_conditions)

Prepare a three-point, three-time-level plane-strain elastic simulation.
Cerjan cells are appended outside the physical model. At material/air faces,
the finite-volume stress flux is zero, which supplies a staircase
traction-free surface. Smooth-normal symbolic conditions remain available
through `famousBoundaryConditions("elasticTractionFree2D")`.
"""
function prepare_elastic_wave_2d(
    model,
    spacing::NTuple{2,<:Real};
    material_mask,
    boundary_conditions,
    config=ElasticThreePointConfig2D(),
    dt=nothing,
    T::Type{<:AbstractFloat}=Float32,
)
    _validate(config)
    padding = isnothing(boundary_conditions.cerjan) ?
              zeros(Int, 2, 2) :
              permutedims(hcat(
                  collect(boundary_conditions.cerjan.lower),
                  collect(boundary_conditions.cerjan.upper),
              ))
    ρ = T.(_extend(model.ρ .* 1_000, padding))
    vp = T.(_extend(model.Vpv .* 1_000, padding))
    vs = T.(_extend(model.Vsv .* 1_000, padding))
    material = _extended_material(BitMatrix(material_mask), padding)
    λ = @. ρ * (vp^2 - 2vs^2)
    μ = @. ρ * vs^2
    λ[.!material] .= 0
    μ[.!material] .= 0
    vmax = maximum(vp[material])
    dx, dz = T.(spacing)
    timestep = isnothing(dt) ?
        T(config.cfl / (vmax * sqrt(inv(dx^2) + inv(dz^2)))) : T(dt)
    timestep > 0 || throw(ArgumentError("dt must be positive"))
    damping = _cerjan(
        size(ρ), padding,
        isnothing(boundary_conditions.cerjan) ? 0.0 :
            boundary_conditions.cerjan.damping,
        T,
    )
    zeroarray = zeros(T, size(ρ))
    ElasticWaveState2D(
        copy(zeroarray), copy(zeroarray), copy(zeroarray), copy(zeroarray),
        copy(zeroarray), copy(zeroarray), ρ, λ, μ,
        copy(zeroarray), copy(zeroarray), copy(zeroarray),
        material, damping, (dx, dz), padding, timestep, zero(T), 0, config,
    )
end

@inline _d1(field, material, i1, k1, i0, k0, h) =
    material[i1, k1] ? (field[i1, k1] - field[i0, k0]) / h : zero(eltype(field))

function step_elastic_wave_2d!(state::ElasticWaveState2D{T}) where T
    (; ux, uz, ux_previous, uz_previous, ux_next, uz_next,
       ρ, λ, μ, σxx, σzz, σxz, material, damping) = state
    dx, dz = state.spacing
    nx, nz = size(ux)
    fill!(σxx, 0); fill!(σzz, 0); fill!(σxz, 0)
    @threads for k in 2:nz-1
        for i in 2:nx-1
            material[i, k] || continue
            ux_x = (ux[i+1,k] - ux[i-1,k]) / (2dx)
            uz_z = (uz[i,k+1] - uz[i,k-1]) / (2dz)
            ux_z = (ux[i,k+1] - ux[i,k-1]) / (2dz)
            uz_x = (uz[i+1,k] - uz[i-1,k]) / (2dx)
            σxx[i,k] = (λ[i,k] + 2μ[i,k]) * ux_x + λ[i,k] * uz_z
            σzz[i,k] = λ[i,k] * ux_x + (λ[i,k] + 2μ[i,k]) * uz_z
            σxz[i,k] = μ[i,k] * (ux_z + uz_x)
        end
    end

    fill!(ux_next, 0); fill!(uz_next, 0)
    dt2 = state.dt^2
    @threads for k in 2:nz-1
        for i in 2:nx-1
            material[i, k] || continue
            # A material/air face has zero stress flux: σ⋅n = 0.
            σxx_r = material[i+1,k] ? (σxx[i,k] + σxx[i+1,k]) / 2 : zero(T)
            σxx_l = material[i-1,k] ? (σxx[i,k] + σxx[i-1,k]) / 2 : zero(T)
            σxz_t = material[i,k+1] ? (σxz[i,k] + σxz[i,k+1]) / 2 : zero(T)
            σxz_b = material[i,k-1] ? (σxz[i,k] + σxz[i,k-1]) / 2 : zero(T)
            σxz_r = material[i+1,k] ? (σxz[i,k] + σxz[i+1,k]) / 2 : zero(T)
            σxz_l = material[i-1,k] ? (σxz[i,k] + σxz[i-1,k]) / 2 : zero(T)
            σzz_t = material[i,k+1] ? (σzz[i,k] + σzz[i,k+1]) / 2 : zero(T)
            σzz_b = material[i,k-1] ? (σzz[i,k] + σzz[i,k-1]) / 2 : zero(T)
            ax = ((σxx_r - σxx_l) / dx + (σxz_t - σxz_b) / dz) / ρ[i,k]
            az = ((σxz_r - σxz_l) / dx + (σzz_t - σzz_b) / dz) / ρ[i,k]
            damp = damping[i,k]
            ux_next[i,k] = damp * (2ux[i,k] - ux_previous[i,k] + dt2 * ax)
            uz_next[i,k] = damp * (2uz[i,k] - uz_previous[i,k] + dt2 * az)
        end
    end
    state.ux_previous, state.ux, state.ux_next =
        state.ux, state.ux_next, state.ux_previous
    state.uz_previous, state.uz, state.uz_next =
        state.uz, state.uz_next, state.uz_previous
    state.step += 1
    state.time += state.dt
    state
end

function add_ricker_source!(state::ElasticWaveState2D, index::CartesianIndex{2};
                            f0, t0=1.5/f0, amplitude=1, component=:z,
                            source_kind::Symbol=:acceleration)
    checkbounds(state.ux, index)
    state.material[index] || throw(ArgumentError("source must be in elastic material"))
    a = π * f0 * (state.time - t0)
    ricker = (1 - 2a^2) * exp(-a^2)
    acceleration = if source_kind === :acceleration
        amplitude * ricker
    elseif source_kind === :force
        # A 2-D simulation represents a unit-thickness slice. `amplitude`
        # is therefore a line force in N/m and rho*dx*dz is the lumped
        # nodal mass per metre out of plane.
        dx, dz = state.spacing
        amplitude * ricker / (state.ρ[index] * dx * dz)
    else
        throw(ArgumentError("source_kind must be :acceleration or :force"))
    end
    value = eltype(state.ux)(state.dt^2 * acceleration)
    field = component === :x ? state.ux :
            component === :z ? state.uz :
            throw(ArgumentError("component must be :x or :z"))
    field[index] += value
    state
end

displacement_magnitude(state::ElasticWaveState2D) = @. hypot(state.ux, state.uz)

function elastic_wave_coordinates(x, z, state::ElasticWaveState2D)
    px, pz = state.padding[1, :]
    dx, dz = state.spacing
    (
        x=range(first(x) - px*dx; step=dx, length=size(state.ux, 1)),
        z=range(first(z) - pz*dz; step=dz, length=size(state.ux, 2)),
    )
end

end
