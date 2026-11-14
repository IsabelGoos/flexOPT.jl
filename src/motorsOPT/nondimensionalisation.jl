"""
Utilities for constructing OPT recipes on the canonical grid `Delta = 1` while
keeping the physical grid spacing in the model scaling.

For

    rho d_t^2 u_i - d_j(C_ijkl d_l u_k) = f_i,

with `x_j = h_j xhat_j`, `t = tau that`, and `u = U uhat`, the dimensionless
elasticity tensor is

    Chat_ijkl = C_ijkl * tau^2 / (rho_ref * h_j * h_l).

The two spatial derivative directions are deliberately retained.  This is
essential for anisotropy and for grids whose spacings differ by direction;
in that case scaled scalar Lame parameters are not sufficient.
"""

function _physical_delta(physicalDelta)
    values = Float64.(Tuple(physicalDelta))
    length(values) >= 2 ||
        throw(ArgumentError("physicalDelta must contain spatial steps followed by dt"))
    all(isfinite, values) && all(>(0), values) ||
        throw(ArgumentError("all physicalDelta values must be finite and positive"))
    return values
end

function _positive_reference_density(rho, reference_density)
    value = if reference_density === :auto
        rho isa Real ? Float64(rho) : maximum(Float64.(rho))
    else
        Float64(reference_density)
    end
    isfinite(value) && value > 0 ||
        throw(ArgumentError("reference_density must be finite and positive"))
    all(x -> isfinite(x) && x > 0, rho isa Real ? (rho,) : rho) ||
        throw(ArgumentError("density must be finite and strictly positive"))
    return value
end

"""
    opt_nondimensionalization(physicalDelta; reference_density, reference_field=1)

Return the complete coordinate and field scaling. `physicalDelta` is ordered
as `(dx1, dx2, ..., dt)`. Recipes use `recipeDelta == (1, ..., 1)`.
"""
function opt_nondimensionalization(
    physicalDelta;
    reference_density::Real,
    reference_field::Real=1.0,
)
    delta = _physical_delta(physicalDelta)
    rho0 = Float64(reference_density)
    u0 = Float64(reference_field)
    rho0 > 0 || throw(ArgumentError("reference_density must be positive"))
    u0 > 0 || throw(ArgumentError("reference_field must be positive"))
    return (
        physicalDelta=delta,
        spatial_scales=delta[1:end-1],
        time_scale=delta[end],
        recipeDelta=ntuple(_ -> 1.0, length(delta)),
        reference_density=rho0,
        reference_field=u0,
        body_force_scale=delta[end]^2 / (rho0 * u0),
        displacement_scale=u0,
        velocity_scale=u0 / delta[end],
        acceleration_scale=u0 / delta[end]^2,
    )
end

"""
    isotropic_elasticity_tensor(lambda, mu, dimension)

Build `C[i,j,k,l,...] = lambda delta_ij delta_kl +
mu (delta_ik delta_jl + delta_il delta_jk)`. Spatial field dimensions, if
present, follow the four tensor dimensions.
"""
function isotropic_elasticity_tensor(lambda, mu, dimension::Integer)
    dimension in (1, 2, 3) || throw(ArgumentError("dimension must be 1, 2, or 3"))
    lambda_is_scalar = lambda isa Real
    mu_is_scalar = mu isa Real
    lambda_is_scalar == mu_is_scalar ||
        throw(ArgumentError("lambda and mu must both be scalars or both be arrays"))
    if !lambda_is_scalar
        size(lambda) == size(mu) || throw(DimensionMismatch("lambda and mu sizes differ"))
    end
    trailing = lambda_is_scalar ? () : size(lambda)
    T = promote_type(lambda_is_scalar ? typeof(lambda) : eltype(lambda),
                     mu_is_scalar ? typeof(mu) : eltype(mu), Float64)
    C = zeros(T, dimension, dimension, dimension, dimension, trailing...)
    tail = ntuple(_ -> Colon(), length(trailing))
    for i in 1:dimension, j in 1:dimension, k in 1:dimension, l in 1:dimension
        a = (i == j && k == l) ? 1.0 : 0.0
        b = ((i == k && j == l) ? 1.0 : 0.0) +
            ((i == l && j == k) ? 1.0 : 0.0)
        if isempty(trailing)
            C[i,j,k,l] = a * lambda + b * mu
        else
            @views C[(i,j,k,l,tail...)...] .= a .* lambda .+ b .* mu
        end
    end
    return C
end

"""
    nondimensionalize_elasticity_tensor(rho, C; physicalDelta,
        reference_density=:auto, reference_field=1)

Directionally scale a general elasticity tensor whose first four dimensions
are `(i,j,k,l)`. The remaining dimensions are arbitrary model-grid dimensions.
"""
function nondimensionalize_elasticity_tensor(
    rho,
    C::AbstractArray;
    physicalDelta,
    reference_density=:auto,
    reference_field::Real=1.0,
)
    delta = _physical_delta(physicalDelta)
    dimension = length(delta) - 1
    ndims(C) >= 4 || throw(DimensionMismatch("C needs four leading tensor dimensions"))
    all(size(C, d) == dimension for d in 1:4) ||
        throw(DimensionMismatch("the four leading C dimensions must equal $dimension"))
    rho0 = _positive_reference_density(rho, reference_density)
    scaling = opt_nondimensionalization(delta;
        reference_density=rho0, reference_field=reference_field)
    # `float` preserves Complex inputs (for frequency-domain attenuation),
    # whereas `Float64.(C)` silently rejects their imaginary parts.
    Chat = float.(C)
    tail = ntuple(_ -> Colon(), ndims(Chat) - 4)
    tau2 = scaling.time_scale^2
    for i in 1:dimension, j in 1:dimension, k in 1:dimension, l in 1:dimension
        factor = tau2 / (rho0 * scaling.spatial_scales[j] * scaling.spatial_scales[l])
        if isempty(tail)
            Chat[i,j,k,l] *= factor
        else
            @views Chat[(i,j,k,l,tail...)...] .*= factor
        end
    end
    return (rho=Float64.(rho) ./ rho0, C=Chat, scaling=scaling)
end

"""Scale a physical body-force density into the dimensionless OPT equation."""
nondimensionalize_body_force(force, scaling) = force .* scaling.body_force_scale

"""
Scale a moment tensor `M[i,j,...]` used through `d_j M_ij`.  Its scaling is
directional because the divergence differentiates along direction `j`.
"""
function nondimensionalize_moment_tensor(moment::AbstractArray, scaling)
    dimension = length(scaling.spatial_scales)
    ndims(moment) >= 2 || throw(DimensionMismatch("moment needs leading (i,j) dimensions"))
    size(moment, 1) == dimension && size(moment, 2) == dimension ||
        throw(DimensionMismatch("moment leading dimensions must equal $dimension"))
    result = Float64.(moment)
    tail = ntuple(_ -> Colon(), ndims(result) - 2)
    for i in 1:dimension, j in 1:dimension
        factor = scaling.time_scale^2 /
            (scaling.reference_density * scaling.reference_field * scaling.spatial_scales[j])
        @views result[(i,j,tail...)...] .*= factor
    end
    return result
end

"""
    prepare_nondimensional_recipe(parameters; physicalDelta)

Copy recipe parameters and replace only their numerical recipe spacing by
ones. The physical spacing remains explicit in the returned metadata and must
be used to scale the model and sources.
"""
function prepare_nondimensional_recipe(parameters::AbstractDict; physicalDelta)
    delta = _physical_delta(physicalDelta)
    recipe_parameters = copy(parameters)
    if keytype(recipe_parameters) <: AbstractString
        recipe_parameters["Delta"] = collect(ntuple(_ -> 1.0, length(delta)))
        recipe_parameters["Δ"] = collect(ntuple(_ -> 1.0, length(delta)))
    else
        recipe_parameters[:Delta] = collect(ntuple(_ -> 1.0, length(delta)))
        recipe_parameters[:Δ] = collect(ntuple(_ -> 1.0, length(delta)))
    end
    return (parameters=recipe_parameters, physicalDelta=delta,
            recipeDelta=ntuple(_ -> 1.0, length(delta)))
end

"""
    prepare_nondimensional_elasticity(rho, lambda, mu; physicalDelta, ...)

Convenience entry point for isotropic physical input. The returned `C` is a
general directionally scaled tensor and is therefore also suitable as the
starting point for anisotropic models.
"""
function prepare_nondimensional_elasticity(
    rho, lambda, mu;
    physicalDelta,
    reference_density=:auto,
    reference_field::Real=1.0,
)
    dimension = length(_physical_delta(physicalDelta)) - 1
    C = isotropic_elasticity_tensor(lambda, mu, dimension)
    result = nondimensionalize_elasticity_tensor(rho, C;
        physicalDelta=physicalDelta,
        reference_density=reference_density,
        reference_field=reference_field)
    return merge(result, (
        models=(result.rho, elasticity_component_models(result.C)...),
        material_layout=(:rho, :C_ijkl),
    ))
end

"""
    elasticity_component_models(C)

Flatten only the four tensor axes in Julia column-major order, returning one
scalar or spatial array per `C[i,j,k,l]`. This is the material ordering used by
the tensor-valued elastic famous equations.
"""
function elasticity_component_models(C::AbstractArray)
    dimension = size(C, 1)
    all(size(C, d) == dimension for d in 1:4) ||
        throw(DimensionMismatch("the first four dimensions of C must agree"))
    tail = ntuple(_ -> Colon(), ndims(C) - 4)
    return Tuple(begin
        component = @view C[(I[1], I[2], I[3], I[4], tail...)...]
        isempty(tail) ? component[] : copy(component)
    end for I in CartesianIndices(ntuple(_ -> dimension, 4)))
end
