using Symbolics
using StaticArrays
using LinearAlgebra

export FreeSurfacePointSet, CerjanBoundarySpec, BoundaryConditionSet
export famousBoundaryConditions, free_surface_points
export cerjan_padding, boundary_geometry

"""
    FreeSurfacePointSet(points, normals)

Material-grid points where the traction-free condition is applied. Normals are
stored separately from `GeoPoint`, in the same order as `points`, so geographic
coordinates remain immutable and reusable.
"""
struct FreeSurfacePointSet{N,T<:Real}
    points::Vector{CartesianIndex{N}}
    normals::Vector{SVector{N,T}}
    function FreeSurfacePointSet(points::Vector{CartesianIndex{N}},
                                 normals::Vector{SVector{N,T}}) where {N,T<:Real}
        length(points) == length(normals) ||
            throw(DimensionMismatch("free-surface points and normals differ"))
        all(normal -> isapprox(norm(normal), 1; atol=1e-6), normals) ||
            throw(ArgumentError("free-surface normals must be unit vectors"))
        new{N,T}(points, normals)
    end
end

"""
    CerjanBoundarySpec(lower, upper; damping=0.0053)

Number of grid points appended on the lower and upper side of every spatial
axis. For example, a 2D x-z model with no absorbing layer above the physical
free surface can use `lower=(12, 12), upper=(12, 0)`.
"""
struct CerjanBoundarySpec{N,T<:Real}
    lower::NTuple{N,Int}
    upper::NTuple{N,Int}
    damping::T
    function CerjanBoundarySpec(
        lower::NTuple{N,Int},
        upper::NTuple{N,Int};
        damping::T=0.0053,
    ) where {N,T<:Real}
        all(>=(0), lower) && all(>=(0), upper) ||
            throw(ArgumentError("Cerjan widths must be non-negative"))
        damping >= 0 || throw(ArgumentError("Cerjan damping must be non-negative"))
        new{N,T}(lower, upper, damping)
    end
end

Base.@kwdef struct BoundaryConditionSet{F,C,M}
    free_surface::F = nothing
    cerjan::C = nothing
    material_mask::M = nothing
    free_surface_mode::Symbol = :traction
end

cerjan_padding(spec::CerjanBoundarySpec{N}) where {N} =
    permutedims(hcat(collect(spec.lower), collect(spec.upper)))

"""
    famousBoundaryConditions("elasticTractionFree2D")
    famousBoundaryConditions("elasticTractionFree3D")

Return symbolic boundary residuals, fields, material/normal variables and local
coordinates, following the same gallery pattern as `famousEquations`.
"""
function famousBoundaryConditions(name::AbstractString)
    return famousBoundaryCondition(Val(Symbol("bc_" * name)))
end

function famousBoundaryCondition(::Val{:bc_elasticTractionFree2D})
    @variables x z t
    @variables λ(x,z) μ(x,z) ux(x,z,t) uz(x,z,t) nx(x,z) nz(x,z)
    divergence = Differential(x)(ux) + Differential(z)(uz)
    σxx = λ * divergence + 2μ * Differential(x)(ux)
    σzz = λ * divergence + 2μ * Differential(z)(uz)
    σxz = μ * (Differential(z)(ux) + Differential(x)(uz))
    residuals = (
        σxx * nx + σxz * nz,
        σxz * nx + σzz * nz,
    )
    return (
        residuals=residuals,
        fields=(ux, uz),
        vars=(λ, μ, nx, nz),
        normals=(nx, nz),
        coordinates=(x, z, t),
        name=:elasticTractionFree2D,
    )
end

function famousBoundaryCondition(::Val{:bc_elasticTractionFree3D})
    @variables x y z t
    @variables λ(x,y,z) μ(x,y,z) u(x,y,z,t)[1:3]
    @variables n(x,y,z)[1:3]
    coordinates = (x, y, z)
    divergence = sum(Differential(coordinates[k])(u[k]) for k in 1:3)
    σ = Matrix{Any}(undef, 3, 3)
    for i in 1:3, j in 1:3
        σ[i,j] = (i == j ? λ * divergence : 0) +
                 μ * (Differential(coordinates[j])(u[i]) +
                      Differential(coordinates[i])(u[j]))
    end
    residuals = ntuple(i -> sum(σ[i,j] * n[j] for j in 1:3), 3)
    return (
        residuals=residuals,
        fields=Tuple(u),
        vars=(λ, μ, Tuple(n)...),
        normals=Tuple(n),
        coordinates=(x, y, z, t),
        name=:elasticTractionFree3D,
    )
end

function _surface_height_indices(material::AbstractArray{Bool,2})
    nx, nz = size(material)
    indices = Vector{CartesianIndex{2}}(undef, nx)
    for i in 1:nx
        k = findlast(@view material[i, :])
        isnothing(k) && throw(ArgumentError("no material in column $i"))
        indices[i] = CartesianIndex(i, k)
    end
    return indices
end

function _surface_height_indices(material::AbstractArray{Bool,3})
    nx, ny, nz = size(material)
    indices = Matrix{CartesianIndex{3}}(undef, nx, ny)
    for j in 1:ny, i in 1:nx
        k = findlast(@view material[i, j, :])
        isnothing(k) && throw(ArgumentError("no material in column $((i,j))"))
        indices[i,j] = CartesianIndex(i, j, k)
    end
    return indices
end

@inline _difference(values, i, spacing) =
    i == firstindex(values) ? (values[i+1] - values[i]) / spacing :
    i == lastindex(values) ? (values[i] - values[i-1]) / spacing :
    (values[i+1] - values[i-1]) / (2spacing)

function free_surface_points(
    material::AbstractArray{Bool,2},
    spacing::NTuple{2,<:Real},
)
    dx, dz = spacing
    surface = _surface_height_indices(material)
    height = [index[2] * dz for index in surface]
    points = vec(surface)
    normals = SVector{2,Float64}[]
    for i in eachindex(height)
        dhdx = _difference(height, i, dx)
        push!(normals, normalize(SVector(-dhdx, 1.0)))
    end
    return FreeSurfacePointSet(points, normals)
end

function free_surface_points(
    material::AbstractArray{Bool,3},
    spacing::NTuple{3,<:Real},
)
    dx, dy, dz = spacing
    surface = _surface_height_indices(material)
    nx, ny = size(surface)
    height = [surface[i,j][3] * dz for i in 1:nx, j in 1:ny]
    points = vec(surface)
    normals = Vector{SVector{3,Float64}}(undef, length(points))
    linear = LinearIndices(surface)
    for j in 1:ny, i in 1:nx
        dhdx = _difference(@view(height[:,j]), i, dx)
        dhdy = _difference(@view(height[i,:]), j, dy)
        normals[linear[i,j]] = normalize(SVector(-dhdx, -dhdy, 1.0))
    end
    return FreeSurfacePointSet(points, normals)
end

"""
    boundary_geometry(material, spacing; cerjan=nothing)

Build the complete geometry object passed to `numericalOperatorConstruction`.
"""
function boundary_geometry(material::AbstractArray{Bool,N},
                           spacing::NTuple{N,<:Real};
                           cerjan=nothing,
                           free_surface_mode::Symbol=:traction) where {N}
    free = free_surface_points(material, spacing)
    if !isnothing(cerjan)
        cerjan isa CerjanBoundarySpec{N} ||
            throw(DimensionMismatch("Cerjan specification dimension differs"))
    end
    free_surface_mode in (:traction, :dietrich, :pinned_void) ||
        throw(ArgumentError(
            "free_surface_mode must be :traction, :dietrich or :pinned_void",
        ))
    return BoundaryConditionSet(
        free_surface=free,
        cerjan=cerjan,
        material_mask=BitArray(material),
        free_surface_mode=free_surface_mode,
    )
end
