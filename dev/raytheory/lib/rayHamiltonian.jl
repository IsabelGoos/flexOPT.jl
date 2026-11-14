# dev/raytheory/lib/rayHamiltonian.jl
#
# Layer 1 of the ray theory: Hamiltonian propagation.
#
# A wave branch is a Hamiltonian H(x, p) = 0 with p = ∇T (slowness). The ray
# is the characteristic
#        dx/dσ = ∂H/∂p ,   dp/dσ = -∂H/∂x ,   dT/dσ = p⋅∂H/∂p .
# We normalise H so that p⋅∂H/∂p = 1 on H = 0, i.e. σ IS the travel time T.
#
# Isotropic branch of velocity v:   H = ½ (v² |p|² - 1)
#   Cartesian (x, z) or (x, y, z):  dx/dT = v² p,       dp/dT = -|p|² v ∇v
#   Spherical1D (r, θ; p_r, p_θ):   |p|² = p_r² + p_θ²/r²
#       dr/dT = v² p_r,  dθ/dT = v² p_θ / r²
#       dp_r/dT = -|p|² v ∂_r v + v² p_θ² / r³,   dp_θ/dT = 0
#   p_θ = r sin(i)/v is the spherical ray parameter (s/rad); in a laterally
#   homogeneous Cartesian model p_x = sin(i)/v (s/km) is conserved instead.
#
# Anisotropic branches (qP, qS1, qS2 from the Christoffel matrix) are
# rayChristoffel.jl; they plug into the same tracer through `rayRHS` and
# `interfaceSlowness`.

abstract type AbstractHamiltonian end
"Isotropic elastic/acoustic branches; the mode (:P or :S) is part of the ray state."
struct IsotropicHamiltonian <: AbstractHamiltonian end

# ---------------------------------------------------------------- medium API
# All media answer, for a region index `k`:
#   velocity(med, k, x, mode), gradVelocity(med, k, x, mode)   (in the geometry's coordinates)
#   boundaryIds(med, k)            ids of the boundaries enclosing region k
#   insideDistance(med, k, b, x)   > 0 inside region k, 0 on boundary b
#   outwardNormal(med, k, b, x)    unit normal leaving region k through b (physical frame)
#   across(med, k, b)              region on the other side of b
#   crossDirection(med, k, b)      :down or :up
#   boundaryName / isFreeSurface / isFluidSolid / isDiscontinuous / matchesBoundary

geometry(m::VerticalModel) = m.geom
geometry(m::PiecewiseMedium) = m.geom

zetaOf(::VerticalModel{Spherical1D}, x) = x[1]
zetaOf(::VerticalModel, x) = -x[end]

velocity(m::VerticalModel, k::Int, x, mode::Symbol) = zoneVelocity(m, k, zetaOf(m, x), mode)
function gradVelocity(m::VerticalModel{Spherical1D}, k::Int, x, mode::Symbol)
    return (zoneVelocityDeriv(m, k, x[1], mode), 0.0)
end
function gradVelocity(m::VerticalModel, k::Int, x, mode::Symbol)
    dvdz = -zoneVelocityDeriv(m, k, -x[end], mode)
    return ntuple(i -> i == length(x) ? dvdz : 0.0, length(x))
end
boundaryIds(m::VerticalModel, k::Int) = (k, k + 1)          # bottom, top (indices into m.boundaries)
function insideDistance(m::VerticalModel, k::Int, b::Int, x)
    ζ = zetaOf(m, x)
    return b == k ? ζ - m.boundaries[b].zeta : m.boundaries[b].zeta - ζ
end
crossDirection(m::VerticalModel, k::Int, b::Int) = b == k ? :down : :up
across(m::VerticalModel, k::Int, b::Int) = b == k ? k - 1 : k + 1
function outwardNormal(m::VerticalModel{Spherical1D}, k::Int, b::Int, x)
    return b == k ? [-1.0, 0.0] : [1.0, 0.0]            # (e_r, e_θ) frame
end
function outwardNormal(m::VerticalModel, k::Int, b::Int, x)
    n = zeros(length(x)); n[end] = b == k ? 1.0 : -1.0  # going down = +z
    return n
end
boundaryName(m::VerticalModel, b::Int) = m.boundaries[b].name
isFreeSurface(m::VerticalModel, b::Int) = b == length(m.boundaries)
isModelBottom(m::VerticalModel, b::Int) = b == 1
isFluidSolid(m::VerticalModel, b::Int) = m.boundaries[b].fluidSolid
isDiscontinuous(m::VerticalModel, b::Int) = m.boundaries[b].discontinuous
matchesBoundary(m::VerticalModel, b::Int, name::Symbol) = boundaryMatches(m, b, name)
function initialRegion(m::VerticalModel, x, goingDown::Bool)
    return zoneAt(m, zetaOf(m, x); side=goingDown ? :below : :above)
end

# PiecewiseMedium: region k between interfaces[k] (above) and interfaces[k+1] (below)
velocity(m::PiecewiseMedium, k::Int, x, mode::Symbol) = mode === :P ? m.regions[k].vp(x) : m.regions[k].vs(x)
gradVelocity(m::PiecewiseMedium, k::Int, x, mode::Symbol) = mode === :P ? m.regions[k].gradvp(x) : m.regions[k].gradvs(x)
boundaryIds(m::PiecewiseMedium, k::Int) = k < length(m.interfaces) ? (k, k + 1) : (k,)
insideDistance(m::PiecewiseMedium, k::Int, b::Int, x) = b == k ? signedDistance(m.interfaces[b], x) : -signedDistance(m.interfaces[b], x)
crossDirection(m::PiecewiseMedium, k::Int, b::Int) = b == k ? :up : :down
across(m::PiecewiseMedium, k::Int, b::Int) = b == k ? k - 1 : k + 1
outwardNormal(m::PiecewiseMedium, k::Int, b::Int, x) = b == k ? -normalAt(m.interfaces[b], x) : normalAt(m.interfaces[b], x)
boundaryName(m::PiecewiseMedium, b::Int) = m.interfaces[b].name
isFreeSurface(m::PiecewiseMedium, b::Int) = b == 1
isModelBottom(m::PiecewiseMedium, b::Int) = false
isFluidSolid(m::PiecewiseMedium, b::Int) = false
isDiscontinuous(m::PiecewiseMedium, b::Int) = true
matchesBoundary(m::PiecewiseMedium, b::Int, name::Symbol) = m.interfaces[b].name === name
initialRegion(m::PiecewiseMedium, x, goingDown::Bool) = regionAt(m, x)

# ---------------------------------------------------------------- frames
# physical slowness q (orthonormal frame) <-> canonical momentum p
physicalSlowness(::Spherical1D, x, p) = [p[1], p[2] / x[1]]
canonicalMomentum(::Spherical1D, x, q) = [q[1], q[2] * x[1]]
physicalSlowness(::RayGeometry, x, p) = collect(p)
canonicalMomentum(::RayGeometry, x, q) = collect(q)

# ---------------------------------------------------------------- isotropic RHS
"""
    rayRHS(H, geom, med, k, x, p, mode) -> (dx/dT, dp/dT)
"""
function rayRHS(::IsotropicHamiltonian, ::Spherical1D, med, k, x, p, mode)
    r = x[1]
    v = velocity(med, k, x, mode)
    dv = gradVelocity(med, k, x, mode)[1]
    pr, pθ = p[1], p[2]
    p2 = pr^2 + pθ^2 / r^2
    return [v^2 * pr, v^2 * pθ / r^2], [-p2 * v * dv + v^2 * pθ^2 / r^3, 0.0]
end
function rayRHS(::IsotropicHamiltonian, ::RayGeometry, med, k, x, p, mode)
    v = velocity(med, k, x, mode)
    g = gradVelocity(med, k, x, mode)
    p2 = sum(abs2, p)
    return v^2 .* p, -(p2 * v) .* collect(g)
end

"H residual v|p| - 1 (should stay ≈ 0 along the ray)"
function hamiltonianResidual(::IsotropicHamiltonian, geom, med, k, x, p, mode)
    q = physicalSlowness(geom, x, p)
    return velocity(med, k, x, mode) * norm(q) - 1
end

"""
    interfaceSlowness(H, geom, med, kOut, x, q, n, mode, action) -> q′ or nothing

Layer 3 (interface matching), isotropic case: generalized Snell law. The
tangential slowness q_t = q - (q⋅n) n is conserved; the outgoing normal
component solves the outgoing branch's dispersion relation
v_out² (|q_t|² + q_n²) = 1. `n` is the outward normal of the region being
left (q⋅n > 0); a transmitted ray keeps going along +n, a reflected one
along -n. `nothing` = evanescent (post-critical, or S into a fluid).
"""
function interfaceSlowness(::IsotropicHamiltonian, geom, med, kOut, x, q, n, mode, action)
    v = velocity(med, kOut, x, mode)
    (v <= 0 || !isfinite(v)) && return nothing
    qn = dot(q, n)
    qt = q .- qn .* n
    s2 = 1 / v^2 - sum(abs2, qt)
    s2 < 0 && return nothing
    return qt .+ (action === :transmit ? 1.0 : -1.0) * sqrt(s2) .* n
end
