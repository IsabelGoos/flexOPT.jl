# dev/raytheory/lib/rayModels.jl
#
# Media for ray theory. Only the Julia standard library is used here, so the
# whole ray core runs without flexOPT's src/ (planet1D is NOT loaded: the
# DSM .psv.model reader below is a 30-line stand-alone copy of its format).
#
# Two families:
#
#  * VerticalModel{G}: velocity depends on one "vertical" coordinate ζ only,
#    with ζ increasing UPWARD.
#       Spherical1D : ζ = r (km), polynomial in r/L with L = planet radius
#                     (exactly the DSM convention: v = c0 + c1 x + c2 x² + c3 x³)
#       Cartesian2D/3D (laterally homogeneous): ζ = -z (z = depth, km)
#    The same object serves the τ-p integrals (taup1D.jl) and the ray tracer.
#
#  * PiecewiseMedium{N}: genuinely 2-D/3-D Cartesian medium made of regions
#    separated by planar (possibly dipping) interfaces, each region with its
#    own v(x) = v0 + g⋅(x - x0) for P and S, or any user function.
#
# Coordinates in Cartesian geometry: x = (x, z) in 2-D and (x, y, z) in 3-D,
# z = depth positive DOWN, free surface at z = 0.

using LinearAlgebra, Printf

# ----------------------------------------------------------------------------
# geometry tags
# ----------------------------------------------------------------------------
abstract type RayGeometry end
"1-D spherical planet, rays in the great-circle plane: coordinates (r, θ)."
struct Spherical1D <: RayGeometry end
"2-D Cartesian (x, z), z depth positive down."
struct Cartesian2D <: RayGeometry end
"3-D Cartesian (x, y, z), z depth positive down."
struct Cartesian3D <: RayGeometry end
ndims_geom(::Spherical1D) = 2
ndims_geom(::Cartesian2D) = 2
ndims_geom(::Cartesian3D) = 3

const WAVE_MODES = (:P, :S)

# ----------------------------------------------------------------------------
# Vertical (1-D) models
# ----------------------------------------------------------------------------
"""
One zone of a VerticalModel: bottom < ζ < top, velocities (km/s) and density
as cubic polynomials in ξ = ζ/L. vs ≡ 0 marks a fluid zone.
"""
struct VZone
    bottom::Float64
    top::Float64
    vp::NTuple{4,Float64}
    vs::NTuple{4,Float64}
    rho::NTuple{4,Float64}
end

"""
Boundary between zone k and zone k+1 (k = 0: bottom of the model, k = nzone:
the free surface). `name` is a Symbol used by phase events (:surface, :cmb,
:icb, :moho, :d410 ...), `discontinuous` false for a mere change of
polynomial, `fluidSolid` true when one side is fluid (such boundaries are
only crossed when a phase event says so, as in TauP: "P" never enters the core).
"""
struct VBoundary
    zeta::Float64
    name::Symbol
    depth::Float64
    discontinuous::Bool
    fluidSolid::Bool
end

struct VerticalModel{G<:RayGeometry}
    geom::G
    name::String
    L::Float64                 # polynomial scale (planet radius, or 1)
    surface::Float64           # ζ of the free surface
    zones::Vector{VZone}       # bottom → top
    boundaries::Vector{VBoundary}  # length nzone+1: boundaries[k+1] is the top of zone k (k=0: bottom)
end

@inline polyval4(c::NTuple{4,Float64}, ξ) = c[1] + ξ * (c[2] + ξ * (c[3] + ξ * c[4]))
@inline polyder4(c::NTuple{4,Float64}, ξ) = c[2] + ξ * (2c[3] + ξ * 3c[4])

nzones(m::VerticalModel) = length(m.zones)

"velocity of `mode` in zone `k` at ζ (no check that ζ is inside the zone: zones are extrapolated smoothly)"
@inline function zoneVelocity(m::VerticalModel, k::Int, ζ, mode::Symbol)
    c = mode === :P ? m.zones[k].vp : m.zones[k].vs
    return polyval4(c, ζ / m.L)
end
"dv/dζ in zone k"
@inline function zoneVelocityDeriv(m::VerticalModel, k::Int, ζ, mode::Symbol)
    c = mode === :P ? m.zones[k].vp : m.zones[k].vs
    return polyder4(c, ζ / m.L) / m.L
end
isFluid(m::VerticalModel, k::Int) = all(iszero, m.zones[k].vs)

"Zone containing ζ; at a boundary `side = :above` or `:below` picks the side."
function zoneAt(m::VerticalModel, ζ; side::Symbol=:below)
    n = nzones(m)
    for k in 1:n
        z = m.zones[k]
        if side === :below
            (z.bottom < ζ <= z.top) && return k
        else
            (z.bottom <= ζ < z.top) && return k
        end
    end
    ζ <= m.zones[1].bottom && return 1
    return n
end

depthOf(m::VerticalModel{Spherical1D}, ζ) = m.surface - ζ
depthOf(m::VerticalModel, ζ) = -ζ           # Cartesian: ζ = -z
zetaOfDepth(m::VerticalModel{Spherical1D}, d) = m.surface - d
zetaOfDepth(m::VerticalModel, d) = -d

"""
    nameBoundaries(...)

Automatic names: :surface (top), :center/:bottom (k=0), :cmb (solid above
fluid), :icb (fluid above solid), :moho (the deepest discontinuity shallower
than 90 km with a P jump > 0.4 km/s), and every discontinuity is also
reachable as :d<depth> (e.g. :d400, :d670) — see `boundaryMatches`.
A zone boundary is a discontinuity only if vp, vs or ρ jumps by more than
`jumpTol` (1e-3 km/s, g/cm³): PREM's 600 and 771 km boundaries have
~4e-3 km/s steps that are mere rounding of the published coefficients and
would otherwise create (physically meaningless) total-reflection branches.
"""
function buildBoundaries(geom, L, surface, zones::Vector{VZone}; mohoMax=90.0, jumpTol=1e-3)
    n = length(zones)
    dep(ζ) = geom isa Spherical1D ? surface - ζ : -ζ
    b = VBoundary[]
    push!(b, VBoundary(zones[1].bottom, geom isa Spherical1D ? :center : :bottom, dep(zones[1].bottom), false, false))
    for k in 1:n-1
        ζ = zones[k].top
        lo, hi = zones[k], zones[k+1]
        ξ = ζ / L
        jump = max(abs(polyval4(lo.vp, ξ) - polyval4(hi.vp, ξ)),
                   abs(polyval4(lo.vs, ξ) - polyval4(hi.vs, ξ)),
                   abs(polyval4(lo.rho, ξ) - polyval4(hi.rho, ξ)))
        fl(z) = all(iszero, z.vs)
        fs = fl(lo) != fl(hi)
        name = fs ? (fl(lo) ? :cmb : :icb) : Symbol("d", round(Int, dep(ζ)))
        push!(b, VBoundary(ζ, name, dep(ζ), jump > jumpTol || fs, fs))
    end
    push!(b, VBoundary(zones[n].top, :surface, dep(zones[n].top), true, false))
    # moho: deepest strong P discontinuity above mohoMax km
    best = 0
    for k in 2:n
        bk = b[k]
        if bk.discontinuous && !bk.fluidSolid && 0 < bk.depth <= mohoMax
            ξ = bk.zeta / L
            jump = polyval4(zones[k-1].vp, ξ) - polyval4(zones[k].vp, ξ)  # below minus above
            (jump > 0.4 && best == 0) && (best = k)   # k increases upward: first hit = deepest
        end
    end
    if best > 0
        bb = b[best]; b[best] = VBoundary(bb.zeta, :moho, bb.depth, true, false)
    end
    return b
end

"""
    boundaryMatches(m, k, name; tol=25.0)

True if boundary index k (1-based into m.boundaries) answers to `name`:
its own name, or :d<depth> within `tol` km of its depth (so :d410 finds
PREM's 400-km discontinuity and :d660 its 670), or :moho/:cmb/:icb/:surface.
"""
function boundaryMatches(m::VerticalModel, k::Int, name::Symbol; tol=25.0)
    b = m.boundaries[k]
    b.name === name && return true
    s = String(name)
    if startswith(s, "d") && b.discontinuous
        d = tryparse(Float64, s[2:end])
        d === nothing && return false
        # pick the CLOSEST discontinuity only
        best = argmin(i -> (m.boundaries[i].discontinuous && i != 1 && i != length(m.boundaries)) ?
                            abs(m.boundaries[i].depth - d) : Inf, eachindex(m.boundaries))
        return best == k && abs(b.depth - d) <= tol
    end
    return false
end

function Base.show(io::IO, m::VerticalModel{G}) where {G}
    println(io, "VerticalModel{$G} \"$(m.name)\": $(nzones(m)) zones, surface ζ = $(m.surface)")
    for (k, b) in enumerate(m.boundaries)
        b.discontinuous || k == 1 || continue
        @printf(io, "   %-10s depth %8.1f km%s\n", b.name, b.depth, b.fluidSolid ? "  (fluid-solid)" : "")
    end
end

# ---------------------------------------------------------------- DSM reader
"""
    readDSMModel(path; isotropic=:voigt, name=basename(path))

Read a DSM `.psv.model` file (flexOPT dataInput/model/*/*.psv.model) into a
VerticalModel{Spherical1D}. Transverse isotropy is reduced to an isotropic
equivalent:
  :voigt      vp² = (vpv² + 4vph²)/5, vs² = (2vsv² + vsh²)/3 (pointwise, then
              re-fitted by a cubic in each zone — exact when vpv == vph)
  :vertical   (vpv, vsv)      :horizontal  (vph, vsh)
The planet radius L is the top of the model (DSM convention: polynomials in r/L).
Zones thinner than `minThickness` km are dropped (zero-thickness zones trap rays).
"""
function readDSMModel(path::AbstractString; isotropic::Symbol=:voigt, name=splitext(splitext(basename(path))[1])[1], minThickness=1e-3)
    nums = Float64[]
    for line in eachline(path)
        s = strip(line)
        (isempty(s) || startswith(s, "#") || startswith(lowercase(s), "end")) && continue
        s = split(s, '#')[1]
        for tok in split(s)
            v = tryparse(Float64, tok)
            v === nothing || push!(nums, v)
        end
    end
    nz = Int(nums[1])
    length(nums) == 1 + 28nz || error("$path: expected DSM format with 28 numbers per zone, got $(length(nums)-1) for $nz zones")
    raw = reshape(nums[2:end], 28, nz)
    L = maximum(raw[2, :])
    zones = VZone[]
    for k in 1:nz
        c = raw[:, k]
        rb, rt = c[1], c[2]
        ρ = Tuple(c[3:6]); vpv = Tuple(c[7:10]); vph = Tuple(c[11:14])
        vsv = Tuple(c[15:18]); vsh = Tuple(c[19:22])
        vp, vs = if isotropic === :vertical
            vpv, vsv
        elseif isotropic === :horizontal
            vph, vsh
        elseif isotropic === :voigt
            (vpv == vph ? vpv : refitCubic(ξ -> sqrt((polyval4(vpv, ξ)^2 + 4polyval4(vph, ξ)^2) / 5), rb / L, rt / L)),
            (vsv == vsh ? vsv : refitCubic(ξ -> sqrt((2polyval4(vsv, ξ)^2 + polyval4(vsh, ξ)^2) / 3), rb / L, rt / L))
        else
            error("isotropic must be :voigt, :vertical or :horizontal")
        end
        rt - rb > minThickness && push!(zones, VZone(rb, rt, vp, vs, ρ))   # drop empty zones (TAYAKr has some)
    end
    for k in 1:length(zones)-1       # close the gaps left by dropped zones
        zones[k].top != zones[k+1].bottom && (zones[k] = VZone(zones[k].bottom, zones[k+1].bottom, zones[k].vp, zones[k].vs, zones[k].rho))
    end
    return VerticalModel(Spherical1D(), String(name), L, L, zones, buildBoundaries(Spherical1D(), L, L, zones))
end

"least-squares cubic in ξ on [a, b] (Chebyshev-spaced samples)"
function refitCubic(f, a, b)
    ξs = [(a + b) / 2 + (b - a) / 2 * cos(π * (j + 0.5) / 40) for j in 0:39]
    A = [ξ^p for ξ in ξs, p in 0:3]
    c = A \ f.(ξs)
    return (c[1], c[2], c[3], c[4])
end

# ---------------------------------------------------------------- Cartesian layered
"""
    layeredCartesian(geom; tops, vp, vs, gp=zeros, gs=zeros, rho=ones, zmax=1e4, names=nothing)

Laterally homogeneous Cartesian model (Cartesian2D() or Cartesian3D()):
layer i spans tops[i] ≤ z < tops[i+1] (tops[1] = 0 is the free surface, the
last layer goes down to zmax) with v(z) = v[i] + g[i] (z - tops[i]).
Interfaces are named :i1, :i2, ... (top of layer 2, 3, ...) unless `names`.
"""
function layeredCartesian(geom::Union{Cartesian2D,Cartesian3D}; tops, vp, vs,
        gp=zeros(length(tops)), gs=zeros(length(tops)), rho=ones(length(tops)),
        zmax=1e4, names=nothing, name="layered")
    n = length(tops)
    @assert tops[1] == 0 "the first layer must start at the free surface z = 0"
    zones = VZone[]
    for i in n:-1:1   # bottom → top in ζ = -z
        zt = tops[i]
        zb = i < n ? tops[i+1] : zmax
        # v = v_i + g_i (z - zt) = (v_i - g_i zt) + g_i z = (v_i - g_i zt) - g_i ζ
        P = (vp[i] - gp[i] * zt, -gp[i], 0.0, 0.0)
        S = (vs[i] - gs[i] * zt, -gs[i], 0.0, 0.0)
        push!(zones, VZone(-zb, -zt, P, S, (Float64(rho[i]), 0.0, 0.0, 0.0)))
    end
    b = buildBoundaries(geom, 1.0, 0.0, zones; mohoMax=-1.0)
    # rename internal boundaries :i1 (shallowest) ...
    nb = length(b)
    for k in 2:nb-1
        j = nb - k                 # 1 for the shallowest internal boundary
        nm = names === nothing ? Symbol("i", j) : names[j]
        bb = b[k]
        b[k] = VBoundary(bb.zeta, nm, bb.depth, bb.discontinuous, bb.fluidSolid)
    end
    return VerticalModel(geom, name, 1.0, 0.0, zones, b)
end

# ----------------------------------------------------------------------------
# Genuinely 2-D / 3-D Cartesian media
# ----------------------------------------------------------------------------
"""
Planar interface {x : n⋅(x - x0) = 0}, n unit normal pointing to the side
whose region index is HIGHER (i.e. downward for a stack ordered from the
surface). `name` is used by phase events.
"""
struct PlaneInterface{N}
    x0::NTuple{N,Float64}
    n::NTuple{N,Float64}
    name::Symbol
    freeSurface::Bool
end
function PlaneInterface(x0, n, name; freeSurface=false)
    N = length(x0)
    nn = collect(Float64, n); nn ./= norm(nn)
    return PlaneInterface{N}(Tuple(Float64.(x0)), Tuple(nn), name, freeSurface)
end
signedDistance(I::PlaneInterface, x) = sum(I.n[i] * (x[i] - I.x0[i]) for i in eachindex(I.n))
normalAt(I::PlaneInterface, x) = collect(I.n)

"""
Region of a PiecewiseMedium: velocities given as functions of the position
vector x (km/s) with their gradients. Use `linearRegion` for v0 + g⋅(x - x0).
"""
struct MRegion{FP,FS,GP,GS}
    vp::FP
    vs::FS
    gradvp::GP
    gradvs::GS
end

"v(x) = v0 + g⋅(x - x0)"
function linearRegion(x0, vp0, gp, vs0, gs)
    x0 = collect(Float64, x0); gp = collect(Float64, gp); gs = collect(Float64, gs)
    MRegion(x -> vp0 + dot(gp, x .- x0), x -> vs0 + dot(gs, x .- x0), x -> gp, x -> gs)
end
"region from functions only; gradients by central differences (step h km)"
function functionRegion(vp, vs; h=1e-3)
    fd(f) = x -> [(f(x .+ h .* (1:length(x) .== i)) - f(x .- h .* (1:length(x) .== i))) / 2h for i in 1:length(x)]
    MRegion(vp, vs, fd(vp), fd(vs))
end

"""
Stack of regions separated by non-crossing planes, ordered from the free
surface down: interfaces[1] is the free surface, region i lies between
interfaces[i] and interfaces[i+1] (the last region is unbounded below).
"""
struct PiecewiseMedium{N,G<:RayGeometry,R}
    geom::G
    interfaces::Vector{PlaneInterface{N}}
    regions::Vector{R}          # MRegion (isotropic) or ARegion (anisotropic, rayChristoffel.jl)
end
function PiecewiseMedium(geom::RayGeometry, interfaces::AbstractVector, regions::AbstractVector)
    N = ndims_geom(geom)
    R = typejoin(map(typeof, regions)...)
    return PiecewiseMedium{N,typeof(geom),R}(geom, Vector{PlaneInterface{N}}(interfaces), Vector{R}(regions))
end

"x = (x, z) or (x, y, z), z depth"
function regionAt(m::PiecewiseMedium, x)
    k = 1
    for i in 2:length(m.interfaces)
        signedDistance(m.interfaces[i], x) > 0 && (k = i)
    end
    return k
end

function Base.show(io::IO, m::PiecewiseMedium{N}) where {N}
    print(io, "PiecewiseMedium{$N}: $(length(m.regions)) regions, interfaces ",
        join(string.(getfield.(m.interfaces, :name)), ", "))
end
