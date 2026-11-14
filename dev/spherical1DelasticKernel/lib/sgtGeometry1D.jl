# dev/spherical1DelasticKernel/lib/sgtGeometry1D.jl
#
# Reconstruction "on the fly" of the SGT at any point Q for a source S and a
# receiver R anywhere on (or in) the planet, from the pole disks of
# lib/sgtDisk1D.jl.  Everything is done with explicit 3×3 rotation matrices
# in Cartesian coordinates (no hand-written azimuth formulas):
#
#   1. Rs rotates the geographic frame so that S is on the north pole
#      (Rs ŝ = ẑ); the azimuth about S is left arbitrary and handled below.
#   2. the source (MT or force), given in the local geographic basis at S,
#      is taken to Cartesian, rotated by Rs and read in the pole-local basis
#      (r̂, θ̂, φ̂) at φ = 0 = (ẑ, x̂, ŷ);
#   3. Q is rotated by Rs → (θ', φ') in the source frame; the disk gives the
#      field in the local basis at (θ', φ') (rotation of the source by -φ'
#      about the pole, `poleField`);
#   4. that local basis is mapped back to Q's local geographic basis by
#      T = E_Q' Rs' E'(θ', φ').
#
# Geographic convention: colatitude θ, longitude φ; local basis (r̂, θ̂, φ̂) =
# (up, south, east).  A moment tensor given in (r, θ, φ) = (Mrr, Mrθ, ...)
# is the Harvard/GCMT (r, t, p) convention.

using LinearAlgebra

sphUnit(θ, φ) = [sin(θ) * cos(φ), sin(θ) * sin(φ), cos(θ)]

"Columns r̂, θ̂, φ̂ (Cartesian) at colatitude θ, longitude φ (radians)."
localBasis(θ, φ) = hcat(sphUnit(θ, φ),
                        [cos(θ) * cos(φ), cos(θ) * sin(φ), -sin(θ)],
                        [-sin(φ), cos(φ), 0.0])

_Rz(a) = [cos(a) -sin(a) 0.0; sin(a) cos(a) 0.0; 0.0 0.0 1.0]
_Ry(a) = [cos(a) 0.0 sin(a); 0.0 1.0 0.0; -sin(a) 0.0 cos(a)]

"Rotation taking the point (θ, φ) to the north pole: Ry(-θ)·Rz(-φ)."
poleRotation(θ, φ) = _Ry(-θ) * _Rz(-φ)

"(θ, φ) (radians) of a Cartesian vector."
function sphAngles(v)
    n = norm(v)
    θ = acos(clamp(v[3] / n, -1.0, 1.0))
    φ = atan(v[2], v[1])
    return θ, φ
end

"""
    PointGeo(r_km, θ, φ) -- a point inside the planet (θ colatitude, φ longitude, radians).
"""
struct PointGeo
    r::Float64
    θ::Float64
    φ::Float64
end
PointGeo(; r, θ, φ) = PointGeo(Float64(r), Float64(θ), Float64(φ))
"PointGeo from depth, latitude and longitude in degrees."
pointFromLatLon(R, depth, latdeg, londeg) = PointGeo(R - depth, deg2rad(90 - latdeg), deg2rad(londeg))

"""
    sourceToPole(family, a, S) -> (Rs, apole)

Rs = poleRotation(S) and the source components `a` (geographic local basis
at S) expressed in the pole-local basis (r̂, θ̂, φ̂) at φ = 0.
"""
function sourceToPole(family::Symbol, a, S::PointGeo)
    Rs = poleRotation(S.θ, S.φ)
    E = localBasis(S.θ, S.φ)
    P0 = localBasis(0.0, 0.0)           # (ẑ, x̂, ŷ)
    B = P0' * Rs * E                    # local(S) → pole-local
    apole = family === :moment ? mtVector(B * mtMatrix(a) * B') : Tuple(B * collect(a))
    return Rs, apole
end

"""
    frameOf(Rs, Q) -> (θ'deg, φ', T)

Q's colatitude (degrees) and azimuth in the frame of Rs, and the matrix T
taking components in the local basis at (θ', φ') of that frame to Q's local
geographic basis.
"""
function frameOf(Rs, Q::PointGeo)
    q = Rs * sphUnit(Q.θ, Q.φ)
    θp, φp = sphAngles(q)
    T = localBasis(Q.θ, Q.φ)' * Rs' * localBasis(θp, φp)
    return rad2deg(θp), φp, T
end

"""
    sgtAt(disk, iω, S, a, Q; which=:strain, order=5)

Field at Q (geographic local basis) of the source `a` (components in the
geographic local basis at S) located at S, from a pole disk computed for
rs = S.r and containing the radius Q.r.  `which = :strain` (3×3) or `:disp`
(3-vector).  For the receiver side, S is the receiver, `a` the receiver
component direction (unit force) and the result is the RSGT (strain at Q
due to a unit force at R ≡ displacement at R due to a moment at Q).
"""
function sgtAt(disk, iω::Integer, S::PointGeo, a, Q::PointGeo; which::Symbol=:strain, order::Integer=5)
    abs(disk.rs_km - S.r) <= 1e-6 * S.r || error("disk computed for rs = $(disk.rs_km), point at $(S.r)")
    Rs, ap = sourceToPole(disk.family, a, S)
    θp, φp, T = frameOf(Rs, Q)
    ir = diskRadiusIndex(disk, Q.r)
    F = poleField(disk, iω, ap, ir, θp, φp; which=which, order=order)
    return which === :strain ? T * F * T' : T * F
end

"Unit vector (local geographic basis at R) of the components :Z, :N, :E, or :R/:T towards/away from S."
function receiverDirection(comp::Symbol, S::PointGeo, Rcv::PointGeo)
    comp === :Z && return (1.0, 0.0, 0.0)
    comp === :N && return (0.0, -1.0, 0.0)
    comp === :E && return (0.0, 0.0, 1.0)
    # radial = along the great circle, away from S; transverse = radial × up... (r̂ × L̂ ... ) see note
    E = localBasis(Rcv.θ, Rcv.φ)
    s = sphUnit(S.θ, S.φ); rr = sphUnit(Rcv.θ, Rcv.φ)
    t = cross(s, rr); t ./= norm(t)               # normal of the S-R plane
    L = cross(t, rr)                              # horizontal, pointing away from S
    if comp === :R
        return Tuple(E' * L)
    elseif comp === :T
        return Tuple(E' * t)                      # = φ̂ when S is on the pole (DSM / p6 convention)
    end
    error("component must be :Z, :N, :E, :R or :T")
end

"""
    sgtAtInterp(itp, S, a, Q; which=:strain, side=:upper)

`sgtAt` at any radius Q.r inside the disk radii, through the interpolants
`itp = diskInterpolants(disk, iω)` (lib/sgtInterp1D.jl).  `side` resolves a
point exactly on a discontinuity.
"""
function sgtAtInterp(itp, S::PointGeo, a, Q::PointGeo; which::Symbol=:strain, side::Symbol=:upper)
    disk = itp.disk
    abs(disk.rs_km - S.r) <= 1e-6 * S.r || error("disk computed for rs = $(disk.rs_km), point at $(S.r)")
    Rs, ap = sourceToPole(disk.family, a, S)
    θp, φp, T = frameOf(Rs, Q)
    F = poleFieldAt(itp, ap, Q.r, θp, φp; which=which, side=side)
    return which === :strain ? T * F * T' : T * F
end
