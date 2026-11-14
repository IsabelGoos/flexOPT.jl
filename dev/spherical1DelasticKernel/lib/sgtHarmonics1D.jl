# dev/spherical1DelasticKernel/lib/sgtHarmonics1D.jl
#
# Angular part of the strain Green's tensors (SGT) for a point source on the
# north pole of a spherically symmetric planet.  Pure functions, no flexOPT.
#
# Conventions (identical to dev/spherical1Delastic/lib/elastic1D.jl):
#   Y_lm = X_lm(θ) e^{imφ}, fully normalised, Condon–Shortley,
#   P = r̂ Y,  B = k⁻¹ ∇₁Y,  C = -k⁻¹ r̂ × ∇₁Y,  k = √(l(l+1)),
#   u = Σ_lm U(r) P + V(r) B + W(r) C .
# Components are always in the local spherical basis (r̂, θ̂, φ̂) and in the
# order
#   displacement  (r, θ, φ)
#   strain        (rr, θθ, φφ, rθ, rφ, θφ)    (tensor components, not Voigt ×2)
#   moment tensor (Mrr, Mrθ, Mrφ, Mθθ, Mθφ, Mφφ)  (pointSource order)
#   force         (fr, fθ, fφ)
#
# Parity.  The meridian plane φ = 0 is a mirror plane of the problem
# (y → -y, i.e. φ̂ → -φ̂).  A source component with an odd number of φ
# indices (Mrφ, Mθφ, fφ) is odd, the others even; on φ = 0 an even source
# only produces the even field components (u_r, u_θ; ε_rr, ε_θθ, ε_φφ, ε_rθ)
# and an odd source only the odd ones (u_φ; ε_rφ, ε_θφ).  This is exactly
# why Zhao & Chevrot (2011) need 20 TSGT and 10 RSGT functions:
#   moment: 4 even × 4 + 2 odd × 2 = 20 strains,  4 × 2 + 2 × 1 = 10 displacements
#   force : 2 even × 4 + 1 odd × 2 = 10 strains,  2 × 2 + 1 × 1 =  5 displacements
# Here the pairs are generated, not hand-written (`sgtPairs`).

const MT_NAMES = ("Mrr", "Mrθ", "Mrφ", "Mθθ", "Mθφ", "Mφφ")
const FORCE_NAMES = ("fr", "fθ", "fφ")
const STRAIN_NAMES = ("rr", "θθ", "φφ", "rθ", "rφ", "θφ")
const DISP_NAMES = ("r", "θ", "φ")
const MT_ODD = (false, false, true, false, true, false)
const FORCE_ODD = (false, false, true)
const STRAIN_ODD = (false, false, false, false, true, true)
const DISP_ODD = (false, false, true)
# strain component index (1..6) of the tensor entry (i, j), i, j ∈ 1:3 = (r, θ, φ)
const STRAIN_INDEX = [1 4 5; 4 2 6; 5 6 3]
# moment-tensor component index of (i, j)
const MT_INDEX = [1 2 3; 2 4 5; 3 5 6]

"""
    sgtPairs(family) -> (strain = [(isrc, icomp), ...], disp = [...])

Non-zero (source component, field component) pairs on the meridian φ = 0
for `family = :moment` (6 sources) or `:force` (3 sources).
"""
function sgtPairs(family::Symbol)
    odd = family === :moment ? MT_ODD : family === :force ? FORCE_ODD :
          error("family must be :moment or :force")
    strain = [(s, c) for s in eachindex(odd) for c in 1:6 if odd[s] == STRAIN_ODD[c]]
    disp = [(s, c) for s in eachindex(odd) for c in 1:3 if odd[s] == DISP_ODD[c]]
    return (strain=strain, disp=disp)
end

"Unit point source of `family` for component index `s` (lib/sources1D.jl `pointSource`)."
function unitSource(family::Symbol, s::Integer)
    v = zeros(family === :moment ? 6 : 3); v[s] = 1.0
    return family === :moment ? pointSource(M=Tuple(v)) : pointSource(f=Tuple(v))
end

# ---------------------------------------------------------------------------
# Legendre table with first and second θ-derivatives
# ---------------------------------------------------------------------------

"""
    legendreTable2(lmax, mmax, θ; θclamp=1e-6)

Fully normalised X_lm(θ) (same recurrence as `normalisedLegendre` of
lib/sphericalTI_SH_step3b.jl), D = ∂θX and DD = ∂²θX from the Legendre
equation DD = -cotθ D - (Λ - m²/sin²θ) X.  θ is clamped to
[θclamp, π - θclamp] (the formulas contain 1/sinθ; the fields themselves are
regular on the axis, so the clamp only moves the point by θclamp·r).
"""
function legendreTable2(lmax::Integer, mmax::Integer, θ::Real; θclamp=1e-6)
    θc = clamp(θ, θclamp, π - θclamp)
    x, s = cos(θc), sin(θc)
    X = zeros(lmax + 1, mmax + 1); D = zeros(lmax + 1, mmax + 1); DD = zeros(lmax + 1, mmax + 1)
    pmm = 1 / sqrt(4π)
    for m in 0:mmax
        m > 0 && (pmm *= -sqrt((2m + 1) / (2m)) * s)
        m <= lmax || break
        X[m+1, m+1] = pmm
        m + 1 <= lmax && (X[m+2, m+1] = sqrt(2m + 3) * x * pmm)
        for l in m+2:lmax
            a = sqrt((4l^2 - 1) / (l^2 - m^2))
            b = sqrt(((l - 1)^2 - m^2) / (4(l - 1)^2 - 1))
            X[l+1, m+1] = a * (x * X[l, m+1] - b * X[l-1, m+1])
        end
        for l in m:lmax
            prev = l - 1 >= m ? X[l, m+1] : 0.0
            c = l > m ? sqrt((2l + 1) / (2l - 1) * (l - m) * (l + m)) : 0.0
            D[l+1, m+1] = (l * x * X[l+1, m+1] - c * prev) / s
            DD[l+1, m+1] = -(x / s) * D[l+1, m+1] - (l * (l + 1) - m^2 / s^2) * X[l+1, m+1]
        end
    end
    return (X=X, D=D, DD=DD, sinθ=s, cosθ=x, θ=θc)
end

"(X, D, DD) of Y_lm for signed m (Y_l,-m = (-1)^m conj Y_lm; X real so only the sign)."
@inline function _xdd(tab, l, m)
    am = abs(m)
    sg = (m < 0 && isodd(am)) ? -1.0 : 1.0
    return sg * tab.X[l+1, am+1], sg * tab.D[l+1, am+1], sg * tab.DD[l+1, am+1]
end

# ---------------------------------------------------------------------------
# one (l, m) term: displacement and strain on (r, θ, φ = 0)
# ---------------------------------------------------------------------------

"""
    modeFields(rad, r, tab, l, m) -> (u::NTuple{3}, ε::NTuple{6})

Displacement (r, θ, φ) and strain (rr, θθ, φφ, rθ, rφ, θφ) at (r, θ, φ = 0)
of the single term U P_lm + V B_lm + W C_lm, `rad = (U, dU, V, dV, W, dW)`
(dU = ∂r U, km⁻¹ if r in km).  Multiply by e^{imφ} for another meridian.
Strain from the spherical-coordinate gradient
  ε_rr = ∂r u_r,  ε_θθ = (∂θu_θ + u_r)/r,  ε_φφ = (∂φu_φ/s + u_r + u_θ c/s)/r,
  ε_rθ = (∂r u_θ - u_θ/r + ∂θu_r/r)/2,  ε_rφ = (∂r u_φ - u_φ/r + ∂φu_r/(r s))/2,
  ε_θφ = (∂θu_φ - u_φ c/s + ∂φu_θ/s)/(2r).
"""
function modeFields(rad, r, tab, l::Integer, m::Integer)
    U, dU, V, dV, W, dW = rad
    X, D, DD = _xdd(tab, l, m)
    s, c = tab.sinθ, tab.cosθ
    im_m = im * m
    ur = U * X; dθur = U * D; dφur = im_m * U * X; drur = dU * X
    if l == 0
        z = zero(ComplexF64)
        return ((ur, z, z),
                (drur, ur / r, ur / r, z, z, z))
    end
    k = sqrt(l * (l + 1.0))
    Xs = X / s                     # X/sinθ (regular for m ≠ 0)
    dθXs = D / s - c * X / s^2     # ∂θ(X/sinθ)
    uθ = (V * D + W * im_m * Xs) / k
    uφ = (V * im_m * Xs - W * D) / k
    druθ = (dV * D + dW * im_m * Xs) / k
    druφ = (dV * im_m * Xs - dW * D) / k
    dθuθ = (V * DD + W * im_m * dθXs) / k
    dφuθ = (V * im_m * D - W * m^2 * Xs) / k
    dθuφ = (V * im_m * dθXs - W * DD) / k
    dφuφ = (-V * m^2 * Xs - W * im_m * D) / k
    err = drur
    eθθ = (dθuθ + ur) / r
    eφφ = (dφuφ / s + ur + uθ * c / s) / r
    erθ = (druθ - uθ / r + dθur / r) / 2
    erφ = (druφ - uφ / r + dφur / (r * s)) / 2
    eθφ = (dθuφ - uφ * c / s + dφuθ / s) / (2r)
    return ((ur, uθ, uφ), (err, eθθ, eφφ, erθ, erφ, eθφ))
end

# ---------------------------------------------------------------------------
# small tensor utilities
# ---------------------------------------------------------------------------

"3×3 symmetric matrix from 6 components in (rr, θθ, φφ, rθ, rφ, θφ) order."
strainMatrix(e) = [e[1] e[4] e[5]; e[4] e[2] e[6]; e[5] e[6] e[3]]
"6 components (rr, θθ, φφ, rθ, rφ, θφ) of a symmetric 3×3 matrix."
strainVector(E) = (E[1, 1], E[2, 2], E[3, 3], E[1, 2], E[1, 3], E[2, 3])
"3×3 symmetric matrix of a moment tensor (Mrr, Mrθ, Mrφ, Mθθ, Mθφ, Mφφ)."
mtMatrix(M) = [M[1] M[2] M[3]; M[2] M[4] M[5]; M[3] M[5] M[6]]
"(Mrr, Mrθ, Mrφ, Mθθ, Mθφ, Mφφ) of a symmetric 3×3 matrix."
mtVector(A) = (A[1, 1], A[1, 2], A[1, 3], A[2, 2], A[2, 3], A[3, 3])
