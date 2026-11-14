# helpers/sphericalTI_PSV.jl
#
# Step P1 of the 1-D spherical TI workflow, P-SV (spheroidal) part
# (notebooks/SphericalTI_PSV_oct2026.ipynb).  Same chain as the SH Step 1:
#
#   ρω²u_i + ∂_j(C_ijkl ∂_l u_k) = -f_i   (Cartesian, TI with axis x/|x|)
#   u = U(r) P_lm + V(r) B_lm,  P_lm = r̂ Y_lm,  B_lm = k⁻¹ ∇₁Y_lm,  k = √(l(l+1))
#   (Y_lm fully normalised; P_lm, B_lm orthonormal)
#   project onto P*_lm and B*_lm  →  two coupled radial equations R_U, R_V
#
# Zonal jet (m = 0) as for SH: q = cot θ, p0 = P, p1 = P', P'' = -qP' - ΛP,
# with Λ = k² so that only k appears.  u = (U p0, V p1/k, 0):
#   E_r = R_U p0 ,  E_θ = R_V p1/k      (checked: no other p0/p1/q dependence)
#
# Fluid regions (Geller & Ohminato 1994, GJI 116, §3): unknown Q = P/ω =
# λ∇·u/ω, strong form ω²Q/λ + ∂_i(Q_,i/ρ) = 0 (source-free); Q = q(r) Y_lm:
#   r²·[ω²q/λ + r⁻²∂r(r² q'/ρ) - Λq/(ρr²)] = ∂r(r²/ρ ∂r q) + (ω² r²/λ - Λ/ρ) q
# Solid–fluid coupling only through the surface terms of G&O eq. (35):
#   solid rows  + ∫ ω w_i Q n_i^(S) dS ,   fluid rows  - ∫ ω w u_i n_i^(F) dS
# (Step P4.)
#
# Requires helpers/sphericalTI_SH.jl (tiElasticTensor, cartesianElasticOperator,
# sphericalTIStressDivergence, zonalJetDθ, _cross, _num).

using Symbolics, LinearAlgebra, Random

# ---------------------------------------------------------------------------
# spheroidal projection (solid, TI)
# ---------------------------------------------------------------------------

"""
    spheroidalProjection(U, V, ρ, A, C, F, L, N; ω, k, q, p0, p1)

Insert u = (U p0, V p1/k, 0) (zonal spheroidal field, normalised P and B)
in ρω²u + ∇·σ (TI, axis r̂) and return the components E and the radial
operators R_U (coefficient of p0 in E_r) and R_V (k × coefficient of p1 in E_θ).
"""
function spheroidalProjection(U, V, ρ, A, C, F, L, N; ω, k, q, p0, p1)
    Λ = k^2
    Dr = e -> expand_derivatives(Differential(x)(Num(e)))
    Dθ = zonalJetDθ(q, p0, p1, Λ)
    Dφ = e -> Num(0)
    u = (U * p0, V * p1 / k, Num(0))
    d = sphericalTIStressDivergence(u; Dr=Dr, Dθ=Dθ, Dφ=Dφ, r=x, q=q,
                                    s⁻¹=Num(1), A=A, C=C, F=F, L=L, N=N)
    E = ntuple(i -> Symbolics.expand(ω^2 * ρ * u[i] + d[i]), 3)
    RU = substitute(E[1], Dict(p0 => 1, p1 => 0, q => 0))
    RV = k * substitute(E[2], Dict(p1 => 1, p0 => 0, q => 0))
    return (E=E, RU=Symbolics.expand(RU), RV=Symbolics.expand(RV))
end

"""
    spheroidalExpectedR2(U, V, ρ, A, C, F, L, N; ω, k)

r² R_U and r² R_V written with the traction fluxes taken out
(r²σ_rr = r²C U' + 2rF U - krF V,  r²σ_rθ = r²L V' - rL V + krL U):

  r²R_U = ∂r(r²σ_rr) - 2rF U' - krL V' + (ω²ρr² - 4A + 4N - k²L) U + k(2A + L - 2N) V
  r²R_V = ∂r(r²σ_rθ) + krF U' + rL V' + (ω²ρr² - k²A - L + 2N) V + k(2A + L - 2N) U

No second derivatives and no material derivatives remain outside the fluxes;
the resulting natural weak form is symmetric (the DSM P-SV ω²T - H).
"""
function spheroidalExpectedR2(U, V, ρ, A, C, F, L, N; ω, k)
    D = Differential(x)
    r = x
    σrr = C * D(U) + F * (2U - k * V) / r
    σrθ = L * (D(V) - V / r + k * U / r)
    eU = D(r^2 * σrr) - 2r * F * D(U) - k * r * L * D(V) +
         (ω^2 * r^2 * ρ - 4A + 4N - k^2 * L) * U + k * (2A + L - 2N) * V
    eV = D(r^2 * σrθ) + k * r * F * D(U) + r * L * D(V) +
         (ω^2 * r^2 * ρ - k^2 * A - L + 2N) * V + k * (2A + L - 2N) * U
    return expand_derivatives(eU), expand_derivatives(eV)
end

"Symbolic derivation for display."
function deriveSpheroidalSymbolic()
    @variables ω k q p0 p1
    @variables U(x) V(x) ρ(x) A(x) C(x) F(x) L(x) N(x)
    proj = spheroidalProjection(U, V, ρ, A, C, F, L, N; ω=ω, k=k, q=q, p0=p0, p1=p1)
    eU, eV = spheroidalExpectedR2(U, V, ρ, A, C, F, L, N; ω=ω, k=k)
    return (E=proj.E, RU=proj.RU, RV=proj.RV, r2RU_expected=eU, r2RV_expected=eV,
            symbols=(; ω, k, q, p0, p1, U, V, ρ, A, C, F, L, N))
end

_psvTestProfiles() = Dict(
    :U => s -> s^2 - s / 3 + 1,
    :V => s -> 1 + s / 2 - s^2 / 5,
    :ρ => s -> 2 + s,
    :A => s -> 5 + s^2,
    :C => s -> 4 + s,
    :F => s -> 1 + s / 2,
    :L => s -> s^2 / 3 + 1,
    :N => s -> 2 - s / 5 + s^2 / 4,
)

"""
    verifySpheroidalProjection(; l=3, ω=1.4, ntrials=5, rng=MersenneTwister(1))

Concrete-profile check at random (r, cot θ, P, P'):
E_φ = 0; E_r has no P' part and E_θ no P part (and no θ dependence);
E_r = (r²R_U)_expected/r² · P and E_θ = (r²R_V)_expected/r² · P'/k.
Returns the maximum discrepancies relative to the size of E.
"""
function verifySpheroidalProjection(; l=3, ω=1.4, ntrials=5, rng=MersenneTwister(1))
    @variables q p0 p1
    k = sqrt(l * (l + 1))
    P = Dict(kk => f(x) for (kk, f) in _psvTestProfiles())
    proj = spheroidalProjection(P[:U], P[:V], P[:ρ], P[:A], P[:C], P[:F], P[:L], P[:N];
                                ω=ω, k=k, q=q, p0=p0, p1=p1)
    eU, eV = spheroidalExpectedR2(P[:U], P[:V], P[:ρ], P[:A], P[:C], P[:F], P[:L], P[:N];
                                  ω=ω, k=k)
    err = (Eφ=0.0, Er=0.0, Eθ=0.0)
    for _ in 1:ntrials
        rv, qv, p0v, p1v = 0.3 + 2rand(rng), randn(rng), randn(rng), 0.5 + rand(rng)
        s = Dict(x => rv, q => qv, p0 => p0v, p1 => p1v)
        Eφ = _num(substitute(proj.E[3], s))
        Er = _num(substitute(proj.E[1], s))
        Eθ = _num(substitute(proj.E[2], s))
        Ur = _num(substitute(eU, Dict(x => rv))) / rv^2
        Vr = _num(substitute(eV, Dict(x => rv))) / rv^2
        scale = max(abs(Er), abs(Eθ), 1.0)
        err = (Eφ=max(err.Eφ, abs(Eφ) / scale),
               Er=max(err.Er, abs(Er - Ur * p0v) / scale),
               Eθ=max(err.Eθ, abs(Eθ - Vr * p1v / k) / scale))
    end
    return err
end

"""
    cartesianSpheroidalCheck(; harmonic=:Y32, ω=1.4, npoints=3, rng=MersenneTwister(1))

Independent m ≠ 0 check: full Cartesian TI tensor (axis x/|x|) applied to
u = U(r) Y r̂ + V(r) ∇₁Y/k for a real Y (l=3, m=2 or l=2, m=1), compared with
R_U Y r̂ + R_V ∇₁Y/k (expected radial operators).  Returns max abs/rel error.
"""
function cartesianSpheroidalCheck(; harmonic=:Y32, ω=1.4, npoints=3, rng=MersenneTwister(1))
    X = (x, y, z)
    rr = sqrt(x^2 + y^2 + z^2)
    prof = _psvTestProfiles()
    P = Dict(kk => f(rr) for (kk, f) in prof)
    n = [x / rr, y / rr, z / rr]
    Cijkl = tiElasticTensor(P[:A], P[:C], P[:F], P[:L], P[:N], n)
    Y, l = harmonic === :Y32 ? ((x^2 - y^2) * z / rr^3, 3) :
           harmonic === :Y21 ? (x * z / rr^2, 2) :
           error("harmonic must be :Y32 or :Y21")
    k = sqrt(l * (l + 1))
    ∇₁Y = [expand_derivatives(Differential(X[i])(Y)) * rr for i in 1:3]
    u = [P[:U] * Y * n[i] + P[:V] * ∇₁Y[i] / k for i in 1:3]
    Lu = expand_derivatives.(cartesianElasticOperator(u, Cijkl, P[:ρ], ω, X))
    eU, eV = spheroidalExpectedR2(prof[:U](x), prof[:V](x), prof[:ρ](x), prof[:A](x),
                                  prof[:C](x), prof[:F](x), prof[:L](x), prof[:N](x);
                                  ω=ω, k=k)
    maxabs, maxref = 0.0, 0.0
    for _ in 1:npoints
        pt = Dict(x => 2rand(rng) - 1, y => 2rand(rng) - 1, z => 0.2 + rand(rng))
        rv = sqrt(pt[x]^2 + pt[y]^2 + pt[z]^2)
        RUv = _num(substitute(eU, Dict(x => rv))) / rv^2
        RVv = _num(substitute(eV, Dict(x => rv))) / rv^2
        Yv = _num(substitute(Y, pt))
        for i in 1:3
            lhs = _num(substitute(Lu[i], pt))
            rhs = RUv * Yv * _num(substitute(n[i], pt)) + RVv * _num(substitute(∇₁Y[i], pt)) / k
            maxabs = max(maxabs, abs(lhs - rhs)); maxref = max(maxref, abs(lhs))
        end
    end
    return (harmonic=harmonic, l=l, maxabs=maxabs, maxrel=maxabs / maxref)
end

# ---------------------------------------------------------------------------
# fluid: Q = P/ω (Geller & Ohminato 1994)
# ---------------------------------------------------------------------------

"""
    fluidProjection(qf, ρ, λ; ω, Λ, qc, p0, p1)

Projection of the G&O fluid equation ω²Q/λ + ∇·(∇Q/ρ) for Q = q(r) P(θ)
(zonal jet), returns the coefficient of p0 (the radial operator) and the
expected r²-form  ∂r(r²/ρ ∂r q) + (ω²r²/λ - Λ/ρ) q.
"""
function fluidProjection(qf, ρ, λ; ω, Λ, qc, p0, p1)
    Dr = e -> expand_derivatives(Differential(x)(Num(e)))
    Dθ = zonalJetDθ(qc, p0, p1, Λ)
    Q = qf * p0
    vr = Dr(Q) / ρ                 # ∇Q/ρ
    vθ = Dθ(Q) / (x * ρ)
    div = Dr(vr) + 2vr / x + Dθ(vθ) / x + qc * vθ / x
    E = Symbolics.expand(ω^2 * Q / λ + div)
    R = substitute(E, Dict(p0 => 1, p1 => 0, qc => 0))
    D = Differential(x)
    expected = expand_derivatives(D(x^2 / ρ * D(qf)) + (ω^2 * x^2 / λ - Λ / ρ) * qf)
    return (E=E, R=R, r2R_expected=expected)
end

"Concrete-profile check of the fluid projection (no P' part, no θ dependence, r²R = expected)."
function verifyFluidProjection(; l=3, ω=1.4, ntrials=5, rng=MersenneTwister(2))
    @variables qc p0 p1
    Λ = l * (l + 1)
    prof = (q=s -> s^2 - s / 3 + 1, ρ=s -> 2 + s, λ=s -> 5 + s^2)
    fp = fluidProjection(prof.q(x), prof.ρ(x), prof.λ(x); ω=ω, Λ=Λ, qc=qc, p0=p0, p1=p1)
    err = 0.0
    for _ in 1:ntrials
        rv, qv, p0v, p1v = 0.3 + 2rand(rng), randn(rng), randn(rng), randn(rng)
        Ev = _num(substitute(fp.E, Dict(x => rv, qc => qv, p0 => p0v, p1 => p1v)))
        ex = _num(substitute(fp.r2R_expected, Dict(x => rv))) / rv^2
        err = max(err, abs(Ev - ex * p0v) / max(1.0, abs(Ev)))
    end
    return err
end

# ---------------------------------------------------------------------------
# famousEquation entries (x ≡ r)
# ---------------------------------------------------------------------------

"""
famousEquation `"1DsphericalPSV_TI_Freq"`: fields (U, V), the r²-weighted
radial P-SV equations with geometric material fields (per l):

  ρr² = ρr²,  Cr² = r²C,  Fr = rF,  Lr² = r²L,  Lr = rL,
  kFr = k r F,  kLr = k r L,
  mUU = -4A + 4N - ΛL,  mVV = -ΛA - L + 2N,  mUV = k(2A + L - 2N)

  U: ∂x(Cr²∂xU) + ∂x(2Fr U) - ∂x(kFr V) - 2Fr ∂xU - kLr ∂xV + (ω²ρr² + mUU) U + mUV V
  V: ∂x(Lr²∂xV) - ∂x(Lr V) + ∂x(kLr U) + kFr ∂xU + Lr ∂xV + (ω²ρr² + mVV) V + mUV U

(= -r² f_U, -r² f_V).  Natural fluxes: r²σ_rr and r²σ_rθ (traction).
"""
function flexOPT.famousEquation(::Val{:eq_1DsphericalPSV_TI_Freq})
    @variables ω
    @variables U(x) V(x) ρr²(x) Cr²(x) Fr(x) Lr²(x) Lr(x) kFr(x) kLr(x)
    @variables mUU(x) mVV(x) mUV(x) fUr²(x) fVr²(x)
    eqU = ∂x(Cr² * ∂x(U)) + ∂x(2 * Fr * U) - ∂x(kFr * V) - 2 * Fr * ∂x(U) -
          kLr * ∂x(V) + ω^2 * ρr² * U + mUU * U + mUV * V
    eqV = ∂x(Lr² * ∂x(V)) - ∂x(Lr * V) + ∂x(kLr * U) + kFr * ∂x(U) + Lr * ∂x(V) +
          ω^2 * ρr² * V + mVV * V + mUV * U
    exprs = (eqU, eqV)
    fields = (U, V)
    vars = (ρr², Cr², Fr, Lr², Lr, kFr, kLr, mUU, mVV, mUV, ω)
    extexprs = (fUr², fVr²)
    extfields = extexprs
    extvars = extexprs
    coordinates = (x,)
    ∂, ∂² = usefulPartials(coordinates)
    return exprs, fields, vars, extexprs, extfields, extvars, coordinates, ∂, ∂²
end

"""
famousEquation `"1DsphericalFluidQ_Freq"`: field q (Q = P/ω, G&O 1994),
  ∂x(Pr² ∂x q) + ω² κr² q - Mₗ q   (= -source)
with Pr² = r²/ρ, κr² = r²/λ, Mₗ = Λ/ρ.  Natural flux r²q'/ρ (∝ normal
displacement); Q itself is continuous (essential) -- the opposite of a solid.
"""
function flexOPT.famousEquation(::Val{:eq_1DsphericalFluidQ_Freq})
    @variables ω
    @variables qf(x) Pr²(x) κr²(x) Mₗ(x) gr²(x)
    expr = ∂x(Pr² * ∂x(qf)) + ω^2 * κr² * qf - Mₗ * qf
    exprs = expr
    fields = qf
    vars = (Pr², κr², Mₗ, ω)
    extexprs = gr²
    extfields = gr²
    extvars = gr²
    coordinates = (x,)
    ∂, ∂² = usefulPartials(coordinates)
    return exprs, fields, vars, extexprs, extfields, extvars, coordinates, ∂, ∂²
end

"""
    verifyPSVEntry(; l=3, ω=1.4, rng=MersenneTwister(4))

The famousEquation entry with its geometric fields replaced by concrete
profiles must equal r²R_U, r²R_V (expected forms).  Returns max rel. error.
"""
function verifyPSVEntry(; l=3, ω=1.4, npoints=5, rng=MersenneTwister(4))
    eq, _ = famousEquations("1DsphericalPSV_TI_Freq")
    prof = _psvTestProfiles()
    k = sqrt(l * (l + 1)); Λ = l * (l + 1)
    p = Dict(kk => f(x) for (kk, f) in prof)
    @variables ω_ U(x) V(x) ρr²(x) Cr²(x) Fr(x) Lr²(x) Lr(x) kFr(x) kLr(x) mUU(x) mVV(x) mUV(x)
    @variables ω
    sub = Dict(U => p[:U], V => p[:V], ρr² => p[:ρ] * x^2, Cr² => x^2 * p[:C], Fr => x * p[:F],
               Lr² => x^2 * p[:L], Lr => x * p[:L], kFr => k * x * p[:F], kLr => k * x * p[:L],
               mUU => -4p[:A] + 4p[:N] - Λ * p[:L], mVV => -Λ * p[:A] - p[:L] + 2p[:N],
               mUV => k * (2p[:A] + p[:L] - 2p[:N]), ω => ω)
    eU, eV = spheroidalExpectedR2(p[:U], p[:V], p[:ρ], p[:A], p[:C], p[:F], p[:L], p[:N]; ω=ω, k=k)
    e1 = expand_derivatives(substitute(eq.exprs[1], sub))
    e2 = expand_derivatives(substitute(eq.exprs[2], sub))
    err = 0.0
    for _ in 1:npoints
        rv = 0.3 + 2rand(rng)
        a1, b1 = _num(substitute(e1, Dict(x => rv))), _num(substitute(eU, Dict(x => rv)))
        a2, b2 = _num(substitute(e2, Dict(x => rv))), _num(substitute(eV, Dict(x => rv)))
        err = max(err, abs(a1 - b1) / max(1.0, abs(b1)), abs(a2 - b2) / max(1.0, abs(b2)))
    end
    return err
end
