# helpers/sphericalTI_SH.jl
#
# Step 1 of the 1-D spherical, frequency-domain, transversely isotropic (TI)
# workflow (SH / toroidal first), following
#   Cummins, Geller, Hatori & Takeuchi (1994a, GRL, 94GL00013)   -- SH, DSM
#   Cummins, Geller & Takeuchi        (1994b, GRL, 94GL01281)   -- P-SV, DSM
#   Takeuchi, Geller & Cummins        (1996,  GRL, 96GL00973)   -- modified ops
#
# Chain implemented here (all symbolic, Symbolics.jl only):
#
#   (1a) ρ ω² u_i + ∂_j( C_ijkl ∂_l u_k ) = -f_i      Cartesian, TI with the
#        symmetry axis n = x/|x| (Love's A, C, F, L, N)
#   (1b) the same operator in spherical physical components (r, θ, φ)
#   (1c) u = W(r) C_lm(θ,φ),  C_lm = -k⁻¹ r̂ × ∇₁Y_lm,  k = √(l(l+1)),
#        Y_lm fully normalised  ->  project onto C*_lm (∫ dΩ)  ->  radial
#        strong form, depending on l only through Λ = l(l+1)
#   (1d) registration of that radial equation as a famousEquation entry
#
# The angular algebra is done with a "jet" trick so that only polynomial
# cancellations are needed (no trig simplification):
#   q = cot θ,  dq/dθ = -(1+q²)
#   p0 = P(θ), p1 = P'(θ),  P'' = -q P' - Λ P   (zonal Legendre ODE, m = 0)
# The operator is rotationally invariant, so the radial equation cannot depend
# on m; (1c) is verified for m ≠ 0 independently by `cartesianToroidalCheck`,
# which applies the full Cartesian TI operator to a real l=3, m=2 (or l=2,
# m=1) toroidal field.
#
# Normalisation (1/k, and the N_lm of Y_lm) is a constant factor common to u
# and to the test harmonic, so it cancels in the operator.  It matters only for
# the source projection  f_W(r) = ∫ f · C*_lm dΩ  (Step 3).
#
# Requirements: `using .commonBatchs, .flexOPT` already done in Main (the
# coordinate `x`, which plays the role of r here, and `usefulPartials` come
# from commonBatchs; `famousEquation` is extended in flexOPT).
#
# Nothing in src/ is modified.

using Symbolics, LinearAlgebra, Random

# ---------------------------------------------------------------------------
# (1a) Cartesian TI operator
# ---------------------------------------------------------------------------

const _δ3 = [i == j ? 1 : 0 for i in 1:3, j in 1:3]

"true iff the (expanded) symbolic expression is the number zero"
_symzero(e) = (v = Symbolics.value(Symbolics.expand(e)); v isa Number && iszero(v))

"""
    tiElasticTensor(A, C, F, L, N, n)

Fourth-order TI elastic tensor with unit symmetry axis `n` (Love notation):

    C_ijkl = (A-2N) δij δkl + N (δik δjl + δil δjk)
           + (F-A+2N)(δij nk nl + ni nj δkl)
           + (L-N)(δik nj nl + δil nj nk + δjk ni nl + δjl ni nk)
           + (A+C-2F-4L) ni nj nk nl

With n = (0,0,1): C₃₃₃₃=C, C₁₁₁₁=A, C₁₁₃₃=F, C₁₃₁₃=L, C₁₂₁₂=N, C₁₁₂₂=A-2N.
"""
function tiElasticTensor(A, C, F, L, N, n)
    δ = _δ3
    T = Array{Num}(undef, 3, 3, 3, 3)
    for i in 1:3, j in 1:3, k in 1:3, l in 1:3
        T[i, j, k, l] =
            (A - 2N) * δ[i, j] * δ[k, l] +
            N * (δ[i, k] * δ[j, l] + δ[i, l] * δ[j, k]) +
            (F - A + 2N) * (δ[i, j] * n[k] * n[l] + n[i] * n[j] * δ[k, l]) +
            (L - N) * (δ[i, k] * n[j] * n[l] + δ[i, l] * n[j] * n[k] +
                       δ[j, k] * n[i] * n[l] + δ[j, l] * n[i] * n[k]) +
            (A + C - 2F - 4L) * n[i] * n[j] * n[k] * n[l]
    end
    return T
end

"""
    cartesianElasticOperator(u, Cijkl, ρ, ω, X)

Returns the three components of ρ ω² u_i + ∂_j(C_ijkl ∂_l u_k) (unexpanded
Differentials; call `expand_derivatives` to evaluate them).
"""
function cartesianElasticOperator(u, Cijkl, ρ, ω, X)
    D = [Differential(X[j]) for j in 1:3]
    return map(1:3) do i
        traction = sum(D[j](sum(Cijkl[i, j, k, l] * D[l](u[k])
                                for k in 1:3, l in 1:3)) for j in 1:3)
        ω^2 * ρ * u[i] + traction
    end
end

"""
    checkTIisotropicLimitAndLoveNotation()

(i) A=C=λ+2μ, F=λ, L=N=μ must give the isotropic tensor exactly, for a
symbolic axis n; (ii) with n = ẑ the tensor entries must be Love's A,C,F,L,N.
Returns a NamedTuple of Bools.
"""
function checkTIisotropicLimitAndLoveNotation()
    @variables λ μ n1 n2 n3 A C F L N
    n = [n1, n2, n3]
    δ = _δ3
    CTI = tiElasticTensor(λ + 2μ, λ + 2μ, λ, μ, μ, n)
    isotropic_ok = all(CartesianIndices(CTI)) do I
        i, j, k, l = Tuple(I)
        Ciso = λ * δ[i, j] * δ[k, l] + μ * (δ[i, k] * δ[j, l] + δ[i, l] * δ[j, k])
        _symzero(CTI[I] - Ciso)
    end
    Cz = tiElasticTensor(A, C, F, L, N, [0, 0, 1])
    love_ok = all(_symzero(a - b) for (a, b) in (
        (Cz[3, 3, 3, 3], C), (Cz[1, 1, 1, 1], A), (Cz[1, 1, 3, 3], F),
        (Cz[1, 3, 1, 3], L), (Cz[1, 2, 1, 2], N), (Cz[1, 1, 2, 2], A - 2N)))
    return (isotropic_limit=isotropic_ok, love_notation=love_ok)
end

"""
famousEquation entry: 3-D frequency-domain elastic equation, TI with a radial
symmetry axis, written exactly as ρ ω² u_i + ∂_j(C_ijkl ∂_l u_k) (= -f_i).
Kept as a reference/display object; it is not meant to be fed to the OPT
recipe builder in 3-D for now.
"""
function flexOPT.famousEquation(::Val{:eq_3DsismoFreqTIradial})
    @variables ω
    @variables ρ(x, y, z) A(x, y, z) C(x, y, z) F(x, y, z) L(x, y, z) N(x, y, z)
    @variables u(x, y, z)[1:3] f(x, y, z)[1:3]
    X = (x, y, z)
    rr = sqrt(x^2 + y^2 + z^2)
    n = [x / rr, y / rr, z / rr]
    Cijkl = tiElasticTensor(A, C, F, L, N, n)
    uu = [u[1], u[2], u[3]]
    ops = cartesianElasticOperator(uu, Cijkl, ρ, ω, X)
    exprs = (ops[1], ops[2], ops[3])
    fields = (u[1], u[2], u[3])
    vars = (ρ, A, C, F, L, N, ω)
    extexprs = (f[1], f[2], f[3])
    extfields = extexprs
    extvars = extexprs
    coordinates = (x, y, z)
    ∂, ∂² = usefulPartials(coordinates)
    return exprs, fields, vars, extexprs, extfields, extvars, coordinates, ∂, ∂²
end

# ---------------------------------------------------------------------------
# (1b) spherical physical components
# ---------------------------------------------------------------------------

"""
    sphericalTIStressDivergence(u; Dr, Dθ, Dφ, r, q, s⁻¹, A, C, F, L, N)

∇·σ in spherical physical components for u = (u_r, u_θ, u_φ), TI with the
axis r̂.  `Dr, Dθ, Dφ` are derivative operators (closures), `q = cot θ`,
`s⁻¹ = 1/sin θ`.  Standard formulas:

    e_rr = ∂r u_r
    e_θθ = (∂θ u_θ + u_r)/r
    e_φφ = (s⁻¹ ∂φ u_φ + u_r + q u_θ)/r
    2e_rθ = ∂θ u_r / r + ∂r u_θ - u_θ/r
    2e_rφ = s⁻¹ ∂φ u_r / r + ∂r u_φ - u_φ/r
    2e_θφ = (∂θ u_φ - q u_φ + s⁻¹ ∂φ u_θ)/r

    σ_rr = C e_rr + F(e_θθ+e_φφ),  σ_θθ = F e_rr + A e_θθ + (A-2N) e_φφ,
    σ_φφ = F e_rr + (A-2N) e_θθ + A e_φφ,
    σ_rθ = 2L e_rθ,  σ_rφ = 2L e_rφ,  σ_θφ = 2N e_θφ
"""
function sphericalTIStressDivergence(u; Dr, Dθ, Dφ, r, q, s⁻¹, A, C, F, L, N)
    ur, uθ, uφ = u
    e_rr = Dr(ur)
    e_θθ = Dθ(uθ) / r + ur / r
    e_φφ = s⁻¹ * Dφ(uφ) / r + ur / r + q * uθ / r
    e_rθ = (Dθ(ur) / r + Dr(uθ) - uθ / r) / 2
    e_rφ = (s⁻¹ * Dφ(ur) / r + Dr(uφ) - uφ / r) / 2
    e_θφ = (Dθ(uφ) / r - q * uφ / r + s⁻¹ * Dφ(uθ) / r) / 2

    σ_rr = C * e_rr + F * (e_θθ + e_φφ)
    σ_θθ = F * e_rr + A * e_θθ + (A - 2N) * e_φφ
    σ_φφ = F * e_rr + (A - 2N) * e_θθ + A * e_φφ
    σ_rθ = 2L * e_rθ
    σ_rφ = 2L * e_rφ
    σ_θφ = 2N * e_θφ

    d_r = Dr(σ_rr) + Dθ(σ_rθ) / r + s⁻¹ * Dφ(σ_rφ) / r +
          (2σ_rr - σ_θθ - σ_φφ + q * σ_rθ) / r
    d_θ = Dr(σ_rθ) + Dθ(σ_θθ) / r + s⁻¹ * Dφ(σ_θφ) / r +
          ((σ_θθ - σ_φφ) * q + 3σ_rθ) / r
    d_φ = Dr(σ_rφ) + Dθ(σ_θφ) / r + s⁻¹ * Dφ(σ_φφ) / r +
          (3σ_rφ + 2q * σ_θφ) / r
    return (d_r, d_θ, d_φ)
end

"θ-derivative on the zonal jet (q = cot θ, p0 = P, p1 = P', P'' = -qP' - ΛP)."
function zonalJetDθ(q, p0, p1, Λ)
    return function (e)
        e = Num(e)
        Symbolics.derivative(e, q) * (-(1 + q^2)) +
            Symbolics.derivative(e, p0) * p1 +
            Symbolics.derivative(e, p1) * (-q * p1 - Λ * p0)
    end
end

# ---------------------------------------------------------------------------
# (1c) toroidal projection
# ---------------------------------------------------------------------------

"""
    toroidalProjection(W, ρ, A, C, F, L, N; ω, Λ, q, p0, p1)

W, ρ, A, … are expressions in `x` (≡ r): either symbolic functions W(x) … or
concrete test profiles.  Inserts the zonal toroidal field

    u = W(r) (r̂ × ∇₁Y)·(-1)  ->  (u_r, u_θ, u_φ) = (0, 0, -W P'(θ))

(the 1/k normalisation cancels) and returns

  * `E`  : the three components of ρω²u + ∇·σ (expanded)
  * `R`  : the radial operator, i.e. E = R(r) · (0, 0, -P')  ⇔  R = -E_φ/p1,
           which is ∫ C*_lm · (ρω²u + ∇·σ) dΩ for normalised C_lm.

`R` is read off at q = p0 = 0, p1 = 1; that this is legitimate (E_φ ∝ p1 with
no q or p0 dependence, E_r = E_θ = 0) is checked by `verifyToroidalProjection`.
"""
function toroidalProjection(W, ρ, A, C, F, L, N; ω, Λ, q, p0, p1)
    Dr = e -> expand_derivatives(Differential(x)(Num(e)))
    Dθ = zonalJetDθ(q, p0, p1, Λ)
    Dφ = e -> Num(0)                       # m = 0
    z0 = Num(0)
    u = (z0, z0, -W * p1)
    d = sphericalTIStressDivergence(u; Dr=Dr, Dθ=Dθ, Dφ=Dφ, r=x, q=q,
                                    s⁻¹=Num(1), A=A, C=C, F=F, L=L, N=N)
    E = ntuple(i -> Symbolics.expand(ω^2 * ρ * u[i] + d[i]), 3)
    R = -substitute(E[3], Dict(p1 => 1, p0 => 0, q => 0))
    return (E=E, R=Symbolics.expand(R))
end

"Expected radial strong form (Takeuchi & Saito toroidal system, TI):
 ω²ρW + r⁻³ d/dr[ r³ L (W' - W/r) ] - (Λ-2) N W / r²"
function toroidalRadialStrongForm(W, ρ, L, N; ω, Λ)
    D = Differential(x)
    return expand_derivatives(ω^2 * ρ * W +
        D(x^3 * L * (D(W) - W / x)) / x^3 - (Λ - 2) * N * W / x^2)
end

"""
Self-adjoint form used for the famousEquation entry (multiplied by r²):

    r² R = ∂r(r²L ∂r W) - ∂r(rL W) + rL ∂r W + (ω² ρ r² - [L + (Λ-2)N]) W

whose natural weak form is exactly the DSM SH operator  ω² T - H  with
    T_ij = ∫ ρ r² X_i X_j dr
    H_ij = ∫ L (r X_i' - X_i)(r X_j' - X_j) dr + (Λ-2) ∫ N X_i X_j dr
and the boundary term  [X r² L (W' - W/r)] = [X r² σ_rφ-type traction], which
vanishes at a free surface.
"""
function toroidalSelfAdjointForm(W, ρ, L, N; ω, Λ)
    D = Differential(x)
    return expand_derivatives(D(x^2 * L * D(W)) - D(x * L * W) + x * L * D(W) +
        (ω^2 * ρ * x^2 - (L + (Λ - 2) * N)) * W)
end

"Symbolic derivation for display: returns symbols, E, R and the expected R."
function deriveToroidalSymbolic()
    @variables ω Λ q p0 p1
    @variables W(x) ρ(x) A(x) C(x) F(x) L(x) N(x)
    proj = toroidalProjection(W, ρ, A, C, F, L, N; ω=ω, Λ=Λ, q=q, p0=p0, p1=p1)
    expected = toroidalRadialStrongForm(W, ρ, L, N; ω=ω, Λ=Λ)
    selfadj = toroidalSelfAdjointForm(W, ρ, L, N; ω=ω, Λ=Λ)
    return (E=proj.E, R=proj.R, expected=expected, selfadjoint=selfadj,
            symbols=(; ω, Λ, q, p0, p1, W, ρ, A, C, F, L, N))
end

# concrete radial test profiles (functions of a scalar s = r)
_testProfiles() = Dict(
    :W => s -> s^2 - s / 3 + 1,
    :ρ => s -> 2 + s,
    :A => s -> 5 + s^2,
    :C => s -> 4 + s,
    :F => s -> 1 + s / 2,
    :L => s -> s^2 / 3 + 1,
    :N => s -> 2 - s / 5 + s^2 / 4,
)

_num(e) = Float64(Symbolics.value(e))

"""
    verifyToroidalProjection(; l=3, ω=1.4, ntrials=5, rng=MersenneTwister(1))

Concrete-profile check of (1c) at random (r, q, p0, p1):
  * E_r = E_θ = 0
  * E_φ(p1 = 0) = 0   (no P term)
  * -E_φ/p1 does not depend on (q, p0) and equals the expected radial form
  * the self-adjoint (flexOPT) form equals r² × expected
Returns the maximum absolute discrepancies.
"""
function verifyToroidalProjection(; l=3, ω=1.4, ntrials=5, rng=MersenneTwister(1))
    @variables q p0 p1
    Λ = l * (l + 1)
    P = Dict(k => f(x) for (k, f) in _testProfiles())
    proj = toroidalProjection(P[:W], P[:ρ], P[:A], P[:C], P[:F], P[:L], P[:N];
                              ω=ω, Λ=Λ, q=q, p0=p0, p1=p1)
    expected = toroidalRadialStrongForm(P[:W], P[:ρ], P[:L], P[:N]; ω=ω, Λ=Λ)
    selfadj = toroidalSelfAdjointForm(P[:W], P[:ρ], P[:L], P[:N]; ω=ω, Λ=Λ)
    err = (Er=0.0, Eθ=0.0, Eφ_p0=0.0, R=0.0, selfadjoint=0.0)
    for _ in 1:ntrials
        rv, qv, p0v, p1v = 0.3 + 2rand(rng), randn(rng), randn(rng), 0.5 + rand(rng)
        s = Dict(x => rv, q => qv, p0 => p0v, p1 => p1v)
        Er = _num(substitute(proj.E[1], s))
        Eθ = _num(substitute(proj.E[2], s))
        Eφ = _num(substitute(proj.E[3], s))
        Eφ0 = _num(substitute(proj.E[3], Dict(x => rv, q => qv, p0 => p0v, p1 => 0)))
        ex = _num(substitute(expected, Dict(x => rv)))
        sa = _num(substitute(selfadj, Dict(x => rv)))
        err = (Er=max(err.Er, abs(Er)), Eθ=max(err.Eθ, abs(Eθ)),
               Eφ_p0=max(err.Eφ_p0, abs(Eφ0)),
               R=max(err.R, abs(-Eφ / p1v - ex)),
               selfadjoint=max(err.selfadjoint, abs(sa - rv^2 * ex)))
    end
    return err
end

_cross(a, b) = [a[2] * b[3] - a[3] * b[2], a[3] * b[1] - a[1] * b[3],
                a[1] * b[2] - a[2] * b[1]]

"""
    cartesianToroidalCheck(; harmonic=:Y32, ω=1.4, npoints=3, rng=MersenneTwister(1))

Independent check of (1a)→(1c) with m ≠ 0: builds the full Cartesian TI
tensor with axis x/|x|, the real toroidal field u = W(r) (-r̂ × ∇₁Y) for
Y ∝ (x²-y²)z/r³ (l=3, m=2) or Y ∝ xz/r² (l=2, m=1), applies
ρω²u_i + ∂_j(C_ijkl ∂_l u_k) and compares, at random points, with
R(r) · (-r̂ × ∇₁Y) where R is the radial strong form.  Returns the maximum
absolute and relative discrepancy.  (Takes a little while: the Cartesian
expression is large.)
"""
function cartesianToroidalCheck(; harmonic=:Y32, ω=1.4, npoints=3,
                                rng=MersenneTwister(1))
    X = (x, y, z)
    rr = sqrt(x^2 + y^2 + z^2)
    prof = _testProfiles()
    P = Dict(k => f(rr) for (k, f) in prof)
    n = [x / rr, y / rr, z / rr]
    Cijkl = tiElasticTensor(P[:A], P[:C], P[:F], P[:L], P[:N], n)
    Y, l = harmonic === :Y32 ? ((x^2 - y^2) * z / rr^3, 3) :
           harmonic === :Y21 ? (x * z / rr^2, 2) :
           error("harmonic must be :Y32 or :Y21")
    ∇₁Y = [expand_derivatives(Differential(X[i])(Y)) * rr for i in 1:3]
    Cvec = -_cross(n, ∇₁Y)
    u = [P[:W] * c for c in Cvec]
    Lu = expand_derivatives.(cartesianElasticOperator(u, Cijkl, P[:ρ], ω, X))

    # radial operator with the same profiles, evaluated in the variable x ≡ r
    Rx = toroidalRadialStrongForm(prof[:W](x), prof[:ρ](x), prof[:L](x),
                                  prof[:N](x); ω=ω, Λ=l * (l + 1))
    maxabs, maxref = 0.0, 0.0
    for _ in 1:npoints
        pt = Dict(x => 2rand(rng) - 1, y => 2rand(rng) - 1, z => 0.2 + rand(rng))
        rv = sqrt(pt[x]^2 + pt[y]^2 + pt[z]^2)
        Rv = _num(substitute(Rx, Dict(x => rv)))
        for i in 1:3
            lhs = _num(substitute(Lu[i], pt))
            rhs = Rv * _num(substitute(Cvec[i], pt))
            maxabs = max(maxabs, abs(lhs - rhs))
            maxref = max(maxref, abs(lhs))
        end
    end
    return (harmonic=harmonic, l=l, maxabs=maxabs, maxrel=maxabs / maxref)
end

# ---------------------------------------------------------------------------
# (1d) radial equation as a famousEquation entry
# ---------------------------------------------------------------------------

"""
famousEquation entry `"1DsphericalSH_TI_Freq"` (coordinate x ≡ r).

    ∂x(Lr² ∂x W) - ∂x(Lr W) + Lr ∂x W + ω² ρr² W - Nₗ W   (= -fr²)

with the *geometric material fields* (precomputed on the grid, per l):

    ρr² = ρ r²,   Lr² = L r²,   Lr = L r,   Nₗ = L + (l(l+1) - 2) N,
    fr² = r² f_W,   f_W(r) = ∫ f · C*_lm dΩ.

Why these fields: PDECoefFinder treats every coefficient as (material field)
× (field derivative), and naturalWeakForm integrates by parts only terms of
the form constant × ∂x(flux).  Absorbing r into the materials keeps both
assumptions true and gives the DSM weak form ω²T - H literally.  Folding Λ
into Nₗ keeps the semi-symbolic recipe independent of l: one recipe, and only
the numerical material arrays change with l.
"""
function flexOPT.famousEquation(::Val{:eq_1DsphericalSH_TI_Freq})
    @variables ω
    @variables W(x) ρr²(x) Lr²(x) Lr(x) Nₗ(x) fr²(x)
    expr = ∂x(Lr² * ∂x(W)) - ∂x(Lr * W) + Lr * ∂x(W) + ω^2 * ρr² * W - Nₗ * W
    exprs = expr
    fields = W
    vars = (ρr², Lr², Lr, Nₗ, ω)
    extexprs = fr²
    extfields = fr²
    extvars = fr²
    coordinates = (x,)
    ∂, ∂² = usefulPartials(coordinates)
    return exprs, fields, vars, extexprs, extfields, extvars, coordinates, ∂, ∂²
end

"""
    sphericalSHMaterialFields(r, ρ, L, N, l)

Numerical geometric material arrays for `eq_1DsphericalSH_TI_Freq` on the
radial grid `r` (L = ρ Vsv², N = ρ Vsh² for a TI model such as PREM).
"""
function sphericalSHMaterialFields(r::AbstractVector, ρ, L, N, l::Integer)
    Λ = l * (l + 1)
    return (ρr²=ρ .* r .^ 2, Lr²=L .* r .^ 2, Lr=L .* r, Nₗ=L .+ (Λ - 2) .* N)
end
