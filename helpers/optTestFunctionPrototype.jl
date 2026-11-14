#!/usr/bin/env julia
#
# helpers/optTestFunctionPrototype.jl   (Claude-authored; independent of src/)
#
# Closed-form 1-D prototype of the natural-weak OPT_P operator for
#     d/dx( κ dT/dx ) = f
# with an ARBITRARY test function W (piecewise polynomial, support inside the
# stencil, W = 0 at both ends):
#
#     -∫ W' κ̃ T̃' dx = ∫ W f̃ dx
#
# T̃, f̃, κ̃ are the Lagrange interpolants on the P stencil nodes. In 1-D this is
# exactly what flexOPT computes (the hierarchical Taylor inverse with exact
# degree P-1 uses up all nodal values, so any μ / Y gives this same
# interpolant); it reproduces flexOPT's OPT3/4/5/6 periodic residuals digit for
# digit (scripts/tmp/poisson_mu_benchmark, wave_fd_vs_opt_weak.jld2).
#
# Moment form of the truncation error (integrate by parts once more):
#     τ = ∫ W g ds / ∫ W,   g(s) = (κ̃ T̃_s')'_s / h² - f̃(s)   (a polynomial)
# so τ depends on W only through M_m = ∫ W s^m / ∫ W. Leading term (symmetric W):
#   P=3:  -[(6M₂-1) κT⁽⁴⁾ + (18M₂-2) κ'T''' + (6M₂-2) κ'''T'] h²/12
#         -> no M₂ cancels all three: O(h²) for any test function once κ varies
#            (hat M₂=1/6 leaves -(κ'T''' - κ'''T')/12, which vanishes for the
#             same-wavelength zero-phase benchmark: the "lucky" O(h⁴)).
#   P=5:  cancelled by M₂ = 8/15, M₄ = 4/5  -> W = 4/15 hat₁ + 11/15 hat₂   (O(h⁶))
#   P=7:  one-parameter family M₂ = M₆/21 + 24/49, M₄ = M₆/5 + 18/35;
#         with hats only: W = 40/231 hat₁ + 191/231 hat₂                    (O(h⁸))
#   (hat_a = hat on [-a,a] with unit integral; derivation:
#    helpers/derivations/optTestFunctionMoments_sympy.py)
#
# Usage:  julia helpers/optTestFunctionPrototype.jl          (no package needed)

module OPTTestFunctionPrototype

using LinearAlgebra
using Printf

export hat, mixture, lagrange_basis, opt_row, periodic_residual, moments,
       MIXTURES, run_prototype

# ---------------------------------------------------------------- polynomials
# coefficient vectors, increasing degree: c[1] + c[2] s + c[3] s² + ...

function padd(a, b)
    c = zeros(max(length(a), length(b)))
    c[1:length(a)] .+= a
    c[1:length(b)] .+= b
    return c
end
function pmul(a, b)
    c = zeros(length(a) + length(b) - 1)
    for i in eachindex(a), j in eachindex(b)
        c[i+j-1] += a[i] * b[j]
    end
    return c
end
pder(c) = length(c) == 1 ? [0.0] : [c[i] * (i - 1) for i in 2:length(c)]
pint(c) = [0.0; [c[i] / i for i in eachindex(c)]]
peval(c, s) = foldr((ci, acc) -> ci + s * acc, c; init=0.0)
pint_def(c, a, b) = (C = pint(c); peval(C, b) - peval(C, a))

function lagrange_basis(nodes)
    basis = Vector{Vector{Float64}}()
    for (j, xj) in enumerate(nodes)
        c = [1.0]
        for (m, xm) in enumerate(nodes)
            m == j && continue
            c = pmul(c, [-xm, 1.0] ./ (xj - xm))
        end
        push!(basis, c)
    end
    return basis
end

# ---------------------------------------------------------------- test functions
# W = list of pieces (left, right, coefficients in s), local coordinate s in units of h

"Hat on [-a, a] with unit integral."
hat(a) = [(-Float64(a), 0.0, [1 / a, 1 / a^2]), (0.0, Float64(a), [1 / a, -1 / a^2])]

"Weighted sum of test functions: mixture((w1, W1), (w2, W2), ...)."
mixture(parts...) = [(l, r, w .* c) for (w, W) in parts for (l, r, c) in W]

"Normalised moments M_m = ∫ W s^m / ∫ W, m = 0..mmax."
function moments(W; mmax=6)
    total = sum(pint_def(c, l, r) for (l, r, c) in W)
    return [sum(pint_def(pmul(c, [zeros(m); 1.0]), l, r) for (l, r, c) in W) / total
            for m in 0:mmax]
end

const MIXTURES = Dict(
    "hat1" => hat(1),
    "hat2" => hat(2),
    "OPT5 mixture 4/15 hat1 + 11/15 hat2" => mixture((4 / 15, hat(1)), (11 / 15, hat(2))),
    "OPT7 mixture 40/231 hat1 + 191/231 hat2" => mixture((40 / 231, hat(1)), (191 / 231, hat(2))),
)

# ---------------------------------------------------------------- operator row

"""
    opt_row(nodes, W) -> (A, G)

A[m, j] = -∫ W' ψ_m φ_j'   (material node m, field node j; multiply by κ_m / h²)
G[j]    =  ∫ W φ_j         (source redistribution Γ)
with φ, ψ the Lagrange bases on `nodes` (field and material on the same nodes).
"""
function opt_row(nodes, W)
    phi = lagrange_basis(nodes)
    P = length(nodes)
    A = zeros(P, P)
    G = zeros(P)
    for m in 1:P, j in 1:P
        g = pmul(phi[m], pder(phi[j]))
        A[m, j] = -sum(pint_def(pmul(pder(c), g), l, r) for (l, r, c) in W)
    end
    for j in 1:P
        G[j] = sum(pint_def(pmul(c, phi[j]), l, r) for (l, r, c) in W)
    end
    return A, G
end

"""
    periodic_residual(P, W; n, q, amp, k=2.0, κ0=2.5)

RMS over the 2π-periodic grid of (A T - Γ f) / Σ|Γ| for T = cos(kx),
κ = κ0 (1 + amp cos(qx)) -- the 1-D scenarios of runWavePeriodicFDvsOPT.jl
(homogeneous amp=0; same_wave q=2; long_material q=1; short_material q=6; amp=0.15).
"""
function periodic_residual(P, W; n, q, amp, k=2.0, κ0=2.5)
    h = 2pi / n
    nodes = collect(0:P-1) .- (P - 1) ÷ 2          # ν at a node (as flexOPT)
    A, G = opt_row(nodes, W)
    κ(x) = κ0 * (1 + amp * cos(q * x))
    dκ(x) = -κ0 * amp * q * sin(q * x)
    u(x) = cos(k * x)
    f(x) = -dκ(x) * k * sin(k * x) - κ(x) * k^2 * cos(k * x)
    r2 = 0.0
    for i in 0:n-1
        x = i * h
        lhs = 0.0
        rhs = 0.0
        for (j, sj) in enumerate(nodes)
            for (m, sm) in enumerate(nodes)
                lhs += A[m, j] * κ(x + sm * h) * u(x + sj * h) / h^2
            end
            rhs += G[j] * f(x + sj * h)
        end
        r2 += ((lhs - rhs) / sum(abs, G))^2
    end
    return sqrt(r2 / n)
end

# ---------------------------------------------------------------- demo

const SCENARIOS = [
    ("homogeneous", 1, 0.0), ("same_wave", 2, 0.15),
    ("long_material", 1, 0.15), ("short_material", 6, 0.15),
]

function run_prototype(; sizes=[16, 24, 32, 48, 64, 96, 128])
    cases = [
        (3, "hat1"), (5, "hat1"), (5, "hat2"),
        (5, "OPT5 mixture 4/15 hat1 + 11/15 hat2"),
        (7, "hat1"), (7, "OPT7 mixture 40/231 hat1 + 191/231 hat2"),
    ]
    for (P, name) in cases
        W = MIXTURES[name]
        M = moments(W)
        @printf("\nP = %d, W = %s   (M2 = %.6f, M4 = %.6f, M6 = %.6f)\n", P, name, M[3], M[5], M[7])
        for (scenario, q, amp) in SCENARIOS
            e = [periodic_residual(P, W; n, q, amp) for n in sizes]
            o = [log(e[i] / e[i+1]) / log(sizes[i+1] / sizes[i]) for i in 1:length(e)-1]
            @printf("  %-15s n=%d: %.4e   n=%d: %.4e   orders %s\n", scenario,
                sizes[1], e[1], sizes[end], e[end], join((@sprintf("%.2f", x) for x in o), " "))
        end
    end
end

end # module

abspath(PROGRAM_FILE) == abspath(@__FILE__) && OPTTestFunctionPrototype.run_prototype()
