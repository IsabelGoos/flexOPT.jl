# src/motorsOPT/testFunctionMoments.jl
#
# Moment theory of the OPT test function W for the 1-D divergence-form operator
#
#     d/dx( κ(x) dT/dx ) = f(x)            (natural weak form, -∫W'κ̃T̃' = ∫W f̃)
#
# on a P-point stencil around the test-function centre ν, with T̃, f̃ (and κ̃)
# the Lagrange interpolants on the stencil nodes. In 1-D this is exactly the
# operator flexOPT builds: the hierarchical Taylor inverse with exact degree
# P-1 uses up all nodal values, so every μ / Y_μ gives the same interpolant.
#
# Integrating by parts once more (W vanishes at the ends of its support),
#
#     τ = ∫ W g ds / ∫ W ,   g(s) = (κ̃ T̃_s')'_s / h² - f̃(s)   (a polynomial in s),
#
# so the truncation error depends on W only through its normalised moments
#
#     M_m = ∫ W s^m ds / ∫ W ds        (s in units of Δ, origin at ν).
#
# τ is bilinear in κ and T. Taking κ = s^a/a!, T = s^b/b! (b ≥ 1) makes the
# contribution exactly homogeneous of order h^(a+b-2):
#
#     τ_(a,b) = Σ_m M_m g_m^(a,b)        at order q = a + b - 2.
#
# A scheme is of order ≥ q₀ iff τ_(a,b) = 0 for every (a,b) with a+b-2 < q₀.
# These are LINEAR conditions on the moments, hence on the weights of a
# mixture of test functions W = Σ w_i W_i (each W_i normalised to ∫W_i = 1).
#
# A mass / reaction term ρ ω² T (SH frequency domain: ρω²T + (μT')' = f) adds
# rows ρ = s^c/c!, T = s^b/b! at order h^(c+b): g = ρ̃ T̃ - I(ρT). For centred
# stencils they start at h^(P+1) for every symmetric W, so they never limit the
# order below the stiffness term (keyword `mass=true`).
#
# Everything is exact (Rational{BigInt}); no Symbolics is needed.
#
# Known results (Poisson, nodes centred on ν):
#   P = 3: the hat is the only W with O(h⁴) for constant κ; no W gives better
#          than O(h²) once κ varies.
#   P = 5: W = 4/15 hat₁ + 11/15 hat₂ gives O(h⁶) for any κ (hat: O(h⁴)).
#   P = 7: W = 40/231 hat₁ + 191/231 hat₂ gives O(h⁸) (hat: O(h⁶)).
# The 2-D tensor product w(x)w(y) of the 1-D mixture has been verified
# numerically (flexOPT n ≤ 32, separable prototype to n = 192), not proven.
#
# Self-contained module (no dependency). Claude-authored, November 2026;
# validation script: helpers/testTestFunctionMoments.jl.

module TestFunctionMoments

export PiecewisePoly, bsplineTestFunction, hatTestFunction, centredBspline,
       mixTestFunctions, testFunctionMoments, defaultStencilNodes,
       truncationRows, predictedOrder, optimalTestFunctionWeights,
       flexOPTRefPoints, evaluateTestFunction

const Q = Rational{BigInt}

_q(x::Integer) = Q(x)
_q(x::Rational) = Q(x)
_q(x::AbstractFloat) = Q(rationalize(BigInt, x))

# ---------------------------------------------------------------- polynomials
# coefficient vectors in increasing degree: c[1] + c[2] s + c[3] s² + ...

function _padd(a::Vector{Q}, b::Vector{Q})
    c = zeros(Q, max(length(a), length(b)))
    c[1:length(a)] .+= a
    c[1:length(b)] .+= b
    return c
end

function _pmul(a::Vector{Q}, b::Vector{Q})
    c = zeros(Q, length(a) + length(b) - 1)
    for i in eachindex(a), j in eachindex(b)
        c[i+j-1] += a[i] * b[j]
    end
    return c
end

_pder(c::Vector{Q}) = length(c) <= 1 ? Q[0] : Q[(i - 1) * c[i] for i in 2:length(c)]
_pint(c::Vector{Q}) = vcat(Q[0], Q[c[i] / i for i in eachindex(c)])

function _peval(c::Vector{Q}, s)
    acc = zero(Q)
    for i in length(c):-1:1
        acc = acc * s + c[i]
    end
    return acc
end

_pint_def(c::Vector{Q}, a, b) = (C = _pint(c); _peval(C, b) - _peval(C, a))
_monomial(m::Integer) = vcat(zeros(Q, m), Q[1])

# ---------------------------------------------------------------- test functions

"""
Piecewise polynomial test function in the local coordinate s (units of Δ,
origin at the test-function centre ν): `pieces[k]` holds on
`[breaks[k], breaks[k+1]]`, zero outside `[breaks[1], breaks[end]]`.
"""
struct PiecewisePoly
    breaks::Vector{Q}
    pieces::Vector{Vector{Q}}
end

"""
    bsplineTestFunction(knots)

B-spline of degree `length(knots) - 2` on the given knots (relative to ν, in
Δ units; may be non-integer or repeated), by the Cox-de Boor recursion in
exact arithmetic. Not normalised (use `testFunctionMoments`, which normalises).
"""
function bsplineTestFunction(knots)
    t = _q.(collect(knots))
    p = length(t) - 2
    p >= 0 || throw(ArgumentError("a B-spline needs at least two knots"))
    issorted(t) || throw(ArgumentError("knots must be sorted"))
    breaks = unique(t)
    pieces = Vector{Vector{Q}}()
    for k in 1:length(breaks)-1
        a, b = breaks[k], breaks[k+1]
        B = [(t[i] <= a && b <= t[i+1] && t[i] < t[i+1]) ? Q[1] : Q[0] for i in 1:p+1]
        for r in 1:p
            Bnext = Vector{Vector{Q}}(undef, p + 1 - r)
            for i in 1:p+1-r
                term = Q[0]
                d1 = t[i+r] - t[i]
                if !iszero(d1)
                    term = _padd(term, _pmul(Q[-t[i] / d1, 1 / d1], B[i]))
                end
                d2 = t[i+r+1] - t[i+1]
                if !iszero(d2)
                    term = _padd(term, _pmul(Q[t[i+r+1] / d2, -1 / d2], B[i+1]))
                end
                Bnext[i] = term
            end
            B = Bnext
        end
        push!(pieces, B[1])
    end
    return PiecewisePoly(breaks, pieces)
end

"Hat on [-a, a] (degree-1 B-spline with knots -a, 0, a)."
hatTestFunction(a) = bsplineTestFunction([-_q(a), zero(Q), _q(a)])

"Centred cardinal B-spline of degree p with knot spacing d."
centredBspline(p::Integer, d=1) =
    bsplineTestFunction([_q(d) * (Q(j) - Q(p + 1) / 2) for j in 0:p+1])

function _piece_on(W::PiecewisePoly, a::Q, b::Q)
    mid = (a + b) / 2
    for k in 1:length(W.breaks)-1
        W.breaks[k] <= mid <= W.breaks[k+1] && return W.pieces[k]
    end
    return Q[0]
end

_integral(W::PiecewisePoly) =
    sum(_pint_def(W.pieces[k], W.breaks[k], W.breaks[k+1]) for k in eachindex(W.pieces))

"""
    mixTestFunctions(weights, Ws)

Σ wᵢ Wᵢ / ∫Wᵢ: a mixture whose components are each normalised to unit
integral, so the weights are the ones returned by `optimalTestFunctionWeights`.
"""
function mixTestFunctions(weights, Ws)
    breaks = sort(unique(vcat((W.breaks for W in Ws)...)))
    norms = [_integral(W) for W in Ws]
    pieces = Vector{Vector{Q}}()
    for k in 1:length(breaks)-1
        c = Q[0]
        for (w, W, n) in zip(weights, Ws, norms)
            c = _padd(c, (_q(w) / n) .* _piece_on(W, breaks[k], breaks[k+1]))
        end
        push!(pieces, c)
    end
    return PiecewisePoly(breaks, pieces)
end

"Value of W at s (Float64 for plotting)."
function evaluateTestFunction(W::PiecewisePoly, s::Real)
    sq = _q(s)
    (sq < W.breaks[1] || sq > W.breaks[end]) && return 0.0
    for k in 1:length(W.breaks)-1
        W.breaks[k] <= sq <= W.breaks[k+1] && return Float64(_peval(W.pieces[k], sq))
    end
    return 0.0
end

"""
    testFunctionMoments(W; mmax=8)

Normalised moments M_m = ∫ W s^m / ∫ W for m = 0..mmax (M[1] = 1 is M₀).
"""
function testFunctionMoments(W::PiecewisePoly; mmax::Integer=8)
    total = _integral(W)
    iszero(total) && throw(ArgumentError("test function has zero integral"))
    return [sum(_pint_def(_pmul(W.pieces[k], _monomial(m)), W.breaks[k], W.breaks[k+1])
                for k in eachindex(W.pieces)) / total for m in 0:mmax]
end

# ---------------------------------------------------------------- truncation rows

"Stencil nodes relative to ν, in flexOPT's convention (ν = node (P-1)÷2 + 1)."
defaultStencilNodes(P::Integer) = collect(0:P-1) .- (P - 1) ÷ 2

"1-based flexOPT node indices (for `trialFunctionRefPoints`) of knots given relative to ν."
function flexOPTRefPoints(knots, P::Integer)
    ν = (P - 1) ÷ 2 + 1
    points = collect(knots) .+ ν
    all(p -> 1 <= p <= P, points) ||
        throw(ArgumentError("knots $(collect(knots)) leave the $(P)-point stencil"))
    return points
end

function _lagrange(nodes::Vector{Q})
    basis = Vector{Vector{Q}}()
    for (j, xj) in enumerate(nodes)
        c = Q[1]
        for (m, xm) in enumerate(nodes)
            m == j && continue
            c = _pmul(c, Q[-xm / (xj - xm), 1 / (xj - xm)])
        end
        push!(basis, c)
    end
    return basis
end

function _interpolate(poly::Vector{Q}, nodes::Vector{Q}, basis)
    c = Q[0]
    for (x, b) in zip(nodes, basis)
        c = _padd(c, _peval(poly, x) .* b)
    end
    return c
end

_taylor_monomial(a::Integer) = vcat(zeros(Q, a), Q[1 // factorial(big(a))])

"""
    truncationRows(nodes, q; mmax, kappa=:interpolated, source=:interpolated,
                   mass=false, rho=:interpolated)

Rows `(a, b) => g` with a + b - 2 = q: the order-hᵠ truncation error for
κ = s^a/a!, T = s^b/b! is Σ_m M_m g[m+1]. `kappa`/`source` = `:exact` use the
exact κ or f instead of their interpolants (to separate the error sources).
With `mass=true`, rows `(:mass, c, b) => g` (c + b = q, c ≥ 1) are added for a
term ρ ω² T with ρ = s^c/c! (multiplied by ω², which does not change the order).
"""
function truncationRows(nodes, q::Integer; mmax::Integer,
                        kappa::Symbol=:interpolated, source::Symbol=:interpolated,
                        mass::Bool=false, rho::Symbol=:interpolated)
    nq = _q.(collect(nodes))
    basis = _lagrange(nq)
    rows = Pair{Any,Vector{Q}}[]
    for a in 0:q+1
        b = q + 2 - a
        b >= 1 || continue
        κ = _taylor_monomial(a)
        T = _taylor_monomial(b)
        f = _pder(_pmul(κ, _pder(T)))
        κt = kappa === :exact ? κ : _interpolate(κ, nq, basis)
        Tt = _interpolate(T, nq, basis)
        ft = source === :exact ? f : _interpolate(f, nq, basis)
        g = _padd(_pder(_pmul(κt, _pder(Tt))), -ft)
        row = [m + 1 <= length(g) ? g[m+1] : zero(Q) for m in 0:mmax]
        any(!iszero, g[mmax+2:end]) &&
            throw(ArgumentError("mmax = $mmax too small for q = $q"))
        push!(rows, (a, b) => row)
    end
    if mass
        for c in 1:q
            b = q - c
            ρ = _taylor_monomial(c)
            T = _taylor_monomial(b)
            ρt = rho === :exact ? ρ : _interpolate(ρ, nq, basis)
            Tt = _interpolate(T, nq, basis)
            ρT = _pmul(ρ, T)
            target = source === :exact ? ρT : _interpolate(ρT, nq, basis)
            g = _padd(_pmul(ρt, Tt), -target)
            any(!iszero, g[mmax+2:end]) &&
                throw(ArgumentError("mmax = $mmax too small for q = $q"))
            push!(rows, (:mass, c, b) => [m + 1 <= length(g) ? g[m+1] : zero(Q) for m in 0:mmax])
        end
    end
    return rows
end

_mmax(nodes, qmax) = 2 * length(nodes) + qmax + 2

"""
    predictedOrder(nodes, W; qmax=12, kwargs...) -> (order, leading)

First order q at which some τ_(a,b) ≠ 0; `leading` lists those (a, b) => τ.
`order = qmax + 1` means exact up to qmax.
"""
function predictedOrder(nodes, W::PiecewisePoly; qmax::Integer=12, kwargs...)
    mmax = _mmax(nodes, qmax)
    M = testFunctionMoments(W; mmax)
    for q in 0:qmax
        leading = [(ab => sum(row .* M)) for (ab, row) in truncationRows(nodes, q; mmax, kwargs...)]
        filter!(p -> !iszero(p.second), leading)
        isempty(leading) || return (order=q, leading=leading)
    end
    return (order=qmax + 1, leading=Pair{Any,Q}[])
end

# ---------------------------------------------------------------- exact linear algebra

function _rref!(A::Matrix{Q})
    nr, nc = size(A)
    pivots = Int[]
    r = 1
    for c in 1:nc
        r > nr && break
        p = findfirst(i -> !iszero(A[i, c]), r:nr)
        p === nothing && continue
        p += r - 1
        A[[r, p], :] = A[[p, r], :]
        A[r, :] = A[r, :] ./ A[r, c]
        for i in 1:nr
            if i != r && !iszero(A[i, c])
                A[i, :] = A[i, :] .- A[i, c] .* A[r, :]
            end
        end
        push!(pivots, c)
        r += 1
    end
    return A, pivots
end

"Solve A x = b exactly: (consistent, particular solution, null-space basis)."
function _solve_exact(A::Matrix{Q}, b::Vector{Q})
    n = size(A, 2)
    R, pivots = _rref!(hcat(A, b))
    (n + 1) in pivots && return (false, nothing, nothing)
    x = zeros(Q, n)
    for (row, c) in enumerate(pivots)
        x[c] = R[row, end]
    end
    free = setdiff(1:n, pivots)
    null = Vector{Vector{Q}}()
    for f in free
        v = zeros(Q, n)
        v[f] = 1
        for (row, c) in enumerate(pivots)
            v[c] = -R[row, f]
        end
        push!(null, v)
    end
    return (true, x, null)
end

"""
    optimalTestFunctionWeights(nodes, candidates; target=nothing, qmax=12, kwargs...)

Weights wᵢ (Σwᵢ = 1, components normalised to unit integral) of
W = Σ wᵢ candidates[i] that make every truncation term vanish below order
`target`. With `target = nothing` the highest achievable order ≤ qmax+1 is
returned. Result: `(order, weights, freeDirections, negative)`; `weights` is
the particular solution with free weights set to zero, `freeDirections` spans
the remaining freedom (usable to cancel the next order or improve the
constant), `negative` flags a non-positive mixture.
"""
function optimalTestFunctionWeights(nodes, candidates::AbstractVector{PiecewisePoly};
                                    target=nothing, qmax::Integer=12, kwargs...)
    mmax = _mmax(nodes, qmax)
    Ms = [testFunctionMoments(W; mmax) for W in candidates]
    n = length(candidates)
    function system(order)
        A = reshape(Q.(ones(Int, n)), 1, n)
        b = Q[1]
        for q in 0:order-1
            for (_, row) in truncationRows(nodes, q; mmax, kwargs...)
                A = vcat(A, reshape([sum(row .* M) for M in Ms], 1, n))
                push!(b, zero(Q))
            end
        end
        return _solve_exact(A, b)
    end
    orders = target === nothing ? (qmax+1:-1:0) : (target:target)
    for order in orders
        ok, x, null = system(order)
        ok || continue
        return (order=order, weights=x, freeDirections=null,
                negative=any(<(0), x))
    end
    return (order=-1, weights=nothing, freeDirections=nothing, negative=false)
end

end # module
