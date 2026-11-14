#!/usr/bin/env python3
# helpers/derivations/optTestFunctionMoments_sympy.py   (Claude-authored derivation; not part of flexOPT)
#
# Symbolic truncation error of the 1-D natural-weak OPT_P operator for (κ T')' = f
# with an arbitrary symmetric test function W, written through its normalised
# moments M_m = ∫ W s^m / ∫ W:
#
#   τ = Σ_m M_m c_m ,  g(s) = (κ̃ T̃_s')'_s / h² - f̃(s) = Σ_m c_m s^m
#
# κ̃, T̃, f̃: Lagrange interpolants on the P stencil nodes (exactly what flexOPT
# builds in 1-D). κ and T are general analytic functions: κ = Σ k_i x^i/i!,
# T = Σ t_i x^i/i! around the test-function centre ν.
#
# Requires sympy (pip install sympy).  Run:  python3 optTestFunctionMoments_sympy.py
import sympy as sp

s, h, X = sp.symbols('s h X')
N = 9
k = sp.symbols('k0:%d' % N)
t = sp.symbols('t0:%d' % N)
kap = sum(k[i] * X**i / sp.factorial(i) for i in range(N))
T = sum(t[i] * X**i / sp.factorial(i) for i in range(N))
f = sp.expand(sp.diff(kap * sp.diff(T, X), X))


def interp(func, nodes):
    return sp.expand(sp.interpolate([(n, func.subs(X, n * h)) for n in nodes], s))


def tau(P, exact_kappa=False, exact_f=False):
    nodes = list(range(-(P - 1) // 2, (P - 1) // 2 + 1))
    K = sp.expand(kap.subs(X, s * h)) if exact_kappa else interp(kap, nodes)
    F = sp.expand(f.subs(X, s * h)) if exact_f else interp(f, nodes)
    g = sp.Poly(sp.expand(sp.diff(K * sp.diff(interp(T, nodes), s), s) / h**2 - F), s)
    M = sp.symbols('M0:%d' % (g.degree() + 1))
    # symmetric W: odd moments vanish, M0 = 1
    expr = sum((M[m] if m > 0 else 1) * g.coeff_monomial(s**m) for m in range(0, g.degree() + 1, 2))
    return sp.expand(sp.series(expr, h, 0, P + 1).removeO()), M


if __name__ == "__main__":
    for P in (3, 5, 7):
        expr, M = tau(P)
        lead = sp.expand(expr.coeff(h, P - 1))
        free = [m for m in M[2::2] if lead.has(m)]
        print(f"\nP = {P}: leading term h^{P-1}:\n  {sp.factor(lead)}")
        print("  moments cancelling it for every κ, T:",
              sp.solve(sp.Poly(lead, *k, *t).coeffs(), free, dict=True) or "none")
    for ek, ef in [(True, False), (True, True)]:
        expr, M = tau(3, ek, ef)
        print(f"\nP = 3, exact κ = {ek}, exact f = {ef}: h^2 term {sp.factor(expr.coeff(h, 2))}")
    # mixtures of hats (hat_a on [-a,a], unit integral: M2 = a²/6, M4 = a⁴/15, M6 = a⁶/28)
    w1, w2, w3 = sp.symbols('w1 w2 w3')
    print("\nP = 5 hat mixture:", sp.solve([w1 + w2 - 1, (w1 + 4 * w2) / 6 - sp.Rational(8, 15),
                                            (w1 + 16 * w2) / 15 - sp.Rational(4, 5)], [w1, w2]))
    M2 = (w1 + 4 * w2 + 9 * w3) / 6; M4 = (w1 + 16 * w2 + 81 * w3) / 15; M6 = (w1 + 64 * w2 + 729 * w3) / 28
    print("P = 7 hat mixture:", sp.solve([w1 + w2 + w3 - 1, M2 - (M6 / 21 + sp.Rational(24, 49)),
                                           M4 - (M6 / 5 + sp.Rational(18, 35))], [w1, w2, w3]))
