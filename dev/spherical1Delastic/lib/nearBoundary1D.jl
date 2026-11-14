# dev/spherical1Delastic/lib/nearBoundary1D.jl
#
# Checks for the near-boundary rows of the OPT schemes (flexOPTZoneBasis,
# nearRows = :fallback | :movingNu | :clippedFamily) on the homogeneous
# traction-free SH shell a ≤ ξ ≤ 1, through the same global assembly as the
# PREM systems (elasticGlobalSystem):
#   * eigenfrequencies vs the exact roots (shellToroidalRoots);
#   * row-by-row residuals of the exact eigenfunction at the exact Ω, which
#     show directly which rows (end row, rows 2, 3, …, interior) carry the
#     error and with which order.
# Requires lib/sphericalTI_SH_step3.jl, lib/elastic1D*.jl, lib/sources1D.jl.

using LinearAlgebra, Printf

_shellModel(; R=6371.0, rcmb=3480.0, ρ=3.3, β=4.5) =
    layeredPSVModel(; R=R, radii=[rcmb, R], ρ=[ρ], Vp=[sqrt(3) * β], Vs=[β])

"Exact toroidal eigenfunction W(ξ) of the shell for wavenumber k (traction-free at a and 1)."
function shellToroidalMode(k, l, a)
    J(x) = sphericalBesselJ(l, k * x);  Jd(x) = k * _dsbj(l, k * x)
    Y(x) = sphericalBesselY(l, k * x);  Yd(x) = k * _dsby(l, k * x)
    cJ = Yd(a) - Y(a) / a; cY = -(Jd(a) - J(a) / a)
    return x -> cJ * J(x) + cY * Y(x)
end

"""
    shellEigenNearRows(schemes, variants; configs, l=10, nIntervals=[50,100,200,400],
                       nModes=8, recipes=Dict(), β=4.5, ...)

Relative eigenfrequency errors of the first `nModes` toroidal modes for
every scheme × near-boundary variant × grid (DSMconv: once).  Generalised
eigenproblem H v = -Ω² T v (dense); complex pairs are counted.
"""
function shellEigenNearRows(schemes, variants; configs, l=10, nIntervals=[50, 100, 200, 400],
                            nModes=8, recipes=Dict{Any,Any}(), R=6371.0, rcmb=3480.0,
                            ρ=3.3, β=4.5)
    model = _shellModel(; R=R, rcmb=rcmb, ρ=ρ, β=β)
    kex = shellToroidalRoots(l, rcmb / R, 1.0; kmax=80.0)
    Ωex = β .* kex[1:nModes]
    @printf("exact Ω (l=%d): %s\n", l, join((@sprintf("%.6f", w) for w in Ωex), " "))
    Λ = l * (l + 1)
    rows = NamedTuple[]
    for n in nIntervals
        grid = elasticRadialGrid(model; h=(R - rcmb) / n, region=:outerSolid, minIntervals=4)
        for sc in schemes, v in (sc == "DSMconv" ? [:fallback] : variants)
            sys = elasticGlobalSystem(:SH, model, grid; scheme=sc, configs=configs,
                                      recipes=recipes, nearRows=v)
            T = Matrix(sys.T); H = Matrix(sys.H0 .+ Λ .* sys.HΛ)
            λ = eigvals(H, -T); λ = λ[isfinite.(λ)]
            good = [abs(imag(z)) <= 1e-8 * abs(real(z)) && real(z) > 0 for z in λ]
            Ω = sort(sqrt.(real.(λ[good])))
            err = [(Ω[argmin(abs.(Ω .- w))] - w) / w for w in Ωex]
            nC = count(!, good)
            push!(rows, (scheme=sc, variant=v, n=n, h_km=(R - rcmb) / n, relErr=err, nComplex=nC))
            @printf("%-8s %-14s n=%4d nC=%2d | %s\n", sc, v, n, nC,
                    join((@sprintf("%+.2e", e) for e in err), " "))
        end
    end
    return rows
end

"""
    shellRowResiduals(schemes, variants; configs, l=10, mode=4, nIntervals=[50,100,200,400],
                      nShow=4, recipes=Dict(), ...)

|row_i · W_exact| / max_i |(Ω²T W_exact)_i| for the exact eigenpair `mode`:
the first and last `nShow` rows individually and the interior maximum.
Rows of a consistent O(h^p) scheme scale like h^p; a row that does not
shows where the order is lost.
"""
function shellRowResiduals(schemes, variants; configs, l=10, mode=4, nIntervals=[50, 100, 200, 400],
                           nShow=4, recipes=Dict{Any,Any}(), R=6371.0, rcmb=3480.0, ρ=3.3, β=4.5)
    model = _shellModel(; R=R, rcmb=rcmb, ρ=ρ, β=β)
    a = rcmb / R
    k = shellToroidalRoots(l, a, 1.0; kmax=80.0)[mode]; Ω = β * k
    W = shellToroidalMode(k, l, a)
    Λ = l * (l + 1)
    rows = NamedTuple[]
    for n in nIntervals
        grid = elasticRadialGrid(model; h=(R - rcmb) / n, region=:outerSolid, minIntervals=4)
        ξ = only(grid.segments).ξ; w = W.(ξ)
        for sc in schemes, v in (sc == "DSMconv" ? [:fallback] : variants)
            sys = elasticGlobalSystem(:SH, model, grid; scheme=sc, configs=configs,
                                      recipes=recipes, nearRows=v)
            dofs = vec(sys.dofmap[1])
            x = zeros(size(sys.T, 1)); x[dofs] .= w
            r = abs.((Ω^2 .* sys.T .+ sys.H0 .+ Λ .* sys.HΛ) * x)[dofs]
            s = maximum(abs.((Ω^2 .* sys.T) * x))
            r ./= s
            nn = length(r)
            bottom = r[1:nShow]; top = r[nn:-1:nn-nShow+1]
            interior = maximum(r[nShow+1:nn-nShow])
            push!(rows, (scheme=sc, variant=v, n=n, bottom=bottom, top=top, interior=interior))
            @printf("%-8s %-14s n=%4d bottom %s | top %s | interior %.1e\n", sc, v, n,
                    join((@sprintf("%.1e", x) for x in bottom), " "),
                    join((@sprintf("%.1e", x) for x in top), " "), interior)
        end
    end
    return rows
end

"Observed orders per (scheme, variant) of a scalar `f(row)` between successive grids."
function nearRowOrders(rows, f; label="")
    for key in unique([(r.scheme, r.variant) for r in rows])
        rr = sort(filter(r -> (r.scheme, r.variant) == key, rows); by=r -> r.n)
        length(rr) < 2 && continue
        e = [f(r) for r in rr]
        @printf("%-8s %-14s %s orders %s\n", key[1], key[2], label,
                join((@sprintf("%5.2f", log2(e[i] / e[i+1])) for i in 1:length(e)-1), " "))
    end
end
