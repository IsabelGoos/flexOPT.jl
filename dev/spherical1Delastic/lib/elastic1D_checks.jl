# dev/spherical1Delastic/lib/elastic1D_checks.jl
#
# Exact checks for the P-SV + fluid machinery of lib/elastic1D.jl:
# homogeneous fluid core (0 ≤ ξ ≤ a) + homogeneous isotropic solid mantle
# (a ≤ ξ ≤ 1), traction-free surface, G&O coupling at the CMB.
#   exact: fluid q ∝ j_l(k_f r); solid = P (j_l, y_l) + S (j_l, y_l) potentials;
#   CMB: U continuous, σ_rr = ω q (= P), σ_rθ = 0; surface σ_rr = σ_rθ = 0
#   → 5×5 determinant D(Ω) = 0 (Python prototype: identical mode list and h²
#   convergence of conventional DSM).
#   discrete: (Ω²T + Ω Cpl + H) v = 0 → companion linearisation.

using LinearAlgebra, Printf

"""
    layeredPSVModel(; R=6371.0, radii, ρ, Vp, Vs, Vph=Vp, Vsh=Vs, η=ones, Qμ, Qκ)

In-memory planet1D model in the full "DSM" format (28 numbers per zone),
so fluid zones (Vs = 0) are recognised by planet1D.  Each property per zone:
number or 4 polynomial coefficients in r/R.
"""
function layeredPSVModel(; R=6371.0, radii, ρ, Vp, Vs, Vph=Vp, Vsh=Vs,
                         η=fill(1.0, length(radii) - 1),
                         Qμ=fill(600.0, length(radii) - 1), Qκ=fill(57823.0, length(radii) - 1))
    nz = length(radii) - 1
    arr = Float64[]
    for i in 1:nz
        append!(arr, [radii[i], radii[i+1]])
        for p in (ρ[i], Vp[i], Vph[i], Vs[i], Vsh[i], η[i])
            append!(arr, _coef4(p))
        end
        append!(arr, [Qμ[i], Qκ[i]])
    end
    return planet1D.planet1DPSVmodel(nz, "DSM", Float64(R), arr)
end

_sbj(l, z) = sphericalBesselJ(l, z)
_sby(l, z) = sphericalBesselY(l, z)
_dsb(f, l, z) = l / z * f(l, z) - f(l + 1, z)                      # f_l'(z)
_ddsb(f, l, z) = -(2 / z) * _dsb(f, l, z) - (1 - l * (l + 1) / z^2) * f(l, z)

"(U, σ_rr, σ_rθ) of the four homogeneous isotropic solid solutions at ξ."
function _solidColumns(Ω, l, ξ, ρ, α, β)
    k = sqrt(l * (l + 1)); Λ = l * (l + 1)
    kp, ks = Ω / α, Ω / β; μ = ρ * β^2; λ = ρ * α^2 - 2μ
    cols = NTuple{3,Float64}[]
    for f in (_sbj, _sby)                       # P
        z = kp * ξ; fv = f(l, z); fd = kp * _dsb(f, l, z); fdd = kp^2 * _ddsb(f, l, z)
        U = fd; V = k * fv / ξ; Ud = fdd; Vd = k * (fd / ξ - fv / ξ^2)
        push!(cols, (U, λ * (-kp^2 * fv) + 2μ * Ud, μ * (Vd - V / ξ + k * U / ξ)))
    end
    for f in (_sbj, _sby)                       # S
        z = ks * ξ; fv = f(l, z); fd = ks * _dsb(f, l, z); fdd = ks^2 * _ddsb(f, l, z)
        U = Λ * fv / ξ; Ud = Λ * (fd / ξ - fv / ξ^2)
        V = k * (fv / ξ + fd); Vd = k * (fd / ξ - fv / ξ^2 + fdd)
        push!(cols, (U, 2μ * Ud, μ * (Vd - V / ξ + k * U / ξ)))
    end
    return cols
end

"5×5 determinant (columns normalised) for the fluid core + solid mantle."
function coreMantleDeterminant(Ω, l; a, ρF, αF, ρS, αS, βS)
    kf = Ω / αF; z = kf * a
    qv = _sbj(l, z); qd = kf * _dsb(_sbj, l, z)
    M = zeros(5, 5)
    M[1, 1] = qd / (ρF * Ω)          # -U_fluid, U_fluid = -q'/(ρω)
    M[2, 1] = -Ω * qv
    cs = _solidColumns(Ω, l, a, ρS, αS, βS); cu = _solidColumns(Ω, l, 1.0, ρS, αS, βS)
    for j in 1:4
        M[1, 1+j] = cs[j][1]; M[2, 1+j] = cs[j][2]; M[3, 1+j] = cs[j][3]
        M[4, 1+j] = cu[j][2]; M[5, 1+j] = cu[j][3]
    end
    for j in 1:5
        s = maximum(abs.(M[:, j])); s > 0 && (M[:, j] ./= s)
    end
    return det(M)
end

"Exact eigenfrequencies Ω (ξ units) below Ωmax: sign changes on a fine grid + bisection."
function coreMantleRoots(l; Ωmax, dΩ=1e-3, Ωmin=0.5, kw...)
    f(Ω) = coreMantleDeterminant(Ω, l; kw...)
    Ωs = collect(Ωmin:dΩ:Ωmax); v = f.(Ωs); roots = Float64[]
    for i in 1:length(Ωs)-1
        sign(v[i]) == sign(v[i+1]) && continue
        lo, hi, flo = Ωs[i], Ωs[i+1], v[i]
        for _ in 1:80
            mid = (lo + hi) / 2; fm = f(mid)
            sign(fm) == sign(flo) ? ((lo, flo) = (mid, fm)) : (hi = mid)
        end
        push!(roots, (lo + hi) / 2)
    end
    return roots
end

"""
    quadraticEigenfrequencies(sys, l)

Real positive Ω of (Ω²T + ΩCpl + H0 + kHk + ΛHΛ) v = 0 (companion
linearisation, dense) and the largest |Im/Re| among them.
"""
function quadraticEigenfrequencies(sys::ElasticGlobalSystem, l)
    Λ = l * (l + 1); k = sqrt(Λ)
    T = Matrix(sys.T); H = Matrix(sys.H0 .+ k .* sys.Hk .+ Λ .* sys.HΛ); C = Matrix(sys.Cpl)
    n = size(T, 1); Z = zeros(n, n); Id = Matrix{Float64}(I, n, n)
    λ = eigvals([Z Id; -H -C], [Id Z; Z T])
    λ = λ[isfinite.(λ)]
    pos = λ[real.(λ) .> 1e-8]
    return (Ω=sort(real.(pos)), maxImRatio=isempty(pos) ? 0.0 : maximum(abs.(imag.(pos)) ./ abs.(real.(pos))))
end

"""
    coreMantleEigenBenchmark(schemes; configs, l=2, nIntervals=[30, 60, 120], nModes=3,
                             a=0.546, ρF=10.0, αF=9.0, ρS=4.5, αS=11.0, βS=6.0)

Relative eigenfrequency errors of the first `nModes` modes for every scheme
("DSMconv" or a config name) and grid; the flexOPT systems use the native
clipped half rows at both ends of both zones and the G&O coupling.
"""
function coreMantleEigenBenchmark(schemes; configs, l=2, nIntervals=[30, 60, 120], nModes=3,
                                  a=0.546, ρF=10.0, αF=9.0, ρS=4.5, αS=11.0, βS=6.0,
                                  recipes=Dict{Any,Any}(), nearRows=:fallback)
    R = 6371.0
    model = layeredPSVModel(; R=R, radii=[0.0, a * R, R], ρ=[ρF, ρS], Vp=[αF, αS], Vs=[0.0, βS])
    @show model.solid_or_fluid
    exact = coreMantleRoots(l; Ωmax=60.0, a=a, ρF=ρF, αF=αF, ρS=ρS, αS=αS, βS=βS)
    @printf("exact Ω (l=%d): %s\n", l, join((@sprintf("%.6f", w) for w in exact), "  "))
    length(exact) >= nModes || (nModes = length(exact))
    out = NamedTuple[]
    for n in nIntervals
        grid = elasticRadialGrid(model; h=(a * R) / n, region=:all)
        for sc in schemes
            sys = elasticGlobalSystem(:PSV, model, grid; scheme=sc, configs=configs, recipes=recipes,
                                      nearRows=nearRows)
            ev = quadraticEigenfrequencies(sys, l)
            err = [(ev.Ω[argmin(abs.(ev.Ω .- w))] - w) / w for w in exact[1:nModes]]
            push!(out, (scheme=sc, n=n, relErr=err, maxImRatio=ev.maxImRatio))
            @printf("%-8s n=%4d  max|Im/Re|=%.1e | %s\n", sc, n, ev.maxImRatio,
                    join((@sprintf("%+.2e", e) for e in err), " "))
        end
    end
    return out
end
