# dev/spherical1DelasticKernel/lib/sgtChecks1D.jl
#
# Checks of the SGT machinery (used by stages/k1_checks.jl and k3_born.jl,
# and run in a standalone harness when the lib was written):
#   checkModeFieldsFD      per-(l,m) strain formulas vs Cartesian finite
#                          differences of the displacement (pure maths);
#   checkReciprocity       u_n(R) of a moment tensor at S two ways: TSGT-side
#                          displacement disk vs RSGT-side strain disk : M, for
#                          generic S and R (tests disks, source coefficients,
#                          rotations); plus the disk vs pointSourceSynthetics;
#   checkInterpolation     on-the-fly (interpolated) fields vs fields computed
#                          exactly at off-grid radii / colatitudes;
#   bornBruteForceCheck    ∫ K δm dV over a shell vs the central difference of
#                          synthetics of the perturbed models.

using LinearAlgebra, Printf, Random

"""
    checkModeFieldsFD(; lmax=6, trials=4, h=1e-5, rng=MersenneTwister(3))

Worst relative difference between `modeFields` strains and centred finite
differences of the Cartesian displacement field (should be ~1e-9).
"""
function checkModeFieldsFD(; lmax=6, trials=4, h=1e-5, rng=MersenneTwister(3), verbose=true)
    fU(r) = (0.3 + 0.2im) * r + 0.05 * r^2;  dfU(r) = (0.3 + 0.2im) + 0.1r
    fV(r) = -0.4 * r^1.5 + 0.1im;            dfV(r) = -0.6 * r^0.5
    fW(r) = 0.7 * sin(r) + 0.2im * r;        dfW(r) = 0.7 * cos(r) + 0.2im
    rad(r) = (fU(r), dfU(r), fV(r), dfV(r), fW(r), dfW(r))
    function uCart(x, l, m)
        r = norm(x); θ, φ = sphAngles(x)
        u, _ = modeFields(rad(r), r, legendreTable2(l, 2, θ; θclamp=0.0), l, m)
        return localBasis(θ, φ) * (collect(u) .* cis(m * φ))
    end
    worst = 0.0
    for l in 0:lmax, m in -min(l, 2):min(l, 2), trial in 1:trials
        r = 1.0 + rand(rng); θ = trial == 1 ? 0.03 : π * rand(rng); φ = 2π * rand(rng)
        x = r * sphUnit(θ, φ)
        G = zeros(ComplexF64, 3, 3)
        for j in 1:3
            e = zeros(3); e[j] = h
            G[:, j] = (uCart(x + e, l, m) - uCart(x - e, l, m)) / (2h)
        end
        E = localBasis(θ, φ)
        εfd = E' * ((G + transpose(G)) / 2) * E
        _, e = modeFields(rad(r), r, legendreTable2(l, 2, θ; θclamp=0.0), l, m)
        εan = strainMatrix(collect(e) .* cis(m * φ))
        worst = max(worst, norm(εan - εfd) / max(norm(εfd), 1e-12))
    end
    verbose && @printf("modeFields vs finite differences (l ≤ %d, |m| ≤ 2): worst relative error %.2e\n", lmax, worst)
    return worst
end

"""
    checkReciprocity(dS, dR, S, M, Rcv; psv=nothing, sh=nothing, p=2)

dS: :moment disk at rs = S.r containing the radius Rcv.r; dR: :force disk at
rs = Rcv.r containing S.r.  Prints and returns the worst relative difference
of u_n(R) = TSGT-side displacement · n vs ε^R(S) : M over components Z, R, T
and all frequencies (exact discrete reciprocity: ~1e-12).  With `psv`, `sh`
also compares the displacement disk with `pointSourceSynthetics` (independent
code, pole frame; differences = θ interpolation of the disk).
"""
function checkReciprocity(dS, dR, S::PointGeo, M, Rcv::PointGeo; psv=nothing, sh=nothing, p=2,
                          omegaiTlen=nothing, verbose=true)
    worst = 0.0
    for comp in (:Z, :R, :T), iω in eachindex(dS.ω)
        n = receiverDirection(comp, S, Rcv)
        u1 = dot(collect(n), sgtAt(dS, iω, S, M, Rcv; which=:disp))
        u2 = sum(sgtAt(dR, iω, Rcv, n, S; which=:strain) .* mtMatrix(M))
        worst = max(worst, abs(u1 - u2) / max(abs(u1), 1e-300))
    end
    verbose && @printf("reciprocity u_n(R): TSGT-side vs RSGT-side, worst relative difference %.2e\n", worst)
    worstRef = NaN
    if !isnothing(psv)
        Rs, Mp = sourceToPole(:moment, M, S); θp, φp, T = frameOf(Rs, Rcv)
        oi = isnothing(omegaiTlen) ? exp(-dS.ωI * dS.tlen) : omegaiTlen
        ref = pointSourceSynthetics(psv, sh, psv.R; rs_km=S.r, sources=Dict("MT" => pointSource(M=Mp)),
            radial=(method=:lagrange, p=p), receiversθdeg=[θp], φr=φp, tlen=dS.tlen,
            fmax=maximum(dS.ω) / 2π + 1e-12, omegaiTlen=oi, verbose=false)["MT"]
        worstRef = 0.0
        for (iω, ω) in enumerate(dS.ω)
            k = findfirst(x -> abs(x - ω) < 1e-9, ref.ω)
            isnothing(k) && continue
            uref = T * [ref.ur[k, 1], ref.uθ[k, 1], ref.uφ[k, 1]]
            worstRef = max(worstRef, norm(sgtAt(dS, iω, S, M, Rcv; which=:disp) - uref) / norm(uref))
        end
        verbose && @printf("disk displacement vs pointSourceSynthetics (rotated): worst relative difference %.2e\n", worstRef)
    end
    return (recip=worst, vsSynthetics=worstRef)
end

"""
    checkInterpolation(disk, exactDisk; iω=length(disk.ω), method=:cubic, rmin=0.0, θmin=2.0)

Relative error of the on-the-fly strain (and displacement) from `disk`
(interpolants) at the radii and colatitudes of `exactDisk` (same source,
same frequencies, computed exactly at off-grid points), excluding r < rmin
and θ < θmin, and points closer than dmin (km) to the source.  Returns (strain, disp) errors =
max |interp - exact| / max |exact| per pair, worst over pairs.
"""
function checkInterpolation(disk, exactDisk; iω=length(disk.ω), method=:cubic, rmin=0.0, θmin=2.0,
                            dmin=0.0, verbose=true)
    itp = diskInterpolants(disk, iω; method=method)
    ok = [i for i in eachindex(exactDisk.r_km) if exactDisk.r_km[i] >= rmin]
    θok = [j for j in eachindex(exactDisk.θdeg) if exactDisk.θdeg[j] >= θmin]
    rs = disk.rs_km
    far(r, θ) = sqrt(r^2 + rs^2 - 2r * rs * cosd(θ)) >= dmin
    errS = 0.0; errD = 0.0
    a(s, ns) = (v = zeros(ns); v[s] = 1.0; Tuple(v))
    ns = disk.family === :moment ? 6 : 3
    for which in (:strain, :disp)
        pairs = which === :strain ? disk.pairs.strain : disk.pairs.disp
        for (ip, (s, c)) in enumerate(pairs)
            ex = which === :strain ? exactDisk.strain[iω, ip, ok, θok] : exactDisk.disp[iω, ip, ok, θok]
            mask = [far(exactDisk.r_km[i], exactDisk.θdeg[j]) for i in ok, j in θok]
            any(mask) || continue
            nrm = maximum(abs.(ex) .* mask); nrm == 0 && continue
            e = 0.0
            for (a_, i) in enumerate(ok), (b_, j) in enumerate(θok)
                mask[a_, b_] || continue
                F = poleFieldAt(itp, a(s, ns), exactDisk.r_km[i], exactDisk.θdeg[j], 0.0;
                                which=which, side=exactDisk.side[i])
                v = which === :strain ? strainVector(F)[c] : F[c]
                e = max(e, abs(v - ex[a_, b_]))
            end
            which === :strain ? (errS = max(errS, e / nrm)) : (errD = max(errD, e / nrm))
        end
    end
    verbose && @printf("interpolation (%s): strain %.2e, displacement %.2e (relative to max per pair)\n",
                       method, errS, errD)
    return (strain=errS, disp=errD)
end

"""
    bornBruteForceCheck(; buildModel, buildSystems, S, M, Rcv, perturbation, ε=1e-3,
        r_nodes, r_weights, δ, tlen, fmax, omegaiTlen, θstep=0.25, nθ, nφ, iωs,
        interp=nothing, p=2)

`buildModel(x)` returns the model for the perturbation amplitude x (x = 0:
reference), `buildSystems(model)` its (psv, sh); `δ(r)` the first-order
perturbation (NamedTuple of ρ, A, C, F, L, N) for x = ε.  The kernel is
integrated over the shell (r_nodes, r_weights; Gauss–Legendre nθ × uniform
nφ) and compared with (u(+ε) - u(-ε))/2 from `pointSourceSynthetics` in the
source frame.  `interp = :cubic | :linear | :lagrange` integrates with
on-the-fly fields from disks stored on the grid (`diskRadii`) instead of
disks computed exactly at r_nodes.
"""
function bornBruteForceCheck(; buildModel, buildSystems, S::PointGeo, M, Rcv::PointGeo, ε=1e-3,
                             r_nodes, r_weights, δ::Function, tlen, fmax, omegaiTlen,
                             θstep=0.25, nθ, nφ, iωs=nothing, interp=nothing, stride=1, p=2,
                             verbose=true)
    model = buildModel(0.0); psv, sh = buildSystems(model)
    R = psv.R
    θg = 0:θstep:180
    radii = isnothing(interp) ? vcat(r_nodes, [Rcv.r, S.r]) : diskRadii(psv; stride=stride)
    dS = computeSGTDisk(psv, sh, model; family=:moment, rs_km=S.r, radii=radii, θdeg=θg, tlen=tlen,
                        fmax=fmax, omegaiTlen=omegaiTlen, p=p, verbose=false)
    dR = computeSGTDisk(psv, sh, model; family=:force, rs_km=Rcv.r, radii=radii, θdeg=θg, tlen=tlen,
                        fmax=fmax, omegaiTlen=omegaiTlen, p=p, verbose=false)
    Rs, Mp = sourceToPole(:moment, M, S); θp, φp, T = frameOf(Rs, Rcv)
    function synth(x)
        mdl = buildModel(x); ps, ss = buildSystems(mdl)
        sp = pointSourceSynthetics(ps, ss, R; rs_km=S.r, sources=Dict("MT" => pointSource(M=Mp)),
            radial=(method=:lagrange, p=p), receiversθdeg=[θp], φr=φp, tlen=tlen, fmax=fmax,
            omegaiTlen=omegaiTlen, verbose=false)["MT"]
        return [T * [sp.ur[i, 1], sp.uθ[i, 1], sp.uφ[i, 1]] for i in eachindex(sp.ω)]
    end
    up = synth(ε); um = synth(-ε); u0 = synth(0.0)
    iωs = isnothing(iωs) ? eachindex(dS.ω) : iωs
    rows = NamedTuple[]
    for iω in iωs
        itps = isnothing(interp) ? nothing :
               (S=diskInterpolants(dS, iω; method=interp), R=diskInterpolants(dR, iω; method=interp))
        for comp in (:Z, :R, :T)
            n = receiverDirection(comp, S, Rcv)
            bf = dot(collect(n), (up[iω] - um[iω]) / 2)
            bo = shellBornIntegral(dS, dR, iω, S, M, Rcv, n; r_nodes=r_nodes, r_weights=r_weights,
                                   δ=δ, nθ=nθ, nφ=nφ, interp=itps)
            rel = abs(bo - bf) / abs(bf)
            push!(rows, (f=dS.ω[iω] / 2π, comp=comp, brute=bf, born=bo, rel=rel,
                         size=abs(bf) / abs(dot(collect(n), u0[iω]))))
            verbose && @printf("  f=%.4f %s  brute %+.5e%+.5ei  Born %+.5e%+.5ei  rel.diff %.2e  (|δu|/|u| %.1e)\n",
                               dS.ω[iω] / 2π, comp, real(bf), imag(bf), real(bo), imag(bo), rel, rows[end].size)
        end
    end
    return rows
end
