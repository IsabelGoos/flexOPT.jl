# dev/spherical1DelasticKernel/stages/k4_traveltime.jl
# K4: first 3-D kernels -- cross-correlation traveltime kernels (ZC2011 eqs 9-12) of the direct
# S wave on the transverse component, PREM, source 600 km deep (vertical strike-slip, Mθφ),
# receiver at Δ = 60° due east, periods ≥ 1/k4fmax.  K = Σ_k Re(c_k K_k(Q)) with the waveform
# kernels K_k and the weights c_k of the reference trace (traveltimeWeights).
#   * check: cross-correlation delay of brute-force synthetics (δln vs = ±0.5 % in PREM zone
#     3630-5600 km) vs the linear prediction Re Σ c_k δu_k;
#   * figures: δln β, δln α, δln ρ kernels on the great-circle section and at 300/1000/2000 km.
# First computed in the cloud harness with DSMconv (same lib), here with the scheme k4Scheme.

k4Scheme = "OPT3"; k4h = 10.0
k4tlen = 2048.0; k4fmax = 0.02; k4Oi = 1e-3
k4Stride = 3; k4θstep = 0.5
k4Window = (940.0, 1060.0)        # direct S window (s), as in the cloud run; check k4_trace.png
prem = planet1D.getSet1Dmodel!(:PREM)
Rp = prem.averagedPlanetRadiusInKilometer
k4S = PointGeo(Rp - 600.0, π / 2, 0.0)
k4R = PointGeo(Rp, π / 2, deg2rad(60.0))
k4M = (0.0, 0.0, 0.0, 0.0, 1.0, 0.0)
k4n = receiverDirection(:T, k4S, k4R)
sys = kerSystems("PREM", k4Scheme, k4h)
psv, sh = sys["psv"], sys["sh"]
k4Disks = Dict{Symbol,Any}()
for (family, rs) in ((:moment, k4S.r), (:force, k4R.r))
    params = Dict{String,Any}("stage" => "k4_disks", "codeVersion" => CODE_VERSION_K["k4_disks"],
        "family" => String(family), "rs" => rs, "scheme" => k4Scheme, "h" => k4h, "tlen" => k4tlen,
        "fmax" => k4fmax, "omegaiTlen" => k4Oi, "stride" => k4Stride, "θstep" => k4θstep, "rmin" => 3480.0)
    devRedo!("disks", params, KER_CACHE_DIR, "k4_disks")
    k4Disks[family] = myProduceOrLoad(params, KER_CACHE_DIR, "k4_disks") do cfg
        Dict("disk" => computeSGTDisk(psv, sh, prem; family=family, rs_km=rs,
                 radii=diskRadii(psv; stride=k4Stride, rmin=3480.0), θdeg=0:k4θstep:180,
                 tlen=k4tlen, fmax=k4fmax, omegaiTlen=k4Oi, p=sourceOrderOf(k4Scheme)))
    end["disk"]
end
dS, dR = k4Disks[:moment], k4Disks[:force]
uω = [dot(collect(k4n), sgtAt(dS, k, k4S, k4M, k4R; which=:disp)) for k in eachindex(dS.ω)]
tw = traveltimeWeights(uω, dS.ω, dS.ωI, dS.tlen; window=k4Window)
let fig = Figure(size=(1100, 300)), ax = Axis(fig[1, 1]; xlabel="time (s)", title="reference T trace and S window")
    lines!(ax, tw.t, tw.u; color=:black); lines!(ax, tw.t, tw.W .* maximum(abs, tw.u); color=:red)
    devSave(fig, "k4_trace")
end

println("\n---- check of the traveltime weights against brute-force cross-correlation")
function k4Scaled(model, ε)
    Vsv = copy(model.C_Vsv); Vsh = copy(model.C_Vsh)
    z = findfirst(i -> model.bottomRadius[i] == 3630.0, 1:model.nzone)
    Vsv[z, :] .*= (1 + ε); Vsh[z, :] .*= (1 + ε)
    m2 = deepcopy(model); m2.C_Vsv .= Vsv; m2.C_Vsh .= Vsh
    return m2
end
Rs, Mp = sourceToPole(:moment, k4M, k4S); θp, φp, Tm = frameOf(Rs, k4R)
function k4SynthT(model)
    p = elasticGlobalSystem(:PSV, model, elasticRadialGrid(model; h=k4h, region=:all); scheme=k4Scheme, configs=shConfigs, recipes=kerRecipes)
    s = elasticGlobalSystem(:SH, model, elasticRadialGrid(model; h=k4h, region=:outerSolid); scheme=k4Scheme, configs=shConfigs, recipes=kerRecipes)
    sp = pointSourceSynthetics(p, s, Rp; rs_km=k4S.r, sources=Dict("MT" => pointSource(M=Mp)),
        radial=(method=:lagrange, p=sourceOrderOf(k4Scheme)), receiversθdeg=[θp], φr=φp, tlen=k4tlen,
        fmax=k4fmax, omegaiTlen=k4Oi, verbose=false)["MT"]
    [dot(collect(k4n), Tm * [sp.ur[i, 1], sp.uθ[i, 1], sp.uφ[i, 1]]) for i in eachindex(sp.ω)]
end
let ε = 0.005
    up = k4SynthT(k4Scaled(prem, ε)); um = k4SynthT(k4Scaled(prem, -ε))
    tr(u) = traveltimeWeights(u, dS.ω, dS.ωI, dS.tlen; window=k4Window).u
    dTp = ccDelay(tw.t, tw.u, tr(up), tw.W); dTm = ccDelay(tw.t, tw.u, tr(um), tw.W)
    @printf("cross-correlation (dT+ - dT-)/2 = %.4f s;  linear Re Σ c δu = %.4f s\n",
            (dTp - dTm) / 2, real(sum(tw.c .* (up .- um) ./ 2)))
end

println("\n---- kernels on the section and maps")
kern(Q) = traveltimeKernelAt(dS, dR, k4S, k4M, k4R, k4n, Q, tw.c; model=prem)
lons = collect(-20.0:0.5:80.0); rs = collect(3490.0:20.0:6370.0)
@time sec = [kern(PointGeo(r, π / 2, deg2rad(lo))) for r in rs, lo in lons]
for (k, lab) in ((:β, "dlnVs"), (:α, "dlnVp"), (:ρ, "dlnRho"))
    K = [getfield(x, k) for x in sec]; lim = quantile_abs(K, 0.99)
    fig = Figure(size=(1100, 450)); ax = Axis(fig[1, 1]; xlabel="longitude along the equator (°)", ylabel="radius (km)",
        title="traveltime kernel $(lab), S on T, Δ = 60°, $(k4Scheme)")
    hm = heatmap!(ax, lons, rs, permutedims(K); colormap=:balance, colorrange=(-lim, lim))
    Colorbar(fig[1, 2], hm); devSave(fig, "k4_section_$(lab)")
end
lats = collect(-45.0:1.0:45.0); lonm = collect(-30.0:1.0:90.0)
for dep in (300.0, 1000.0, 2000.0)
    A = [kern(PointGeo(Rp - dep, deg2rad(90 - la), deg2rad(lo))) for la in lats, lo in lonm]
    K = [x.β for x in A]; lim = quantile_abs(K, 0.995)
    fig = Figure(size=(800, 500)); ax = Axis(fig[1, 1]; aspect=DataAspect(), xlabel="longitude (°)", ylabel="latitude (°)",
        title="dlnVs traveltime kernel at $(Int(dep)) km")
    hm = heatmap!(ax, lonm, lats, permutedims(K); colormap=:balance, colorrange=(-lim, lim))
    Colorbar(fig[1, 2], hm); devSave(fig, "k4_map_dlnVs_$(Int(dep))km")
end
