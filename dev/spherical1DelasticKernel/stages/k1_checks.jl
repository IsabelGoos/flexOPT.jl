# dev/spherical1DelasticKernel/stages/k1_checks.jl
# K1: checks of the SGT machinery on PREM (lib/sgtChecks1D.jl, lib/sgtZC1D.jl):
#   (a) per-(l, m) strain formulas vs Cartesian finite differences (pure maths, ~1e-9);
#   (b) reciprocity: u_n(R) of a moment tensor at S from the :moment disk (displacement)
#       and from the :force disk (ε^R(S) : M), generic S and R -- exact (~1e-12);
#       and the :moment disk vs pointSourceSynthetics (independent code, ~1e-7 = θ interp.);
#   (c) the ZC route (h1-h10, H1-H20 + B10/B11/B14, corrected extraction) vs the
#       rotation route (~1e-13);
#   (d) accuracy of the on-the-fly fields (disk stored on the grid, every `k1Stride`-th
#       node, θ step k1θstep) vs fields computed exactly off the grid.
# Standalone run of the same checks (DSMconv, layered model) when the lib was written:
# (a) 3e-9, (b) 1.3e-13 / 1e-7, (c) 6e-15 / 4e-14.

k1Schemes = ["DSMconv", "OPT3"]
k1h = 10.0
k1Depth = 600.0
k1tlen = 1638.4; k1fmax = 0.01; k1Oi = 1e-4
k1Stride = 2; k1θstep = 0.25
prem = planet1D.getSet1Dmodel!(:PREM)
Rp = prem.averagedPlanetRadiusInKilometer
k1S = pointFromLatLon(Rp, k1Depth, 40.0, 20.0)          # source, 600 km deep
k1R = pointFromLatLon(Rp, 0.0, -15.0, 75.0)              # receiver at the surface
k1M = (1.0, 0.6, -0.8, -0.3, 0.4, -0.7)                  # (Mrr, Mrθ, Mrφ, Mθθ, Mθφ, Mφφ), 1e25 dyn·cm
k1rQ = 5000.0                                            # radius of the ZC-route check
@printf("source-receiver distance %.2f°\n", epicentralDeg(k1S, k1R))

println("\n---- (a) strain formulas vs finite differences")
checkModeFieldsFD()

for scheme in k1Schemes
    println("\n================ scheme ", scheme, " (h = $(k1h) km)")
    sys = kerSystems("PREM", scheme, k1h)
    psv, sh = sys["psv"], sys["sh"]
    params = Dict{String,Any}("stage" => "k1_disks", "codeVersion" => CODE_VERSION_K["k1_disks"],
        "scheme" => scheme, "h" => k1h, "systemsVersion" => CODE_VERSION_K["k_systems"],
        "S" => (k1S.r, k1S.θ, k1S.φ), "R" => (k1R.r, k1R.θ, k1R.φ), "tlen" => k1tlen, "fmax" => k1fmax,
        "omegaiTlen" => k1Oi, "stride" => k1Stride, "θstep" => k1θstep, "rQ" => k1rQ)
    devRedo!("disks", params, KER_CACHE_DIR, "k1_disks")
    disks = myProduceOrLoad(params, KER_CACHE_DIR, "k1_disks") do cfg
        p = sourceOrderOf(scheme)
        θg = 0:0.1:180
        dS = computeSGTDisk(psv, sh, prem; family=:moment, rs_km=k1S.r, radii=[k1R.r, k1rQ], θdeg=θg,
                            tlen=k1tlen, fmax=k1fmax, omegaiTlen=k1Oi, p=p)
        dR = computeSGTDisk(psv, sh, prem; family=:force, rs_km=k1R.r, radii=[k1S.r, k1rQ], θdeg=θg,
                            tlen=k1tlen, fmax=k1fmax, omegaiTlen=k1Oi, p=p)
        # (d): one frequency (the highest), grid-stored disk and exact off-grid fields
        ωd = [2π * floor(Int, k1fmax * k1tlen) / k1tlen]
        dG = computeSGTDisk(psv, sh, prem; family=:moment, rs_km=k1S.r, radii=diskRadii(psv; stride=k1Stride),
                            θdeg=0:k1θstep:180, tlen=k1tlen, fmax=k1fmax, ωlist=ωd, omegaiTlen=k1Oi, p=p)
        rx = [1000.0 + 37.3, 2711.1, 3480.0 + 12.7, 4123.4, 4777.7, 5300.3, 6000.1, 6300.0]
        dX = computeSGTDisk(psv, sh, prem; family=:moment, rs_km=k1S.r, radii=rx, θdeg=collect(2.13:3.31:178),
                            tlen=k1tlen, fmax=k1fmax, ωlist=ωd, omegaiTlen=k1Oi, p=p)
        Dict("dS" => dS, "dR" => dR, "dG" => dG, "dX" => dX)
    end
    dS, dR = disks["dS"], disks["dR"]
    println("\n---- (b) reciprocity and vs pointSourceSynthetics")
    checkReciprocity(dS, dR, k1S, k1M, k1R; psv=psv, sh=sh, p=sourceOrderOf(scheme), omegaiTlen=k1Oi)
    println("\n---- (c) ZC route (10 + 20 functions) vs rotation route")
    checkZCRoute(dS, dR, length(dS.ω))
    println("\n---- (d) on-the-fly interpolation (stride $(k1Stride), θ step $(k1θstep)°), > 500 km from the source")
    for method in (:lagrange, :cubic, :linear)
        checkInterpolation(disks["dG"], disks["dX"]; iω=1, method=method, dmin=500.0)
    end
    @printf("truncation tail (last 10 degrees / max) at r = R, highest frequency: %.1e\n", dS.tail[end, 1])
end
