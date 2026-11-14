# dev/raytheory/stages/s2_spherical.jl
# 1-D spherical planets (cases: PREM [default], ak135, iasp, DWThot, TAYAKr, vprem):
#   * tracer (Hamilton's equations in (r, θ)) vs τ-p integrals at the same p
#   * PREM: τ-p arrivals vs TauP (ObsPy) on the same model (reference/taup_premDSM.csv)
#   * figures: travel-time curves (+ TauP points for PREM), ray paths in the disk
RT_S2_PHASES = ["P", "PP", "PcP", "PKP", "PKIKP", "PKiKP", "S", "SS", "ScS", "SKS"]
RT_S2_RAYS = [("P", 30.0), ("P", 60.0), ("P", 90.0), ("PcP", 40.0), ("PP", 120.0), ("PKP", 150.0),
    ("PKIKP", 175.0), ("PKiKP", 120.0), ("SKS", 110.0), ("ScS", 60.0), ("P410s", 50.0)]

for case in (isempty(DEV.cases) ? ["PREM"] : DEV.cases)
    m = rtModel(case)
    show(stdout, m)
    printChecks(checkSphericalTracer(m); title="$case: tracer vs τ-p (same ray parameter)")
    case == "PREM" && printChecks(checkTauP(m, RT_TAUP_CSV); title="PREM: τ-p vs TauP (same model, reference/taup_premDSM.csv)")

    # arrivals table at a few distances, surface source
    println("\n$case: arrivals for a surface source")
    for d in (30.0, 60.0, 90.0, 120.0, 150.0)
        for ph in RT_S2_PHASES, a in taupArrivals(m, ph, 0.0, d)
            println("  ", a)
        end
    end

    fig = Figure(size=(1500, 750))
    ax = Axis(fig[1, 1]; title="$case: travel times, surface source (lines: τ-p, dots: TauP)", xlabel="Δ (°)", ylabel="T (min)")
    for (i, ph) in enumerate(RT_S2_PHASES)
        first = true
        for (Δs, Ts, _) in travelTimeCurve(m, ph, 0.0)
            lines!(ax, Δs, Ts ./ 60; color=rtColor(i), linewidth=1.5, label=first ? ph : nothing)
            first = false
        end
    end
    if case == "PREM"
        rows = [split(l, ',') for l in readlines(RT_TAUP_CSV)[2:end]]
        sel = [r for r in rows if r[1] in RT_S2_PHASES && parse(Float64, r[3]) == 0.0]
        scatter!(ax, [parse(Float64, r[4]) for r in sel], [parse(Float64, r[5]) / 60 for r in sel];
            color=:black, markersize=5)
    end
    xlims!(ax, 0, 180)
    axislegend(ax; position=:lt, framevisible=false)

    ax = Axis(fig[1, 2]; title="$case: rays from a surface source (blue P / K, orange S)", aspect=DataAspect())
    hidedecorations!(ax); hidespines!(ax)
    plotPlanetDisk!(ax, m)
    for (ph, d) in RT_S2_RAYS
        for a in taupArrivals(m, ph, 0.0, d)
            x0, p0 = sphericalStartP(m, 0.0, a.p, asPhase(ph).mode, a.startDir)
            r = traceRay(m, asPhase(ph), x0, p0; dT=0.5)
            r.ok && plotSphericalRay!(ax, r)
        end
    end
    devSave(fig, "spherical_$(case)")
end
