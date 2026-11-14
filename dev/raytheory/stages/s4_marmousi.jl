# dev/raytheory/stages/s4_marmousi.jl
# Acoustic rays in Marmousi (vp only), from notebooks/tmp/seismicModelMarmousi.jld2
# (the model the Marmousi notebooks build from dataInput/model/random/marmousi.png
# with src/imageReader; 802 × 676 samples spanning 9.2 km × 3.0 km, last
# column = surface).
#   * the grid is smoothed (Gaussian σ, cases s100 [default], s50, s200 = σ in m)
#     and interpolated by a C² cubic B-spline (lib/rayGrid.jl);
#   * checks: gridded gradient vs closed form, eikonal convergence, and the
#     earliest Marmousi ray vs the eikonal first arrival;
#   * figures: ray fans + eikonal wavefronts (surface and buried source), and
#     T(x) at the surface: every ray branch (multipathing, triplications) vs the
#     eikonal first arrival (which also contains non-ray energy: the direct
#     wave grazing the free surface, waves along the model bottom).
RT_MARMOUSI_SIGMA = Dict("s50" => 0.05, "s100" => 0.1, "s200" => 0.2)
isfile(RT_MARMOUSI_FILE) || error("$RT_MARMOUSI_FILE not found: run notebooks/marmousi_elastic_OPT_debug.ipynb (or src/imageReader on dataInput/model/random/marmousi.png) once")
Vm, dxm, dzm = loadMarmousiVp()
@printf("Marmousi vp: %d × %d, dx = %.2f m, dz = %.2f m, vp ∈ [%.2f, %.2f] km/s\n", size(Vm)..., 1e3dxm, 1e3dzm, extrema(Vm)...)

"fan of rays, keeping the paths; returns RayResult vector"
marmousiFan(med, bbox, xs, takeoffs; dT=0.005, Tmax=6.0) =
    [traceRay(med, Phase("P", :P, :any, PhaseEvent[]), cartesianStart(med, xs, i, :P)...; dT=dT, Tmax=Tmax, bbox=bbox) for i in takeoffs]

"positions of all rays at time t (NaN-separated where neighbours are far apart: wavefront pieces)"
function rayFront(rays, t; gap=0.15)
    xs = Float64[]; zs = Float64[]
    last = nothing
    for r in rays
        i = findfirst(>=(t), r.pathT)
        if i === nothing
            last = nothing; push!(xs, NaN); push!(zs, NaN); continue
        end
        p = r.path[i]
        if last !== nothing && hypot(p[1] - last[1], p[2] - last[2]) > gap
            push!(xs, NaN); push!(zs, NaN)
        end
        push!(xs, p[1]); push!(zs, p[2]); last = p
    end
    return xs, zs
end

for case in (isempty(DEV.cases) ? ["s100"] : DEV.cases)
    σ = RT_MARMOUSI_SIGMA[case]
    gp = smoothGrid(Vm, 0.0, 0.0, dxm, dzm; σ=σ)
    med, bbox = gridMedium2D(gp)
    if case == "s100"
        printChecks(checkGrid(gp); title="Marmousi, σ = $(1e3σ) m: grid checks")
    else
        # other σ: the earliest ray is not always the first arrival (rougher model:
        # shadow zones, non-ray energy), so print the comparison, no pass/fail
        Eg = eikonalOnGrid(gp, 3.0, 0.0; h=0.005); Ef = eikonalOnGrid(gp, 3.0, 0.0; h=0.0025)
        println("\nMarmousi σ = $(round(Int, 1e3σ)) m: earliest ray vs eikonal first arrival (Richardson 5/2.5 m)")
        for X in 6.5:0.5:9.0
            s = shootCartesian2D(med, "P", [3.0, 0.0], X; nfan=2881, dT=0.005, Tmax=10.0, bbox=bbox)
            e = 2 * Ef[1][argmin(abs.(Ef[2] .- X)), 1] - Eg[1][argmin(abs.(Eg[2] .- X)), 1]
            @printf("  offset %.1f km: %2d rays, earliest %s, eikonal %.4f s\n", X - 3.0, length(s),
                isempty(s) ? "  —   " : @sprintf("%.4f", s[1].time), e)
        end
    end

    xsS = [3.0, 0.0]; xsB = [6.0, 2.2]
    raysS = marmousiFan(med, bbox, xsS, range(deg2rad(-89), deg2rad(89); length=900))
    raysB = marmousiFan(med, bbox, xsB, range(-π, π; length=1441)[1:end-1])
    TS, xgS, zgS = eikonalOnGrid(gp, xsS[1], xsS[2]; h=0.01)
    TB, xgB, zgB = eikonalOnGrid(gp, xsB[1], xsB[2]; h=0.01)

    fig = Figure(size=(1700, 1250))
    xr = extentX(gp); zr = extentZ(gp)
    ax = Axis(fig[1, 1]; title="Marmousi vp (raw grid)", xlabel="x (km)", ylabel="z (km)", yreversed=true, aspect=DataAspect())
    hm = heatmap!(ax, gridX(gp), gridZ(gp), Vm; colormap=:viridis)
    Colorbar(fig[1, 2], hm; label="vp (km/s)")
    ax = Axis(fig[1, 3]; title="smoothed σ = $(round(Int, 1e3σ)) m (what the rays see)", xlabel="x (km)", yreversed=true, aspect=DataAspect())
    heatmap!(ax, gridX(gp), gridZ(gp), gp.values; colormap=:viridis)

    for (row, rays, xs, T, xg, zg, lab) in ((2, raysS, xsS, TS, xgS, zgS, "surface source"),
                                           (3, raysB, xsB, TB, xgB, zgB, "buried source"))
        ax = Axis(fig[row, 1]; title="$lab: rays (every 4th) and eikonal first-arrival wavefronts (0.2 s)",
            xlabel="x (km)", ylabel="z (km)", yreversed=true, aspect=DataAspect())
        heatmap!(ax, gridX(gp), gridZ(gp), gp.values; colormap=(:grays, 0.35))
        for r in rays[1:4:end]
            lines!(ax, [p[1] for p in r.path], [p[2] for p in r.path]; color=(RT_COLORS[1], 0.5), linewidth=0.5)
        end
        contour!(ax, xg, zg, T; levels=0.2:0.2:6, color=:black, linewidth=0.8)
        scatter!(ax, [xs[1]], [xs[2]]; marker=:star5, markersize=18, color=:red)
        xlims!(ax, xr...); ylims!(ax, zr...)
        ax = Axis(fig[row, 3]; title="$lab: ray wavefronts (dots, 0.2 s) — folds = caustics/triplications",
            xlabel="x (km)", yreversed=true, aspect=DataAspect())
        heatmap!(ax, gridX(gp), gridZ(gp), gp.values; colormap=(:grays, 0.35))
        for t in 0.2:0.2:6.0
            fx, fz = rayFront(rays, t)
            lines!(ax, fx, fz; color=RT_COLORS[2], linewidth=1)
        end
        contour!(ax, xg, zg, T; levels=0.2:0.2:6, color=(:black, 0.5), linewidth=0.5, linestyle=:dash)
        scatter!(ax, [xs[1]], [xs[2]]; marker=:star5, markersize=18, color=:red)
        xlims!(ax, xr...); ylims!(ax, zr...)
    end
    # T(x) at the surface
    ax = Axis(fig[4, 1:3]; title="surface source: travel times at the surface — all ray branches (dots) vs eikonal first arrival (line)",
        xlabel="x (km)", ylabel="T (s)")
    arr = [(r.xEnd[1], r.T) for r in raysS if r.ok]
    scatter!(ax, first.(arr), last.(arr); color=RT_COLORS[1], markersize=4, label="rays (arrive at z = 0)")
    lines!(ax, xgS, TS[:, 1]; color=:black, linewidth=1.5, label="eikonal (first arrival, h = 10 m)")
    axislegend(ax; position=:lt, framevisible=false)
    xlims!(ax, xr...); ylims!(ax, 0, 4.5)
    rowsize!(fig.layout, 4, Relative(0.18))
    devSave(fig, "marmousi_$(case)")
end
