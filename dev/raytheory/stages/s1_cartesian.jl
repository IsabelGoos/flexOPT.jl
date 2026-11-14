# dev/raytheory/stages/s1_cartesian.jl
# Isotropic Cartesian 2-D/3-D: analytic checks + ray figures.
#   gradient medium (circular rays), layer over half-space (PmP, PmS: P-S
#   conversion = an event, not a new equation), dipping reflector (2-D) and
#   an oblique dipping reflector in 3-D.
printChecks(checkCartesian(); title="Cartesian, isotropic")

M = cartesianModels()
fig = Figure(size=(1500, 950))

# (a) gradient: fan + wavefronts
ax = Axis(fig[1, 1]; title="v = $(M.v0) + $(M.g) z km/s: rays and wavefronts (T = 1, 2, … s)", xlabel="x (km)", ylabel="z (km)",
    yreversed=true, aspect=DataAspect())
fan = [traceRay(M.grad, parsePhase("P"), cartesianStart(M.grad, [0.0, 0.0], i, :P)...; dT=0.02, Tmax=40.0)
       for i in range(deg2rad(1), deg2rad(89); length=45)]
foreach(r -> plotRayPath!(ax, r; linewidth=0.8), fan)
for Tk in 1:9
    pts = [r.path[findfirst(t -> t >= Tk - 1e-9, r.pathT)] for r in fan if r.pathT[end] >= Tk]
    lines!(ax, [p[1] for p in pts], [p[2] for p in pts]; color=:black, linewidth=0.6)
end
xlims!(ax, -1, 60); ylims!(ax, -1, 30)

# (b) layer over half-space: P (direct), PmP, PmS to X = 30 km
ax = Axis(fig[1, 2]; title="layer over half-space: PmP (blue) and PmS (P→S at the Moho)", xlabel="x (km)", ylabel="z (km)",
    yreversed=true, aspect=DataAspect())
hlines!(ax, [M.H]; color=:gray40, linewidth=1)
for (ph, X) in (("PmP", 30.0), ("PmS", 30.0), ("PmP", 12.0), ("PmS", 12.0))
    s = shootCartesian2D(M.layer, ph, [0.0, 0.0], X; nfan=361, dT=0.02, Tmax=100.0)
    isempty(s) || plotRayPath!(ax, s[1].ray; linewidth=2)
    isempty(s) || @printf("  %-4s X = %4.1f km: T = %.6f s\n", ph, X, s[1].time)
end
xlims!(ax, -1, 32); ylims!(ax, -1, 12)

# (c) dipping reflector (2-D)
ax = Axis(fig[2, 1]; title="15° dipping reflector: Pv{refl}P to receivers every 4 km", xlabel="x (km)", ylabel="z (km)",
    yreversed=true, aspect=DataAspect())
Irefl = M.dipping.interfaces[2]
xx = [-20.0, 30.0]
lines!(ax, xx, [Irefl.x0[2] - Irefl.n[1] / Irefl.n[2] * (x - Irefl.x0[1]) for x in xx]; color=:gray40, linewidth=1.5)
for X in -16:4:24
    s = shootCartesian2D(M.dipping, "Pv{refl}P", [0.0, 0.0], Float64(X); nfan=361, dT=0.02, Tmax=60.0,
        bbox=[(-200.0, 200.0), (-1.0, 200.0)])
    isempty(s) || plotRayPath!(ax, s[1].ray; linewidth=1.2)
end
scatter!(ax, [0.0], [0.0]; marker=:star5, markersize=16, color=:red)

# (d) 3-D oblique dipping reflector: map view of the reflection points
ax = Axis(fig[2, 2]; title="3-D oblique dipping reflector: rays (map view), reflection points", xlabel="x (km)", ylabel="y (km)",
    aspect=DataAspect())
for (a, b) in ((10.0, 5.0), (-6.0, 12.0), (15.0, -10.0), (-12.0, -8.0), (0.0, 18.0), (20.0, 15.0))
    s = shootCartesian3D(M.dipping3, "Pv{refl}P", [0.0, 0, 0], [a, b, 0.0]; dT=0.02, Tmax=60.0,
        bbox=[(-200.0, 200.0), (-200.0, 200.0), (-1.0, 200.0)])
    s === nothing && continue
    lines!(ax, [p[1] for p in s.ray.path], [p[2] for p in s.ray.path]; color=RT_COLORS[1])
    h = s.ray.hits[1]
    scatter!(ax, [h.x[1]], [h.x[2]]; color=RT_COLORS[2], markersize=10)
    scatter!(ax, [a], [b]; marker=:utriangle, color=:black, markersize=12)
    @printf("  3-D receiver (%5.1f, %5.1f): T = %.6f s, reflection point (%.2f, %.2f, %.2f), misfit %.1e\n",
        a, b, s.time, h.x..., s.misfit)
end
scatter!(ax, [0.0], [0.0]; marker=:star5, markersize=16, color=:red)
devSave(fig, "cartesian_rays")
