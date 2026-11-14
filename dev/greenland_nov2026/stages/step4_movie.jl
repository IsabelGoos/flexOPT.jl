# dev/greenland_nov2026/stages/step4_movie.jl
# Propagation movie (MP4, CairoMakie `record`) from the surface snapshots of
# step2: left, u_z on the top material cell (map view, coast/ice outline,
# stations, source); right, the vertical traces of every station with a moving
# time cursor. Frame rate follows the snapshot interval (config: snapshotEvery).

snaps = optResult["snapshots"]
snapTimes = optResult["snapshotTimes"]
isempty(snaps) && error("no snapshots in the step2 result")
nFrames = length(snaps)
uzFrames = [asMap(s[3, :, :]) for s in snaps]
# colour scale: robust (98th percentile over all frames), symmetric
allAbs = sort!(abs.(reduce(vcat, vec.(uzFrames))))
clim = max(allAbs[clamp(round(Int, 0.98 * length(allAbs)), 1, length(allAbs))], eps(Float32))

frame = Observable(1)
figMovie = Figure(size=(1600, 820))
axMap = Axis(figMovie[1, 1]; aspect=DataAspect(), xlabel="local y ≈ east (km)", ylabel="local -x ≈ north (km)",
    title=@lift(@sprintf("u_z at the surface, t = %4.0f s", snapTimes[$frame])))
hMap = heatmap!(axMap, yKm, northKm, @lift(uzFrames[$frame]); colormap=:balance, colorrange=(-clim, clim))
landMask = asMap(Float64.((iceThickness .> 0) .| (waterDepth .<= 0)))
contour!(axMap, yKm, northKm, landMask; levels=[0.5], color=:black, linewidth=1)
contour!(axMap, yKm, northKm, asMap(iceThickness); levels=[500.0], color=(:gray40, 0.8), linewidth=0.8,
    linestyle=:dash)
scatter!(axMap, [s.y / 1e3 for s in stations], [-s.x / 1e3 for s in stations];
    marker=:dtriangle, markersize=11, color=:white, strokecolor=:black, strokewidth=1)
text!(axMap, [s.y / 1e3 for s in stations], [-s.x / 1e3 for s in stations];
    text=[s.code for s in stations], fontsize=9, offset=(4, 3))
scatter!(axMap, [sourceSite.y / 1e3], [-sourceSite.x / 1e3]; marker=:star5, markersize=18, color=:yellow,
    strokecolor=:black, strokewidth=1)
Colorbar(figMovie[1, 2], hMap; label="u_z (m), clipped at the 98th percentile")

distancesMovie = [greatCircleKm(s, SOURCE) for s in stations]
orderMovie = sortperm(distancesMovie)
axTr = Axis(figMovie[1, 3]; xlabel="time (s)", ylabel="distance from the source (km)",
    title="Z traces (unfiltered, trace-normalised)")
for is in orderMovie
    tr = tracesZNE[is, 1, :]
    a = max(maximum(abs, tr), eps())
    lines!(axTr, traceTimes, distancesMovie[is] .+ 70 .* tr ./ a; color=:black, linewidth=0.7)
    text!(axTr, last(traceTimes), distancesMovie[is]; text=stations[is].code, fontsize=9, align=(:left, :center))
end
vlines!(axTr, @lift([snapTimes[$frame]]); color=:red, linewidth=2)
xlims!(axTr, 0, last(traceTimes) * 1.1)
colsize!(figMovie.layout, 1, Relative(0.55))
Label(figMovie[0, :], "Greenland OPT3 3-D ($(greenlandCase), dx = $(dx/1e3) km): ice slide near Nuuk, " *
    "single force, peak $(SOURCE.peakForce) N, $(SOURCE.duration) s", fontsize=15)

movieDir = joinpath(DEV_OUTPUT_DIR[], DEV_STEP[])
mkpath(movieDir)
moviePath = joinpath(movieDir, "propagation_$(greenlandCase).mp4")
framerate = clamp(round(Int, 50 / snapshotEvery), 2, 24)     # ~50 s of simulation per second of video
record(figMovie, moviePath, 1:nFrames; framerate=framerate) do n
    frame[] = n
end
println("movie -> ", moviePath, "  (", nFrames, " frames, ", framerate, " fps)")
