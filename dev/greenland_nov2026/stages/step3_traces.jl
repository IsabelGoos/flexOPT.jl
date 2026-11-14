# dev/greenland_nov2026/stages/step3_traces.jl
# Station traces (Z, N, E) in the requested band, record section, surface
# snapshots, P-C diagnostics, and one CSV per station. Cheap: rerun freely.

fsTrace = 1 / traceSampling
bandShort, bandLong = filterBand
bandShortEffective = max(bandShort, minResolvedPeriod)
bandShortEffective > bandShort && @warn(
    "the $(greenlandCase) grid resolves periods ≳ $(minResolvedPeriod) s; traces are filtered from $(bandShortEffective) s instead of $(bandShort) s")
function bandpassTrace(trace, shortPeriod, longPeriod)
    tr = trace .- mean(trace)
    taper = ones(length(tr)); nt = max(2, length(tr) ÷ 20)
    taper[1:nt] .= 0.5 .* (1 .- cos.(π .* (0:nt-1) ./ nt)); taper[end-nt+1:end] .= reverse(taper[1:nt])
    flt = digitalfilter(Bandpass(1 / longPeriod, 1 / shortPeriod; fs=fsTrace), Butterworth(2))
    return filtfilt(flt, tr .* taper)
end
filtered = similar(tracesZNE)
for is in axes(tracesZNE, 1), c in 1:3
    filtered[is, c, :] .= bandpassTrace(tracesZNE[is, c, :], bandShortEffective, bandLong)
end

distances = [greatCircleKm(s, SOURCE) for s in stations]
order = sortperm(distances)
componentNames = ("Z", "N", "E")

# ---- record section, one panel per component, each trace normalised
figSection = Figure(size=(1500, 900))
for c in 1:3
    ax = Axis(figSection[1, c]; title="$(componentNames[c]), $(bandShortEffective)–$(bandLong) s, trace-normalised",
        xlabel="time (s)", ylabel="distance from the source (km)")
    for is in order
        tr = filtered[is, c, :]
        a = maximum(abs, tr)
        lines!(ax, traceTimes, distances[is] .+ (a > 0 ? 70 .* tr ./ a : tr); color=:black, linewidth=0.8)
        text!(ax, last(traceTimes), distances[is]; text=stations[is].code, fontsize=9, align=(:left, :center))
    end
    for v in (3.0, 4.0)       # reference group velocities
        lines!(ax, SOURCE.onset .+ SOURCE.duration / 2 .+ [0, 3000] ./ v, [0, 3000];
            color=(:red, 0.5), linestyle=:dash)
    end
    xlims!(ax, 0, last(traceTimes) * 1.08)
end
Label(figSection[0, :], "Greenland OPT3 3-D ($(greenlandCase), dx = $(dx/1e3) km): single force near Nuuk " *
    "(peak $(SOURCE.peakForce) N, $(SOURCE.duration) s, slide az. $(SOURCE.slideAzimuth)°); dashed: 3 and 4 km/s",
    fontsize=15)
devSave(figSection, "record_section_$(greenlandCase)")

# ---- absolute amplitudes (µm), Z and the transverse picture
figTraces = Figure(size=(1300, 140 * length(stations) + 100))
for (row, is) in enumerate(order)
    ax = Axis(figTraces[row, 1]; ylabel=stations[is].code, yticklabelsize=8)
    for c in 1:3
        lines!(ax, traceTimes, 1e6 .* filtered[is, c, :]; label=componentNames[c], linewidth=0.8)
    end
    text!(ax, 0.01, 0.95; text=@sprintf("Δ = %.0f km", distances[is]), space=:relative,
        align=(:left, :top), fontsize=9)
    row < length(stations) && hidexdecorations!(ax; grid=false)
    row == 1 && axislegend(ax; position=:rt, orientation=:horizontal, framevisible=false, labelsize=9)
end
Label(figTraces[0, 1], "displacement (µm), $(bandShortEffective)–$(bandLong) s", fontsize=15)
devSave(figTraces, "traces_ZNE_$(greenlandCase)")

# ---- surface snapshots (box-frame vertical displacement on the top material cell)
snaps = optResult["snapshots"]
if !isempty(snaps)
    picks = unique(round.(Int, range(1, length(snaps); length=min(8, length(snaps)))))
    figSnap = Figure(size=(1600, 820))
    vmax = maximum(maximum(abs, s[3, :, :]) for s in snaps[picks])
    for (n, ip) in enumerate(picks)
        ax = Axis(figSnap[(n - 1) ÷ 4 + 1, (n - 1) % 4 + 1]; aspect=DataAspect(),
            title=@sprintf("u_z, t = %.0f s", optResult["snapshotTimes"][ip]))
        heatmap!(ax, yKm, northKm, asMap(snaps[ip][3, :, :]); colormap=:balance, colorrange=(-vmax, vmax) .* 0.3)
        scatter!(ax, [s.y / 1e3 for s in stations], [-s.x / 1e3 for s in stations];
            marker=:dtriangle, markersize=7, color=:black)
        hidedecorations!(ax)
    end
    devSave(figSnap, "snapshots_$(greenlandCase)")
end

# ---- 3-D view: u_z draped on the true curved surface (vertical exaggeration of the relief)
if !isempty(snaps)
    zSurf = [only(GeoPoints_to_local([GeoPoint(box.allGridsInGeoPoints[i, j, end].lat,
                box.allGridsInGeoPoints[i, j, end].lon; alt=max(modelResult["surfaceElevation"][i, j], 0.0))], box))[3] / 1e3
             for i in 1:Nx, j in 1:Ny]
    E3 = [yKm[j] for i in 1:Nx, j in 1:Ny]; N3 = [-xKm[i] for i in 1:Nx, j in 1:Ny]
    picks3 = unique(round.(Int, range(1, length(snaps); length=min(4, length(snaps)))))
    fig3d = Figure(size=(1600, 1250))
    for (n, ip) in enumerate(picks3)
        ax = Axis3(fig3d[(n - 1) ÷ 2 + 1, (n - 1) % 2 + 1]; aspect=(1, 1, 0.3),
            azimuth=-0.55π, elevation=0.15π, xlabel="≈ east (km)", ylabel="≈ north (km)", zlabel="z (km)",
            title=@sprintf("u_z on the surface, t = %.0f s", optResult["snapshotTimes"][ip]))
        uz = snaps[ip][3, :, :]
        a = max(maximum(abs, uz), eps(Float32))
        surface!(ax, E3, N3, zSurf; color=uz, colormap=:balance, colorrange=(-a, a) .* 0.5)
        scatter!(ax, [s.y / 1e3 for s in stations], [-s.x / 1e3 for s in stations],
            [zKm[s.index[3]] + 15 for s in stations]; marker=:dtriangle, markersize=10, color=:black)
    end
    devSave(fig3d, "snapshots3d_$(greenlandCase)")
end

# ---- predictor-corrector diagnostics
figPC = Figure(size=(1000, 380))
axPC = Axis(figPC[1, 1]; yscale=log10, xlabel="step", ylabel="accepted relative residual")
lines!(axPC, 1:length(pc.residuals), max.(pc.residuals, 1e-16))
axIt = Axis(figPC[1, 2]; xlabel="step", ylabel="corrector iterations")
lines!(axIt, 1:length(pc.iterations), pc.iterations)
devSave(figPC, "predictor_corrector_$(greenlandCase)")

# ---- CSV per station: time, Z, N, E (raw, unfiltered, metres)
csvDir = joinpath(DEV_OUTPUT_DIR[], DEV_STEP[], "traces")
mkpath(csvDir)
for (is, s) in enumerate(stations)
    open(joinpath(csvDir, "$(s.code)_$(greenlandCase).csv"), "w") do io
        println(io, "# $(s.code) lat=$(s.lat) lon=$(s.lon) distance_km=$(round(distances[is]; digits=1)) displacement_m unfiltered")
        println(io, "time_s,Z,N,E")
        writedlm(io, hcat(traceTimes, tracesZNE[is, 1, :], tracesZNE[is, 2, :], tracesZNE[is, 3, :]), ',')
    end
end
println("traces -> ", csvDir)
