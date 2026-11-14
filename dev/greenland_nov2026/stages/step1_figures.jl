# dev/greenland_nov2026/stages/step1_figures.jl
# Stop-and-review figure of the model. Map panels are drawn with the box
# y axis (~east) to the right and -x (~north) up, so they look like Fig. 1.

# xKm, yKm, northKm, asMap: defined in step1_sites.jl
topCell(field) = [field[i, j, topIndex[i, j]] for i in 1:Nx, j in 1:Ny]

function addSites!(ax)
    scatter!(ax, [s.y / 1e3 for s in stations], [-s.x / 1e3 for s in stations];
        marker=:dtriangle, markersize=13, color=:red, strokecolor=:black, strokewidth=1)
    text!(ax, [s.y / 1e3 for s in stations], [-s.x / 1e3 for s in stations];
        text=[s.code for s in stations], fontsize=9, offset=(4, 4))
    scatter!(ax, [sourceSite.y / 1e3], [-sourceSite.x / 1e3];
        marker=:star5, markersize=20, color=:yellow, strokecolor=:black, strokewidth=1)
    arrows!(ax, [sourceSite.y / 1e3], [-sourceSite.x / 1e3],
        [250 * slideDirection[2]], [-250 * slideDirection[1]]; color=:black, linewidth=2)
end

figModel = Figure(size=(1500, 1250))
axIce = Axis(figModel[1, 1]; title="ice thickness (ETOPO 2022 surface - bed), m",
    xlabel="local y ≈ east (km)", ylabel="local -x ≈ north (km)", aspect=DataAspect())
hIce = heatmap!(axIce, yKm, northKm, asMap(replace(iceThickness, 0.0 => NaN)); colormap=:turbo, colorrange=(0, 3500))
Colorbar(figModel[1, 2], hIce)
addSites!(axIce)

axWater = Axis(figModel[1, 3]; title="ocean depth, m", xlabel="local y ≈ east (km)",
    ylabel="local -x ≈ north (km)", aspect=DataAspect())
hWater = heatmap!(axWater, yKm, northKm, asMap(replace(waterDepth, 0.0 => NaN)); colormap=:deep, colorrange=(0, 4000))
Colorbar(figModel[1, 4], hWater)
addSites!(axWater)

axTop = Axis(figModel[2, 1]; title="Vs of the top material cell (mixed ice/water/rock), km/s",
    xlabel="local y ≈ east (km)", ylabel="local -x ≈ north (km)", aspect=DataAspect())
hTop = heatmap!(axTop, yKm, northKm, asMap(topCell(VsModel)); colormap=:viridis)
Colorbar(figModel[2, 2], hTop)
addSites!(axTop)

jSource = sourceIndex[2]
section = copy(VsModel[:, jSource, :]); section[.!materialMask[:, jSource, :]] .= NaN
axSec = Axis(figModel[2, 3]; title="Vs section through the source (y = $(round(Int, yKm[jSource])) km), km/s",
    xlabel="local x ≈ south (km)", ylabel="local z (km)")
hSec = heatmap!(axSec, xKm, zKm, section; colormap=:viridis)
scatter!(axSec, [sourceSite.x / 1e3], [zKm[sourceIndex[3]]]; marker=:star5, markersize=18, color=:yellow,
    strokecolor=:black, strokewidth=1)
Colorbar(figModel[2, 4], hSec)
Label(figModel[0, :], "Greenland box ($(greenlandCase), dx = $(dx/1e3) km): $(Nx)×$(Ny)×$(Nz) points, " *
    "PREM + ETOPO 2022 ice and ocean, arrow = slide direction", fontsize=16)
devSave(figModel, "model_$(greenlandCase)")
