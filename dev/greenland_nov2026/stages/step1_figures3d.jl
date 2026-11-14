# dev/greenland_nov2026/stages/step1_figures3d.jl
# 3-D views of the model (CairoMakie Axis3), in the local box frame (km):
#   model3d_icesheet  -- ETOPO 2022 ice surface over the bedrock / sea floor,
#                        vertical exaggeration, stations and the source;
#   model3d_curved    -- the true geometry of the box: curved Earth surface
#                        (ice/ocean coloured), stations, source, slide arrow;
#   model3d_fence     -- Vs on two vertical slices (through the source and
#                        through SUMG) under the semi-transparent surface.

# local (x, y) on the tangent plane -> GeoPoint below it (lat/lon of that column)
localToGeo(x, y) = GeoPoint(p_local_to_ECEF(x, y, 0.0, box.pOriginECEF, Rbox))
localZ(lat, lon, alt) = only(GeoPoints_to_local([GeoPoint(lat, lon; alt=alt)], box))[3]

nFine = 220
xFine = range(first(xBox), last(xBox); length=nFine)
yFine = range(first(yBox), last(yBox); length=nFine)
geoFine = [localToGeo(x, y) for x in xFine, y in yFine]
surfFine = [surfItp(g.lon, g.lat) for g in geoFine]
bedFine = [bedItp(g.lon, g.lat) for g in geoFine]
iceFine = max.(surfFine .- bedFine, 0.0) .* (surfFine .> 0)
iceMask = iceFine .> ICE_THRESHOLD
# map orientation as in the 2-D figures: east = local y, north = -local x
E = [y / 1e3 for x in xFine, y in yFine]
N = [-x / 1e3 for x in xFine, y in yFine]
stationE = [s.y / 1e3 for s in stations]; stationN = [-s.x / 1e3 for s in stations]

# ---------------------------------------------------------------- 1. ice sheet
exaggeration = 40.0
figIce3d = Figure(size=(1500, 1100))
axIce3d = Axis3(figIce3d[1, 1]; aspect=(1, 1, 0.3), azimuth=-0.55π, elevation=0.25π,
    xlabel="≈ east (km)", ylabel="≈ north (km)", zlabel="elevation × $(Int(exaggeration)) (km)",
    title="ETOPO 2022: ice surface (colour = thickness) over bedrock / sea floor")
hBed = surface!(axIce3d, E, N, exaggeration .* bedFine ./ 1e3; color=bedFine, colormap=:oleron,
    colorrange=(-4000, 4000))
iceTop = [m ? exaggeration * s / 1e3 : NaN for (m, s) in zip(iceMask, surfFine)]
hIce3d = surface!(axIce3d, E, N, iceTop; color=[m ? t : NaN for (m, t) in zip(iceMask, iceFine)],
    colormap=:ice, colorrange=(0, 3500), nan_color=:transparent)
stationTop = [exaggeration * max(surfItp(s.lon, s.lat), 0.0) / 1e3 + 15 for s in stations]
scatter!(axIce3d, stationE, stationN, stationTop; marker=:dtriangle, markersize=14, color=:red,
    strokecolor=:black, strokewidth=1)
text!(axIce3d, stationE, stationN, stationTop .+ 25; text=[s.code for s in stations], fontsize=11)
scatter!(axIce3d, [sourceSite.y / 1e3], [-sourceSite.x / 1e3], [40.0]; marker=:star5, markersize=24,
    color=:yellow, strokecolor=:black, strokewidth=1)
Colorbar(figIce3d[1, 2], hIce3d; label="ice thickness (m)")
Colorbar(figIce3d[1, 3], hBed; label="bedrock / sea-floor elevation (m)")
devSave(figIce3d, "model3d_icesheet_$(greenlandCase)")

# ---------------------------------------------------------------- 2. true geometry
# Earth's surface (ice/ocean/land top) in the box frame: curvature drops it
# ~230 km below the tangent plane at the corners.
zTrue = [localZ(g.lat, g.lon, max(s, 0.0)) / 1e3 for (g, s) in zip(geoFine, surfFine)]
surfaceKind = [m ? 2.0 : (s <= 0 ? 0.0 : 1.0) for (m, s) in zip(iceMask, surfFine)]   # 0 ocean, 1 land, 2 ice
figCurved = Figure(size=(1500, 1000))
axCurved = Axis3(figCurved[1, 1]; aspect=(1, 1, 0.3), azimuth=-0.55π, elevation=0.12π,
    xlabel="≈ east (km)", ylabel="≈ north (km)", zlabel="local z (km)",
    title="Box geometry ($(greenlandCase)): true curved surface, ocean / land / ice, stations and source")
surface!(axCurved, E, N, zTrue; color=surfaceKind, colormap=[:steelblue, :tan, :white],
    colorrange=(0, 2))
stationZ = [zKm[s.index[3]] + dx / 2e3 for s in stations]
scatter!(axCurved, stationE, stationN, stationZ .+ 5; marker=:dtriangle, markersize=13, color=:red,
    strokecolor=:black, strokewidth=1)
text!(axCurved, stationE, stationN, stationZ .+ 20; text=[s.code for s in stations], fontsize=10)
sourceZ = zKm[sourceIndex[3]] + dx / 2e3
scatter!(axCurved, [sourceSite.y / 1e3], [-sourceSite.x / 1e3], [sourceZ + 5]; marker=:star5,
    markersize=24, color=:yellow, strokecolor=:black, strokewidth=1)
arrows!(axCurved, Point3f[(sourceSite.y / 1e3, -sourceSite.x / 1e3, sourceZ + 10)],
    Vec3f[(300 * slideDirection[2], -300 * slideDirection[1], 300 * slideDirection[3])];
    color=:black, linewidth=4, arrowsize=Vec3f(40, 40, 30))
# outline of the physical box
bx = [first(yKm), last(yKm)]; by = [-last(xKm), -first(xKm)]; bz = [first(zKm), last(zKm)]
for a in bx, b in by
    lines!(axCurved, [a, a], [b, b], bz; color=(:gray30, 0.6))
end
for zz in bz
    lines!(axCurved, [bx[1], bx[2], bx[2], bx[1], bx[1]], [by[1], by[1], by[2], by[2], by[1]],
        fill(zz, 5); color=(:gray30, 0.6))
end
devSave(figCurved, "model3d_curved_$(greenlandCase)")

# ---------------------------------------------------------------- 3. fence diagram
sumg = stations[findfirst(s -> s.code == "SUMG", stations)]
iSlice, jSlice = sumg.index[1], sourceIndex[2]
vsShow = copy(VsModel); vsShow[.!materialMask] .= NaN
figFence = Figure(size=(1500, 1000))
axFence = Axis3(figFence[1, 1]; aspect=(1, 1, 0.45), azimuth=-0.62π, elevation=0.20π,
    xlabel="≈ east (km)", ylabel="≈ north (km)", zlabel="local z (km)",
    title="Vs (km/s) on vertical slices through the source and through SUMG ($(greenlandCase), dx = $(dx/1e3) km)")
crange = (1.5, 6.0)
# slice 1: constant local y through the source (runs ≈ north-south)
X1 = [yKm[jSlice] for i in 1:Nx, k in 1:Nz]; Y1 = [-xKm[i] for i in 1:Nx, k in 1:Nz]
Z1 = [zKm[k] for i in 1:Nx, k in 1:Nz]
hFence = surface!(axFence, X1, Y1, Z1; color=vsShow[:, jSlice, :], colormap=:viridis, colorrange=crange,
    nan_color=:transparent, shading=NoShading)
# slice 2: constant local x through SUMG (runs ≈ east-west)
X2 = [yKm[j] for j in 1:Ny, k in 1:Nz]; Y2 = [-xKm[iSlice] for j in 1:Ny, k in 1:Nz]
Z2 = [zKm[k] for j in 1:Ny, k in 1:Nz]
surface!(axFence, X2, Y2, Z2; color=vsShow[iSlice, :, :], colormap=:viridis, colorrange=crange,
    nan_color=:transparent, shading=NoShading)
# the Earth surface as a light wireframe, so the slices stay visible
sub = 1:10:nFine
wireframe!(axFence, E[sub, sub], N[sub, sub], zTrue[sub, sub]; color=(:black, 0.35), linewidth=0.6)
scatter!(axFence, stationE, stationN, stationZ .+ 5; marker=:dtriangle, markersize=11, color=:red)
scatter!(axFence, [sourceSite.y / 1e3], [-sourceSite.x / 1e3], [sourceZ + 5]; marker=:star5,
    markersize=22, color=:yellow, strokecolor=:black, strokewidth=1)
Colorbar(figFence[1, 2], hFence; label="Vs (km/s)")
devSave(figFence, "model3d_fence_$(greenlandCase)")
