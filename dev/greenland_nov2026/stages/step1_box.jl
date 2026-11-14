# dev/greenland_nov2026/stages/step1_box.jl
# Local Cartesian box (GeoPoints.constructLocalBox) containing the true,
# curved Earth surface: at the box corners the surface lies ~230 km below the
# tangent plane at the centre, so the box goes deep enough to keep
# DEPTH_BELOW_SURFACE of Earth there, and just high enough to hold the
# summit of the ice sheet (~3.2 km) plus one air layer.

p0 = GeoPoint(BOX_CENTRE.lat, BOX_CENTRE.lon)
# Size the box from where the stations and the source really fall in the
# GeoPoints (WGS84) frame: a tiny symmetric box has exactly the same frame
# (origin p0, same axes) as the final symmetric one.
frameBox = constructLocalBox(p0, dx, dx, dx, -dx, dx, -dx, dx, -dx, dx;
    axis_angle_deg=BOX_AXIS_AZIMUTH)
siteXY = GeoPoints_to_local(vcat([GeoPoint(s.lat, s.lon) for s in GLISN_STATIONS],
    [GeoPoint(SOURCE.lat, SOURCE.lon)]), frameBox)
farthestSite = maximum(max(abs(v[1]), abs(v[2])) for v in siteXY)
boxHalfWidth = max(BOX_HALF_WIDTH_MIN,
    ceil((farthestSite + MIN_STATION_EDGE_DISTANCE + dx) / dx) * dx)
@printf("farthest site %.0f km from the centre -> box half-width %.0f km\n",
    farthestSite / 1e3, boxHalfWidth / 1e3)
cornerDistance = sqrt(2) * boxHalfWidth
earthRadius = planet1D.my1DDSMmodel.averagedPlanetRadiusInKilometer * 1e3
cornerDrop = earthRadius - sqrt(earthRadius^2 - cornerDistance^2) # surface below the tangent plane
nBelow = ceil(Int, (cornerDrop + DEPTH_BELOW_SURFACE) / dx)
nAbove = ceil(Int, (4.0e3 + dx) / dx)                              # ice summit + >= 1 air layer
altMin = -nBelow * dx
altMax = nAbove * dx
box = constructLocalBox(
    p0, dx, dx, dx,
    -boxHalfWidth, boxHalfWidth,
    -boxHalfWidth, boxHalfWidth,
    altMin, altMax;
    axis_angle_deg=BOX_AXIS_AZIMUTH,
)
Nx, Ny, Nz = size(box.allGridsInGeoPoints)
xBox = [p.xyz[1] for p in box.allGridsInCartesian[:, 1, 1]]
yBox = [p.xyz[2] for p in box.allGridsInCartesian[1, :, 1]]
zBox = [p.xyz[3] for p in box.allGridsInCartesian[1, 1, :]]
@assert isapprox(xBox[2] - xBox[1], dx) && isapprox(zBox[2] - zBox[1], dx)
nWhole = (Nx + 2cerjanCells) * (Ny + 2cerjanCells) * (Nz + cerjanCells)
@printf("box: %d x %d x %d physical points (%.2f M), %.2f M with Cerjan layers\n",
    Nx, Ny, Nz, Nx * Ny * Nz / 1e6, nWhole / 1e6)
@printf("     x, y in [%.0f, %.0f] km, z in [%.0f, %.0f] km; surface drop at corners %.0f km\n",
    first(xBox) / 1e3, last(xBox) / 1e3, first(zBox) / 1e3, last(zBox) / 1e3, cornerDrop / 1e3)
@printf("     OPT operators will need roughly %.1f GB of memory\n",
    estimate_operator_memory_gb(nWhole))
