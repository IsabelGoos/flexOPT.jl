# dev/greenland_nov2026/stages/step1_sites.jl
# Stations and source on the grid: nearest column, top material node of that
# column (the receiver/force sits right below the free surface), local ZNE
# bases, and the force direction in the box frame.

Rbox = box.rotationMatrix
function siteOnGrid(lat, lon; depth=0.0)
    g = GeoPoint(lat, lon)
    xyz = only(GeoPoints_to_local([g], box))
    i = argmin(abs.(xBox .- xyz[1]))
    j = argmin(abs.(yBox .- xyz[2]))
    k = topIndex[i, j]
    # sensor depth (BORG: 95 m) is far below one cell: go down only if > dx/2
    k = max(1, k - round(Int, depth / dx))
    edge = minimum((xyz[1] - first(xBox), last(xBox) - xyz[1], xyz[2] - first(yBox), last(yBox) - xyz[2]))
    return (geo=g, x=xyz[1], y=xyz[2], index=CartesianIndex(i, j, k), edgeDistance=edge)
end

stations = map(GLISN_STATIONS) do s
    site = siteOnGrid(s.lat, s.lon; depth=s.depth)
    z, n, e = zne_basis_local(site.geo, Rbox)
    merge(s, site, (zLocal=z, nLocal=n, eLocal=e,
        iceInCell=modelResult["fIce"][site.index], waterInCell=modelResult["fWater"][site.index]))
end
for s in stations
    @assert s.edgeDistance >= MIN_STATION_EDGE_DISTANCE "station $(s.code) is only $(round(s.edgeDistance/1e3)) km from the box edge"
end

sourceSite = siteOnGrid(SOURCE.lat, SOURCE.lon)
sourceIndex = sourceSite.index
@assert materialMask[sourceIndex] && VsModel[sourceIndex] > 0 "the force must sit in solid material"
# unit force ON THE EARTH while the mass accelerates (-slide direction); the
# time function landslide_force_shape changes its sign during deceleration
slideDirection = local_direction(sourceSite.geo, Rbox; azimuth=SOURCE.slideAzimuth, plunge=SOURCE.slidePlunge)
forceDirection = -slideDirection
greatCircleKm(a, b) = 6371.0 * acos(clamp(
    sind(a.lat) * sind(b.lat) + cosd(a.lat) * cosd(b.lat) * cosd(a.lon - b.lon), -1, 1))
for s in stations
    s_dist = greatCircleKm(s, SOURCE)
    @printf("%-6s x=%6.0f y=%6.0f km  node %-16s ice %.2f water %.2f  Δ=%5.0f km\n",
        s.code, s.x / 1e3, s.y / 1e3, string(Tuple(s.index)), s.iceInCell, s.waterInCell, s_dist)
end
@printf("source: node %s (x=%.0f, y=%.0f km), Vs=%.2f km/s, force dir (box) = %s\n",
    string(Tuple(sourceIndex)), sourceSite.x / 1e3, sourceSite.y / 1e3, VsModel[sourceIndex],
    string(round.(forceDirection; digits=3)))

# map coordinates for figures (y ≈ east to the right, -x ≈ north up)
xKm, yKm, zKm = xBox ./ 1e3, yBox ./ 1e3, zBox ./ 1e3
northKm = reverse(-xKm)
asMap(field) = permutedims(field[end:-1:1, :])        # (y, north) for heatmap!
