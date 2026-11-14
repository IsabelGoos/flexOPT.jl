# dev/kirishima_oct2026/stages/step1_box.jl
# Ported from notebooks/Elastic2DtopoBenchmark_oct2026.ipynb, cells step1-box-model + step1-topography-plots.
# The models are cached by name with lazyProduceOrLoad in KIRISHIMA_MODEL_CACHE
# (= notebooks/tmp, where the notebook already saved them: same names, reused).
# Kirishima summit -- same point and same box/model construction as
# SimuKirishima.ipynb's "making an initial model" section.

# Keep the 3D model moderate; the 2D slice is resolved independently
# (its OPT recipe is nondimensional, so it is not rebuilt when these
# physical spacings change).
Δx3D = 200.0 # m
Δy3D = 200.0 # m
Δz3D = 200.0 # m
Δx2D = 100.0 # m
Δz2D = 100.0 # m
altMax = 2.e3 # m
altMin = -40.e3 # m
horizontalHalfWidth = 40.e3 # m; the model spans 80 km horizontally

boxGrids3D = constructLocalBox(
    p0, Δx3D, Δy3D, Δz3D,
    -horizontalHalfWidth, horizontalHalfWidth,
    -horizontalHalfWidth, horizontalHalfWidth,
    altMin, altMax,
)
boxGrids2D = constructLocalBox(
    p0, Δx2D, Δz2D,
    -horizontalHalfWidth, horizontalHalfWidth,
    altMin, altMax,
)


gridCacheTag3D = "$(join(size(boxGrids3D.allGridsInGeoPoints), 'x'))_" *
    "d$(boxGrids3D.Δx)_$(boxGrids3D.Δy)_$(boxGrids3D.Δz)"
modelCacheName3D = "seismicModel3D_Kirishima_$(gridCacheTag3D)_NIED_$(velocitySourceKey)"
seismicModel3D = lazyProduceOrLoad(
    modelCacheName3D,
    getParamsAndTopo,
    boxGrids3D.allGridsInGeoPoints,
    boxGrids3D.effectiveRadii,
    0.1;
    folder=KIRISHIMA_MODEL_CACHE, velocity_model=:NIED,
    nied_source=niedVelocitySource,
    nied_confidence_max=0.8,
    nied_outside=:planet1D,
    nied_low_confidence=:planet1D,
)
@assert size(seismicModel3D.ρ) == size(boxGrids3D.allGridsInGeoPoints)

gridCacheTag2D = "$(join(size(boxGrids2D.allGridsInGeoPoints), 'x'))_" *
    "d$(boxGrids2D.Δx)_$(boxGrids2D.Δz)"
modelCacheName2D = "seismicModel2D_Kirishima_$(gridCacheTag2D)_NIED_$(velocitySourceKey)"
seismicModel2D = lazyProduceOrLoad(
    modelCacheName2D,
    getParamsAndTopo,
    boxGrids2D.allGridsInGeoPoints,
    boxGrids2D.effectiveRadii,
    0.1;
    folder=KIRISHIMA_MODEL_CACHE, velocity_model=:NIED,
    nied_source=niedVelocitySource,
    nied_confidence_max=0.8,
    nied_outside=:planet1D,
    nied_low_confidence=:planet1D,
)
@assert size(seismicModel2D.ρ) == size(boxGrids2D.allGridsInGeoPoints)

@show seismicModel3D.velocity_model seismicModel3D.nied_source
@show count(seismicModel3D.nied_mask) / length(seismicModel3D.nied_mask)
@show seismicModel2D.velocity_model seismicModel2D.nied_source
# getParamsAndTopo uses ρAir=0.001 by default; 0.01 cleanly separates
# rock from air/void (same threshold as SimuKirishima.ipynb). No
# free-surface case, material-filling or boundary-condition logic yet --
# that is Steps 2-4.
# air_density_cutoff: see shared_site.jl
material3D = (seismicModel3D.ρ .> air_density_cutoff) .& (seismicModel3D.Vsv .> 0)
material2D = (seismicModel2D.ρ .> air_density_cutoff) .& (seismicModel2D.Vsv .> 0)

# --- 2D topography (x, z in metres; what Steps 2-4 will build directly on) ---
xCoordinates2D = [p.xz[1] for p in boxGrids2D.allGridsInCartesian[:, 1]]
zCoordinates2D = [p.xz[2] for p in boxGrids2D.allGridsInCartesian[1, :]]
topographicSurfaceIndices2D = [
    findlast(@view material2D[ix, :]) for ix in axes(material2D, 1)]
@assert all(!isnothing, topographicSurfaceIndices2D) "a column is all air/void: widen altMax, or all solid: widen altMin"
topographicSurfaceIndices2D = Int.(topographicSurfaceIndices2D)
surfaceZ2D = [zCoordinates2D[k] for k in topographicSurfaceIndices2D]

# --- 3D topography (local x, y, z in km, matching SimuKirishima.ipynb's
#     own plotting convention) ---
x3DKm = [p.xyz[1] for p in boxGrids3D.allGridsInCartesian[:, 1, 1]] .* 1e-3
y3DKm = [p.xyz[2] for p in boxGrids3D.allGridsInCartesian[1, :, 1]] .* 1e-3
z3DKm = [p.xyz[3] for p in boxGrids3D.allGridsInCartesian[1, 1, :]] .* 1e-3
topography3D = Matrix{Float64}(undef, length(x3DKm), length(y3DKm))
for j in eachindex(y3DKm), i in eachindex(x3DKm)
    k = findlast(@view material3D[i, j, :])
    isnothing(k) && error("No material in column (i, j) = $((i, j))")
    k == length(z3DKm) && error("Topography reaches altMax in column (i, j) = $((i, j)); increase altMax")
    all(@view material3D[i, j, 1:k]) || error("Material mask has a cavity at $((i, j))")
    topography3D[i, j] = (z3DKm[k] + z3DKm[k + 1]) / 2
end

@show size(material2D) count(material2D) extrema(surfaceZ2D)
@show size(material3D) count(material3D) extrema(topography3D)
@show extrema(seismicModel2D.Vpv[material2D]) extrema(seismicModel2D.Vsv[material2D])

# 2D Vp/Vs cross-section, same as the figure the later FD3/OPT3/OPT4/OPT5
# steps will use to sanity-check their own heterogeneous model.
step1Figure2D = Figure(size=(1100, 420))
vpAxis2D = Axis(step1Figure2D[1, 1]; xlabel="x (km)", ylabel="z (km)",
    aspect=DataAspect(), title="Kirishima 2-D slice — Vp")
vsAxis2D = Axis(step1Figure2D[1, 2]; xlabel="x (km)", ylabel="z (km)",
    aspect=DataAspect(), title="Kirishima 2-D slice — Vs")
vpRange2D = extrema(seismicModel2D.Vpv[material2D])
vsRange2D = extrema(seismicModel2D.Vsv[material2D])
vpPlot2D = heatmap!(vpAxis2D, xCoordinates2D ./ 1e3, zCoordinates2D ./ 1e3,
    seismicModel2D.Vpv; colormap=:turbo, colorrange=vpRange2D)
vsPlot2D = heatmap!(vsAxis2D, xCoordinates2D ./ 1e3, zCoordinates2D ./ 1e3,
    seismicModel2D.Vsv; colormap=:turbo, colorrange=vsRange2D)
lines!(vpAxis2D, xCoordinates2D ./ 1e3, surfaceZ2D ./ 1e3; color=:white, linewidth=2)
lines!(vsAxis2D, xCoordinates2D ./ 1e3, surfaceZ2D ./ 1e3; color=:white, linewidth=2)
Colorbar(step1Figure2D[2, 1], vpPlot2D; vertical=false, label="Vp (km/s)")
Colorbar(step1Figure2D[2, 2], vsPlot2D; vertical=false, label="Vs (km/s)")
devSave(step1Figure2D, "kirishima_2D_slice_Vp_Vs")

# 3D shaded-relief check, same recipe as SimuKirishima.ipynb.
step1Figure3D = Figure(size=(1000, 800))
axis3D = Axis3(step1Figure3D[1, 1]; xlabel="local x (km)", ylabel="local y (km)",
    zlabel="elevation (km)", title="Kirishima 3-D box topography", aspect=:data)
terrain3D = surface!(axis3D, x3DKm, y3DKm, topography3D; color=topography3D,
    colormap=:terrain, colorrange=extrema(topography3D))
contour3d!(axis3D, x3DKm, y3DKm, topography3D;
    levels=range(extrema(topography3D)...; length=15), color=(:black, 0.35), linewidth=1)
Colorbar(step1Figure3D[1, 2], terrain3D; label="elevation (km)")
devSave(step1Figure3D, "kirishima_3D_topography")
nothing
