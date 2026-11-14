# dev/kirishima_oct2026/stages/step2_box.jl
# Ported from notebooks/Elastic2DtopoBenchmark_oct2026.ipynb, cell step2-box-model.

# A Kirishima sub-volume that stays entirely below the real topography, so
# there is no free surface to model yet. Reuses `p0`/`niedVelocitySource`/
# `velocitySourceKey` from Step 1 as the geographic reference, but the box
# itself is centred ON THE SOURCE (offset step2SourceOffsetX from p0,
# depth step2CenterDepth), exactly like the existing
# Homogeneous/HeterogeneousElastic2DBenchmark notebooks -- that guarantees
# the full step2HalfWidth margin in every direction from the source to the
# absorbing boundary, which the receiver-grid safety assertion below needs.
# (A box centred on p0 instead left only ~15 km on the near side, since the
# source sits 10 km off-centre -- not enough margin at quickRun's 20 km
# threshold.)
referenceDx = 500.0        # same physical scale as the existing benchmarks
fdOptSpatialRefinement = 2 # dx = 250 m, the proven stable FD3/OPT3 setup
dx2D = referenceDx / fdOptSpatialRefinement
step2HalfWidth = 50.e3     # enlarged from 45 km: the real Kirishima vpMaximum
                           # (possibly uppermost-mantle, not just crustal) needs
                           # more margin to clear the edge-arrival assertion below
step2CenterDepth = -55.e3  # keeps the shallow bound at -5 km (step2CenterDepth +
                           # step2HalfWidth == -5e3), same safety margin as before
step2SourceOffsetX = -10.e3 # same source offset used by the existing benchmarks
boxGrids2DStep2 = constructLocalBox(
    p0, dx2D, dx2D,
    step2SourceOffsetX - step2HalfWidth, step2SourceOffsetX + step2HalfWidth,
    step2CenterDepth - step2HalfWidth, step2CenterDepth + step2HalfWidth,
)
gridCacheTagStep2 = "$(join(size(boxGrids2DStep2.allGridsInGeoPoints), 'x'))_" *
    "d$(boxGrids2DStep2.Δx)_$(boxGrids2DStep2.Δz)"
modelCacheNameStep2 = "seismicModel2DStep2_Kirishima_$(gridCacheTagStep2)_NIED_$(velocitySourceKey)"
seismicModel2DStep2 = lazyProduceOrLoad(
    modelCacheNameStep2, getParamsAndTopo,
    boxGrids2DStep2.allGridsInGeoPoints, boxGrids2DStep2.effectiveRadii, 0.1;
    folder=KIRISHIMA_MODEL_CACHE, velocity_model=:NIED, nied_source=niedVelocitySource,
    nied_confidence_max=0.8, nied_outside=:planet1D, nied_low_confidence=:planet1D,
)
@assert size(seismicModel2DStep2.ρ) == size(boxGrids2DStep2.allGridsInGeoPoints)
# Sanity check: the whole box must already be solid material (no void/air),
# since Step 2 has no free surface by construction, not by masking.
@assert all(seismicModel2DStep2.ρ .> air_density_cutoff) "step2CenterDepth/step2HalfWidth let the box reach the surface -- push it deeper"
xCoordinatesStep2 = [p.xz[1] for p in boxGrids2DStep2.allGridsInCartesian[:, 1]]
zCoordinatesStep2 = [p.xz[2] for p in boxGrids2DStep2.allGridsInCartesian[1, :]]
nxStep2, nzStep2 = length(xCoordinatesStep2), length(zCoordinatesStep2)
sourcePositionStep2 = (x=step2SourceOffsetX, z=step2CenterDepth)
@show size(seismicModel2DStep2.ρ) extrema(seismicModel2DStep2.Vpv) extrema(seismicModel2DStep2.Vsv)
@show nxStep2 nzStep2 dx2D step2HalfWidth step2CenterDepth sourcePositionStep2
