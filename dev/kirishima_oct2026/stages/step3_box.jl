# dev/kirishima_oct2026/stages/step3_box.jl
# Step 3: flat free surface over the real Kirishima subsurface structure.
# Reuses Step 1's p0/niedVelocitySource/velocitySourceKey and Step 2's
# box-centred-on-the-source convention (offset step3SourceOffsetX from p0),
# but the surface is flattened to z=0 (sea level) by construction, not
# masked from the real NIED topography -- that is Step 4. This matches the
# existing applyFreeSurface branches already wired into
# step2_opt_build.jl/step2_specfem_run.jl (both assume a flat z<=0 split
# whenever applyFreeSurface=true -- see notebooks/
# HomogeneousElastic2DBenchmark_FreeSurface.ipynb, which uses the same
# z<=0.0 convention).
#
# Lateral/below-source margins reuse Step 2's own empirically-fixed 50 km
# clearance: the SAME real Kirishima field's vpMaximum (likely
# uppermost-mantle, sampled whenever the box goes deep enough) needed that
# much room to clear the edge-arrival safety assertion there, and nothing
# about this shallower box rules out sampling the same material at its
# own, shallower bottom.
referenceDx = 500.0        # same physical scale as Step 2 / the existing benchmarks
fdOptSpatialRefinement = 2 # dx = 250 m, the proven stable FD3/OPT3 setup
dx2D = referenceDx / fdOptSpatialRefinement
step3SourceOffsetX = -10.e3    # same x convention as Step 2
step3SourceDepth = -8.e3       # shallow/crustal: a meaningful free-surface test
step3LateralHalfWidth = 50.e3  # same clearance that fixed Step 2's heterogeneous case
step3BelowMargin = 50.e3       # idem, measured below the source
step3AboveSurfaceVoid = 5.e3   # void headroom above the flat z=0 surface
boxGrids2DStep3 = constructLocalBox(
    p0, dx2D, dx2D,
    step3SourceOffsetX - step3LateralHalfWidth, step3SourceOffsetX + step3LateralHalfWidth,
    step3SourceDepth - step3BelowMargin, step3AboveSurfaceVoid,
)
gridCacheTagStep3 = "$(join(size(boxGrids2DStep3.allGridsInGeoPoints), 'x'))_" *
    "d$(boxGrids2DStep3.Δx)_$(boxGrids2DStep3.Δz)"
modelCacheNameStep3 = "seismicModel2DStep3_Kirishima_$(gridCacheTagStep3)_NIED_$(velocitySourceKey)"
seismicModel2DStep3 = lazyProduceOrLoad(
    modelCacheNameStep3, getParamsAndTopo,
    boxGrids2DStep3.allGridsInGeoPoints, boxGrids2DStep3.effectiveRadii, 0.1;
    folder=KIRISHIMA_MODEL_CACHE, velocity_model=:NIED, nied_source=niedVelocitySource,
    nied_confidence_max=0.8, nied_outside=:planet1D, nied_low_confidence=:planet1D,
)
@assert size(seismicModel2DStep3.ρ) == size(boxGrids2DStep3.allGridsInGeoPoints)
xCoordinatesStep3 = [p.xz[1] for p in boxGrids2DStep3.allGridsInCartesian[:, 1]]
zCoordinatesStep3 = [p.xz[2] for p in boxGrids2DStep3.allGridsInCartesian[1, :]]
nxStep3, nzStep3 = length(xCoordinatesStep3), length(zCoordinatesStep3)
sourcePositionStep3 = (x=step3SourceOffsetX, z=step3SourceDepth)
# The material mask below is flat (z<=0), not the real NIED air/rock split:
# Step 3 deliberately idealises the local topography as flat sea level, so
# only the SOLID (z<=0) half of this fetched field needs to be real rock --
# the void half's fetched values are never used (masked out downstream).
solidMaskStep3 = zCoordinatesStep3 .<= 0.0
@assert any(solidMaskStep3) && !all(solidMaskStep3) "step3AboveSurfaceVoid/step3BelowMargin must leave both a solid and a void part of the box"
# The flat z<=0 cut is an idealisation: a few columns can still land on
# real coastline/low-lying terrain (true air or water below the chosen
# sea-level datum) even though most of the solid half is genuine rock.
# Patch those cells to the solid region's own median instead of hard
# failing the whole box over what should be a small edge effect -- but
# still fail loudly if the patched fraction is too large to trust (the
# box would then be genuinely badly placed, not just brushing the coast).
solidRhoStep3 = @view seismicModel2DStep3.ρ[:, solidMaskStep3]
badSolidStep3 = solidRhoStep3 .<= air_density_cutoff
nBadSolidStep3, nSolidTotalStep3 = count(badSolidStep3), length(solidRhoStep3)
badSolidFractionStep3 = nBadSolidStep3 / nSolidTotalStep3
if nBadSolidStep3 > 0
    solidVpStep3 = @view seismicModel2DStep3.Vpv[:, solidMaskStep3]
    solidVsStep3 = @view seismicModel2DStep3.Vsv[:, solidMaskStep3]
    goodMedianVpStep3 = median(solidVpStep3[.!badSolidStep3])
    goodMedianVsStep3 = median(solidVsStep3[.!badSolidStep3])
    goodMedianRhoStep3 = median(solidRhoStep3[.!badSolidStep3])
    solidVpStep3[badSolidStep3] .= goodMedianVpStep3
    solidVsStep3[badSolidStep3] .= goodMedianVsStep3
    solidRhoStep3[badSolidStep3] .= goodMedianRhoStep3
    @warn "step3_box: patched real air/water cells inside the flat z<=0 solid region with the solid median (likely coastline/low terrain inside this box's footprint)" nBadSolidStep3 nSolidTotalStep3 badSolidFractionStep3
end
@assert badSolidFractionStep3 < 0.05 "too much of the nominally-solid (z<=0) part of this box is real air/water ($(round(100*badSolidFractionStep3; digits=1))%) -- the box reaches real coastline/sea; reduce step3LateralHalfWidth or move step3SourceOffsetX"
@show size(seismicModel2DStep3.ρ) extrema(seismicModel2DStep3.Vpv[:, solidMaskStep3]) extrema(seismicModel2DStep3.Vsv[:, solidMaskStep3])
@show nxStep3 nzStep3 dx2D step3LateralHalfWidth step3SourceDepth sourcePositionStep3
