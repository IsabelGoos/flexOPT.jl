# dev/kirishima_oct2026/stages/step3_config.jl
# Shared by both sub-cases (homog/heter). Same structure as
# step2_config.jl, adapted for a flat (z<=0) free surface: FD3 gets its
# own zero-traction-flux closure via boundary_geometry's default mode;
# OPT3/OPT5/convFD3 get the :pinned_void volume + the additive
# elasticTractionFree2D closure already wired into step2_opt_build.jl
# (its applyFreeSurface branches were written for exactly this case).
quickRun = !DEV.full   # run.jl --full  =>  quickRun = false
optSpatialRefinement = 1   # FD3, convFD3 and OPT3 share dx (no OPT refinement yet)
buildHigherOrderRecipes = true
runOPT4 = false # 4-point stencil is inherently asymmetric -- dropped, nothing to fix
runOPT5 = true  # ElasticLHSOperatorAudit's coefficient gate has already passed
runOPTWithoutGamma = true
runConvFDBeyondBoundarySafeTime = false
waveformQuantity = :displacement
runSPECFEM2D = true
temporalRefinement = 2
optTimeMarchingScheme = :predictor_corrector
optPredictorMode = :lumped_future
optCorrectorIterations = 3
optCorrectorTolerance = 1e-8
optCorrectorRelaxation = 1.0
# Keeps Step 3's SPECFEM2D case directory separate from Step 2's
# "..._no_free_surface" one -- see the caseDirectoryTag fallback added to
# step2_specfem_run.jl (shared, unchanged otherwise).
caseDirectoryTag = "flat_free_surface"

dx = dx2D
x = xCoordinatesStep3
z = zCoordinatesStep3
nx, nz = length(x), length(z)
sourcePosition = sourcePositionStep3
applyFreeSurface = true
applyCerjan = !quickRun
material = repeat(reshape(solidMaskStep3, 1, :), nx, 1)
vp0, vs0, rho0 = (
    median(seismicModel2DStep3.Vpv[:, solidMaskStep3]) * 1e3,
    median(seismicModel2DStep3.Vsv[:, solidMaskStep3]) * 1e3,
    median(seismicModel2DStep3.ρ[:, solidMaskStep3]) * 1e3,
)
cerjanCells = round(Int, 24 * referenceDx / dx)
cerjan = CerjanBoundarySpec(
    (cerjanCells, cerjanCells),
    (cerjanCells, 0); # no absorption at the top: it is the physical free surface
    damping=0.0053 / fdOptSpatialRefinement^2,
)
activeCerjan = applyCerjan ? cerjan : nothing
# FD3's own zero-traction-flux closure (boundary_geometry's default free
# surface mode for a solid/void material mask), matching the existing
# HomogeneousElastic2DBenchmark_FreeSurface.ipynb convention.
bcFD = boundary_geometry(material, (dx, dx); cerjan=activeCerjan)
receiverZ = sourcePosition.z + 5e3
receiverX = quickRun ?
    (sourcePosition.x .+ [-20e3, -15e3, -10e3, -5e3, 0.0]) :
    (sourcePosition.x .+ [-20e3, -10e3, 0.0, 10e3, 20e3])
receiverGridX = quickRun ? collect((sourcePosition.x-20e3):5e3:(sourcePosition.x)) :
    collect((sourcePosition.x-30e3):10e3:(sourcePosition.x+30e3))
# Clamped at -dx (not 0.0): a receiver belongs to the solid interior, not
# sitting exactly on the free-surface boundary node itself.
receiverGridZ = quickRun ? collect((sourcePosition.z-15e3):5e3:(sourcePosition.z)) :
    collect((sourcePosition.z-30e3):10e3:min(sourcePosition.z+10e3, -dx))
receiverGrid = [(x=xr, z=zr) for zr in receiverGridZ for xr in receiverGridX
    if hypot(xr - sourcePosition.x, zr - sourcePosition.z) >= 5e3]
# Left/right/bottom margins only -- like distanceToNearestPhysicalEdge
# above, the top is the physical free surface (Cerjan's upper-z damping
# is explicitly 0 cells there, see `cerjan` above), not an artificial
# absorbing edge that a receiver needs 20-30 km of clearance from.
receiverBoundaryMargins = [minimum((
    station.x - first(x), last(x) - station.x,
    station.z - first(z),
)) for station in receiverGrid]
@assert minimum(receiverBoundaryMargins) >= (quickRun ? 20e3 : 30e3)
duration = quickRun ? 6.4 : 12.0
outputSampling = 0.10
referencePointsPerSWavelength = 16
rickerSharpness = 1.5
minimumAcceptedPointsPerSWavelength = 10
coarsestSpacing = referenceDx
sourceFrequency = rickerSharpness * vs0 /
    (referencePointsPerSWavelength * coarsestSpacing)
sourceDelay = 1.5 / sourceFrequency
sourceForce = 1.0e10
sourceMode = :gaussian
sourceSpatialSigma = referenceDx
sourceCutoffSigmas = 2.0
function spatial_source_specs(xaxis, zaxis; center=sourcePosition,
    mode=sourceMode, sigma=sourceSpatialSigma, cutoff=sourceCutoffSigmas)
    if mode === :point
        return [(x=xaxis[argmin(abs.(xaxis .- center.x))],
            z=zaxis[argmin(abs.(zaxis .- center.z))], weight=1.0)]
    end
    specs = [(x=Float64(xv), z=Float64(zv),
        weight=exp(-((xv-center.x)^2 + (zv-center.z)^2) / (2sigma^2)))
        for xv in xaxis for zv in zaxis
        if hypot(xv-center.x, zv-center.z) <= cutoff*sigma]
    total = sum(s.weight for s in specs)
    [(x=s.x, z=s.z, weight=s.weight / total) for s in specs]
end
sourceSpecsReference = spatial_source_specs(x, z)
@assert isapprox(sum(s.weight for s in sourceSpecsReference), 1.0; atol=1e-14)
rickerSource(time) = begin
    a = π * sourceFrequency * (time - sourceDelay)
    (1 - 2a^2) * exp(-a^2)
end
# Artificial (absorbing) edges only -- the top is the physical free
# surface, not an edge to worry about reflections from (matches
# HomogeneousElastic2DBenchmark_FreeSurface.ipynb's own
# distanceToNearestArtificialEdge).
distanceToNearestPhysicalEdge = minimum((
    sourcePosition.x - first(x), last(x) - sourcePosition.x,
    sourcePosition.z - first(z),
))
significantSourceStart = sourceDelay - 1 / sourceFrequency
@show applyFreeSurface applyCerjan (nx, nz) sourcePosition receiverZ dx
@show sourceFrequency sourceDelay rickerSharpness sourceForce sourceMode
@show vp0 vs0 rho0 cerjanCells distanceToNearestPhysicalEdge
receiverGeometryFigure = Figure(size=(720, 620))
receiverGeometryAxis = Axis(receiverGeometryFigure[1, 1];
    xlabel="x (km)", ylabel="z (km)", aspect=DataAspect(),
    title="Step 3 receiver grid around the Kirishima-structure source (flat free surface)")
scatter!(receiverGeometryAxis, getproperty.(receiverGrid, :x) ./ 1e3,
    getproperty.(receiverGrid, :z) ./ 1e3; label="receiver grid")
hlines!(receiverGeometryAxis, [0.0]; color=:black, linestyle=:dash, label="flat free surface")
scatter!(receiverGeometryAxis, [sourcePosition.x / 1e3], [sourcePosition.z / 1e3];
    marker=:star5, markersize=20, color=:gold, strokecolor=:black, label="source")
axislegend(receiverGeometryAxis)
devSave(receiverGeometryFigure, "receiver_geometry")
