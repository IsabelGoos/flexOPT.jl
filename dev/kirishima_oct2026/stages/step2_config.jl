# dev/kirishima_oct2026/stages/step2_config.jl
# Ported from notebooks/Elastic2DtopoBenchmark_oct2026.ipynb, cell step2-shared-config.
# Shared by both sub-cases (homog / heter).
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
optTimeMarchingScheme = :predictor_corrector # exercise the P-C scheme now
optPredictorMode = :lumped_future
optCorrectorIterations = 3
optCorrectorTolerance = 1e-8
optCorrectorRelaxation = 1.0

dx = dx2D
x = xCoordinatesStep2
z = zCoordinatesStep2
nx, nz = length(x), length(z)
sourcePosition = sourcePositionStep2
applyFreeSurface = false
applyCerjan = !quickRun
material = trues(nx, nz)
vp0, vs0, rho0 = (
    median(seismicModel2DStep2.Vpv) * 1e3,
    median(seismicModel2DStep2.Vsv) * 1e3,
    median(seismicModel2DStep2.ρ) * 1e3,
)
cerjanCells = round(Int, 24 * referenceDx / dx)
cerjan = CerjanBoundarySpec(
    (cerjanCells, cerjanCells),
    (cerjanCells, cerjanCells); # no free surface: absorb on every side
    damping=0.0053 / fdOptSpatialRefinement^2,
)
activeCerjan = applyCerjan ? cerjan : nothing
bcFD = BoundaryConditionSet(free_surface=nothing, cerjan=activeCerjan,
    material_mask=BitArray(material), free_surface_mode=:pinned_void)
receiverZ = sourcePosition.z + 5e3
receiverX = quickRun ?
    (sourcePosition.x .+ [-20e3, -15e3, -10e3, -5e3, 0.0]) :
    (sourcePosition.x .+ [-20e3, -10e3, 0.0, 10e3, 20e3])
receiverGridX = quickRun ? collect((sourcePosition.x-20e3):5e3:(sourcePosition.x)) :
    collect((sourcePosition.x-30e3):10e3:(sourcePosition.x+30e3))
receiverGridZ = quickRun ? collect((sourcePosition.z-15e3):5e3:(sourcePosition.z)) :
    collect((sourcePosition.z-30e3):10e3:(sourcePosition.z+10e3))
receiverGrid = [(x=xr, z=zr) for zr in receiverGridZ for xr in receiverGridX
    if hypot(xr - sourcePosition.x, zr - sourcePosition.z) >= 5e3]
receiverBoundaryMargins = [minimum((
    station.x - first(x), last(x) - station.x,
    station.z - first(z), last(z) - station.z,
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
distanceToNearestPhysicalEdge = minimum((
    sourcePosition.x - first(x), last(x) - sourcePosition.x,
    sourcePosition.z - first(z), last(z) - sourcePosition.z,
))
significantSourceStart = sourceDelay - 1 / sourceFrequency
@show applyFreeSurface applyCerjan (nx, nz) sourcePosition receiverZ dx
@show sourceFrequency sourceDelay rickerSharpness sourceForce sourceMode
@show vp0 vs0 rho0 cerjanCells distanceToNearestPhysicalEdge
receiverGeometryFigure = Figure(size=(720, 620))
receiverGeometryAxis = Axis(receiverGeometryFigure[1, 1];
    xlabel="x (km)", ylabel="z (km)", aspect=DataAspect(),
    title="Step 2 receiver grid around the Kirishima-structure source")
scatter!(receiverGeometryAxis, getproperty.(receiverGrid, :x) ./ 1e3,
    getproperty.(receiverGrid, :z) ./ 1e3; label="receiver grid")
scatter!(receiverGeometryAxis, [sourcePosition.x / 1e3], [sourcePosition.z / 1e3];
    marker=:star5, markersize=20, color=:gold, strokecolor=:black, label="source")
axislegend(receiverGeometryAxis)
devSave(receiverGeometryFigure, "receiver_geometry")
