# dev/kirishima_oct2026/stages/step2_specfem.jl
#
# SPECFEM2D reference for the current sub-case, cached with myProduceOrLoad;
# on a miss step2_specfem_run.jl (notebook code) prepares and runs SPECFEM2D
# (into data/specfem2d_benchmarks/<case>_kirishima_no_free_surface) and reads
# the traces. Below: the notebook's post-processing, unchanged.

specfemCacheParams = Dict{String,Any}(
    "stage" => "step2_specfem",
    "codeVersion" => CODE_VERSION["step2_specfem"],
    "case" => caseLabel,
    "modelHash" => modelHash,
    "grid" => (x0=first(x), z0=first(z), nx=nx, nz=nz, dx=dx, referenceDx=referenceDx),
    "dt" => dtOPT,
    "duration" => duration,
    "outputSampling" => outputSampling,
    "receivers" => (x=receiverX, z=receiverZ),
    "source" => (position=sourcePosition, frequency=sourceFrequency, delay=sourceDelay,
        force=sourceForce, mode=sourceMode, sigma=sourceSpatialSigma, cutoff=sourceCutoffSigmas),
    "applyFreeSurface" => applyFreeSurface,
    "runSPECFEM2D" => runSPECFEM2D,
)
specfemCachePrefix = "specfem_$(caseLabel)"
devRedo!("specfem", specfemCacheParams, KIRISHIMA_CACHE_DIR, specfemCachePrefix)
specfemCacheFile = devCacheFile(specfemCacheParams, KIRISHIMA_CACHE_DIR, specfemCachePrefix)
println("SPECFEM2D cache: ", specfemCacheFile, isfile(specfemCacheFile) ? "  (hit)" : "  (miss -> computing)")

specfemResult = myProduceOrLoad(specfemCacheParams, KIRISHIMA_CACHE_DIR, specfemCachePrefix) do config
    include(joinpath(KIRISHIMA_STAGES, "step2_specfem_run.jl"))  # last expression = Dict
end

dtSPECFEM = dtOPT
specfemTracesX, specfemTracesZ = specfemResult["specfemTracesX"], specfemResult["specfemTracesZ"]
specfemCase = specfemResult["specfemCase"]
specfemCFL = specfemResult["specfemCFL"]
specfemSnapshotVideoPath = specfemResult["specfemSnapshotVideoPath"]
zSPECFEM = specfemResult["zSPECFEM"]
@assert length(specfemTracesZ) == length(receiverX)

function integrate_velocity_trace(trace)
    displacement = zeros(Float64, length(trace.values))
    displacement[2:end] .= cumsum(
        ((trace.values[1:end-1] .+ trace.values[2:end]) ./ 2) .*
        diff(trace.time))
    (time=Float64.(trace.time), values=displacement)
end
specfemWaveformsX = waveformQuantity === :velocity ? specfemTracesX :
    integrate_velocity_trace.(specfemTracesX)
specfemWaveformsZ = waveformQuantity === :velocity ? specfemTracesZ :
    integrate_velocity_trace.(specfemTracesZ)
specfemTraces = specfemWaveformsZ # selected-quantity compatibility alias
specfemTrace = specfemWaveformsZ[receiver]

cflSummary = (
    FD3_axis=vpMaximum * fdDt / dx,
    FD3_2D=sqrt(2) * vpMaximum * fdDt / dx,
    OPT3_axis=vpMaximum * dtOPT / dxOPT,
    OPT3_2D=sqrt(2) * vpMaximum * dtOPT / dxOPT,
    SPECFEM2D_reported=specfemCFL,
)
# SPECFEM uses 5 GLL points per element in this build. Its default number
# of elements makes the mean GLL-node interval comparable to dx.
specfemElementSize = (
    x=(last(x)-first(x)) / specfemCase.nx_elements,
    z=(last(zSPECFEM)-first(zSPECFEM)) / specfemCase.nz_elements,
)
discretizationSummary = (
    OPT3=(dx=dxOPT, dt=dtOPT),
    SPECFEM2D=(element_size=specfemElementSize,
        mean_GLL_interval=(x=specfemElementSize.x/4,
                           z=specfemElementSize.z/4), dt=dtSPECFEM),
)
@show cflSummary discretizationSummary specfemCase.time_axis_shift
@show specfemCase.case_directory dtSPECFEM specfemSnapshotVideoPath
