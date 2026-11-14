# dev/kirishima_oct2026/stages/step2_fd.jl
#
# Explicit FD3 (elasticWave2D) for the current sub-case, cached with
# myProduceOrLoad; on a miss step2_fd_propagate.jl (notebook code) runs.

fdCacheParams = Dict{String,Any}(
    "stage" => "step2_fd",
    "codeVersion" => CODE_VERSION["step2_fd"],
    "case" => caseLabel,
    "modelHash" => modelHash,
    "grid" => (x0=first(x), z0=first(z), nx=nx, nz=nz, dx=dx),
    "duration" => duration,
    "outputSampling" => outputSampling,
    "source" => (position=sourcePosition, frequency=sourceFrequency, delay=sourceDelay,
        force=sourceForce, mode=sourceMode, sigma=sourceSpatialSigma, cutoff=sourceCutoffSigmas),
    "boundaries" => (applyCerjan=applyCerjan, cerjanCells=cerjanCells, damping=cerjan.damping),
)
fdCachePrefix = "fd_$(caseLabel)"
devRedo!("fd", fdCacheParams, KIRISHIMA_CACHE_DIR, fdCachePrefix)
fdCacheFile = devCacheFile(fdCacheParams, KIRISHIMA_CACHE_DIR, fdCachePrefix)
println("FD3 cache: ", fdCacheFile, isfile(fdCacheFile) ? "  (hit)" : "  (miss -> computing)")

fdResult = myProduceOrLoad(fdCacheParams, KIRISHIMA_CACHE_DIR, fdCachePrefix) do config
    include(joinpath(KIRISHIMA_STAGES, "step2_fd_propagate.jl"))  # last expression = Dict
end

uxFD, uzFD, fdTimes = fdResult["uxFD"], fdResult["uzFD"], fdResult["fdTimes"]
fdCoordinates, fdPadding, fdDt = fdResult["fdCoordinates"], fdResult["fdPadding"], fdResult["fdDt"]
@show fdDt size(uzFD) maximum(abs, uzFD)
