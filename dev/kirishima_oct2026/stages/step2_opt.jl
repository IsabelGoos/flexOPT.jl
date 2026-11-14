# dev/kirishima_oct2026/stages/step2_opt.jl
#
# convFD3 / OPT3 / OPT3 without Γ / OPT4 / OPT5 for the current sub-case,
# cached with myProduceOrLoad. On a cache hit nothing is assembled at all
# (no recipe loading, no operator construction): the wavefields are loaded.
# On a miss: step2_opt_build.jl (recipes + operators) then
# step2_opt_propagate.jl (time marching) run, exactly the notebook code.

optCacheParams = Dict{String,Any}(
    "stage" => "step2_opt",
    "codeVersion" => CODE_VERSION["step2_opt"],
    "case" => caseLabel,
    "modelHash" => modelHash,
    "grid" => (x0=first(x), z0=first(z), nx=nx, nz=nz, dx=dx),
    "optSpatialRefinement" => optSpatialRefinement,
    "temporalRefinement" => temporalRefinement,
    "duration" => duration,
    "outputSampling" => outputSampling,
    "source" => (position=sourcePosition, frequency=sourceFrequency, delay=sourceDelay,
        force=sourceForce, mode=sourceMode, sigma=sourceSpatialSigma, cutoff=sourceCutoffSigmas),
    "boundaries" => (applyCerjan=applyCerjan, applyFreeSurface=applyFreeSurface,
        cerjanCells=cerjanCells, damping=cerjan.damping),
    "runs" => (OPT4=runOPT4, OPT5=runOPT5, withoutGamma=runOPTWithoutGamma,
        convFDBeyondSafeTime=runConvFDBeyondBoundarySafeTime,
        higherOrderRecipes=buildHigherOrderRecipes),
    "marching" => (scheme=optTimeMarchingScheme, predictor=optPredictorMode,
        iterations=optCorrectorIterations, tolerance=optCorrectorTolerance,
        relaxation=optCorrectorRelaxation),
)
optCachePrefix = "opt_$(caseLabel)"
devRedo!("opt", optCacheParams, KIRISHIMA_CACHE_DIR, optCachePrefix)
optCacheFile = devCacheFile(optCacheParams, KIRISHIMA_CACHE_DIR, optCachePrefix)
println("OPT family cache: ", optCacheFile, isfile(optCacheFile) ? "  (hit)" : "  (miss -> computing)")

optResult = myProduceOrLoad(optCacheParams, KIRISHIMA_CACHE_DIR, optCachePrefix) do config
    include(joinpath(KIRISHIMA_STAGES, "step2_opt_build.jl"))
    include(joinpath(KIRISHIMA_STAGES, "step2_opt_propagate.jl"))  # last expression = Dict
end

xOPT, zOPT = optResult["xOPT"], optResult["zOPT"]
dxOPT, dtOPT, optSteps = optResult["dxOPT"], optResult["dtOPT"], optResult["optSteps"]
uxOPT, uzOPT, optTimes = optResult["uxOPT"], optResult["uzOPT"], optResult["optTimes"]
uxOPTWithoutGamma, uzOPTWithoutGamma = optResult["uxOPTWithoutGamma"], optResult["uzOPTWithoutGamma"]
uxOPT4, uzOPT4, opt4Times = optResult["uxOPT4"], optResult["uzOPT4"], optResult["opt4Times"]
uxOPT5, uzOPT5, opt5Times = optResult["uxOPT5"], optResult["uzOPT5"], optResult["opt5Times"]
uxConvFD, uzConvFD, convFDTimes = optResult["uxConvFD"], optResult["uzConvFD"], optResult["convFDTimes"]
convFDBoundarySafeTime, convFDDuration = optResult["convFDBoundarySafeTime"], optResult["convFDDuration"]
operatorTimings = optResult["operatorTimings"]

# predictor-corrector diagnostics are printed on every run, cached or not
if !isnothing(optResult["pcDiagnostics"])
    pc = optResult["pcDiagnostics"]
    @show extrema(pc.iterations_per_step) extrema(pc.relative_residuals)
end
@show optResult["optTiming"] operatorTimings
@show dtOPT size(uzOPT) maximum(abs, uzOPT) convFDDuration
