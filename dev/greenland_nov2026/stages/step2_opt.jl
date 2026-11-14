# dev/greenland_nov2026/stages/step2_opt.jl
#
# OPT3 3-D elastic propagation, cached with myProduceOrLoad. On a cache hit
# nothing is assembled: the station traces and surface snapshots are loaded.
# On a miss:
#   step2_build.jl  -- recipe (makeOPTsemiSymbolic, itself cached in
#                      data/semiSymbolic), numerical operators, prepared
#                      linear system, 3-D free-surface closure;
#   step2_march.jl  -- predictor-corrector time marching, receivers only.

dtOPT = courant * dx / (sqrt(3) * vpMaximum * 1e3)
optSteps = ceil(Int, duration / dtOPT)
optCacheParams = Dict{String,Any}(
    "stage" => "opt",
    "codeVersion" => CODE_VERSION["opt"],
    "modelHash" => modelHash,
    "case" => greenlandCase,
    "dt" => dtOPT, "duration" => duration,
    "cerjan" => (cells=cerjanCells, damping=cerjanDamping),
    "source" => SOURCE, "sourceIndex" => Tuple(sourceIndex),
    "stations" => [(s.code, Tuple(s.index)) for s in stations],
    "recipe" => (supplementaryOrder=optSupplementaryOrder, gpu=recipeOnGPU),
    "freeSurface" => freeSurfaceClosure,
    "marching" => (scheme=optTimeMarchingScheme, predictor=optPredictorMode,
        iterations=optCorrectorIterations, tolerance=optCorrectorTolerance,
        relaxation=optCorrectorRelaxation),
    "snapshotEvery" => snapshotEvery,
)
devRedo!("opt", optCacheParams, GREENLAND_CACHE_DIR, "opt_$(greenlandCase)")
optCacheFile = devCacheFile(optCacheParams, GREENLAND_CACHE_DIR, "opt_$(greenlandCase)")
println("OPT cache: ", optCacheFile, isfile(optCacheFile) ? "  (hit)" : "  (miss -> computing)")
@printf("dt = %.3f s, %d steps for %.0f s\n", dtOPT, optSteps, duration)

optResult = myProduceOrLoad(optCacheParams, GREENLAND_CACHE_DIR, "opt_$(greenlandCase)") do cfg
    include(joinpath(GREENLAND_STAGES, "step2_build.jl"))
    include(joinpath(GREENLAND_STAGES, "step2_march.jl"))   # last expression = Dict
end

traceTimes = optResult["times"]
tracesBox = optResult["tracesBox"]          # (station, x|y|z of the box, time), metres
tracesZNE = optResult["tracesZNE"]          # (station, Z|N|E, time), metres
pc = optResult["pcDiagnostics"]
@printf("P-C: iterations per step %s, accepted relative residual in [%.1e, %.1e]\n",
    string(extrema(pc.iterations)), extrema(pc.residuals)...)
@show optResult["timings"]
optResult["stoppedEarly"] && @warn "the run stopped early (blow-up guard); traces are partial"
