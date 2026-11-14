# dev/kirishima_oct2026/run.jl  -- Kirishima elastic 2-D benchmark, script version
#
# Same code as notebooks/Elastic2DtopoBenchmark_oct2026.ipynb (Steps 1-2),
# split into stages/*.jl; the expensive stages are cached with
# myProduceOrLoad (see stages/step2_opt.jl, step2_fd.jl, step2_specfem.jl).
#
# From a terminal, in the repository root:
#   julia --project=. dev/kirishima_oct2026/run.jl step2              # both cases
#   julia --project=. dev/kirishima_oct2026/run.jl step2 heter        # one case
#   julia --project=. dev/kirishima_oct2026/run.jl step2 heter --redo=opt
#   julia --project=. dev/kirishima_oct2026/run.jl step3              # flat free surface, both cases
#   julia --project=. dev/kirishima_oct2026/run.jl step1
# From a REPL (fastest: everything stays compiled and loaded):
#   include("dev/kirishima_oct2026/run.jl")      # first time (runs the default: step2, both)
#   kirishima("step2", "heter", "--redo=opt")    # afterwards
#   kirishima("step3", "homog")
# Figures and metrics: dev/kirishima_oct2026/outputs/<step>[_<case>]/
#
# Step 3 (flat free surface) reuses Step 2's opt/fd/traces/specfem/compare
# stages unchanged (they already branch on applyFreeSurface); only the box
# (step3_box.jl), shared config (step3_config.jl) and per-case model
# (step3_model.jl) are new -- see those files' own comments.

include(joinpath(@__DIR__, "..", "common", "devBootstrap.jl"))
KIRISHIMA_STAGES = joinpath(@__DIR__, "stages")
kirishima(args...) = devRun(joinpath(@__DIR__, "run.jl"), args...)

DEV = parseDevArgs(@isdefined(DEV_ARGS) ? DEV_ARGS : ARGS;
    steps=["step1", "step2", "step3"], defaultStep="step2",
    cases=["homog", "heter"])
DEV_OUTPUT_DIR[] = joinpath(@__DIR__, "outputs")
println("Kirishima dev run: ", DEV)

include(joinpath(KIRISHIMA_STAGES, "env_step1.jl"))
include(joinpath(@__DIR__, "codeVersions.jl"))
include(joinpath(KIRISHIMA_STAGES, "shared_site.jl"))

STEP2_CASE_STAGES = (
    "step2_model.jl",     # homog / heter fields, safety checks, modelHash
    "step2_opt.jl",       # cached: recipes + operators + convFD3/OPT3/(OPT4)/OPT5 propagation
    "step2_fd.jl",        # cached: explicit FD3
    "step2_traces.jl",    # receiver traces, pre-SPECFEM figure
    "step2_specfem.jl",   # cached: SPECFEM2D reference
    "step2_compare.jl",   # metrics vs SPECFEM2D, figure + metrics file
)
STEP3_CASE_STAGES = (
    "step3_model.jl",     # homog / heter fields restricted to the solid z<=0 part
    "step2_opt.jl",       # same stage file: applyFreeSurface=true now closes the
                           # OPT/convFD volume with the additive elasticTractionFree2D recipe
    "step2_fd.jl",        # same stage file: bcFD (from step3_config.jl) already has FD3's
                           # own zero-traction-flux closure
    "step2_traces.jl",
    "step2_specfem.jl",   # case directory tagged "flat_free_surface", separate from Step 2's
    "step2_compare.jl",
)

for step in DEV.steps
    DEV_STEP[] = step
    if step == "step1"
        @time "step1" include(joinpath(KIRISHIMA_STAGES, "step1_box.jl"))
    elseif step == "step2"
        include(joinpath(KIRISHIMA_STAGES, "step2_env.jl"))
        include(joinpath(KIRISHIMA_STAGES, "step2_box.jl"))
        include(joinpath(KIRISHIMA_STAGES, "step2_config.jl"))
        for case in DEV.cases
            global step2Case = Symbol(case)
            DEV_STEP[] = "step2_" * case
            println("\n================ Step 2 -- ", case, " ================")
            for stage in STEP2_CASE_STAGES
                println("---- ", stage)
                @time stage include(joinpath(KIRISHIMA_STAGES, stage))
            end
        end
    elseif step == "step3"
        include(joinpath(KIRISHIMA_STAGES, "step2_env.jl"))
        include(joinpath(KIRISHIMA_STAGES, "step3_box.jl"))
        include(joinpath(KIRISHIMA_STAGES, "step3_config.jl"))
        for case in DEV.cases
            global step3Case = Symbol(case)
            DEV_STEP[] = "step3_" * case
            println("\n================ Step 3 -- ", case, " (flat free surface) ================")
            for stage in STEP3_CASE_STAGES
                println("---- ", stage)
                @time stage include(joinpath(KIRISHIMA_STAGES, stage))
            end
        end
    end
end
