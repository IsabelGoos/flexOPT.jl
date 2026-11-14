# dev/greenland_nov2026/run.jl -- 3-D elastic waves through Greenland (OPT3)
#
# 16 GLISN stations of Toyokuni et al. (2018, PEPI), PREM + ETOPO 2022 ice
# sheet and ocean, true curved surface in a GeoPoints local box, single force
# (ice slide) on the coast near Nuuk.
#
# From a terminal, in the repository root:
#   (step1: model + figure; step2: recipe, operators, propagation;
#    step3: traces, record section, snapshots; 2nd argument = grid case)
#   julia --project=. -t auto dev/greenland_nov2026/run.jl step1
#   julia --project=. -t auto dev/greenland_nov2026/run.jl step2
#   julia --project=. -t auto dev/greenland_nov2026/run.jl step3
#   julia --project=. -t auto dev/greenland_nov2026/run.jl step4      (propagation movie, MP4)
#   julia --project=. -t auto dev/greenland_nov2026/run.jl all standard
#   julia --project=. -t auto dev/greenland_nov2026/run.jl step2 full --redo=opt
# From a REPL (fastest):
#   include("dev/greenland_nov2026/run.jl")        # default: step1 quick
#   greenland("step2")                             # afterwards
# Grids (one case at a time): quick (50 km, periods ≳ 60 s), standard (30 km, ≳ 40 s),
# full (20 km, ≳ 25 s). Outputs: dev/greenland_nov2026/outputs/<step>_<case>/
# Cached stages (--redo=...): etopo, model, opt.  Recipes: data/semiSymbolic.

include(joinpath(@__DIR__, "..", "common", "devBootstrap.jl"))
GREENLAND_STAGES = joinpath(@__DIR__, "stages")
greenland(args...) = devRun(joinpath(@__DIR__, "run.jl"), args...)

# zsh passes "# comment" text as arguments (interactive_comments is off by
# default): drop everything from the first argument starting with '#'.
greenlandArgs = let a = collect(String, @isdefined(DEV_ARGS) ? DEV_ARGS : ARGS)
    k = findfirst(x -> startswith(x, "#"), a)
    isnothing(k) ? a : a[1:k-1]
end
DEV = parseDevArgs(greenlandArgs;
    steps=["step1", "step2", "step3", "step4"], defaultStep="step1",
    cases=["quick", "standard", "full"], defaultCases=["quick"])
length(DEV.cases) == 1 || error("give one grid case at a time (quick | standard | full)")
DEV_OUTPUT_DIR[] = joinpath(@__DIR__, "outputs")
println("Greenland dev run: ", DEV)

include(joinpath(GREENLAND_STAGES, "env.jl"))
include(joinpath(@__DIR__, "codeVersions.jl"))
greenlandCase = only(DEV.cases)
include(joinpath(GREENLAND_STAGES, "config.jl"))

# every step needs the model (cache hit after the first step1)
for stage in ("step1_box.jl", "step1_model.jl", "step1_sites.jl")
    println("---- ", stage)
    @time stage include(joinpath(GREENLAND_STAGES, stage))
end
for step in DEV.steps
    DEV_STEP[] = step * "_" * greenlandCase
    println("\n================ Greenland ", step, " -- ", greenlandCase, " ================")
    if step == "step1"
        include(joinpath(GREENLAND_STAGES, "step1_figures.jl"))
        include(joinpath(GREENLAND_STAGES, "step1_figures3d.jl"))   # 3-D views
    elseif step == "step2"
        @time "step2" include(joinpath(GREENLAND_STAGES, "step2_opt.jl"))
    elseif step == "step3"
        @isdefined(optResult) || include(joinpath(GREENLAND_STAGES, "step2_opt.jl"))  # cache hit
        include(joinpath(GREENLAND_STAGES, "step3_traces.jl"))
    elseif step == "step4"
        @isdefined(optResult) || include(joinpath(GREENLAND_STAGES, "step2_opt.jl"))  # cache hit
        include(joinpath(GREENLAND_STAGES, "step4_movie.jl"))
    end
end
