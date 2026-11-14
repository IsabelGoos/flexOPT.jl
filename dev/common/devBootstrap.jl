# dev/common/devBootstrap.jl
#
# Shared plumbing for the dev/<sub_project>/ scripts (Claude-authored, see
# dev/README.md). Nothing here touches src/: it only
#   * finds the repository root and activates its environment,
#   * parses the command-line options of every dev/<sub>/run.jl,
#   * gives two tiny utilities: `devRedo!` (force one cached stage to be
#     recomputed) and `devSave` (save a figure under dev/<sub>/outputs/).
#
# Caching itself is ALWAYS done explicitly at the call site with
# `myProduceOrLoad(f, params, cacheDir, prefix)` from
# src/batchFiles/batchDrWatson.jl, so you can read in each stage file what
# is cached, where, and under which parameters.

import Pkg

function find_flexopt_root(start_dir=@__DIR__)
    candidates = String[]
    haskey(ENV, "FLEXOPT_ROOT") && push!(candidates, abspath(expanduser(ENV["FLEXOPT_ROOT"])))
    directory = abspath(start_dir)
    while true
        push!(candidates, directory)
        parent = dirname(directory)
        parent == directory && break
        directory = parent
    end
    for candidate in unique(candidates)
        isfile(joinpath(candidate, "src", "flexOPT.jl")) && return candidate
    end
    error("Cannot locate flexOPT; set ENV[\"FLEXOPT_ROOT\"]")
end

if !@isdefined(flexopt_root)
    flexopt_root = find_flexopt_root()
    Pkg.activate(flexopt_root)   # => DrWatson's datadir() == <repo>/data
end
using DrWatson, JLD2

"""
    parseDevArgs(args; steps, defaultStep, cases=String[], defaultCases=cases)

Command line of every dev/<sub>/run.jl:

    run.jl <step> [case ...] [--redo=stageA,stageB | --redo=all] [--full] [--reload]

* `<step>`    one of `steps` (also `all`)
* `case`      optional sub-cases (Kirishima: homog / heter), default all
* `--redo=`   delete those stages' cache files before running, so they are recomputed
* `--full`    the long (non-quickRun) setting where the sub-project has one
* `--reload`  re-read helpers/ in a running REPL (after editing src/: restart Julia;
              stages/*.jl are re-read on every run anyway)
"""
function parseDevArgs(args; steps, defaultStep, cases=String[], defaultCases=cases)
    positional = filter(a -> !startswith(a, "--"), args)
    redo = Set{String}()
    full = false
    reload = false
    for a in args
        if startswith(a, "--redo=")
            union!(redo, split(a[8:end], ','))
        elseif a == "--full"
            full = true
        elseif a == "--reload"
            reload = true
        elseif startswith(a, "--")
            error("unknown option $a")
        end
    end
    step = isempty(positional) ? defaultStep : positional[1]
    step == "all" || step in steps || error("unknown step '$step'; choose among $(join(steps, ", ")) or all")
    chosenCases = length(positional) > 1 ? positional[2:end] : defaultCases
    for c in chosenCases
        c in cases || error("unknown case '$c'; choose among $(join(cases, ", "))")
    end
    return (step=step, steps=(step == "all" ? collect(steps) : [step]),
        cases=collect(chosenCases), redo=redo, full=full, reload=reload)
end

"""
    devRedo!(stage, params, cacheDir, prefix)

If `stage` (or "all") was given with `--redo=`, delete the file that
`myProduceOrLoad(f, params, cacheDir, prefix)` would load, so the very next
`myProduceOrLoad` call recomputes and re-saves it. Uses the same hash as
src/batchFiles/batchDrWatson.jl. Call it right before the myProduceOrLoad.
"""
function devRedo!(stage::AbstractString, params, cacheDir::AbstractString, prefix::AbstractString)
    (stage in DEV.redo || "all" in DEV.redo) || return nothing
    file = datadir(cacheDir, "$(prefix)_$(commonBatchs.hash_parameters(params)).jld2")
    if isfile(file)
        rm(file)
        println("--redo=$stage: removed ", file)
    else
        println("--redo=$stage: nothing cached yet at ", file)
    end
    return nothing
end

"Cache file a `myProduceOrLoad(f, params, cacheDir, prefix)` call uses (for printing)."
devCacheFile(params, cacheDir, prefix) =
    datadir(cacheDir, "$(prefix)_$(commonBatchs.hash_parameters(params)).jld2")

"Save a Makie figure as dev/<sub>/outputs/<subdir>/<name>.png and print the path."
function devSave(fig, name::AbstractString; subdir=DEV_STEP[])
    clean = replace(String(name), r"[^A-Za-z0-9_.=+-]+" => "_")
    clean = strip(clean, '_')
    dir = joinpath(DEV_OUTPUT_DIR[], subdir)
    mkpath(dir)
    path = joinpath(dir, clean * ".png")
    save(path, fig)
    println("figure -> ", path)
    return path
end

# set by each run.jl (guarded so the file can be re-included in a REPL)
if !@isdefined(DEV_OUTPUT_DIR)
    const DEV_OUTPUT_DIR = Ref("")
    const DEV_STEP = Ref("")
end

"""
    devRun(runfile, args...)

REPL entry point used by the `kirishima(...)` / `spherical(...)` functions:
sets `DEV_ARGS` and includes the sub-project's run.jl in Main, in the same
session (so packages, src/ modules and cached results stay loaded).
"""
function devRun(runfile::AbstractString, args...)
    Core.eval(Main, :(DEV_ARGS = $(collect(String, string.(args)))))
    Base.include(Main, runfile)
end
