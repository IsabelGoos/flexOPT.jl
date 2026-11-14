# application/common/appBootstrap.jl
#
# Shared loader for every application/<Name>/run.jl (Claude-authored, see
# application/README.md). It does NOT modify src/: it only
#   * finds the repository root and activates its environment,
#   * includes the src/ modules an application asks for, in the right order,
#     at most once per Julia session (safe to re-include in a REPL),
#   * gives a few tiny utilities: `appArgs` (command-line parsing),
#     `appOutputDir`, `appSave`.
#
# Usage, at the top of a run.jl:
#
#     FLEXOPT_MODULES = [:flexOPT]               # what this application needs
#     include(joinpath(@__DIR__, "..", "common", "appBootstrap.jl"))
#
# Recognised module names (dependencies are pulled in automatically):
#   :commonBatchs  :planet1D  :GeoPoints  :flexOPT  :fieldOutput
#   :sonification  :neutrinoOscillation  :neutrinoFlux  :neutrinoDetector
#   :neutrinoCrossSections  :seismicStations  :elasticWave2D  :elasticWave3D
#   :elasticGreens2D  :specfemBenchmark
#
# GPU: the default is the CPU (KernelAbstractions.CPU()), which is what the
# recipes are built on anyway. To use the Apple GPU like the notebooks, start
# Julia with FLEXOPT_GPU=metal in the environment.

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
    error("Cannot locate flexOPT; set ENV[\"FLEXOPT_ROOT\"] to the clone's absolute path")
end

if !@isdefined(flexopt_root)
    flexopt_root = find_flexopt_root()
    Base.active_project() == joinpath(flexopt_root, "Project.toml") || Pkg.activate(flexopt_root)
end

# ---------------------------------------------------------------- src modules
# name => (file under src/, modules that must be loaded first)
_APP_SRC_MODULES = Dict{Symbol,Tuple{String,Vector{Symbol}}}(
    :commonBatchs          => ("commonBatchs.jl", Symbol[]),
    :planet1D              => ("planet1D.jl", Symbol[]),
    :GeoPoints             => ("GeoPoints.jl", [:commonBatchs, :planet1D]),
    :flexOPT               => ("flexOPT.jl", [:commonBatchs]),
    :fieldOutput           => ("fieldOutput.jl", Symbol[]),
    :sonification          => ("sonification.jl", Symbol[]),
    :neutrinoOscillation   => ("neutrinoOscillation.jl", Symbol[]),
    :neutrinoFlux          => ("neutrinoFlux.jl", [:commonBatchs]),
    :neutrinoDetector      => ("neutrinoDetector.jl", Symbol[]),
    :neutrinoCrossSections => ("neutrinoCrossSections.jl", Symbol[]),
    :seismicStations       => ("seismicStations.jl", [:GeoPoints]),
    :elasticWave2D         => ("elasticWave2D.jl", Symbol[]),
    :elasticWave3D         => ("elasticWave3D.jl", Symbol[]),
    :elasticGreens2D       => ("elasticGreens2D.jl", Symbol[]),
    :specfemBenchmark      => ("specfemBenchmark.jl", Symbol[]),
)
# the Julia module name inside the file, when it differs from the key
_APP_MODULE_NAME = Dict(:neutrinoCrossSections => :Neutrino_Cross_Sections)

function _appLoadOrder(wanted)
    order = Symbol[]
    visit(m) = begin
        haskey(_APP_SRC_MODULES, m) || error("unknown flexOPT module :$m; known: $(sort(collect(keys(_APP_SRC_MODULES))))")
        m in order && return
        foreach(visit, _APP_SRC_MODULES[m][2])
        push!(order, m)
    end
    foreach(visit, wanted)
    return order
end

# GPU backend (optional): batchGPU.jl must be included before flexOPT.jl
if get(ENV, "FLEXOPT_GPU", "none") == "metal" && !isdefined(Main, :backend)
    Core.eval(Main, :(using Metal))
    Base.include(Main, joinpath(flexopt_root, "src", "batchFiles", "batchGPU.jl"))
end

for _m in _appLoadOrder(@isdefined(FLEXOPT_MODULES) ? FLEXOPT_MODULES : [:flexOPT])
    _modname = get(_APP_MODULE_NAME, _m, _m)
    isdefined(Main, _modname) && continue           # already loaded in this session
    Base.include(Main, joinpath(flexopt_root, "src", _APP_SRC_MODULES[_m][1]))
    Core.eval(Main, Expr(:using, Expr(:., :., _modname)))
    _m == :planet1D && Main.planet1D.configure_input!()
end

# ---------------------------------------------------------------- utilities
"""
    appArgs(args; defaults...) -> NamedTuple

Minimal parser shared by the application run.jl files.
`key=value` pairs override `defaults`; bare words are returned in `positional`.
Values are parsed as Int, Float64, comma lists, or kept as strings.

    julia --project=. application/Poisson_1D/run.jl convergence schemes=convFD3,OPT3 case=sameλ
"""
function appArgs(args; defaults...)
    values = Dict{Symbol,Any}(defaults)
    positional = String[]
    for a in args
        a = replace(a, r"^--" => "")
        if occursin('=', a)
            k, v = split(a, '='; limit=2)
            values[Symbol(k)] = _appParseValue(v)
        else
            push!(positional, a)
        end
    end
    return (; positional, values...)
end
function _appParseValue(v)
    occursin(',', v) && return [_appParseValue(s) for s in split(v, ',') if !isempty(s)]
    i = tryparse(Int, v);     i === nothing || return i
    f = tryparse(Float64, v); f === nothing || return f
    v in ("true", "false") && return v == "true"
    return String(v)
end

"application/<Name>/outputs/<sub>/, created if needed (git-ignored)."
appOutputDir(appdir::AbstractString, sub::AbstractString="") =
    mkpath(joinpath(appdir, "outputs", sub))

"Save a Makie figure as PNG in `dir` and print the path."
function appSave(fig, dir::AbstractString, name::AbstractString)
    path = joinpath(dir, replace(name, r"[^A-Za-z0-9_.=+-]+" => "_") * ".png")
    Main.save(path, fig)
    println("figure -> ", path)
    return path
end
