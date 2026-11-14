# dev/spherical1DelasticKernel/stages/env.jl
# Same session environment as dev/spherical1Delastic (src/ modules, Metal,
# lib/*.jl, shConfigs), then the kernel lib (re-read with --reload).
SPH_DIR = joinpath(flexopt_root, "dev", "spherical1Delastic")
include(joinpath(SPH_DIR, "stages", "env.jl"))
include(joinpath(SPH_DIR, "stages", "configs.jl"))          # shConfigs (convFD3, OPT3, OPT4, OPT5)
using Interpolations
kerLoadLib = !@isdefined(KERNEL_ENV) || DEV.reload
KER_LIB = joinpath(flexopt_root, "dev", "spherical1DelasticKernel", "lib")
kerLoadLib && include(joinpath(KER_LIB, "sgtHarmonics1D.jl"))
kerLoadLib && include(joinpath(KER_LIB, "sgtDisk1D.jl"))
kerLoadLib && include(joinpath(KER_LIB, "sgtGeometry1D.jl"))
kerLoadLib && include(joinpath(KER_LIB, "sgtInterp1D.jl"))
kerLoadLib && include(joinpath(KER_LIB, "sgtKernels1D.jl"))
kerLoadLib && include(joinpath(KER_LIB, "sgtChecks1D.jl"))
kerLoadLib && include(joinpath(KER_LIB, "sgtZC1D.jl"))
KERNEL_ENV = true
kerRecipes = @isdefined(kerRecipes) ? kerRecipes : Dict{Any,Any}()

"PREM (or any planet1D model) P-SV + SH systems of `scheme` at uniform h (km), cached."
function kerSystems(modelName::AbstractString, scheme::AbstractString, h::Real; stage="systems")
    params = Dict{String,Any}("stage" => "k_systems", "codeVersion" => CODE_VERSION_K["k_systems"],
        "model" => modelName, "scheme" => scheme, "h" => Float64(h), "configs" => shConfigs,
        "equationRevisions" => EQUATION_REVISIONS)
    devRedo!(stage, params, KER_CACHE_DIR, "k_systems")
    return myProduceOrLoad(params, KER_CACHE_DIR, "k_systems") do cfg
        model = planet1D.getSet1Dmodel!(Symbol(cfg["model"]))
        gP = elasticRadialGrid(model; h=cfg["h"], region=:all)
        gS = elasticRadialGrid(model; h=cfg["h"], region=:outerSolid)
        kw = (scheme=cfg["scheme"], configs=cfg["configs"], recipes=kerRecipes, label=cfg["scheme"])
        Dict("psv" => elasticGlobalSystem(:PSV, model, gP; kw...), "sh" => elasticGlobalSystem(:SH, model, gS; kw...))
    end
end
