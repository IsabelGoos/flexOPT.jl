# dev/spherical1Delastic/stages/step2_operators.jl
# Ported from notebooks/SphericalTI_SH_oct2026.ipynb, Step 2 (cells sh-step2-03 … sh-step2-13).
# Recipes: cached by shRecipe itself (myProduceOrLoad, data/semiSymbolic).
# Convergence tables, stencil exactness and zone comparisons: cached here.
# (On a cache hit the tables printed by the helper functions are not
#  re-printed; the figures and ratio tables below are always redone.)

# ---- 2a: configurations (stages/configs.jl) and recipes
shRecipes = Dict{String,Any}()
for c in shConfigs
    @time shRecipes[c.name] = shRecipe(c)
    r = shRecipes[c.name]["recette"]
    println(c.name, ": form=", r.form, "  varM rows=", size(r.lhs.varM, 1),
            "  stencil nodes=", length(r.nodes[1]), "  Ajiννᶜ size=", size(r.lhs.Ajiννᶜ))
end

# ---- 2b-i: a single grid first (homogeneous isotropic shell) -- cheap, not cached
homShell = layeredSHModel(; R=6371.0, radii=[3480.0, 6371.0], ρ=[3.3], Vsv=[4.5], Vsh=[4.5])
@show homShell.nzone homShell.solid_or_fluid outerSolidZones(homShell)
quick = shConvergenceHomogeneous(shConfigs; l=10, nIntervals=[100])

# ---- 2b-ii: convergence, uniform h (default grid): cached
convUParams = Dict{String,Any}(
    "stage" => "step2_convergence",
    "codeVersion" => CODE_VERSION["step2_convergence"],
    "equationRevision" => SH_EQUATION_REVISION,
    "configs" => shConfigs, "l" => 10, "grid" => :uniform,
    "nIntervals" => [50, 100, 200, 400, 800],
)
devRedo!("convergence", convUParams, SPH_CACHE_DIR, "step2_convU")
convU = myProduceOrLoad(convUParams, SPH_CACHE_DIR, "step2_convU") do cfg
    Dict{String,Any}("rows" => shConvergenceHomogeneous(cfg["configs"]; l=cfg["l"],
        nIntervals=cfg["nIntervals"]))
end["rows"]

# ---- 2b-iii: convergence, planet1D DSM grid (re): cached
convDParams = Dict{String,Any}(
    "stage" => "step2_convergence",
    "codeVersion" => CODE_VERSION["step2_convergence"],
    "equationRevision" => SH_EQUATION_REVISION,
    "configs" => shConfigs, "l" => 10, "grid" => :dsm,
    "reList" => [1e-1, 3e-2, 1e-2, 3e-3, 1e-3],
)
devRedo!("convergence", convDParams, SPH_CACHE_DIR, "step2_convD")
convD = myProduceOrLoad(convDParams, SPH_CACHE_DIR, "step2_convD") do cfg
    Dict{String,Any}("rows" => shConvergenceHomogeneous(cfg["configs"]; l=cfg["l"],
        grid=cfg["grid"], reList=cfg["reList"]))
end["rows"]

plotSHConvergence(convU; title="(uniform h)")
plotSHConvergence(convD; title="(planet1D re-grid)")
println("uniform h:");        shErrorRatios(convU)
println("planet1D re-grid:"); shErrorRatios(convD)

# ---- 2b-ii': is the small-h floor an inexactness of the Δ=1 recipes?: cached
exactnessParams = Dict{String,Any}(
    "stage" => "step2_exactness",
    "codeVersion" => CODE_VERSION["step2_exactness"],
    "equationRevision" => SH_EQUATION_REVISION,
    "configs" => shConfigs,
)
devRedo!("exactness", exactnessParams, SPH_CACHE_DIR, "step2_exactness")
exactness = myProduceOrLoad(exactnessParams, SPH_CACHE_DIR, "step2_exactness") do cfg
    Dict{String,Any}("result" => shStencilExactness(cfg["configs"]))
end["result"]
@show exactness

# ---- 2b-iv: real models, zone by zone: cached
ΩT100 = 2π * 6371.0 / 100.0
twoZoneSpec = (R=6371.0, radii=[3480.0, 5701.0, 6371.0], ρ=[5.0, 3.6], Vsv=[6.5, 4.6], Vsh=[6.5, 4.9])
twoZoneTI = layeredSHModel(; twoZoneSpec...)
gTwo = shRadialGrid(twoZoneTI; mode=:uniform, h=40.0)
for s in gTwo.segments
    @printf("zone %d: %7.1f – %7.1f km, %4d intervals, h = %.3f km\n", s.zone, s.r[1], s.r[end], length(s.r) - 1, s.h_km)
end
cmpTwoParams = Dict{String,Any}(
    "stage" => "step2_zones",
    "codeVersion" => CODE_VERSION["step2_zones"],
    "equationRevision" => SH_EQUATION_REVISION,
    "configs" => shConfigs, "model" => twoZoneSpec, "l" => 10, "Ω" => ΩT100, "h" => 40.0,
)
devRedo!("zones", cmpTwoParams, SPH_CACHE_DIR, "step2_zones_twoZoneTI")
cmpTwo = myProduceOrLoad(cmpTwoParams, SPH_CACHE_DIR, "step2_zones_twoZoneTI") do cfg
    Dict{String,Any}("rows" => shZoneComparison(layeredSHModel(; cfg["model"]...), cfg["configs"];
        l=cfg["l"], Ω=cfg["Ω"], h=cfg["h"]))
end["rows"]

prem = planet1D.getSet1Dmodel!(:PREM)
gPrem = shRadialGrid(prem; mode=:uniform, h=25.0)
for s in gPrem.segments
    @printf("zone %2d: %7.1f – %7.1f km, %4d intervals, h = %7.3f km\n", s.zone, s.r[1], s.r[end], length(s.r) - 1, s.h_km)
end
cmpPremParams = Dict{String,Any}(
    "stage" => "step2_zones",
    "codeVersion" => CODE_VERSION["step2_zones"],
    "equationRevision" => SH_EQUATION_REVISION,
    "configs" => shConfigs, "model" => "PREM", "l" => 10, "Ω" => ΩT100, "h" => 25.0,
)
devRedo!("zones", cmpPremParams, SPH_CACHE_DIR, "step2_zones_PREM")
cmpPrem = myProduceOrLoad(cmpPremParams, SPH_CACHE_DIR, "step2_zones_PREM") do cfg
    Dict{String,Any}("rows" => shZoneComparison(planet1D.getSet1Dmodel!(Symbol(cfg["model"])),
        cfg["configs"]; l=cfg["l"], Ω=cfg["Ω"], h=cfg["h"]))
end["rows"]
