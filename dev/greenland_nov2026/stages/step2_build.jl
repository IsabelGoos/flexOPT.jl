# dev/greenland_nov2026/stages/step2_build.jl
# Runs ONLY on a cache miss of step2_opt.jl.
#   1. OPT3 recipe of famousEquations("3DsismoTimeIsoHeteroForce")
#      (ρ ∂²u/∂t² - ∇·(C:∇u) = F, isotropic λ, μ, body force F), in scaled
#      coordinates x' = x/Δx, t' = t/Δt (Δ = (1,1,1,1)): the same recipe serves
#      every dx/dt and every model. Cached in data/semiSymbolic.
#   2. numericalOperatorConstruction + prepareLinearSystem on the Greenland
#      model, Cerjan layers on the four sides and the bottom (none on top).
#   3. free surface: apply_free_surface_3d (helpers), since prepareLinearSystem
#      applies its pinned-void / Dietrich closure in 2-D only.

# ---- 1. recipe (cached; the 4-D symbolic construction is the slow part)
recipeBackend = recipeOnGPU ? backend : KernelAbstractions.CPU()
optParameters = Dict{String,Any}(
    "famousEquationType" => "3DsismoTimeIsoHeteroForce",
    "Δ" => (1.0, 1.0, 1.0, 1.0),
    "orderBtime" => 1, "orderBspace" => 1,
    "pointsInSpace" => 3, "pointsInTime" => 3,
    "supplementaryOrder" => optSupplementaryOrder,
    "taylorInverseMode" => :hierarchical_constrained,
    "fieldItpl" => (ptsSpace=1, ptsTime=1, offsetSpace=1.0,
        offsetTime=1, YorderBspace=-1, YorderBtime=-1),
    "materItpl" => (ptsSpace=1, ptsTime=1, offsetSpace=1.0,
        offsetTime=1, YorderBspace=-1, YorderBtime=-1),
    "recipe_backend" => recipeBackend,
)
recipeCacheParameters = Dict{String,Any}(k => v for (k, v) in optParameters if k != "recipe_backend")
recipeCacheParameters["recipe_cache_version"] = 3                  # hierarchical Taylor inverse
recipeCacheParameters["recipe_backend_name"] = string(nameof(typeof(recipeBackend)))
recipeFile = devCacheFile(recipeCacheParameters, "semiSymbolic", "elastic3DForce_OPT3")
println("OPT3 3-D recipe: ", recipeFile, isfile(recipeFile) ? "  (hit)" : "  (miss -> makeOPTsemiSymbolic, slow)")
timedRecipe = @timed myProduceOrLoad(recipeCacheParameters, "semiSymbolic", "elastic3DForce_OPT3") do cfg
    runtime = Dict{String,Any}(cfg)
    for key in ("hash_id", "recipe_cache_version", "recipe_backend_name")
        pop!(runtime, key, nothing)
    end
    runtime["recipe_backend"] = recipeBackend
    makeOPTsemiSymbolic(runtime)
end
optRecipe = timedRecipe.value

# ---- 2. operators. Models in SI, scaled exactly as the Kirishima benchmark:
# ρ' = ρ/ρ0, λ' = λ/ρ0 (Δt/Δx)², μ' = μ/ρ0 (Δt/Δx)²
ρSI = ρModel .* 1e3
vpSI = VpModel .* 1e3
vsSI = VsModel .* 1e3
rho0 = median(ρSI[materialMask])
scale2 = (dtOPT / dx)^2
rhoOPT = max.(ρSI ./ rho0, 1e-3)
muOPT = (ρSI .* vsSI .^ 2 ./ rho0) .* scale2
lambdaOPT = (ρSI .* (vpSI .^ 2 .- 2 .* vsSI .^ 2) ./ rho0) .* scale2
if freeSurfaceClosure === :dietrich3d
    muOPT[.!materialMask] .= 0.0           # air rows are pinned anyway
    lambdaOPT[.!materialMask] .= 0.0
end
@assert all(isfinite, lambdaOPT) && all(isfinite, muOPT)
@assert minimum(lambdaOPT[materialMask] .+ (2 / 3) .* muOPT[materialMask]) > 0 "negative bulk modulus"

cerjanSpec = CerjanBoundarySpec((cerjanCells, cerjanCells, cerjanCells),
    (cerjanCells, cerjanCells, 0); damping=cerjanDamping)
modelsOPT = [rhoOPT, lambdaOPT, muOPT]
pointsOPT = getModelPoints(modelsOPT[1], 3,
    optRecipe["recette"].numbersOfTheSystem.numbersOfTheSystemL.timeMarching)
familyOPT = (models=modelsOPT, modelPoints=pointsOPT, Δ=(1.0, 1.0, 1.0, 1.0),
    modelName="greenland_$(greenlandCase)_$(modelHash)")
timedAssembly = @timed numericalOperatorConstruction(Dict{String,Any}(
    "optRec" => optRecipe, "modelFam" => familyOPT,
    "absorbingBoundaries" => nothing,
    "maskedRegionInSpace" => nothing,
    "boundaryConditions" => BoundaryConditionSet(free_surface=nothing, cerjan=cerjanSpec,
        material_mask=nothing),
    "representation" => "matrixfree",
))["numOperators"]
timedPreparation = @timed prepareLinearSystem(timedAssembly.value)
preparedVolume = timedPreparation.value
optPadding = cerjan_padding(cerjanSpec)                    # [lower; upper] × (x, y, z)
@assert preparedVolume.spaceShape == (Nx, Ny, Nz) .+ Tuple(optPadding[1, :]) .+ Tuple(optPadding[2, :])
@assert preparedVolume.NField == 3 && preparedVolume.NForceField == 3
@assert preparedVolume.timePointsUsedForOneStep == 3

# ---- 3. free surface
timedClosure = @timed (freeSurfaceClosure === :dietrich3d ?
    apply_free_surface_3d(preparedVolume, materialMask;
        prepared_with_matrices=flexOPT._prepared_with_matrices) :
    preparedVolume)
preparedOPT = timedClosure.value
operatorTimings = (recipe=timedRecipe.time, assembly=timedAssembly.time,
    preparation=timedPreparation.time, freeSurface=timedClosure.time)
preparedVolume = nothing; timedAssembly = nothing; timedPreparation = nothing; timedClosure = nothing
GC.gc()
@show operatorTimings
@printf("OPT3 3-D: whole grid %s, %d unknowns, nnz(A) = %.1f M, nnz(L) = %.1f M, nnz(R) = %.1f M\n",
    string(preparedOPT.spaceShape), size(preparedOPT.A_unknown, 1),
    nnz(preparedOPT.A_unknown) / 1e6, nnz(preparedOPT.L_known) / 1e6, nnz(preparedOPT.R_force) / 1e6)
