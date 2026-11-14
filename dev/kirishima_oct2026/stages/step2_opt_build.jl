# dev/kirishima_oct2026/stages/step2_opt_build.jl
# Ported from notebooks/Elastic2DtopoBenchmark_oct2026.ipynb, cells step2-*-opt-recipes + step2-*-opt-operator.
# Runs ONLY on a cache miss of step2_opt.jl. Recipes are themselves cached
# (myProduceOrLoad, data/semiSymbolic); operators are rebuilt (too large to cache).
dxOPT = dx / optSpatialRefinement
xOPTPhysical = collect(first(x):dxOPT:last(x))
zOPTPhysical = collect(first(z):dxOPT:last(z))
solidOPT = applyFreeSurface ?
    repeat(reshape(zOPTPhysical .<= 0.0, 1, :),
        length(xOPTPhysical), 1) :
    trues(length(xOPTPhysical), length(zOPTPhysical))
dtOPT = 0.20 * dxOPT / (sqrt(2) * vpMaximum * temporalRefinement)
@assert optSpatialRefinement == 1 "OPT refinement != 1 needs field resampling, not implemented for the Kirishima case yet"
vpFieldOPT, vsFieldOPT, rhoFieldOPT = vpField, vsField, rhoField
# Stable reference used by the previously working OPT benchmark. This is
# an exact change of coordinates x'=x/Δx, z'=z/Δz, t'=t/Δt, not a
# change of the physical model. The direct physical-recipe path remains
# available for conditioning experiments, but is not the benchmark default.
optFormulation = :scaled_coordinates # or :physical_recipe_experimental
if optFormulation === :scaled_coordinates
    recipeSpacing = (1.0, 1.0, 1.0)
    rhoOPT = rhoFieldOPT ./ rho0
    muOPT = (rhoFieldOPT .* vsFieldOPT.^2 ./ rho0) .* (dtOPT/dxOPT)^2
    lambdaOPT = (rhoFieldOPT .* (vpFieldOPT.^2 .- 2vsFieldOPT.^2) ./ rho0) .*
        (dtOPT/dxOPT)^2
elseif optFormulation === :physical_recipe_experimental
    recipeSpacing = (dxOPT, dxOPT, dtOPT)
    rhoOPT = rhoFieldOPT
    muOPT = rhoFieldOPT .* vsFieldOPT.^2
    lambdaOPT = rhoFieldOPT .* (vpFieldOPT.^2 .- 2vsFieldOPT.^2)
else
    error("unknown optFormulation=$optFormulation")
end
muOPT[.!solidOPT] .= 0.0
lambdaOPT[.!solidOPT] .= 0.0
optParameters = Dict{String,Any}(
    "famousEquationType" => "2DsismoTimeIsoHeteroSingleForce",
    "Δ" => recipeSpacing,
    "orderBtime" => 1, "orderBspace" => 1,
    "pointsInSpace" => 3, "pointsInTime" => 3,
    "supplementaryOrder" => 2,
    # Preserve all total-degree <= 2 Taylor conditions exactly; use
    # supplementaryOrder only to optimise the remaining null space.
    "taylorInverseMode" => :hierarchical_constrained,
    "fieldItpl" => (ptsSpace=1, ptsTime=1, offsetSpace=1.0,
        offsetTime=1, YorderBspace=-1, YorderBtime=-1),
    "materItpl" => (ptsSpace=1, ptsTime=1, offsetSpace=1.0,
        offsetTime=1, YorderBspace=-1, YorderBtime=-1),
    # The expensive semi-symbolic coefficient construction runs on the
    # GPU selected in the bootstrap cell (Metal here).
    "recipe_backend" => backend,
)
# Cache complete semi-symbolic recipes, not the very large prepared
# operators. The runtime Metal/CPU object is deliberately excluded from
# the cache key and is injected only when a recipe must be produced.
recipeCacheDirectory = "semiSymbolic"
recipeCacheVersion = 3 # hierarchical Taylor inverse; reject legacy recipes
function cachedOPTRecipe(parameters, prefix; existing=nothing)
    cacheParameters = Dict{String,Any}(
        key => value for (key, value) in parameters
        if key != "recipe_backend"
    )
    cacheParameters["recipe_cache_version"] = recipeCacheVersion
    function produceRecipe(config)
        runtimeParameters = Dict{String,Any}(config)
        pop!(runtimeParameters, "hash_id", nothing)
        pop!(runtimeParameters, "recipe_cache_version", nothing)
        runtimeParameters["recipe_backend"] = backend
        isnothing(existing) ?
            makeOPTsemiSymbolic(runtimeParameters) : existing
    end
    myProduceOrLoad(produceRecipe, cacheParameters,
        recipeCacheDirectory, prefix)
end
existingOPT3Recipe = @isdefined(optRecipe) ? optRecipe : nothing
optRecipe = cachedOPTRecipe(optParameters, "elastic2D_OPT3";
    existing=existingOPT3Recipe)
# P-C predictor: same three-point support and physical model, but only
# the supplementary Taylor optimisation is removed. The full OPT3
# operator remains the residual operator used by every correction.
optPredictorParameters = copy(optParameters)
optPredictorParameters["supplementaryOrder"] = 0
optPredictorRecipe = optTimeMarchingScheme === :predictor_corrector &&
    optPredictorMode === :prepared ?
    cachedOPTRecipe(optPredictorParameters, "elastic2D_OPT3_predictor_s0") : nothing

# Four/five-point spatial recipes; time remains a three-point stencil.
opt4Parameters = copy(optParameters)
opt4Parameters["pointsInSpace"] = 4
opt4Parameters["fieldItpl"] = merge(optParameters["fieldItpl"],
    (offsetSpace=1.5,))
opt4Parameters["materItpl"] = merge(optParameters["materItpl"],
    (offsetSpace=1.5,))
opt4Recipe = buildHigherOrderRecipes ? cachedOPTRecipe(
    opt4Parameters, "elastic2D_OPT4"; existing=nothing) : nothing

# Five-point OPT recipe: the time stencil remains three points so that
# only the spatial approximation changes relative to OPT3.
opt5Parameters = copy(optParameters)
# Validated OPT5 family used by seismo1Dbenchmark: five trial points
# retain the linear B-spline test order. orderBspace=3 is not OPT5 and
# generated a nearly singular first-step solve in this elastic problem.
opt5Parameters["orderBspace"] = 1
opt5Parameters["pointsInSpace"] = 5
opt5Parameters["fieldItpl"] = merge(optParameters["fieldItpl"],
    (offsetSpace=2.0,))
opt5Parameters["materItpl"] = merge(optParameters["materItpl"],
    (offsetSpace=2.0,))
# Keep OPT5 available, but do not even load/build it in the quick run.
# When enabled, never migrate the former invalid orderBspace=3 object.
opt5Recipe = buildHigherOrderRecipes ? cachedOPTRecipe(
    opt5Parameters, "elastic2D_OPT5"; existing=nothing) : nothing

# Collocation/delta-test recipe corresponding to the conventional FD3
# definition used in seismo1Dbenchmark.ipynb. Keep it beside OPT3: it is
# a coefficient reference, not an alias for the explicit flux solver.
conventionalFDParameters = copy(optParameters)
conventionalFDParameters["orderBspace"] = -1
conventionalFDParameters["orderBtime"] = -1
conventionalFDParameters["supplementaryOrder"] = 0
existingConvFD3Recipe = @isdefined(conventionalFDRecipe) ?
    conventionalFDRecipe : nothing
conventionalFDRecipe = cachedOPTRecipe(
    conventionalFDParameters, "elastic2D_convFD3";
    existing=existingConvFD3Recipe)
recipeComparison = (
    conventionalFD3=(orderBspace=-1, orderBtime=-1,
        supplementaryOrder=0, pointsInSpace=3, pointsInTime=3),
    OPT3=(orderBspace=optParameters["orderBspace"],
        orderBtime=optParameters["orderBtime"],
        supplementaryOrder=optParameters["supplementaryOrder"],
        pointsInSpace=optParameters["pointsInSpace"],
        pointsInTime=optParameters["pointsInTime"]),
    OPT4=(orderBspace=opt4Parameters["orderBspace"],
        orderBtime=opt4Parameters["orderBtime"],
        supplementaryOrder=opt4Parameters["supplementaryOrder"],
        pointsInSpace=opt4Parameters["pointsInSpace"],
        pointsInTime=opt4Parameters["pointsInTime"]),
    OPT5=(orderBspace=opt5Parameters["orderBspace"],
        orderBtime=opt5Parameters["orderBtime"],
        supplementaryOrder=opt5Parameters["supplementaryOrder"],
        pointsInSpace=opt5Parameters["pointsInSpace"],
        pointsInTime=opt5Parameters["pointsInTime"]),
)
# Public handles for inspecting the complete symbolic recipes in IJulia,
# e.g. `optRecipes.OPT5["recette"]`.
optRecipes = (convFD3=conventionalFDRecipe, OPT3=optRecipe,
    OPT4=opt4Recipe, OPT5=opt5Recipe)
@show recipeComparison
if runOPT4 || runOPT5
    coefficientGatePath = joinpath(flexopt_root, "data",
        "elastic_lhs_coefficient_gate.jld2")
    isfile(coefficientGatePath) || error(
        "Run ElasticLHSOperatorAudit.ipynb before OPT4/OPT5 propagation")
    coefficientGate = load(coefficientGatePath)
    @assert all(result.passed for result in
        values(coefficientGate["coefficient_health"]))
end


optCerjan = CerjanBoundarySpec(
    (cerjan.lower[1] * optSpatialRefinement,
        cerjan.lower[2] * optSpatialRefinement),
    (cerjan.upper[1] * optSpatialRefinement, 0);
    damping=cerjan.damping / optSpatialRefinement^2,
)
activeOptCerjan = applyCerjan ? optCerjan : nothing
bcOPT = applyFreeSurface ?
    boundary_geometry(solidOPT, (dxOPT, dxOPT);
        free_surface_mode=:pinned_void, cerjan=activeOptCerjan) :
    BoundaryConditionSet(free_surface=nothing, cerjan=activeOptCerjan,
        material_mask=BitArray(solidOPT), free_surface_mode=:pinned_void)
modelsOPT = [rhoOPT, lambdaOPT, muOPT]
pointsOPT = getModelPoints(modelsOPT[1], 3,
    optRecipe["recette"].numbersOfTheSystem.numbersOfTheSystemL.timeMarching)
familyOPT = (models=modelsOPT, modelPoints=pointsOPT,
    Δ=recipeSpacing, modelName="$(caseLabel)_OPT3_$(optFormulation)")
operatorTimings = Dict{Symbol,Float64}()
timedOPT3Assembly = @timed numericalOperatorConstruction(Dict{String,Any}(
    "optRec" => optRecipe, "modelFam" => familyOPT,
    "absorbingBoundaries" => nothing,
    "maskedRegionInSpace" => nothing,
    "boundaryConditions" => bcOPT,
    "representation" => "matrixfree",
))["numOperators"]
numericalVolume = timedOPT3Assembly.value
operatorTimings[:OPT3_assembly] = timedOPT3Assembly.time
timedOPT3Preparation = @timed prepareLinearSystem(numericalVolume;
    free_surface_spacing=recipeSpacing[1:2])
preparedVolume = timedOPT3Preparation.value
operatorTimings[:OPT3_preparation] = timedOPT3Preparation.time
preparedOPT4Volume = if runOPT4
    pointsOPT4 = getModelPoints(modelsOPT[1], 4,
        opt4Recipe["recette"].numbersOfTheSystem.numbersOfTheSystemL.timeMarching)
    familyOPT4 = (models=modelsOPT, modelPoints=pointsOPT4,
        Δ=recipeSpacing, modelName="$(caseLabel)_OPT4_$(optFormulation)")
    numericalOPT4 = numericalOperatorConstruction(Dict{String,Any}(
        "optRec" => opt4Recipe, "modelFam" => familyOPT4,
        "absorbingBoundaries" => nothing, "maskedRegionInSpace" => nothing,
        "boundaryConditions" => bcOPT, "representation" => "matrixfree",
    ))["numOperators"]
    prepareLinearSystem(numericalOPT4; free_surface_spacing=recipeSpacing[1:2])
else
    nothing
end
preparedOPT5Volume = if runOPT5
    pointsOPT5 = getModelPoints(modelsOPT[1], 5,
        opt5Recipe["recette"].numbersOfTheSystem.numbersOfTheSystemL.timeMarching)
    familyOPT5 = (models=modelsOPT, modelPoints=pointsOPT5,
        Δ=recipeSpacing, modelName="$(caseLabel)_OPT5_$(optFormulation)")
    numericalOPT5 = numericalOperatorConstruction(Dict{String,Any}(
        "optRec" => opt5Recipe, "modelFam" => familyOPT5,
        "absorbingBoundaries" => nothing,
        "maskedRegionInSpace" => nothing,
        "boundaryConditions" => bcOPT,
        "representation" => "matrixfree",
    ))["numOperators"]
    prepareLinearSystem(numericalOPT5;
        free_surface_spacing=recipeSpacing[1:2])
else
    nothing
end
timedConvFDAssembly = @timed numericalOperatorConstruction(Dict{String,Any}(
    "optRec" => conventionalFDRecipe, "modelFam" => familyOPT,
    "absorbingBoundaries" => nothing,
    "maskedRegionInSpace" => nothing,
    "boundaryConditions" => bcOPT,
    "representation" => "matrixfree",
))["numOperators"]
numericalConvFD = timedConvFDAssembly.value
operatorTimings[:convFD3_assembly] = timedConvFDAssembly.time
timedConvFDPreparation = @timed prepareLinearSystem(numericalConvFD;
    free_surface_spacing=recipeSpacing[1:2])
preparedConvFDVolume = timedConvFDPreparation.value
operatorTimings[:convFD3_preparation] = timedConvFDPreparation.time

# Optional traction recipe σn=0, assembled additively on surface rows.
if applyFreeSurface
boundaryParameters = copy(optParameters)
boundaryParameters["famousEquationType"] = "elasticTractionFree2D"
boundaryRecipe = cachedOPTRecipe(
    boundaryParameters, "elasticTractionFree2D_OPT3")
normalX = zeros(size(solidOPT)); normalZ = zeros(size(solidOPT))
for (point, normal) in zip(bcOPT.free_surface.points, bcOPT.free_surface.normals)
    normalX[point], normalZ[point] = normal
end
boundaryModels = [lambdaOPT, muOPT, normalX, normalZ]
boundaryPoints = getModelPoints(boundaryModels[1], 3,
    boundaryRecipe["recette"].numbersOfTheSystem.numbersOfTheSystemL.timeMarching)
boundaryFamily = (models=boundaryModels, modelPoints=boundaryPoints,
    Δ=recipeSpacing, modelName="$(caseLabel)_free_surface")
numericalBoundary = numericalOperatorConstruction(Dict{String,Any}(
    "optRec" => boundaryRecipe, "modelFam" => boundaryFamily,
    "absorbingBoundaries" => (isnothing(bcOPT.cerjan) ? nothing : cerjan_padding(bcOPT.cerjan)),
    "maskedRegionInSpace" => bcOPT.free_surface.points,
    "representation" => "matrixfree",
))["numOperators"]
preparedBoundary = prepareLinearSystem(numericalBoundary)
surfaceWhole = numericalVolume.numericalOperators.left.geometry.freeSurfaceBoundary.points
    preparedOPT = overlapBoundaryLinearSystem(
    preparedVolume, preparedBoundary, surfaceWhole;
    mode=:additive, boundary_weight=1.0,
    )
    preparedConvFD = overlapBoundaryLinearSystem(
        preparedConvFDVolume, preparedBoundary, surfaceWhole;
        mode=:additive, boundary_weight=1.0,
    )
    preparedOPT4 = runOPT4 ? overlapBoundaryLinearSystem(
        preparedOPT4Volume, preparedBoundary, surfaceWhole;
        mode=:additive, boundary_weight=1.0) : nothing
    preparedOPT5 = runOPT5 ? overlapBoundaryLinearSystem(
        preparedOPT5Volume, preparedBoundary, surfaceWhole;
        mode=:additive, boundary_weight=1.0) : nothing
else
    preparedOPT = preparedVolume
    preparedConvFD = preparedConvFDVolume
    preparedOPT4 = runOPT4 ? preparedOPT4Volume : nothing
    preparedOPT5 = runOPT5 ? preparedOPT5Volume : nothing
end
@assert preparedOPT.NForceField == 2
@assert preparedConvFD.NForceField == 2
runOPT4 && @assert preparedOPT4.NForceField == 2
runOPT5 && @assert preparedOPT5.NForceField == 2
optOperatorUnitConvention = optFormulation
boundaryOverlapRows = hasproperty(preparedOPT, :boundary_overlap_rows) ?
    preparedOPT.boundary_overlap_rows : Int[]
@show quickRun optFormulation optSpatialRefinement dxOPT dtOPT recipeSpacing
@show operatorTimings
@show preparedOPT.spaceShape preparedConvFD.spaceShape runOPT5 length(boundaryOverlapRows)
runOPT5 && @show preparedOPT5.spaceShape
