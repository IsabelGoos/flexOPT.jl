# dev/greenland_nov2026/stages/step1_model.jl
# Material model on the box, cached with myProduceOrLoad:
#   1. planet1D (PREM, no ocean layer) as a function of depth below the local
#      rock surface, averaged over the rock part of each cell (point sampling
#      by absolute radius, as getParamsAndTopo does, makes the 24 km PREM crust
#      appear or vanish cell by cell along the curved surface);
#   2. ETOPO 2022 ice surface + bedrock (cached separately) -> where the ice,
#      the ocean and the air are;
#   3. greenland_material (helpers): volume fractions of air / ice / water /
#      rock in every cell, Voigt average of the non-air part; cells with
#      >= 50 % air are air (pinned by the free-surface closure in step 2).

# ---- 2. ETOPO grids (cached on their own: the slowest download)
etopoParams = Dict{String,Any}(
    "stage" => "etopo2022_60s",
    "codeVersion" => CODE_VERSION["etopo"],
    "region" => etopo_region(box.allGridsInGeoPoints; margin=1.0),
    "strideArcmin" => max(1, min(4, floor(Int, dx / 10e3))),   # 1-4 arcmin, finer than dx/5
    "localDir" => ETOPO_LOCAL_DIR,
)
devRedo!("etopo", etopoParams, GREENLAND_CACHE_DIR, "etopo")
println("ETOPO cache: ", devCacheFile(etopoParams, GREENLAND_CACHE_DIR, "etopo"))
etopo = myProduceOrLoad(fetch_etopo, etopoParams, GREENLAND_CACHE_DIR, "etopo")
surfItp, bedItp = etopo_interpolators(etopo)

# ---- 1 + 3. model
phaseTuple(p) = (ρ=p.ρ, vp=p.vp, vs=p.vs)   # plain values in cache keys / files
modelParams = Dict{String,Any}(
    "stage" => "model",
    "codeVersion" => CODE_VERSION["model"],
    "boxCentre" => BOX_CENTRE, "axis" => BOX_AXIS_AZIMUTH,
    "halfWidth" => boxHalfWidth, "dx" => dx, "altMin" => altMin, "altMax" => altMax,
    "planet1D" => string(planet1D.DEFAULT_1D_MODEL[]), "rock" => "depth-below-bed average",
    "ice" => phaseTuple(ICE), "water" => phaseTuple(WATER), "air" => phaseTuple(AIR),
    "airThreshold" => AIR_THRESHOLD, "iceThreshold" => ICE_THRESHOLD,
    "etopoHash" => commonBatchs.hash_parameters(etopoParams),
)
devRedo!("model", modelParams, GREENLAND_CACHE_DIR, "model")
modelCacheFile = devCacheFile(modelParams, GREENLAND_CACHE_DIR, "model")
println("model cache: ", modelCacheFile, isfile(modelCacheFile) ? "  (hit)" : "  (miss -> computing)")
modelHash = commonBatchs.hash_parameters(modelParams)

modelResult = myProduceOrLoad(modelParams, GREENLAND_CACHE_DIR, "model") do cfg
    # planet1D (PREM) as a function of depth below the local rock surface,
    # 0-1200 km every 0.5 km, discontinuities kept ("below" convention)
    R0km = planet1D.my1DDSMmodel.averagedPlanetRadiusInKilometer
    radiiKm = collect(R0km .- (0.0:0.5:1200.0))
    newRadii, prof = planet1D.compute1DseismicParamtersFromPolynomialCoefficientsWithGivenRadiiArray(
        planet1D.my1DDSMmodel, radiiKm, "below")
    rockProfile = depth_profile((R0km .- newRadii) .* 1e3, prof.ρ, prof.Vpv, prof.Vsv)
    rockTop = average_rock(rockProfile, 0.0, 1000.0)
    @show rockTop average_rock(rockProfile, 0.0, 50e3)
    m = greenland_material(rockProfile, box.allGridsInGeoPoints, dx, surfItp, bedItp;
        ice=ICE, water=WATER, air=AIR, airThreshold=AIR_THRESHOLD, iceThreshold=ICE_THRESHOLD)
    Dict{String,Any}(
        "ρ" => m.ρ, "Vp" => m.Vp, "Vs" => m.Vs, "material" => Array(m.material),
        "fAir" => Float32.(m.fAir), "fIce" => Float32.(m.fIce),
        "fWater" => Float32.(m.fWater), "fRock" => Float32.(m.fRock),
        "iceThickness" => m.iceThickness, "waterDepth" => m.waterDepth,
        "surfaceElevation" => m.surfaceElevation, "rockTop" => phaseTuple(rockTop),
    )
end
ρModel, VpModel, VsModel = modelResult["ρ"], modelResult["Vp"], modelResult["Vs"]
materialMask = BitArray(modelResult["material"])
iceThickness, waterDepth = modelResult["iceThickness"], modelResult["waterDepth"]

# ---- sanity checks
@assert size(ρModel) == (Nx, Ny, Nz)
@assert all(materialMask[:, :, 1]) "the bottom layer must be all rock"
@assert !any(materialMask[:, :, end]) "the top layer must be all air (raise altMax)"
for j in 1:Ny, i in 1:Nx   # no air pocket below material in any column
    k = findlast(@view materialMask[i, j, :])
    @assert all(@view materialMask[i, j, 1:k]) "air below material in column ($i, $j)"
end
topIndex = [findlast(@view materialMask[i, j, :]) for i in 1:Nx, j in 1:Ny]
vpMaximum = maximum(VpModel[materialMask])          # km/s
vsMinimumSolid = minimum(VsModel[materialMask .& (VsModel .> 0)])
nFluid = count(materialMask .& (VsModel .<= 0))
@printf("model: %.1f %% material, max ice %.0f m, max water %.0f m, Vp max %.2f km/s, Vs min %.2f km/s, %d pure-fluid cells\n",
    100count(materialMask) / length(materialMask), maximum(iceThickness), maximum(waterDepth),
    vpMaximum, vsMinimumSolid, nFluid)
nFluid > 0 && @warn("pure-water cells (Vs = 0) exist; the elastic OPT operator treats them with μ = 0, untested in flexOPT so far")
