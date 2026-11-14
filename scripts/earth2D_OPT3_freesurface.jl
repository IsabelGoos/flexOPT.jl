# earth2D_OPT3_freesurface.jl
#
# OPT3 (3 points in space and time, natural weak form) on the whole 2-D Earth
# disk with a closed traction-free surface, for the same 1-D / 2-D StagYY models
# as scripts/earth2D_StagYY_neutrino_seismo.jl (FD3), and comparison with FD3.
#
# Nothing in src/ is modified. The surface treatment re-uses the validated
# `clipped_available` recipe of HomogeneousElastic2DBenchmark_FreeSurface.ipynb
# (3×2 stencil, natural weak flux σ·n = 0, test integral clipped at the
# surface). That recipe is written for a surface whose void lies at +z. On a
# circle every surface node is given one of four orientations (+z, −z, +x, −x),
# chosen as the void neighbour best aligned with the outward radial normal.
# For the other three orientations the same recipe is assembled on a flipped /
# transposed copy of the model and the rows are mapped back with a signed field
# permutation (e.g. ux ↔ uz, uz → −uz). The isotropic elastic equations are
# invariant under these maps, so each mapped row is exactly the clipped weak
# equation for that orientation (checked on a model operator: mapping error 0).
#
# Workflow
#   1. run the FD3 script once with the SAME environment variables; it writes
#      data/earth2D_stagyy/<case>/earth2D_products.jld2 (model + FD3 traces);
#   2. run this script; it reads that file, builds OPT3, marches, compares.
#
# Coarse smoke test first (a few minutes; direct LU is affordable):
#   FLEXOPT_EARTH2D_DX=50e3 FLEXOPT_EARTH2D_PERIOD=100 FLEXOPT_EARTH2D_DURATION=1500 \
#     julia --project=. --threads=auto scripts/earth2D_StagYY_neutrino_seismo.jl
#   FLEXOPT_EARTH2D_DX=50e3 FLEXOPT_EARTH2D_PERIOD=100 FLEXOPT_EARTH2D_DURATION=1500 \
#     julia --project=. --threads=auto scripts/earth2D_OPT3_freesurface.jl
# Production (10 km, T = 30 s) then uses the lumped predictor–corrector:
#   FLEXOPT_OPT3_SCHEME=predictor_corrector julia --project=. --threads=auto scripts/earth2D_OPT3_freesurface.jl

import Pkg
function find_flexopt_root(start=@__DIR__)
    candidates = haskey(ENV, "FLEXOPT_ROOT") ? [ENV["FLEXOPT_ROOT"]] : String[]
    directory = abspath(start)
    while true
        push!(candidates, directory)
        parent = dirname(directory)
        parent == directory && break
        directory = parent
    end
    for candidate in unique(abspath.(expanduser.(candidates)))
        isfile(joinpath(candidate, "src", "flexOPT.jl")) && return candidate
    end
    error("Cannot locate flexOPT; set ENV[\"FLEXOPT_ROOT\"]")
end
flexopt_root = find_flexopt_root()
Pkg.activate(flexopt_root)
# Same backend bootstrap as HomogeneousElastic2DBenchmark_FreeSurface.ipynb.
using Metal
include(joinpath(flexopt_root, "src", "batchFiles", "batchGPU.jl"))
@show backend
include(joinpath(flexopt_root, "src", "commonBatchs.jl"))
include(joinpath(flexopt_root, "src", "flexOPT.jl"))
using .commonBatchs, .flexOPT
set_default_form!(:weak)
@assert DEFAULT_FORM[] === :weak
using LinearAlgebra, SparseArrays, Statistics, JLD2, CairoMakie, Printf
using Base.Threads
CairoMakie.activate!(type="png")
@show VERSION Threads.nthreads()

# ─────────────────────────────────────────────────────────────────────────────
# 0. Parameters (the case keys must match the FD3 run)
# ─────────────────────────────────────────────────────────────────────────────
stagyySnapshot = parse(Int, get(ENV, "FLEXOPT_STAGYY_SNAPSHOT", "199"))
dx = parse(Float64, get(ENV, "FLEXOPT_EARTH2D_DX", "10e3"))
dominantPeriod = parse(Float64, get(ENV, "FLEXOPT_EARTH2D_PERIOD", "30.0"))
caseDirectory = joinpath(flexopt_root, "data", "earth2D_stagyy",
    @sprintf("snap%05d_dx%dkm_T%ds", stagyySnapshot, round(Int, dx / 1e3), round(Int, dominantPeriod)))
productsFile = joinpath(caseDirectory, "earth2D_products.jld2")
isfile(productsFile) || error("run scripts/earth2D_StagYY_neutrino_seismo.jl first with the same " *
    "FLEXOPT_EARTH2D_DX / FLEXOPT_EARTH2D_PERIOD; missing $productsFile")
products = load(productsFile)
xs, zs = products["xs"], products["zs"]
material = BitMatrix(products["material"])
duration = parse(Float64, get(ENV, "FLEXOPT_EARTH2D_DURATION", string(last(products["time"]))))

modelsToRun = Symbol.(split(get(ENV, "FLEXOPT_OPT3_MODELS", "1D,2D"), ","))
# :direct factorises the full future operator once (fine at 50 km, heavy at
# 10 km). :predictor_corrector never factorises it: lumped 2×2 future blocks
# + Richardson corrections (the notebook's optTimeMarchingScheme option).
marchingScheme = Symbol(get(ENV, "FLEXOPT_OPT3_SCHEME", dx >= 25e3 ? "direct" : "predictor_corrector"))
correctorIterations = parse(Int, get(ENV, "FLEXOPT_OPT3_CORRECTOR_ITERATIONS", "3"))
correctorTolerance = 1e-8
# dt = courant · dx / Vp,max. The flat benchmark used 0.2/(√2·2) ≈ 0.07;
# 0.2 is ~3× cheaper and is the first thing to check on the coarse test
# (rerun with FLEXOPT_OPT3_COURANT=0.07 and compare).
courant = parse(Float64, get(ENV, "FLEXOPT_OPT3_COURANT", "0.2"))
outputSampling = 1.0
snapshotTimes = (300.0, 450.0, 600.0)
sourceForce = 1.0e15
recipeCacheVersion = 9   # same cache as the free-surface notebook
checkFactorization = length(material) <= 400_000

N = size(material, 1)
@assert size(material) == (N, N) "the rotated-recipe trick needs a square grid"
@assert xs == zs

# ─────────────────────────────────────────────────────────────────────────────
# 1. OPT3 recipes (identical parameters to the free-surface notebook → cached)
# ─────────────────────────────────────────────────────────────────────────────
recipeSpacing = (1.0, 1.0, 1.0)   # :scaled_coordinates formulation
optParameters = Dict{String,Any}(
    "famousEquationType" => "2DsismoTimeIsoHeteroSingleForce",
    "Δ" => recipeSpacing,
    "orderBtime" => 1, "orderBspace" => 1,
    "pointsInSpace" => 3, "pointsInTime" => 3,
    "supplementaryOrder" => 2,
    "variationalForm" => DEFAULT_FORM[],
    "taylor_inverse_mode" => :weak_operator_optimized,
    "fieldItpl" => (ptsSpace=1, ptsTime=1, offsetSpace=1.0,
        offsetTime=1, YorderBspace=-1, YorderBtime=-1),
    "materItpl" => (ptsSpace=1, ptsTime=1, offsetSpace=1.0,
        offsetTime=1, YorderBspace=-1, YorderBtime=-1),
    "recipe_backend" => backend,
)
clippedSurfaceParameters = copy(optParameters)
clippedSurfaceParameters["pointsInSpace"] = (3, 2)
clippedSurfaceParameters["nuCentre"] = (2, 2)
clippedSurfaceParameters["fieldItpl"] = (ptsSpace=(3, 2), ptsTime=1,
    offsetSpace=(0.0, 0.0), offsetTime=1, YorderBspace=1, YorderBtime=-1)
clippedSurfaceParameters["trialFunctionRefPoints"] = ((1, 2, 3), (1, 2, 3), (1, 2, 3))
clippedSurfaceParameters["testIntegrationBounds"] = ((1.0, 3.0), (1.0, 2.0), (1.0, 3.0))
clippedSurfaceParameters["exactTaylorTotalDegree"] = 2

function cachedOPTRecipe(parameters, prefix)
    cacheParameters = Dict{String,Any}(k => v for (k, v) in parameters if k != "recipe_backend")
    cacheParameters["recipe_cache_version"] = recipeCacheVersion
    function produceRecipe(config)
        runtime = Dict{String,Any}(config)
        pop!(runtime, "hash_id", nothing); pop!(runtime, "recipe_cache_version", nothing)
        runtime["recipe_backend"] = backend
        makeOPTsemiSymbolic(runtime)
    end
    myProduceOrLoad(produceRecipe, cacheParameters, "semiSymbolic", prefix)
end
optRecipe = cachedOPTRecipe(optParameters, "elastic2D_OPT3")
clippedSurfaceRecipe = cachedOPTRecipe(clippedSurfaceParameters, "elastic2D_OPT3_surface_clipped_available")
@assert optRecipe["recette"].variationalForm === :natural_weak
timeMarching = optRecipe["recette"].numbersOfTheSystem.numbersOfTheSystemL.timeMarching

# ─────────────────────────────────────────────────────────────────────────────
# 2. Closed free surface: classification and orientations
# ─────────────────────────────────────────────────────────────────────────────
# Orientation o ↦ (direction of the void neighbour, model map, point map T,
# field permutation π, field signs σ). In the transformed frame the void is
# always at +z, where the clipped recipe is valid:  u'_a = σ_a u_{π(a)}.
const ORIENTATIONS = (
    top=(dir=(0, 1), model=identity, point=(i, k) -> (i, k), π=(1, 2), σ=(1, 1)),
    bottom=(dir=(0, -1), model=a -> a[:, end:-1:1], point=(i, k) -> (i, N + 1 - k), π=(1, 2), σ=(1, -1)),
    right=(dir=(1, 0), model=a -> permutedims(a), point=(i, k) -> (k, i), π=(2, 1), σ=(1, 1)),
    left=(dir=(-1, 0), model=a -> permutedims(a)[:, end:-1:1], point=(i, k) -> (k, N + 1 - i), π=(2, 1), σ=(1, -1)),
)

surfaceGroups = Dict(o => CartesianIndex{2}[] for o in keys(ORIENTATIONS))
surfacePoints = CartesianIndex{2}[]
surfaceNormals = SVector{2,Float64}[]
cornerCount = 0
for k in 2:N-1, i in 2:N-1
    material[i, k] || continue
    voids = [o for (o, spec) in pairs(ORIENTATIONS) if !material[i+spec.dir[1], k+spec.dir[2]]]
    isempty(voids) && continue
    normal = SVector(xs[i], zs[k]) / hypot(xs[i], zs[k])
    best = voids[argmax([dot(normal, SVector(Float64.(ORIENTATIONS[o].dir))) for o in voids])]
    push!(surfaceGroups[best], CartesianIndex(i, k))
    push!(surfacePoints, CartesianIndex(i, k)); push!(surfaceNormals, normal)
    length(voids) > 1 && (global cornerCount += 1)
end
# Material nodes whose 3×3 volume stencil touches a void diagonal only
# (concave staircase corners) keep the volume row: one pinned void node.
concaveCount = count(CartesianIndices((2:N-1, 2:N-1))) do I
    material[I] && all(material[I+CartesianIndex(d)] for d in ((1,0),(-1,0),(0,1),(0,-1))) &&
        !all(material[I+CartesianIndex(d)] for d in ((1,1),(1,-1),(-1,1),(-1,-1)))
end
@show N length(surfacePoints) Dict(o => length(v) for (o, v) in surfaceGroups) cornerCount concaveCount

# ─────────────────────────────────────────────────────────────────────────────
# 3. Operators for one model
# ─────────────────────────────────────────────────────────────────────────────
function scaled_models(ρ, vp, vs, dt; ρref)
    ρSI, vpSI, vsSI = ρ .* 1e3, vp .* 1e3, vs .* 1e3
    rhoOPT = ifelse.(material, ρSI ./ ρref, 1.0)        # void: neutral, rows pinned
    muOPT = ifelse.(material, ρSI .* vsSI .^ 2 ./ ρref .* (dt / dx)^2, 0.0)
    lambdaOPT = ifelse.(material, ρSI .* (vpSI .^ 2 .- 2vsSI .^ 2) ./ ρref .* (dt / dx)^2, 0.0)
    [rhoOPT, lambdaOPT, muOPT]
end

function assemble(recipe, models, name; masked=nothing, bc=nothing)
    points = getModelPoints(models[1], 3,
        recipe["recette"].numbersOfTheSystem.numbersOfTheSystemL.timeMarching)
    family = (models=models, modelPoints=points, Δ=recipeSpacing, modelName=name)
    numerical = numericalOperatorConstruction(Dict{String,Any}(
        "optRec" => recipe, "modelFam" => family,
        "absorbingBoundaries" => nothing, "maskedRegionInSpace" => masked,
        "boundaryConditions" => bc, "representation" => "matrixfree",
    ))["numOperators"]
    isnothing(bc) ? prepareLinearSystem(numerical) :
        prepareLinearSystem(numerical; free_surface_spacing=recipeSpacing[1:2])
end

# Map an operator assembled in a transformed frame back to the original one.
# Rows: (field, point) with field fastest; columns: (point, field, time) with
# point fastest — the layouts of prepareNumericalLinearSystem.
function map_back(prepared_t, spec)
    Np, Nf = prepared_t.NpointsSpace, prepared_t.NField
    LIp = LinearIndices((N, N))
    pointMap = Vector{Int}(undef, Np)
    for k in 1:N, i in 1:N
        pointMap[LIp[i, k]] = LIp[spec.point(i, k)...]
    end
    rowPerm = Vector{Int}(undef, Nf * Np); rowSign = similar(rowPerm, Float64)
    for p in 1:Np, f in 1:Nf
        a = spec.π[f]
        rowPerm[f + Nf * (p - 1)] = a + Nf * (pointMap[p] - 1)
        rowSign[f + Nf * (p - 1)] = spec.σ[a]
    end
    function column_map(ncol)
        nslots = ncol ÷ (Np * Nf)
        perm = Vector{Int}(undef, ncol); sgn = Vector{Float64}(undef, ncol)
        for t in 1:nslots, f in 1:Nf, p in 1:Np
            a = spec.π[f]
            c = p + Np * (f - 1) + Np * Nf * (t - 1)
            perm[c] = pointMap[p] + Np * (a - 1) + Np * Nf * (t - 1)
            sgn[c] = spec.σ[a]
        end
        perm, sgn
    end
    Dr = Diagonal(rowSign)
    transform(M) = begin
        cp, cs = column_map(size(M, 2))
        Dr * M[rowPerm, cp] * Diagonal(cs)
    end
    flexOPT._prepared_with_matrices(prepared_t,
        transform(prepared_t.A_unknown), transform(prepared_t.L_known), transform(prepared_t.R_force))
end

function build_operator(ρ, vp, vs, dt; label)
    ρref = median(ρ[material]) * 1e3
    models = scaled_models(ρ, vp, vs, dt; ρref)
    bc = BoundaryConditionSet(
        free_surface=flexOPT.FreeSurfacePointSet(copy(surfacePoints), copy(surfaceNormals)),
        cerjan=nothing, material_mask=BitArray(material), free_surface_mode=:pinned_void)
    t = @elapsed volume = assemble(optRecipe, models, "earth2D_$(label)_OPT3_volume"; bc=bc)
    @info "$label volume operator" seconds = t size = size(volume.A_unknown)
    prepared = volume
    for (o, spec) in pairs(ORIENTATIONS)
        points = surfaceGroups[o]
        isempty(points) && continue
        modelsT = [spec.model(m) for m in models]
        pointsT = [CartesianIndex(spec.point(Tuple(p)...)) for p in points]
        t = @elapsed surfaceT = assemble(clippedSurfaceRecipe, modelsT,
            "earth2D_$(label)_OPT3_surface_$(o)"; masked=pointsT)
        surface = o === :top ? surfaceT : map_back(surfaceT, spec)
        prepared = overlapBoundaryLinearSystem(prepared, surface, points; mode=:replace)
        @info "$label surface rows" orientation = o points = length(points) seconds = t
    end
    @assert prepared.NForceField == 2
    return prepared, ρref
end

# ─────────────────────────────────────────────────────────────────────────────
# 4. Time marching with receiver recording (propagateLinearSystem keeps every
#    snapshot and a dense source array, which does not fit at 1.6 M nodes)
# ─────────────────────────────────────────────────────────────────────────────
sourceIndex = products["sourceIndex"]
receiverIndices = products["receiverIndices"]
epicentralDistances = products["epicentralDistances"]
sourceFrequency = 1 / dominantPeriod
sourceDelay = 1.5 / sourceFrequency
ricker(t) = (a = π * sourceFrequency * (t - sourceDelay); (1 - 2a^2) * exp(-a^2))

# Threaded y = α·M·x + β·y using Mᵀ stored as CSC (= M in CSR): one row per task.
struct RowMajor{T}
    Mt::SparseMatrixCSC{T,Int}
end
rowmajor(M::SparseMatrixCSC) = RowMajor(sparse(transpose(M)))
function rowmul!(y, R::RowMajor, x, α=1.0, β=0.0)
    Mt = R.Mt; rv = rowvals(Mt); nz = nonzeros(Mt)
    @threads for r in 1:size(Mt, 2)
        acc = 0.0
        @inbounds for idx in nzrange(Mt, r)
            acc += nz[idx] * x[rv[idx]]
        end
        @inbounds y[r] = iszero(β) ? α * acc : α * acc + β * y[r]
    end
    y
end

function march(prepared, dt, ρref; label)
    Np, Nf = prepared.NpointsSpace, prepared.NField
    nT = prepared.timePointsUsedForOneStep
    @assert nT == 3
    steps = ceil(Int, duration / dt)
    stride = max(1, round(Int, outputSampling / dt))
    LIp = LinearIndices((N, N))
    srcLinear = LIp[sourceIndex]
    # force column (point, field=z, time slot) — layout of known_rhs_template
    srcColumns = [srcLinear + Np * (2 - 1) + Np * prepared.NForceField * (s - 1) for s in 1:nT]
    Rsrc = prepared.R_force[:, srcColumns]
    # scaled-coordinate body force: f' = F Δt² / (ρref Δx Δz)
    scale = sourceForce * dt^2 / (ρref * dx^2)
    L = rowmajor(sparse(prepared.L_known))
    b = zeros(Np * Nf); u = zeros(Np * Nf)
    known = zeros(Np * Nf * (nT - 1))      # [u(n−1); u(n)] in (point, field) blocks
    A = flexOPT.materializeConstantMatrix(prepared; sparse_output=true)
    setup = @elapsed if marchingScheme === :direct
        factor = lu(A)
    else
        predictor = prepareLumpedFuturePredictor(A, prepared)
        Arow = rowmajor(A)
        residual = similar(b); correction = similar(b)
    end
    @info "$label matrix setup" scheme = marchingScheme seconds = setup unknowns = size(A, 1)
    recLinear = [LIp[I] for I in receiverIndices]
    nstore = steps ÷ stride + 1
    ux = zeros(nstore, length(recLinear)); uz = similar(ux); times = zeros(nstore)
    snapshots = Dict{Float64,NTuple{2,Matrix{Float32}}}()
    stored = 1; maxResidual = 0.0
    blk = Np * Nf
    elapsed = @elapsed for step in 1:steps
        # source levels t_{n−1}, t_n, t_{n+1} for the unknown u(t_{n+1}), t_{n+1} = step·dt
        f3 = scale .* ricker.(((step - 2):step) .* dt)
        rowmul!(b, L, known, -1.0, 0.0)
        mul!(b, Rsrc, f3, 1.0, 1.0)
        if marchingScheme === :direct
            ldiv!(u, factor, b)
        else
            applyLumpedFuturePredictor!(u, predictor, b)
            for _ in 1:correctorIterations
                rowmul!(residual, Arow, u); residual .= b .- residual
                norm(residual) <= correctorTolerance * max(norm(b), eps()) && break
                applyLumpedFuturePredictor!(correction, predictor, residual)
                u .+= correction
            end
            if step % 500 == 0
                rowmul!(residual, Arow, u); residual .= b .- residual
                maxResidual = max(maxResidual, norm(residual) / max(norm(b), eps()))
            end
        end
        @views known[1:blk] .= known[blk+1:2blk]
        @views known[blk+1:2blk] .= u
        time = step * dt
        if step % stride == 0 && stored < nstore
            stored += 1; times[stored] = time
            for (j, p) in enumerate(recLinear)
                ux[stored, j] = u[p]; uz[stored, j] = u[p+Np]
            end
        end
        for ts in snapshotTimes
            abs(time - ts) < dt / 2 && (snapshots[ts] = (Float32.(reshape(u[1:Np], N, N)),
                Float32.(reshape(u[Np+1:2Np], N, N))))
        end
        if step % 1000 == 0
            amplitude = maximum(abs, u)
            isfinite(amplitude) && amplitude < 1e8 || error("$label OPT3 diverged at t=$(time) s")
            @info "$label" t = round(time) max_u = amplitude
        end
    end
    ur = similar(ux)
    for (j, I) in enumerate(receiverIndices)
        ex, ez = (xs[I[1]], zs[I[2]]) ./ hypot(xs[I[1]], zs[I[2]])
        ur[:, j] .= ux[:, j] .* ex .+ uz[:, j] .* ez
    end
    @info "$label OPT3 done" dt steps elapsed maxResidual
    (time=times[1:stored], ur=ur[1:stored, :], snapshots, dt, elapsed, maxResidual)
end

# ─────────────────────────────────────────────────────────────────────────────
# 5. Run
# ─────────────────────────────────────────────────────────────────────────────
vpMaximum = 1e3 * max(maximum(products["vp1D"][material]), maximum(products["vp2D"][material]))
dtOPT = courant * dx / vpMaximum
@show dx dtOPT marchingScheme courant

results = Dict{Symbol,Any}()
for label in modelsToRun
    ρ, vp, vs = products["ρ$(label)"], products["vp$(label)"], products["vs$(label)"]
    prepared, ρref = build_operator(ρ, vp, vs, dtOPT; label=String(label))
    if checkFactorization
        status = try
            lu(flexOPT.materializeConstantMatrix(prepared)); :regular
        catch err
            err isa SingularException ? :singular : rethrow()
        end
        @show label status
        @assert status === :regular "future operator singular: check the rotated surface rows"
    end
    results[label] = march(prepared, dtOPT, ρref; label=String(label))
end

# ─────────────────────────────────────────────────────────────────────────────
# 6. Comparison with FD3 (same model, same source, same receivers)
# ─────────────────────────────────────────────────────────────────────────────
fdTime = products["time"]
function resample(t, y, tq)
    out = similar(tq, Float64)
    for (n, τ) in enumerate(tq)
        j = clamp(searchsortedlast(t, τ), 1, length(t) - 1)
        w = (τ - t[j]) / (t[j+1] - t[j])
        out[n] = (1 - w) * y[j] + w * y[j+1]
    end
    out
end
comparisonWindow = min(duration, 1500.0)
metrics = Dict{Symbol,Any}()
for label in keys(results)
    fd = products["ur$(label)"]; opt = results[label]
    tq = collect(0.0:outputSampling:min(comparisonWindow, last(opt.time), last(fdTime)))
    corr = Float64[]; relerr = Float64[]
    for j in eachindex(epicentralDistances)
        a = resample(fdTime, fd[:, j], tq); b = resample(opt.time, opt.ur[:, j], tq)
        push!(corr, cor(a, b)); push!(relerr, norm(b - a) / max(norm(a), eps()))
    end
    metrics[label] = (correlation=corr, relative_error=relerr)
    @info "OPT3 vs FD3 ($label, 0–$(round(Int, comparisonWindow)) s)" median_correlation = median(corr) median_relative_error = median(relerr)
end

opt3File = joinpath(caseDirectory, @sprintf("earth2D_OPT3_%s_c%.2f.jld2", marchingScheme, courant))
jldsave(opt3File; dtOPT, courant, marchingScheme, epicentralDistances,
    surfaceCounts=Dict(String(k) => length(v) for (k, v) in surfaceGroups), cornerCount, concaveCount,
    time=Dict(String(k) => v.time for (k, v) in results),
    ur=Dict(String(k) => v.ur for (k, v) in results),
    snapshots=Dict(String(k) => v.snapshots for (k, v) in results),
    metrics=Dict(String(k) => v for (k, v) in metrics))

figure = Figure(size=(1300, 900))
for (column, label) in enumerate(sort(collect(keys(results))))
    ax = Axis(figure[1, column]; xlabel="time (s)", ylabel="epicentral distance (°)",
        title="$(label): FD3 (black) vs OPT3 closed free surface (red)")
    fd = products["ur$(label)"]; opt = results[label]
    for j in eachindex(epicentralDistances)
        reference = maximum(abs, fd[fdTime .< comparisonWindow, j])
        reference > 0 || continue
        lines!(ax, fdTime, epicentralDistances[j] .+ 2.4 .* clamp.(fd[:, j] ./ reference, -1.6, 1.6);
            color=:black, linewidth=0.7)
        lines!(ax, opt.time, epicentralDistances[j] .+ 2.4 .* clamp.(opt.ur[:, j] ./ reference, -1.6, 1.6);
            color=(:firebrick, 0.85), linewidth=0.7)
    end
    xlims!(ax, 0, comparisonWindow)
    axm = Axis(figure[2, column]; xlabel="epicentral distance (°)", ylabel="correlation",
        title="waveform correlation OPT3–FD3 (0–$(round(Int, comparisonWindow)) s)")
    scatterlines!(axm, epicentralDistances, metrics[label].correlation)
    ylims!(axm, -0.1, 1.05)
end
rowsize!(figure.layout, 2, Relative(0.25))
figurePath = joinpath(caseDirectory, @sprintf("earth2D_OPT3_vs_FD3_%s_c%.2f.png", marchingScheme, courant))
save(figurePath, figure)
@info "written" opt3File figurePath
figure
