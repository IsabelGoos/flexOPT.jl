#!/usr/bin/env julia
#
# helpers/runPoissonMixedTestFunction.jl   (Claude-authored; does not modify src/)
#
# flexOPT test of mixed test functions W for
#   1DpoissonHetero / 2DpoissonHetero   (κ T')' = f
#   1DsismoFreqHetero                   ρ ω² u + (μ u')' = f     (SH, frequency domain)
#
# Theory: src/motorsOPT/testFunctionMoments.jl. The truncation error depends on
# W only through its moments M_m = ∫W s^m / ∫W; the weights of
# W = w₁ hat₁ + w₂ hat₂ are computed here with optimalTestFunctionWeights
# (mass term included, it never limits the order):
#   P = 5:  4/15, 11/15     -> predicted O(h⁶) instead of O(h⁴) for the hat
#   P = 7:  40/231, 191/231 -> predicted O(h⁸) instead of O(h⁶)
#   P = 3:  no improvement possible (O(h²) for any W once the medium varies).
# hat_a is the hat on [ν-a, ν+a]; hat₂ is obtained in flexOPT through the knots
# of W: trialFunctionRefPoints = flexOPTRefPoints([-2, 0, 2], P).
#
# The weak form is linear in W, so the mixture is assembled from the recipes
# of its components, each normalised by ΣΓ = ∫W:
#     A = Σ_i w_i A_i / ΣΓ_i,    Γ = Σ_i w_i Γ_i / ΣΓ_i.
# In 2-D the tested W is the tensor product w(x) w(y) of the 1-D mixture.
#
# Measure: boundary-free periodic residual RMS[(A T - Γ f)/Σ|Γ|], same material
# scenarios and machinery as scripts/runWavePeriodicFDvsOPT.jl.
#
# Usage (repository root):
#   FLEXOPT_TF_QUICK=1 julia --project=. helpers/runPoissonMixedTestFunction.jl
#   FLEXOPT_TF_EQUATIONS=1DsismoFreqHetero julia --project=. helpers/runPoissonMixedTestFunction.jl
#
# Environment: FLEXOPT_TF_EQUATIONS (1DpoissonHetero,2DpoissonHetero,1DsismoFreqHetero)
#              FLEXOPT_TF_POINTS_1D (3,5,7)  FLEXOPT_TF_POINTS_2D (5)
#              FLEXOPT_TF_MAX_N_2D (32)  FLEXOPT_TF_RESUME (1)
# Outputs: scripts/tmp/poisson_mixed_test_function/  (Poisson rows of the v1 run are kept)

include(joinpath(@__DIR__, "..", "scripts", "runWavePeriodicFDvsOPT.jl"))

using CairoMakie
CairoMakie.activate!()
# all exports of the moment module, including defaultStencilNodes (not re-exported by flexOPT)
using .flexOPT.TestFunctionMoments

const TF_SCHEMA_VERSION = "mixed_test_function_v2"

tf_list(name, default) = String.(filter(!isempty, strip.(split(get(ENV, name, default), ","))))
tf_ints(name, default) = parse.(Int, tf_list(name, default))

const TF_QUICK = get(ENV, "FLEXOPT_TF_QUICK", "0") == "1"
const TF_RESUME = get(ENV, "FLEXOPT_TF_RESUME", "1") == "1"
const TF_EQUATIONS = TF_QUICK ? ["1DsismoFreqHetero"] :
    tf_list("FLEXOPT_TF_EQUATIONS", "1DpoissonHetero,2DpoissonHetero,1DsismoFreqHetero")
const TF_POINTS = Dict(
    1 => TF_QUICK ? [5] : tf_ints("FLEXOPT_TF_POINTS_1D", "3,5,7"),
    2 => tf_ints("FLEXOPT_TF_POINTS_2D", "5"),
)
const TF_SIZES = Dict(
    1 => TF_QUICK ? [16, 32] : [16, 24, 32, 48, 64, 96, 128],
    2 => filter(<=(parse(Int, get(ENV, "FLEXOPT_TF_MAX_N_2D", "32"))), [8, 12, 16, 24, 32, 48]),
)
# scenarios not defined for an equation are skipped (material_scenarios)
const TF_SCENARIOS = TF_QUICK ? ["homogeneous", "long_material_phase0"] : [
    "homogeneous", "same_wave_phase0", "same_wave_phase_pi2",
    "long_material_phase0", "short_material_phase0",
    "density_only_phase0", "stiffness_only_phase0",
]
const TF_OUTPUT_DIR = joinpath(FLEXOPT_ROOT, "scripts", "tmp", "poisson_mixed_test_function")

# ---------------------------------------------------------------- weights (moment theory)

"Mixture weights (hat₁, hat₂) and predicted order for a P-point stencil, or nothing."
function tf_weights(P)
    P >= 5 || return nothing
    result = optimalTestFunctionWeights(defaultStencilNodes(P),
        [hatTestFunction(1), hatTestFunction(2)]; mass=true)
    result.weights === nothing && return nothing
    result.negative && @warn "negative mixture weight" P weights=result.weights
    return (w=Tuple(Float64.(result.weights)), order=result.order,
        hat_order=predictedOrder(defaultStencilNodes(P), hatTestFunction(1); mass=true).order)
end

const TF_WEIGHTS = let computed = [(P, tf_weights(P)) for P in unique(vcat(values(TF_POINTS)...))]
    Dict(P => w for (P, w) in computed if w !== nothing)
end

# ---------------------------------------------------------------- W components

"Knots (1-based flexOPT node indices) of hat_a centred on ν; `nothing` = default (all nodes)."
hat_knots(P, a) = a == 1 ? nothing : flexOPTRefPoints([-a, 0, a], P)

"Per-axis knot specification for trialFunctionRefPoints."
function refpoints(dimension, P, spacings)
    all(==(1), spacings) && return nothing
    knots = [something(hat_knots(P, a), collect(1:P)) for a in spacings]
    return dimension == 1 ? knots[1] : Tuple(knots)
end

"Test-function variants: name => list of (per-axis hat spacings, weight)."
function tf_variants(dimension, P)
    variants = Pair{String,Vector{Any}}[]
    push!(variants, "hat1" => Any[(ntuple(_ -> 1, dimension), 1.0)])
    P >= 5 || return variants
    push!(variants, "hat2" => Any[(ntuple(_ -> 2, dimension), 1.0)])
    if haskey(TF_WEIGHTS, P)
        w = TF_WEIGHTS[P].w
        parts = Any[]
        for combo in Iterators.product(ntuple(_ -> (1, 2), dimension)...)
            push!(parts, (Tuple(combo), prod(w[a] for a in combo)))
        end
        label = @sprintf("mixture %.4g hat1 + %.4g hat2", w[1], w[2])
        push!(variants, (dimension == 1 ? label : label * " (tensor)") => parts)
    end
    return variants
end

# ---------------------------------------------------------------- recipes

central_itpl(P) = (ptsSpace=1, ptsTime=1, offsetSpace=(P - 1) / 2, offsetTime=0.0,
    YorderBspace=-1, YorderBtime=-1)

function tf_component_recipe(equation_id, dimension, P, spacings, dx)
    params = Dict{String,Any}(
        "famousEquationType" => equation_id,
        "Δ" => dimension == 1 ? Float64(dx) : Tuple(fill(Float64(dx), dimension)),
        "orderBtime" => 0,
        "orderBspace" => 1,
        "pointsInSpace" => P,
        "pointsInTime" => 1,
        "supplementaryOrder" => 2,
        "fieldItpl" => central_itpl(P),
        "materItpl" => central_itpl(P),
        "nuGeometryMode" => :middle,
        "hierarchicalTestFunctions" => false,
        "evenOrderHalfShiftMode" => :none,
        "taylorInverseMode" => :hierarchical_constrained,
        "recipe_backend" => CPU(),
        "variationalForm" => :weak,
    )
    knots = refpoints(dimension, P, spacings)
    knots === nothing || (params["trialFunctionRefPoints"] = knots)
    return makeOPTsemiSymbolic(params)["recette"]
end

"ΣΓ of a recipe (= ∫W up to the common h scaling); placeholder force symbols are set to one."
function gamma_sum(recipe)
    compiled = flexOPT._compile_coefficient_recipes(recipe.rhs.Γjiννᶜ)
    return real(sum(evaluate_coefficient(compiled[j, 1, 1, 1], Dict{Any,Any}())
                    for j in axes(compiled, 1)))
end

function combine_recipes(recipes, weights)
    reference = first(recipes)
    for r in recipes
        size(r.lhs.Ajiννᶜ) == size(reference.lhs.Ajiννᶜ) ||
            error("component recipes have different stencils")
    end
    scales = [w / gamma_sum(r) for (r, w) in zip(recipes, weights)]
    A = sum(scale .* r.lhs.Ajiννᶜ for (r, scale) in zip(recipes, scales))
    G = sum(scale .* r.rhs.Γjiννᶜ for (r, scale) in zip(recipes, scales))
    return merge(reference, (
        lhs=merge(reference.lhs, (; Ajiννᶜ=A)),
        rhs=merge(reference.rhs, (; Γjiννᶜ=G)),
    ))
end

short_err(err) = first(replace(sprint(showerror, err), '\n' => ' '), 300)

# ---------------------------------------------------------------- run

function tf_run!(rows, done, equation_id, checkpoint)
    equation = only(filter(e -> e.equation == equation_id, EQUATIONS))
    dimension = equation.space_dimension
    scenarios = filter(s -> s.name in TF_SCENARIOS, material_scenarios(equation))
    for P in TF_POINTS[dimension], n in TF_SIZES[dimension]
        dx = 2pi / n
        components = Dict{Any,Any}()          # spacings => recipe (shared by the variants)
        for (name, parts) in tf_variants(dimension, P)
            key = (equation_id, P, name, n)
            key in done && continue
            @info "mixed test function" equation_id P n W=name
            status = "ok"
            numeric = nothing
            t0 = time()
            try
                recipes = map(parts) do (spacings, _)
                    get!(() -> tf_component_recipe(equation_id, dimension, P, spacings, dx),
                        components, spacings)
                end
                numeric = compile_recipe(combine_recipes(recipes, [p[2] for p in parts]))
            catch err
                err isa InterruptException && rethrow()
                status = "recipe failed: " * short_err(err)
                @warn "recipe failed" equation_id P n W=name status
            end
            seconds = time() - t0
            for scenario in scenarios
                absolute = NaN
                relative = NaN
                row_status = status
                if numeric !== nothing
                    try
                        absolute, relative, _ = residual_error(numeric, equation, scenario, n, dx, nothing)
                    catch err
                        err isa InterruptException && rethrow()
                        row_status = "residual failed: " * short_err(err)
                    end
                end
                push!(rows, (; equation=equation_id, dimension, P, W=name,
                    scenario=scenario.name, n, dx,
                    absolute=Float64(absolute), relative=Float64(relative),
                    seconds=Float64(seconds), status=row_status))
            end
            push!(done, key)
            checkpoint()
        end
    end
end

function tf_orders(rows)
    out = NamedTuple[]
    for key in unique((r.equation, r.P, r.W, r.scenario) for r in rows)
        rs = sort(filter(r -> (r.equation, r.P, r.W, r.scenario) == key &&
            isfinite(r.absolute) && r.absolute > 0, rows); by=r -> r.n)
        length(rs) < 2 && continue
        orders = [log(rs[i].absolute / rs[i+1].absolute) / log(rs[i+1].n / rs[i].n)
                  for i in 1:length(rs)-1]
        predicted = startswith(key[3], "mixture") && haskey(TF_WEIGHTS, key[2]) ?
            TF_WEIGHTS[key[2]].order : NaN
        push!(out, (; equation=key[1], dimension=rs[1].dimension, P=key[2], W=key[3],
            scenario=key[4], n_finest=rs[end].n, error_finest=rs[end].absolute,
            order_last=orders[end], order_median=median(orders), predicted_order=predicted))
    end
    return out
end

function tf_print(summary)
    for (equation, P) in unique((s.equation, s.P) for s in summary)
        println("\n== $(equation), P = $P   (order at the last pair of grids; error at the finest grid)")
        haskey(TF_WEIGHTS, P) && println("   moment theory: hat₁ order $(TF_WEIGHTS[P].hat_order), mixture order $(TF_WEIGHTS[P].order)")
        sel = filter(s -> s.equation == equation && s.P == P, summary)
        scenarios = unique(s.scenario for s in sel)
        @printf("%-40s", "W")
        foreach(sc -> @printf("%24s", first(sc, 22)), scenarios)
        println()
        for W in unique(s.W for s in sel)
            @printf("%-40s", W)
            for sc in scenarios
                m = filter(s -> s.W == W && s.scenario == sc, sel)
                isempty(m) ? @printf("%24s", "-") :
                    @printf("%24s", @sprintf("%.2f  (%.2e)", m[1].order_last, m[1].error_finest))
            end
            println()
        end
    end
end

const TF_COLORS = Dict(3 => :gray40, 4 => :purple, 5 => :royalblue, 6 => :teal, 7 => :crimson)

function tf_plot(rows, equation)
    sel = filter(r -> r.equation == equation && isfinite(r.absolute) && r.absolute > 0, rows)
    isempty(sel) && return nothing
    scenarios = unique(r.scenario for r in sel)
    fig = Figure(size=(300 * length(scenarios) + 300, 340))
    local ax_last
    for (j, scenario) in enumerate(scenarios)
        ax = Axis(fig[1, j]; title=scenario, xlabel="n", ylabel=j == 1 ? "RMS[(AT-Γf)/Σ|Γ|]" : "",
            xscale=log10, yscale=log10)
        for (P, W) in unique((r.P, r.W) for r in sel)
            rs = sort(filter(r -> r.P == P && r.W == W && r.scenario == scenario, sel); by=r -> r.n)
            isempty(rs) && continue
            style = startswith(W, "mixture") ? :solid : W == "hat2" ? :dot : :dash
            scatterlines!(ax, [r.n for r in rs], [r.absolute for r in rs];
                color=get(TF_COLORS, P, :black), linestyle=style, markersize=5,
                label="P=$P $(W)")
        end
        ax_last = ax
    end
    Legend(fig[1, length(scenarios) + 1], ax_last; labelsize=9, framevisible=false)
    Label(fig[0, 1:length(scenarios)], equation; fontsize=14)
    file = joinpath(TF_OUTPUT_DIR, "mixed_test_function_$(equation)" * (TF_QUICK ? "_quick" : "") * ".png")
    save(file, fig)
    return file
end

"Rows of the v1 (Poisson-only) run, given an equation field so they are not recomputed."
function tf_migrate_v1(previous_rows)
    return [merge((; equation=r.dimension == 1 ? "1DpoissonHetero" : "2DpoissonHetero"), r)
            for r in previous_rows]
end

function tf_main()
    mkpath(TF_OUTPUT_DIR)
    tag = TF_QUICK ? "_quick" : ""
    output_file = joinpath(TF_OUTPUT_DIR, "mixed_test_function$(tag).jld2")
    v1_file = joinpath(TF_OUTPUT_DIR, "poisson_mixed_test_function$(tag).jld2")
    rows = NamedTuple[]
    if TF_RESUME
        if isfile(output_file)
            previous = load(output_file)
            get(previous, "schema_version", "") == TF_SCHEMA_VERSION && append!(rows, previous["rows"])
        elseif isfile(v1_file)
            previous = load(v1_file)
            if get(previous, "schema_version", "") == "poisson_mixed_test_function_v1"
                append!(rows, tf_migrate_v1(previous["rows"]))
                @info "Poisson rows of the v1 run reused" count=length(rows)
            end
        end
    end
    done = Set{Any}((r.equation, r.P, r.W, r.n) for r in rows)
    checkpoint() = jldsave(output_file; schema_version=TF_SCHEMA_VERSION, rows)
    for P in sort(collect(keys(TF_WEIGHTS)))
        @info "moment theory" P weights=TF_WEIGHTS[P].w mixture_order=TF_WEIGHTS[P].order hat_order=TF_WEIGHTS[P].hat_order
    end
    for equation_id in TF_EQUATIONS
        tf_run!(rows, done, equation_id, checkpoint)
    end
    summary = tf_orders(rows)
    jldsave(output_file; schema_version=TF_SCHEMA_VERSION, rows, summary)
    for (file, data) in (("rows", rows), ("summary", summary))
        isempty(data) && continue
        open(joinpath(TF_OUTPUT_DIR, "mixed_test_function$(tag)_$(file).csv"), "w") do io
            names = collect(keys(first(data)))
            println(io, join(string.(names), ","))
            for r in data
                println(io, join((v isa AbstractString ? "\"$v\"" : string(v)
                                  for v in (getproperty(r, k) for k in names)), ","))
            end
        end
    end
    tf_print(summary)
    figures = filter(!isnothing, [tf_plot(rows, e) for e in unique(r.equation for r in rows)])
    failures = filter(r -> r.status != "ok", rows)
    isempty(failures) || @warn "Some rows failed (status column)" count=length(failures)
    @info "mixed test-function benchmark complete" output_file figures
end

abspath(PROGRAM_FILE) == abspath(@__FILE__) && tf_main()
