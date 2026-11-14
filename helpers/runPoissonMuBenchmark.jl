#!/usr/bin/env julia
#
# helpers/runPoissonMuBenchmark.jl   (Claude-authored; does not modify src/)
#
# Poisson 1-D / 2-D benchmark of the interpolated Taylor expansion of the
# manuscript (section 2.1.5, eq. compactTaylorWeighted):
#
#     ψ⁽ⁿ⁾(y) = Σ_μ Σ_l  Y_μ(y) ψ⁽ˡ⁾[μ] K_{l-n}(y - y_μ)
#
# Question: does the choice of the Taylor centres μ (field T) and μᶜ (material
# κ) -- in particular non-integer (half-, quarter-grid) positions -- control the
# accuracy of the weak-form integral ∫ ∂W Y_μᶜ Y_μ K K that builds A?
#
# Fixed: trial stencil (P points), hat test function W, supplementary order,
# CPU / Float64 recipes, one recipe per grid spacing (never rescaled).
# Varied: the μ set and the order of Y, applied to the field only
# ("field"), to the material only ("material") or to both ("both").
# References in every panel: FD_P (strong form) and the one-centre OPT_P
# ("OPT<P>-central", i.e. the recipe used so far).
#
# μ sets, written relative to the stencil centre in units of Δ:
#   central   {0}, indicator Y                  (the usual OPT_P)
#   nodes     every node, hat Y                 (integer)
#   nodesY0   every node, box-car Y
#   nodesY2   every node, quadratic B-spline Y
#   midpoints the P-1 half-grid points, hat Y   (non-integer for odd P)
#   dense     nodes + midpoints, hat Y
#   inner025  {-1/4, 0, 1/4}, hat Y
#   inner05   {-1/2, 0, 1/2}, hat Y
#   pair05    {-1/2, 1/2}, hat Y
#   outer     P+1 half-grid points spanning [-(P/2), P/2] (EGU "staggered")
# Sets that coincide for a given P (e.g. pair05 = midpoints for P=3) are run once.
#
# Two measures:
#   periodic   boundary-free truncation residual RMS[(A T - Γ f) / Σ|Γ|] on a
#              2π-periodic domain, 1-D and 2-D; same machinery and material
#              scenarios as scripts/runWavePeriodicFDvsOPT.jl (manuscript
#              section "Numerical validation").
#   dirichlet  1-D solution error ‖T - cos x‖ / √N on [0, 10π] with exact
#              Dirichlet data; cases of scripts/runPoisson1DBenchmark_EGU.jl.
#
# Usage (repository root):
#   FLEXOPT_MU_QUICK=1 julia --project=. helpers/runPoissonMuBenchmark.jl   # smoke test
#   julia --project=. helpers/runPoissonMuBenchmark.jl                      # full campaign
#   FLEXOPT_MU_PARTS=periodic FLEXOPT_MU_DIMS=2 julia --project=. helpers/runPoissonMuBenchmark.jl
#
# Environment options (comma-separated lists):
#   FLEXOPT_MU_PARTS      periodic,dirichlet
#   FLEXOPT_MU_DIMS       1,2                   (periodic part only)
#   FLEXOPT_MU_POINTS_1D  3,4,5
#   FLEXOPT_MU_POINTS_2D  3
#   FLEXOPT_MU_KINDS_1D   central,nodes,nodesY0,nodesY2,midpoints,dense,inner025,inner05,pair05,outer
#   FLEXOPT_MU_KINDS_2D   central,nodes,midpoints,inner05
#   FLEXOPT_MU_ROLES      field,material,both
#   FLEXOPT_MU_FORM       weak | strong         (OPT recipes; FD is always strong)
#   FLEXOPT_MU_SUPP       2                     (supplementary order of the OPT recipes)
#   FLEXOPT_MU_MAX_N_2D   48
#   FLEXOPT_MU_RESUME     1                     (skip rows already in the checkpoint)
#
# Outputs in scripts/tmp/poisson_mu_benchmark/ (tag = <form>_supp<S>[_quick]):
#   poisson_mu_<tag>.jld2                rows + summary (checkpointed after every scheme)
#   poisson_mu_<tag>_rows.csv / _summary.csv
#   <part>_<d>D_P<P>_<tag>.png           error vs Δ, one row per role
#   <part>_<d>D_P<P>_<tag>_ratio.png     log10(error / error of OPT<P>-central)
#                                        at the finest grid, scheme × scenario

include(joinpath(@__DIR__, "..", "scripts", "runWavePeriodicFDvsOPT.jl"))
# ^ activates the environment, loads src/ (flexOPT, commonBatchs, ...), and
#   defines compile_recipe, residual_error, EQUATIONS, material_scenarios.
#   Its main() is guarded and does not run.

using CairoMakie
using SparseArrays

CairoMakie.activate!()

const MU_SCHEMA_VERSION = "poisson_mu_benchmark_v1"

env_list(name, default) =
    String.(filter(!isempty, strip.(split(get(ENV, name, default), ","))))
env_ints(name, default) = parse.(Int, env_list(name, default))
env_syms(name, default) = Symbol.(env_list(name, default))

const MU_QUICK = get(ENV, "FLEXOPT_MU_QUICK", "0") == "1"
const MU_FORM = Symbol(lowercase(get(ENV, "FLEXOPT_MU_FORM", "weak")))
MU_FORM in (:weak, :strong) || error("FLEXOPT_MU_FORM must be weak or strong")
const MU_SUPP = parse(Int, get(ENV, "FLEXOPT_MU_SUPP", "2"))
const MU_RESUME = get(ENV, "FLEXOPT_MU_RESUME", "1") == "1"
const MU_TAG = "$(MU_FORM)_supp$(MU_SUPP)" * (MU_QUICK ? "_quick" : "")

const ALL_KINDS = [:central, :nodes, :nodesY0, :nodesY2, :midpoints, :dense,
    :inner025, :inner05, :pair05, :outer]

const MU_PARTS = MU_QUICK ? ["periodic", "dirichlet"] :
    env_list("FLEXOPT_MU_PARTS", "periodic,dirichlet")
const MU_DIMS = MU_QUICK ? [1] : env_ints("FLEXOPT_MU_DIMS", "1,2")
const MU_POINTS = Dict(
    1 => MU_QUICK ? [3] : env_ints("FLEXOPT_MU_POINTS_1D", "3,4,5"),
    2 => MU_QUICK ? [3] : env_ints("FLEXOPT_MU_POINTS_2D", "3"),
)
const MU_KINDS = Dict(
    1 => MU_QUICK ? [:central, :nodes, :midpoints] :
        env_syms("FLEXOPT_MU_KINDS_1D", join(string.(ALL_KINDS), ",")),
    2 => MU_QUICK ? [:central, :nodes, :midpoints] :
        env_syms("FLEXOPT_MU_KINDS_2D", "central,nodes,midpoints,inner05"),
)
const MU_ROLES = MU_QUICK ? [:field, :material] :
    env_syms("FLEXOPT_MU_ROLES", "field,material,both")

const PERIODIC_SCENARIOS = MU_QUICK ? ["same_wave_phase0"] : [
    "homogeneous", "same_wave_phase0", "same_wave_phase_pi2",
    "long_material_phase0", "short_material_phase0",
]
const PERIODIC_SIZES = Dict(
    1 => MU_QUICK ? [16, 32] : [16, 24, 32, 48, 64, 96, 128, 192],
    2 => MU_QUICK ? [8, 16] :
        filter(<=(parse(Int, get(ENV, "FLEXOPT_MU_MAX_N_2D", "48"))),
            [8, 12, 16, 24, 32, 48, 64]),
)
const DIRICHLET_LOG_H_INVERSE = MU_QUICK ? [1.0, 1.5] : collect(0.0:0.5:3.0)
const DIRICHLET_L = 10pi
const MU_OUTPUT_DIR = joinpath(FLEXOPT_ROOT, "scripts", "tmp", "poisson_mu_benchmark")

# ------------------------------------------------------------------ μ sets

"(points, offset, Y order) in the convention of investigateDependencies (src/motorsOPT/others.jl)."
function mu_geometry(P::Int, kind::Symbol)
    c = (P - 1) / 2
    kind === :central   && return (points=1, offset=c, order=-1)
    kind === :nodes     && return (points=P, offset=0.0, order=1)
    kind === :nodesY0   && return (points=P, offset=0.0, order=0)
    kind === :nodesY2   && return (points=P, offset=0.0, order=2)
    kind === :midpoints && return (points=P - 1, offset=0.5, order=1)
    kind === :dense     && return (points=2P - 1, offset=0.0, order=1)
    kind === :inner025  && return (points=3, offset=c - 0.25, order=1)
    kind === :inner05   && return (points=3, offset=c - 0.5, order=1)
    kind === :pair05    && return (points=2, offset=c - 0.5, order=1)
    kind === :outer     && return (points=P + 1, offset=-0.5, order=1)
    throw(ArgumentError("unknown μ kind $kind"))
end

"μ positions relative to the stencil centre, in Δ (same formula as investigateDependencies)."
function mu_positions(P, g)
    c = (P - 1) / 2
    g.points == 1 && return [g.offset - c]
    step = (P - 1 - 2 * g.offset) / (g.points - 1)
    return [g.offset + (j - 1) * step - c for j in 1:g.points]
end

_isint(p) = isapprox(p, round(p); atol=1e-10)
format_position(p) = _isint(p) ? string(Int(round(p))) : string(round(p; digits=3))
format_positions(ps) = "{" * join(format_position.(ps), ",") * "}"
function mu_class(ps)
    all(_isint, ps) && return "integer"
    any(_isint, ps) && return "mixed"
    return "non-integer"
end

mu_itpl(g) = (
    ptsSpace=g.points, ptsTime=1,
    offsetSpace=Float64(g.offset), offsetTime=0.0,
    YorderBspace=g.order, YorderBtime=-1,
)

function mu_scheme(name, family, P, role, kind, field, material)
    pf = mu_positions(P, field)
    pm = mu_positions(P, material)
    shown, shown_order = role == "material" ? (pm, material.order) : (pf, field.order)
    label = role == "reference" ? name :
        "$(kind) μ=$(format_positions(shown)) Y$(shown_order)"
    return (; name, family, P, role, kind=String(kind), field, material,
        mu_field=format_positions(pf), mu_material=format_positions(pm),
        # integer / non-integer with respect to the grid nodes (not to the centre)
        mu_class=mu_class(shown .+ (P - 1) / 2), label,
        form=family === :convFD ? :strong : MU_FORM)
end

function mu_schemes(P, kinds, roles)
    central = mu_geometry(P, :central)
    schemes = NamedTuple[
        mu_scheme("FD$P", :convFD, P, "reference", :central, central, central),
        mu_scheme("OPT$P-central", :OPT, P, "reference", :central, central, central),
    ]
    seen = Set{Any}()
    for kind in kinds
        kind === :central && continue
        g = mu_geometry(P, kind)
        key = (round.(mu_positions(P, g); digits=8), g.order)
        if key in seen
            @info "μ set coincides with an earlier one for this P; skipped" P kind
            continue
        end
        push!(seen, key)
        for role in roles
            field = role in (:field, :both) ? g : central
            material = role in (:material, :both) ? g : central
            push!(schemes, mu_scheme("OPT$P-$(role)-$(kind)", :OPT, P,
                String(role), kind, field, material))
        end
    end
    return schemes
end

# ------------------------------------------------------------------ recipes

function mu_recipe(equation_id, dimension, scheme, dx; geometry_mode=:middle)
    is_fd = scheme.family === :convFD
    params = Dict{String,Any}(
        "famousEquationType" => equation_id,
        "Δ" => dimension == 1 ? Float64(dx) : Tuple(fill(Float64(dx), dimension)),
        "orderBtime" => 0,
        "orderBspace" => is_fd ? -1 : 1,
        "pointsInSpace" => scheme.P,
        "pointsInTime" => 1,
        "supplementaryOrder" => is_fd ? 0 : MU_SUPP,
        "fieldItpl" => mu_itpl(scheme.field),
        "materItpl" => mu_itpl(scheme.material),
        "nuGeometryMode" => geometry_mode,
        "hierarchicalTestFunctions" => false,
        "evenOrderHalfShiftMode" => :none,
        "taylorInverseMode" => :hierarchical_constrained,
        "recipe_backend" => CPU(),
        # never depend silently on the mutable package default
        "variationalForm" => scheme.form,
    )
    return makeOPTsemiSymbolic(params)
end

short_error(err) = first(replace(sprint(showerror, err), '\n' => ' '), 300)

mu_row(part, dimension, scheme, scenario, n, dx, absolute, relative, seconds, status) = (;
    part, dimension, P=scheme.P, scheme=scheme.name, family=String(scheme.family),
    role=scheme.role, kind=scheme.kind, label=scheme.label,
    mu_field=scheme.mu_field, mu_material=scheme.mu_material, mu_class=scheme.mu_class,
    form=String(scheme.form), scenario=String(scenario), n, dx=Float64(dx),
    absolute=Float64(absolute), relative=Float64(relative),
    recipe_seconds=Float64(seconds), status,
)

# ------------------------------------------------------------------ periodic part

function run_periodic!(rows, done, dimension, checkpoint)
    equation_id = dimension == 1 ? "1DpoissonHetero" : "2DpoissonHetero"
    equation = only(filter(e -> e.equation == equation_id, EQUATIONS))
    scenarios = filter(s -> s.name in PERIODIC_SCENARIOS, material_scenarios(equation))
    for P in MU_POINTS[dimension]
        for scheme in mu_schemes(P, MU_KINDS[dimension], MU_ROLES)
            for n in PERIODIC_SIZES[dimension]
                key = ("periodic", dimension, scheme.name, n)
                key in done && continue
                dx = 2pi / n
                @info "μ benchmark (periodic)" dimension n scheme=scheme.name mu_field=scheme.mu_field mu_material=scheme.mu_material
                status = "ok"
                numeric = nothing
                seconds = NaN
                try
                    t0 = time()
                    recipe = mu_recipe(equation_id, dimension, scheme, dx)["recette"]
                    seconds = time() - t0
                    numeric = compile_recipe(recipe)
                catch err
                    err isa InterruptException && rethrow()
                    status = "recipe failed: " * short_error(err)
                    @warn "recipe failed" scheme=scheme.name n status
                end
                for scenario in scenarios
                    absolute = NaN
                    relative = NaN
                    row_status = status
                    if numeric !== nothing
                        try
                            absolute, relative, _ = residual_error(
                                numeric, equation, scenario, n, dx, nothing)
                        catch err
                            err isa InterruptException && rethrow()
                            row_status = "residual failed: " * short_error(err)
                            @warn "residual failed" scheme=scheme.name n scenario=scenario.name row_status
                        end
                    end
                    push!(rows, mu_row("periodic", dimension, scheme, scenario.name,
                        n, dx, absolute, relative, seconds, row_status))
                end
                push!(done, key)
            end
            checkpoint()
        end
    end
end

# ------------------------------------------------------------------ Dirichlet part (EGU)

const DIRICHLET_CASE_SPECS = [
    ("homogeneous", ξ -> 1.0 + 0 * ξ),
    ("same_lambda", ξ -> sin(ξ) + 2),
    ("twice_lambda", ξ -> sin(ξ / 2) + 2),
    ("shifted_pi_3", ξ -> sin(ξ + pi / 3) + 2),
    ("lambda_2", ξ -> cos(ξ)^2 + 1),
    ("quadratic", ξ -> ξ^2 + 1),
]

function dirichlet_cases()
    @variables ξ
    D = Differential(ξ)
    specs = MU_QUICK ? filter(s -> s[1] == "same_lambda", DIRICHLET_CASE_SPECS) :
        DIRICHLET_CASE_SPECS
    return [(name=name, ξ=ξ, kappa=kappa(ξ), u=cos(ξ),
             force=expand_derivatives(D(kappa(ξ) * D(cos(ξ)))))
            for (name, kappa) in specs]
end

evaluate_on_grid(expr, ξ, grid) =
    Float64[Float64(Symbolics.value(substitute(expr, Dict(ξ => x)))) for x in grid]

function dirichlet_grid(case, log_h_inverse)
    requested_dx = exp(-log_h_inverse)
    nx = max(3, floor(Int, DIRICHLET_L / requested_dx) + 1)
    dx = DIRICHLET_L / (nx - 1)
    grid = collect(range(0.0, DIRICHLET_L; length=nx))
    return (; nx, dx, grid,
        beta=evaluate_on_grid(case.kappa, case.ξ, grid),
        exact=evaluate_on_grid(case.u, case.ξ, grid),
        force=evaluate_on_grid(case.force, case.ξ, grid))
end

"Solve A u = Γ f with exact Dirichlet values on the two end nodes (as in the EGU script)."
function dirichlet_solve(recipe, data, scheme)
    model_family = (
        models=[data.beta],
        modelPoints=(data.nx,),
        Δ=data.dx,
        modelName="poisson_mu_$(scheme.name)_$(data.nx)",
    )
    numerical = numericalOperatorConstruction(Dict{String,Any}(
        "optRec" => recipe,
        "modelFam" => model_family,
        "absorbingBoundaries" => nothing,
        "maskedRegionInSpace" => nothing,
        "backend" => :cpu,
        "representation" => "matrixfree",
        "compatibility_outputs" => false,
    ))
    prepared = prepareLinearSystem(numerical["numOperators"])
    prepared.NField == 1 || error("expected one Poisson field, got $(prepared.NField)")
    if scheme.family === :OPT
        offdiag = copy(prepared.R_force) - spdiagm(0 => diag(prepared.R_force))
        dropzeros!(offdiag)
        nnz(offdiag) > 0 || @warn "OPT source operator Γ is diagonal" scheme=scheme.name
    end
    known_inputs = vcat(vec(prepared.known_lhs_template), vec(reshape(data.force, :, 1, 1)))
    A = sparse(prepared.A_template)
    b = copy(prepared.b_template)
    prepared.b_fun!(b, known_inputs)
    A[1, :] .= 0
    A[1, 1] = 1
    b[1] = first(data.exact)
    A[end, :] .= 0
    A[end, end] = 1
    b[end] = last(data.exact)
    return real.(lu(A) \ b)
end

mu_rms(v) = norm(v) / sqrt(length(v))

function run_dirichlet!(rows, done, checkpoint)
    cases = dirichlet_cases()
    for P in MU_POINTS[1]
        for scheme in mu_schemes(P, MU_KINDS[1], MU_ROLES)
            for log_h_inverse in DIRICHLET_LOG_H_INVERSE
                representative = dirichlet_grid(first(cases), log_h_inverse)
                key = ("dirichlet", 1, scheme.name, representative.nx)
                key in done && continue
                @info "μ benchmark (Dirichlet)" nx=representative.nx scheme=scheme.name mu_field=scheme.mu_field mu_material=scheme.mu_material
                status = "ok"
                recipe = nothing
                seconds = NaN
                try
                    t0 = time()
                    recipe = mu_recipe("1DpoissonHetero", 1, scheme, representative.dx;
                        geometry_mode=:all)
                    seconds = time() - t0
                catch err
                    err isa InterruptException && rethrow()
                    status = "recipe failed: " * short_error(err)
                    @warn "recipe failed" scheme=scheme.name status
                end
                for case in cases
                    data = dirichlet_grid(case, log_h_inverse)
                    absolute = NaN
                    relative = NaN
                    row_status = status
                    if recipe !== nothing
                        try
                            solution = dirichlet_solve(recipe, data, scheme)
                            absolute = mu_rms(solution - data.exact)
                            relative = absolute / mu_rms(data.exact)
                        catch err
                            err isa InterruptException && rethrow()
                            row_status = "solve failed: " * short_error(err)
                            @warn "solve failed" scheme=scheme.name case=case.name row_status
                        end
                    end
                    push!(rows, mu_row("dirichlet", 1, scheme, case.name, data.nx, data.dx,
                        absolute, relative, seconds, row_status))
                end
                push!(done, key)
            end
            checkpoint()
        end
    end
end

# ------------------------------------------------------------------ summary

function mu_summary(rows)
    ok = filter(r -> isfinite(r.absolute) && r.absolute > 0, rows)
    central = Dict{Tuple,Float64}()
    for r in ok
        r.scheme == "OPT$(r.P)-central" &&
            (central[(r.part, r.dimension, r.P, r.scenario, r.n)] = r.absolute)
    end
    groups = Dict{Tuple,Vector{Any}}()
    for r in ok
        push!(get!(groups, (r.part, r.dimension, r.P, r.scheme, r.scenario), Any[]), r)
    end
    out = NamedTuple[]
    for (_, rs) in groups
        sort!(rs; by=r -> -r.dx)                       # coarse -> fine
        finest = rs[end]
        orders = [log(rs[i-1].absolute / rs[i].absolute) / log(rs[i-1].dx / rs[i].dx)
                  for i in 2:length(rs)]
        ratios = Float64[]
        for r in rs
            ref = get(central, (r.part, r.dimension, r.P, r.scenario, r.n), NaN)
            isfinite(ref) && ref > 0 && push!(ratios, r.absolute / ref)
        end
        ref_finest = get(central,
            (finest.part, finest.dimension, finest.P, finest.scenario, finest.n), NaN)
        push!(out, (;
            part=finest.part, dimension=finest.dimension, P=finest.P,
            scheme=finest.scheme, role=finest.role, kind=finest.kind,
            label=finest.label, mu_field=finest.mu_field,
            mu_material=finest.mu_material, mu_class=finest.mu_class,
            scenario=finest.scenario, n_finest=finest.n,
            error_finest=finest.absolute,
            ratio_vs_central_finest=finest.absolute / ref_finest,
            ratio_vs_central_geomean=isempty(ratios) ? NaN : exp(mean(log.(ratios))),
            order_last=isempty(orders) ? NaN : orders[end],
            order_median=isempty(orders) ? NaN : median(orders),
        ))
    end
    sort!(out; by=s -> (s.part, s.dimension, s.P, s.role, s.scheme, s.scenario))
    return out
end

csv_cell(x::AbstractString) = "\"" * replace(x, "\"" => "\"\"") * "\""
csv_cell(x) = string(x)

function write_csv(path, rows)
    isempty(rows) && return nothing
    names = collect(keys(first(rows)))
    open(path, "w") do io
        println(io, join(string.(names), ","))
        for r in rows
            println(io, join((csv_cell(getproperty(r, k)) for k in names), ","))
        end
    end
    return path
end

function print_summary(summary)
    for key in unique((s.part, s.dimension, s.P) for s in summary)
        part, dimension, P = key
        sel = filter(s -> (s.part, s.dimension, s.P) == key, summary)
        println("\n== $(part) $(dimension)-D, P = $(P)   (ratio < 1: better than OPT$(P)-central)")
        @printf("%-34s %-30s %-12s %12s %12s %8s\n",
            "scheme", "μ field / μ material", "class", "ratio(fine)", "ratio(geo)", "order")
        for name in unique(s.scheme for s in sel)
            ss = filter(s -> s.scheme == name, sel)
            finite_or_nan(v) = (w = filter(isfinite, v); isempty(w) ? NaN : median(w))
            @printf("%-34s %-30s %-12s %12.3g %12.3g %8.2f\n", name,
                first(ss).mu_field * " / " * first(ss).mu_material, first(ss).mu_class,
                finite_or_nan([s.ratio_vs_central_finest for s in ss]),
                finite_or_nan([s.ratio_vs_central_geomean for s in ss]),
                finite_or_nan([s.order_median for s in ss]))
        end
    end
end

# ------------------------------------------------------------------ figures

const KIND_COLORS = Dict(
    "nodes" => :royalblue, "nodesY0" => :deepskyblue, "nodesY2" => :navy,
    "midpoints" => :crimson, "dense" => :darkorange, "inner025" => :forestgreen,
    "inner05" => :limegreen, "pair05" => :purple, "outer" => :sienna,
)

function plot_errors(rows, part, dimension, P, output_dir)
    sel = filter(r -> r.part == part && r.dimension == dimension && r.P == P &&
        isfinite(r.absolute) && r.absolute > 0, rows)
    isempty(sel) && return nothing
    scenarios = unique(r.scenario for r in sel)
    roles = filter(role -> any(r -> r.role == role, sel), ["field", "material", "both"])
    isempty(roles) && return nothing
    ylabel = part == "periodic" ? "RMS[(AT-Γf)/Σ|Γ|]" : "RMS(T - exact)"
    fig = Figure(size=(300 * length(scenarios) + 330, 270 * length(roles) + 40))
    for (i, role) in enumerate(roles)
        local ax_last
        for (j, scenario) in enumerate(scenarios)
            ax = Axis(fig[i, j];
                title=i == 1 ? scenario : "",
                xlabel=i == length(roles) ? "Δ" : "",
                ylabel=j == 1 ? "μ on $(role)\n$(ylabel)" : "",
                xscale=log10, yscale=log10)
            for name in unique(r.scheme for r in sel)
                rs = sort(filter(r -> r.scheme == name && r.scenario == scenario, sel);
                    by=r -> r.dx)
                isempty(rs) && continue
                r0 = first(rs)
                xs = [r.dx for r in rs]
                ys = [r.absolute for r in rs]
                if r0.role == "reference"
                    lines!(ax, xs, ys; linewidth=2.5,
                        color=r0.family == "convFD" ? :gray60 : :black,
                        linestyle=r0.family == "convFD" ? :dash : :solid,
                        label=r0.label)
                elseif r0.role == role
                    scatterlines!(ax, xs, ys; markersize=6,
                        color=get(KIND_COLORS, r0.kind, :magenta),
                        marker=r0.mu_class == "integer" ? :circle :
                            r0.mu_class == "mixed" ? :diamond : :utriangle,
                        label=r0.label)
                end
            end
            ax_last = ax
        end
        Legend(fig[i, length(scenarios) + 1], ax_last; labelsize=9, framevisible=false)
    end
    Label(fig[0, 1:length(scenarios)],
        "Poisson $(dimension)-D, $(part), P = $(P), OPT form = $(MU_FORM)   " *
        "(markers: circle = μ on nodes, diamond = mixed, triangle = off-node μ)"; fontsize=14)
    file = joinpath(output_dir, "$(part)_$(dimension)D_P$(P)_$(MU_TAG).png")
    save(file, fig)
    return file
end

function plot_ratio(summary, part, dimension, P, output_dir)
    sel = filter(s -> s.part == part && s.dimension == dimension && s.P == P &&
        s.role != "reference" || (s.part == part && s.dimension == dimension &&
        s.P == P && s.scheme == "FD$(P)"), summary)
    isempty(sel) && return nothing
    names = unique(s.scheme for s in sel)
    scenarios = unique(s.scenario for s in sel)
    M = fill(NaN, length(scenarios), length(names))
    for s in sel
        v = s.ratio_vs_central_finest
        isfinite(v) && v > 0 &&
            (M[findfirst(==(s.scenario), scenarios), findfirst(==(s.scheme), names)] = log10(v))
    end
    ylabels = [begin
        s = first(filter(x -> x.scheme == name, sel))
        s.role == "reference" ? name : "$(s.role): $(s.label)"
    end for name in names]
    fig = Figure(size=(330 + 95 * length(scenarios), 140 + 20 * length(names)))
    ax = Axis(fig[1, 1];
        xticks=(1:length(scenarios), scenarios), xticklabelrotation=pi / 6,
        yticks=(1:length(names), ylabels), yticklabelsize=9,
        title="log10(error / error OPT$(P)-central), finest grid -- $(part) $(dimension)-D, P=$(P), $(MU_FORM)")
    finite = filter(isfinite, vec(M))
    lim = isempty(finite) ? 1.0 : max(0.5, maximum(abs, finite))
    hm = heatmap!(ax, 1:length(scenarios), 1:length(names), M;
        colormap=:balance, colorrange=(-lim, lim), nan_color=:gray90)
    for i in eachindex(scenarios), j in eachindex(names)
        isfinite(M[i, j]) && text!(ax, i, j; text=@sprintf("%.2f", M[i, j]),
            align=(:center, :center), fontsize=8)
    end
    Colorbar(fig[1, 2], hm; label="< 0 : better than the one-centre OPT")
    file = joinpath(output_dir, "$(part)_$(dimension)D_P$(P)_$(MU_TAG)_ratio.png")
    save(file, fig)
    return file
end

# ------------------------------------------------------------------ main

function mu_main()
    mkpath(MU_OUTPUT_DIR)
    output_file = joinpath(MU_OUTPUT_DIR, "poisson_mu_$(MU_TAG).jld2")
    rows = NamedTuple[]
    if MU_RESUME && isfile(output_file)
        previous = load(output_file)
        if get(previous, "schema_version", "") == MU_SCHEMA_VERSION
            append!(rows, previous["rows"])
            @info "Resuming from checkpoint" output_file rows=length(rows)
        end
    end
    done = Set{Any}((r.part, r.dimension, r.scheme, r.n) for r in rows)
    checkpoint() = jldsave(output_file;
        schema_version=MU_SCHEMA_VERSION, form=String(MU_FORM),
        supplementary_order=MU_SUPP, rows)

    for part in MU_PARTS
        if part == "periodic"
            for dimension in MU_DIMS
                run_periodic!(rows, done, dimension, checkpoint)
            end
        elseif part == "dirichlet"
            run_dirichlet!(rows, done, checkpoint)
        else
            @warn "unknown part (expected periodic or dirichlet)" part
        end
    end

    summary = mu_summary(rows)
    jldsave(output_file; schema_version=MU_SCHEMA_VERSION, form=String(MU_FORM),
        supplementary_order=MU_SUPP, rows, summary)
    write_csv(joinpath(MU_OUTPUT_DIR, "poisson_mu_$(MU_TAG)_rows.csv"), rows)
    write_csv(joinpath(MU_OUTPUT_DIR, "poisson_mu_$(MU_TAG)_summary.csv"), summary)
    print_summary(summary)

    figures = String[]
    for key in unique((r.part, r.dimension, r.P) for r in rows)
        for f in (plot_errors(rows, key..., MU_OUTPUT_DIR),
                  plot_ratio(summary, key..., MU_OUTPUT_DIR))
            f === nothing || push!(figures, f)
        end
    end
    failures = filter(r -> r.status != "ok", rows)
    isempty(failures) || @warn "Some rows failed (see status column of the CSV)" count=length(failures) schemes=unique(r.scheme for r in failures)
    @info "μ benchmark complete" output_file figures
    return output_file
end

abspath(PROGRAM_FILE) == abspath(@__FILE__) && mu_main()
