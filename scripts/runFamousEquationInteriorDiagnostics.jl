#!/usr/bin/env julia
#
# Boundary-free, local-recipe diagnostics for the PDE gallery in
# src/motorsOPT/famousEquations.jl.  This is the inexpensive first stage:
# recipes that fail here should not be sent to a periodic convergence sweep.

import Pkg

const FLEXOPT_ROOT = get(ENV, "FLEXOPT_ROOT", normpath(joinpath(@__DIR__, "..")))
Pkg.activate(FLEXOPT_ROOT)

using JLD2
using KernelAbstractions: CPU
using Printf
using Statistics
using Symbolics

include(joinpath(FLEXOPT_ROOT, "src", "batchFiles", "batchGPU.jl"))
include(joinpath(FLEXOPT_ROOT, "src", "commonBatchs.jl"))
include(joinpath(FLEXOPT_ROOT, "src", "planet1D.jl"))
include(joinpath(FLEXOPT_ROOT, "src", "GeoPoints.jl"))
include(joinpath(FLEXOPT_ROOT, "src", "flexOPT.jl"))
include(joinpath(@__DIR__, "famousEquationBenchmarkMatrix.jl"))

using .commonBatchs
using .flexOPT
using .FamousEquationBenchmarkMatrix

const RECIPES = SCHEMES

function make_parameters(equation, recipe; delta=1.0)
    dimensions = equation.space_dimension + equation.has_time
    Δ = Tuple(fill(Float64(delta), dimensions))
    points_time = equation.has_time ? recipe.points_time : 1
    order_time = equation.has_time ? recipe.order_b_time : 0
    field_itpl = (
        ptsSpace=recipe.field_points,
        ptsTime=1,
        offsetSpace=recipe.field_offset,
        offsetTime=equation.has_time ? 1.0 : 0.0,
        YorderBspace=recipe.field_order,
        YorderBtime=-1,
    )
    material_itpl = (
        ptsSpace=recipe.material_points,
        ptsTime=1,
        offsetSpace=recipe.material_offset,
        offsetTime=equation.has_time ? 1.0 : 0.0,
        YorderBspace=recipe.material_order,
        YorderBtime=-1,
    )
    return Dict{String,Any}(
        "famousEquationType" => equation.equation,
        "Δ" => Δ,
        "orderBtime" => order_time,
        "orderBspace" => recipe.order_b_space,
        "pointsInSpace" => recipe.points_space,
        "pointsInTime" => points_time,
        "supplementaryOrder" => recipe.supplementary_order,
        "fieldItpl" => field_itpl,
        "materItpl" => material_itpl,
        "nuGeometryMode" => :middle,
        "hierarchicalTestFunctions" => false,
        "evenOrderHalfShiftMode" => :none,
        "taylorInverseMode" => recipe.taylor_inverse_mode,
        "recipe_backend" => CPU(),
    )
end

function recipe_summary(equation, recipe)
    elapsed = @elapsed opt = makeOPTsemiSymbolic(
        make_parameters(equation, recipe))
    r = opt["recette"]
    lhs = r.lhs.Ajiννᶜ
    rhs = r.rhs.Γjiννᶜ
    lhs_numbers = r.numbersOfTheSystem.numbersOfTheSystemL
    rhs_numbers = r.numbersOfTheSystem.numbersOfTheSystemR
    blocks = ndims(lhs) >= 5 ? size(lhs, 5) : 1
    finite_lhs = all(x -> try
        isfinite(Float64(Symbolics.value(x)))
    catch
        true
    end, lhs)
    return (
        equation_label=equation.label,
        equation=equation.equation,
        recipe=recipe.name,
        status="constructed",
        build_seconds=elapsed,
        coordinates=equation.space_dimension + equation.has_time,
        fields=lhs_numbers.NtypeofFields,
        equations=lhs_numbers.NtypeofExpr,
        force_fields=rhs_numbers.NtypeofFields,
        lhs_shape=size(lhs),
        gamma_shape=size(rhs),
        test_blocks=blocks,
        gamma_nonzero=count(!iszero, rhs),
        lhs_symbolically_finite=finite_lhs,
        error="",
    )
end

function failed_summary(equation, recipe, exception, elapsed)
    return (
        equation_label=equation.label,
        equation=equation.equation,
        recipe=recipe.name,
        status="failed",
        build_seconds=elapsed,
        coordinates=equation.space_dimension + equation.has_time,
        fields=0,
        equations=0,
        force_fields=0,
        lhs_shape=(),
        gamma_shape=(),
        test_blocks=0,
        gamma_nonzero=0,
        lhs_symbolically_finite=false,
        error=sprint(showerror, exception),
    )
end

function main()
    quick = get(ENV, "FLEXOPT_DIAGNOSTIC_QUICK", "0") == "1"
    recipes = quick ? filter(r -> r.name in ("FD3", "OPT3-center"), RECIPES) : RECIPES
    equations = quick ? EQUATIONS[1:3] : EQUATIONS
    if get(ENV, "FLEXOPT_DIAGNOSTIC_MAX_2D", "0") == "1"
        equations = filter(equation -> equation.space_dimension <= 2, equations)
    end
    if get(ENV, "FLEXOPT_DIAGNOSTIC_WAVES_ONLY", "0") == "1"
        equations = filter(equation -> equation.has_time, equations)
    end
    if haskey(ENV, "FLEXOPT_DIAGNOSTIC_EQUATIONS")
        selected_equations = Set(split(ENV["FLEXOPT_DIAGNOSTIC_EQUATIONS"], ","))
        equations = filter(
            equation -> equation.equation in selected_equations, equations)
    end
    if haskey(ENV, "FLEXOPT_DIAGNOSTIC_SCHEMES")
        selected_recipes = Set(split(ENV["FLEXOPT_DIAGNOSTIC_SCHEMES"], ","))
        recipes = filter(recipe -> recipe.name in selected_recipes, recipes)
    end
    rows = NamedTuple[]
    output_dir = joinpath(FLEXOPT_ROOT, "scripts", "tmp",
        "famous_equation_interior_diagnostics")
    mkpath(output_dir)
    mode_label = quick ? "quick" :
        (get(ENV, "FLEXOPT_DIAGNOSTIC_WAVES_ONLY", "0") == "1" ?
            "waves" : "all")
    dimension_label =
        get(ENV, "FLEXOPT_DIAGNOSTIC_MAX_2D", "0") == "1" ? "_max2d" : ""
    time_label = "_time3_hierarchical"
    output_file = joinpath(output_dir,
        "recipe_construction_$(mode_label)$(dimension_label)$(time_label).jld2")

    for equation in equations, recipe in recipes
        @info "Interior recipe diagnostic" equation=equation.equation recipe=recipe.name
        start = time()
        row = try
            recipe_summary(equation, recipe)
        catch exception
            failed_summary(equation, recipe, exception, time() - start)
        end
        push!(rows, row)
        @printf("%-28s %-28s %-11s %8.3f s  Γ nnz=%d\n",
            equation.label, recipe.name, row.status,
            row.build_seconds, row.gamma_nonzero)
        row.status == "failed" && println("  ", row.error)
        # Preserve completed rows even when a later 3-D construction is
        # interrupted or exhausts the practical time/memory budget.
        jldsave(output_file;
            rows,
            equations=EQUATIONS,
            recipes=RECIPES,
            boundary_mode="none; local interior recipe only",
            diagnostic_stage="construction and Gamma-presence smoke test",
        )
    end

    println("Saved: ", output_file)
    return output_file
end

abspath(PROGRAM_FILE) == abspath(@__FILE__) && main()
