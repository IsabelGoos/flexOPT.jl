#!/usr/bin/env julia
#
# Boundary-free, local-recipe diagnostics for the PDE gallery in
# src/motorsOPT/famousEquations.jl.  This is the inexpensive first stage:
# recipes that fail here should not be sent to a periodic convergence sweep.

import Pkg

const FLEXOPT_ROOT = normpath(joinpath(@__DIR__, ".."))
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

using .commonBatchs
using .flexOPT

const EQUATIONS = [
    (
        label="Poisson 1-D",
        equation="1DpoissonHetero",
        space_dimension=1,
        has_time=false,
        branches=1,
    ),
    (
        label="Poisson 2-D",
        equation="2DpoissonHetero",
        space_dimension=2,
        has_time=false,
        branches=1,
    ),
    (
        label="SH frequency 1-D",
        equation="1DsismoFreqHetero",
        space_dimension=1,
        has_time=false,
        branches=1,
    ),
    (
        label="Acoustic time 2-D",
        equation="2DacousticTime",
        space_dimension=2,
        has_time=true,
        branches=1,
    ),
    (
        label="Elastic time 2-D",
        equation="2DsismoTimeIsoHeteroForce",
        space_dimension=2,
        has_time=true,
        branches=2,
    ),
    (
        label="Elastic time 3-D",
        equation="3DsismoTimeIso",
        space_dimension=3,
        has_time=true,
        branches=3,
    ),
]

interpolation(points, offset, order) = (
    ptsSpace=points,
    ptsTime=points,
    offsetSpace=offset,
    offsetTime=offset,
    YorderBspace=order,
    YorderBtime=order,
)

const RECIPES = [
    (
        name="OPT3",
        points=3,
        order_b=1,
        supplementary_order=2,
        interpolation=interpolation(1, 1.0, -1),
        hierarchical=false,
        half_shift_mode=:none,
    ),
    (
        name="OPT5-ordinary-hat-supp0",
        points=5,
        order_b=1,
        supplementary_order=0,
        interpolation=interpolation(1, 2.0, -1),
        hierarchical=false,
        half_shift_mode=:none,
    ),
]

function make_parameters(equation, recipe; delta=1.0)
    dimensions = equation.space_dimension + equation.has_time
    Δ = Tuple(fill(Float64(delta), dimensions))
    points_time = equation.has_time ? recipe.points : 1
    order_time = equation.has_time ? recipe.order_b : 0
    field_itpl = merge(recipe.interpolation, (
        ptsTime=equation.has_time ? recipe.interpolation.ptsTime : 1,
        offsetTime=equation.has_time ? recipe.interpolation.offsetTime : 0.0,
        YorderBtime=equation.has_time ? recipe.interpolation.YorderBtime : -1,
    ))
    return Dict{String,Any}(
        "famousEquationType" => equation.equation,
        "Δ" => Δ,
        "orderBtime" => order_time,
        "orderBspace" => recipe.order_b,
        "pointsInSpace" => recipe.points,
        "pointsInTime" => points_time,
        "supplementaryOrder" => recipe.supplementary_order,
        "fieldItpl" => field_itpl,
        "materItpl" => field_itpl,
        "nuGeometryMode" => :middle,
        "hierarchicalTestFunctions" => recipe.hierarchical,
        "evenOrderHalfShiftMode" => recipe.half_shift_mode,
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
    recipes = quick ? RECIPES[1:1] : RECIPES
    equations = quick ? EQUATIONS[1:3] : EQUATIONS
    rows = NamedTuple[]

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
    end

    output_dir = joinpath(FLEXOPT_ROOT, "scripts", "tmp",
        "famous_equation_interior_diagnostics")
    mkpath(output_dir)
    output_file = joinpath(output_dir, "recipe_construction.jld2")
    jldsave(output_file;
        rows,
        equations=EQUATIONS,
        recipes=RECIPES,
        boundary_mode="none; local interior recipe only",
        diagnostic_stage="construction and Gamma-presence smoke test",
    )
    println("Saved: ", output_file)
    return output_file
end

abspath(PROGRAM_FILE) == abspath(@__FILE__) && main()
