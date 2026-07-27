#!/usr/bin/env julia
#
# Interior consistency diagnostic for the one-dimensional Poisson recipes.
# It deliberately uses only the central local geometry: there is no global
# assembly and therefore no boundary closure in this test.

import Pkg

const FLEXOPT_ROOT = normpath(joinpath(@__DIR__, ".."))
Pkg.activate(FLEXOPT_ROOT)

using KernelAbstractions: CPU
using LinearAlgebra
using Printf
using Symbolics

include(joinpath(FLEXOPT_ROOT, "src", "batchFiles", "batchGPU.jl"))
include(joinpath(FLEXOPT_ROOT, "src", "commonBatchs.jl"))
include(joinpath(FLEXOPT_ROOT, "src", "planet1D.jl"))
include(joinpath(FLEXOPT_ROOT, "src", "GeoPoints.jl"))
include(joinpath(FLEXOPT_ROOT, "src", "flexOPT.jl"))

using .commonBatchs
using .flexOPT

const EQUATION = "1DpoissonHetero"

interpolation(points, offset, order) = (
    ptsSpace=points,
    ptsTime=1,
    offsetSpace=offset,
    offsetTime=1,
    YorderBspace=order,
    YorderBtime=1,
)

function configuration(
    name,
    points;
    order_b=-1,
    supplementary_order=0,
    field_points=1,
    field_offset=(points - 1) / 2,
    material_points=field_points,
    material_offset=field_offset,
    interpolation_order=-1,
    field_interpolation_order=interpolation_order,
    material_interpolation_order=interpolation_order,
    taylor_inverse_mode=:scaled_svd,
    trial_function_ref_points=nothing,
)
    return (
        name=name,
        orderBtime=1,
        orderBspace=order_b,
        pointsInSpace=points,
        pointsInTime=1,
        supplementaryOrder=supplementary_order,
        fieldItpl=interpolation(field_points, field_offset, field_interpolation_order),
        materItpl=interpolation(material_points, material_offset, material_interpolation_order),
        taylorInverseMode=taylor_inverse_mode,
        trialFunctionRefPoints=trial_function_ref_points,
    )
end

function configurations()
    return [
        configuration("FD3", 3),
        configuration("FD4-half", 4),
        configuration("FD5", 5),
        configuration(
            "OPT5-box", 5;
            order_b=-1, supplementary_order=2,
            field_offset=2.0, material_offset=2.0,
        ),
        configuration(
            "OPT5-ordinary-hat-supp0", 5;
            order_b=1, supplementary_order=0,
            field_offset=2.0, material_offset=2.0,
        ),
        configuration(
            "OPT5-wide-hat-1-to-5-supp0", 5;
            order_b=1, supplementary_order=0,
            field_offset=2.0, material_offset=2.0,
            trial_function_ref_points=[1, 3, 5],
        ),
        configuration(
            "OPT5-cubic-supp0", 5;
            order_b=3, supplementary_order=0,
            field_offset=2.0, material_offset=2.0,
        ),
        configuration(
            "OPT5-wide-hat-1-to-5-supp2", 5;
            order_b=1, supplementary_order=2,
            field_offset=2.0, material_offset=2.0,
            trial_function_ref_points=[1, 3, 5],
        ),
        configuration(
            "OPT5-ordinary-hat-supp2-requested", 5;
            order_b=1, supplementary_order=2,
            field_offset=2.0, material_offset=2.0,
        ),
        configuration(
            "OPT5-cubic-supp2-requested", 5;
            order_b=3, supplementary_order=2,
            field_offset=2.0, material_offset=2.0,
        ),
        configuration("OPT3-normal", 3; order_b=1, supplementary_order=2),
        configuration(
            "OPT3-staggered", 3;
            order_b=1, supplementary_order=2,
            field_points=3, field_offset=0.0,
            material_points=4, material_offset=-0.5,
            interpolation_order=1,
        ),
        configuration(
            "OPT4-central-mu", 4;
            order_b=1, supplementary_order=2,
            field_points=1, field_offset=1.5,
            material_points=1, material_offset=1.5,
            interpolation_order=-1,
        ),
        configuration(
            "OPT4-field23-material-centre", 4;
            order_b=1, supplementary_order=2,
            field_points=2, field_offset=1.0,
            material_points=1, material_offset=1.5,
            field_interpolation_order=1,
            material_interpolation_order=-1,
        ),
        configuration(
            "OPT4-field-centre-material23", 4;
            order_b=1, supplementary_order=2,
            field_points=1, field_offset=1.5,
            material_points=2, material_offset=1.0,
            field_interpolation_order=-1,
            material_interpolation_order=1,
        ),
        configuration(
            "OPT4-field23-material23", 4;
            order_b=1, supplementary_order=2,
            field_points=2, field_offset=1.0,
            material_points=2, material_offset=1.0,
            interpolation_order=1,
        ),
        configuration(
            "OPT4-B2-field23-material23", 4;
            order_b=2, supplementary_order=2,
            field_points=2, field_offset=1.0,
            material_points=2, material_offset=1.0,
            interpolation_order=1,
        ),
        configuration(
            "OPT4-B2-field23-material-all4", 4;
            order_b=2, supplementary_order=2,
            field_points=2, field_offset=1.0,
            material_points=4, material_offset=0.0,
            interpolation_order=1,
        ),
        configuration(
            "OPT4-field23-material-all4", 4;
            order_b=1, supplementary_order=2,
            field_points=2, field_offset=1.0,
            material_points=4, material_offset=0.0,
            interpolation_order=1,
        ),
        configuration(
            "OPT4-field23-material-stagger5", 4;
            order_b=1, supplementary_order=2,
            field_points=2, field_offset=1.0,
            material_points=5, material_offset=-0.5,
            interpolation_order=1,
        ),
        configuration(
            "OPT4-collocated", 4;
            order_b=1, supplementary_order=2,
            field_points=4, field_offset=0.0,
            material_points=4, material_offset=0.0,
            interpolation_order=1,
        ),
        configuration(
            "OPT4-staggered", 4;
            order_b=1, supplementary_order=2,
            field_points=4, field_offset=0.0,
            material_points=5, material_offset=-0.5,
            interpolation_order=1,
        ),
        configuration(
            "OPT5-central-mu", 5;
            order_b=1, supplementary_order=2,
            field_points=1, field_offset=2.0,
            material_points=1, material_offset=2.0,
            interpolation_order=-1,
        ),
        configuration(
            "OPT5-central-mu-moore-penrose", 5;
            order_b=1, supplementary_order=2,
            field_points=1, field_offset=2.0,
            material_points=1, material_offset=2.0,
            interpolation_order=-1,
            taylor_inverse_mode=:moore_penrose_svd,
        ),
        configuration(
            "OPT5-B3-central-mu", 5;
            order_b=3, supplementary_order=2,
            field_points=1, field_offset=2.0,
            material_points=1, material_offset=2.0,
            interpolation_order=-1,
        ),
        configuration(
            "OPT5-B3-central-mu-moore-penrose", 5;
            order_b=3, supplementary_order=2,
            field_points=1, field_offset=2.0,
            material_points=1, material_offset=2.0,
            interpolation_order=-1,
            taylor_inverse_mode=:moore_penrose_svd,
        ),
        configuration(
            "OPT5-field234-material-centre", 5;
            order_b=1, supplementary_order=2,
            field_points=3, field_offset=1.0,
            material_points=1, material_offset=2.0,
            field_interpolation_order=1,
            material_interpolation_order=-1,
        ),
        configuration(
            "OPT5-field-centre-material234", 5;
            order_b=1, supplementary_order=2,
            field_points=1, field_offset=2.0,
            material_points=3, material_offset=1.0,
            field_interpolation_order=-1,
            material_interpolation_order=1,
        ),
        configuration(
            "OPT5-field234-material234", 5;
            order_b=1, supplementary_order=2,
            field_points=3, field_offset=1.0,
            material_points=3, material_offset=1.0,
            interpolation_order=1,
        ),
        configuration(
            "OPT5-collocated", 5;
            order_b=1, supplementary_order=2,
            field_points=5, field_offset=0.0,
            material_points=5, material_offset=0.0,
            interpolation_order=1,
        ),
        configuration(
            "OPT5-staggered", 5;
            order_b=1, supplementary_order=2,
            field_points=5, field_offset=0.0,
            material_points=6, material_offset=-0.5,
            interpolation_order=1,
        ),
    ]
end

function make_recipe(config; delta=1.0)
    return makeOPTsemiSymbolic(Dict{String,Any}(
        "famousEquationType" => EQUATION,
        "Δ" => delta,
        "orderBtime" => config.orderBtime,
        "orderBspace" => config.orderBspace,
        "pointsInSpace" => config.pointsInSpace,
        "pointsInTime" => config.pointsInTime,
        "supplementaryOrder" => config.supplementaryOrder,
        "fieldItpl" => config.fieldItpl,
        "materItpl" => config.materItpl,
        "nuGeometryMode" => :middle,
        "taylorInverseMode" => config.taylorInverseMode,
        "trialFunctionRefPoints" => config.trialFunctionRefPoints,
        "recipe_backend" => CPU(),
    ))["recette"]
end

function homogeneous_stencil(symbolic_stencil, material_symbols)
    mapping = Dict{Any,Any}()
    for entry in vec(material_symbols)
        mapping[entry[]] = 1.0
    end
    return Float64[
        Num2Float64(Symbolics.substitute(value, mapping))
        for value in symbolic_stencil
    ]
end

function symbol_series(stencil, offsets, maximum_order)
    return ComplexF64[
        sum(stencil[j] * (im * offsets[j])^order for j in eachindex(stencil)) /
        factorial(order)
        for order in 0:maximum_order
    ]
end

function analyse(config; maximum_order=12, tolerance=1e-10)
    recipe = make_recipe(config)
    geometry = 1
    nodes = recipe.nodes[geometry]
    centre = nodes[recipe.centresIndices[geometry]][1]
    offsets = Float64[node[1] - centre for node in nodes]

    lhs = recipe.lhs
    rhs = recipe.rhs
    a = homogeneous_stencil(lhs.Ajiννᶜ[:, 1, 1, geometry], lhs.varM)
    gamma = homogeneous_stencil(rhs.Γjiννᶜ[:, 1, 1, geometry], rhs.varF)

    a_series = symbol_series(a, offsets, maximum_order)
    gamma_series = symbol_series(gamma, offsets, maximum_order)
    residual = copy(a_series)
    for order in 2:maximum_order
        residual[order + 1] += gamma_series[order - 1]
    end

    scale = max(maximum(abs, a_series), maximum(abs, gamma_series), 1.0)
    first_failure = findfirst(order ->
        abs(residual[order + 1]) > tolerance * scale,
        0:maximum_order)

    println("\n", config.name)
    println("  offsets = ", offsets)
    println("  A       = ", a)
    println("  Gamma   = ", gamma)
    if first_failure === nothing
        println("  no residual coefficient above tolerance through order ", maximum_order)
    else
        println("  first nonzero residual order = ", first_failure - 1)
    end
    println("  residual Taylor coefficients:")
    for order in 0:maximum_order
        value = residual[order + 1]
        @printf("    theta^%-2d  % .8e %+.8ei\n", order, real(value), imag(value))
    end

    return (; name=config.name, offsets, a, gamma, residual)
end

function main()
    selected = filter(config -> !occursin("wide-hat", config.name), configurations())
    results = [analyse(config) for config in selected]
    return results
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
