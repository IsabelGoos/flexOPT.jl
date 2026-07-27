#!/usr/bin/env julia
#
# Operator-level validation against:
#   Geller & Takeuchi (1998), GJI 135, Eq. (18)
#   Takeuchi & Geller (2000), PEPI 119, Eqs. (23)-(25)
#
# Usage:
#   julia --project=.. scripts/compareTakeuchiGellerRecette.jl
#   FLEXOPT_RECETTE_QUICK=1 julia --project=.. scripts/compareTakeuchiGellerRecette.jl
#
# This is deliberately not a wave-propagation benchmark. It asks flexOPT to
# generate each local recette and compares its numerical coefficients with the
# published operator after substituting homogeneous material properties.

import Pkg

const FLEXOPT_ROOT = normpath(joinpath(@__DIR__, ".."))
Pkg.activate(FLEXOPT_ROOT)

using KernelAbstractions: CPU
using Printf
using Symbolics

include(joinpath(FLEXOPT_ROOT, "src", "batchFiles", "batchGPU.jl"))
include(joinpath(FLEXOPT_ROOT, "src", "commonBatchs.jl"))
include(joinpath(FLEXOPT_ROOT, "src", "planet1D.jl"))
include(joinpath(FLEXOPT_ROOT, "src", "flexOPT.jl"))

using .commonBatchs: Num2Float64
using .flexOPT: makeOPTsemiSymbolic

const RTOL = 2e-11
const ATOL = 2e-11
const RECIPE_CACHE = Dict{Tuple{String,Tuple},Any}()
const QUICK = get(ENV, "FLEXOPT_RECETTE_QUICK", "0") == "1"

interpolation(offset) = (
    ptsSpace=1,
    ptsTime=1,
    offsetSpace=offset,
    offsetTime=1,
    YorderBspace=-1,
    YorderBtime=-1,
)

function opt3_parameters(equation, spacings, supplementary_order)
    return Dict{String,Any}(
        "famousEquationType" => equation,
        "Δ" => spacings,
        "orderBtime" => 1,
        "orderBspace" => 1,
        "pointsInSpace" => 3,
        "pointsInTime" => 3,
        "supplementaryOrder" => supplementary_order,
        "fieldItpl" => interpolation(1.0),
        "materItpl" => interpolation(1.0),
        "nuGeometryMode" => :middle,
        "recipe_backend" => CPU(),
    )
end

function generated_recipe(equation, spacings; supplementary_order=2)
    key = (equation, (Tuple(spacings)..., supplementary_order))
    return get!(RECIPE_CACHE, key) do
        makeOPTsemiSymbolic(
            opt3_parameters(equation, spacings, supplementary_order),
        )["recette"]
    end
end

function material_value(symbol, values)
    name = string(symbol)
    occursin("ρ", name) && return values.rho
    occursin("μ", name) && return values.mu
    occursin("v", name) && return values.velocity
    error("No material substitution was supplied for flexOPT symbol $name")
end

function numerical_coefficients(recipe, values; scale=1.0)
    replacements = Dict{Any,Any}()
    for symbol in vec(recipe.lhs.varM)
        raw_symbol = symbol isa Num ? symbol[] : symbol
        replacements[raw_symbol] = material_value(symbol, values)
    end

    geometry = 1
    symbolic = recipe.lhs.Ajiννᶜ[:, 1, 1, geometry]
    coefficients = Float64[
        scale * Num2Float64(Symbolics.substitute(value, replacements))
        for value in symbolic
    ]

    nodes = recipe.nodes[geometry]
    centre = nodes[recipe.centresIndices[geometry]]
    return Dict(
        Tuple(round.(Int, node .- centre)) => coefficient
        for (node, coefficient) in zip(nodes, coefficients)
    )
end

function paper_1998(delta_x, delta_t, rho, mu)
    result = Dict{Tuple{Int,Int},Float64}()
    smoothing = (1.0, 10.0, 1.0)
    second = (1.0, -2.0, 1.0)
    for (ix, x) in enumerate(-1:1), (it, t) in enumerate(-1:1)
        mass = rho / delta_t^2 * smoothing[ix] * second[it] / 12
        stiffness = mu / delta_x^2 * second[ix] * smoothing[it] / 12
        result[(x, t)] = mass - stiffness
    end
    return result
end

function paper_2000(delta_x, delta_z, delta_t, rho, mu)
    result = Dict{Tuple{Int,Int,Int},Float64}()
    smoothing = (1.0, 10.0, 1.0)
    second = (1.0, -2.0, 1.0)
    for (ix, x) in enumerate(-1:1),
        (iz, z) in enumerate(-1:1),
        (it, t) in enumerate(-1:1)

        mass = rho / delta_t^2 *
               smoothing[ix] * smoothing[iz] * second[it] / 144
        stiffness_x = mu / delta_x^2 *
                      second[ix] * smoothing[iz] * smoothing[it] / 144
        stiffness_z = mu / delta_z^2 *
                      smoothing[ix] * second[iz] * smoothing[it] / 144
        result[(x, z, t)] = mass - stiffness_x - stiffness_z
    end
    return result
end

function compare_coefficients(label, generated, published)
    generated_keys = Set(keys(generated))
    published_keys = Set(keys(published))
    generated_keys == published_keys || error(
        "$label has different supports.\n" *
        "  only flexOPT: $(sort!(collect(setdiff(generated_keys, published_keys))))\n" *
        "  only paper:   $(sort!(collect(setdiff(published_keys, generated_keys))))",
    )

    maximum_error = 0.0
    maximum_scale = 0.0
    worst_offset = first(keys(published))
    for offset in keys(published)
        error_at_offset = abs(generated[offset] - published[offset])
        if error_at_offset > maximum_error
            maximum_error = error_at_offset
            worst_offset = offset
        end
        maximum_scale = max(
            maximum_scale,
            abs(generated[offset]),
            abs(published[offset]),
        )
    end

    tolerance = ATOL + RTOL * maximum_scale
    passed = maximum_error <= tolerance
    @printf(
        "%-30s  max|Δc| = %.4e at %-12s  %s\n",
        label,
        maximum_error,
        string(worst_offset),
        passed ? "PASS" : "FAIL",
    )
    passed || error(
        "$label disagrees with the published recette: " *
        "maximum error $maximum_error exceeds $tolerance; " *
        "flexOPT=$(generated[worst_offset]), " *
        "paper=$(published[worst_offset])",
    )
    return maximum_error
end

function validate_1998(delta_x, delta_t, rho, mu)
    delta_x == delta_t || error(
        "The current dimensionless recette supports one common grid scale; " *
        "got Δx=$delta_x and Δt=$delta_t",
    )
    recipe = generated_recipe("1DsismoTimeHomo", (delta_x, delta_t))
    values = (rho=rho, mu=mu, velocity=sqrt(mu / rho))
    generated = numerical_coefficients(recipe, values; scale=inv(delta_x^2))
    published = paper_1998(delta_x, delta_t, rho, mu)
    return compare_coefficients(
        "1998 Eq. (18), Δ=($delta_x,$delta_t)",
        generated,
        published,
    )
end

function validate_2000(delta_x, delta_z, delta_t, rho, mu)
    # flexOPT's homogeneous acoustic equation is the scalar SH equation divided
    # by rho: u_tt - v^2 (u_xx + u_zz) = 0. With v^2=mu/rho,
    # multiplying its generated coefficients by rho gives A'-Kx'-Kz'.
    delta_x == delta_z == delta_t || error(
        "The current dimensionless recette supports one common grid scale; " *
        "got Δx=$delta_x, Δz=$delta_z, and Δt=$delta_t",
    )
    recipe = generated_recipe(
        "2DacousticHomoTime",
        (delta_x, delta_z, delta_t),
        supplementary_order=0,
    )
    values = (rho=rho, mu=mu, velocity=sqrt(mu / rho))
    generated = numerical_coefficients(
        recipe,
        values;
        scale=rho / delta_x^2,
    )
    published = paper_2000(delta_x, delta_z, delta_t, rho, mu)
    return compare_coefficients(
        "2000 Eqs. (23)-(25), Δ=($delta_x,$delta_z,$delta_t)",
        generated,
        published,
    )
end

function main()
    println("Takeuchi-Geller local-recette validation")
    println("Comparing combined interior coefficients; no time marching is used.\n")

    material_pairs = QUICK ? (
        (rho=1.0, mu=1.0),
    ) : (
        (rho=1.0, mu=1.0),
        (rho=2.7, mu=31.0),
    )
    spacing_pairs_1d = QUICK ? (
        (delta_x=1.0, delta_t=1.0),
    ) : (
        (delta_x=1.0, delta_t=1.0),
        (delta_x=0.1, delta_t=0.1),
    )
    spacing_pairs_2d = QUICK ? (
        (delta_x=1.0, delta_z=1.0, delta_t=1.0),
    ) : (
        (delta_x=1.0, delta_z=1.0, delta_t=1.0),
        (delta_x=0.1, delta_z=0.1, delta_t=0.1),
    )

    for material in material_pairs
        @printf("ρ = %.6g, μ = %.6g\n", material.rho, material.mu)
        for spacing in spacing_pairs_1d
            validate_1998(
                spacing.delta_x,
                spacing.delta_t,
                material.rho,
                material.mu,
            )
        end
        for spacing in spacing_pairs_2d
            validate_2000(
                spacing.delta_x,
                spacing.delta_z,
                spacing.delta_t,
                material.rho,
                material.mu,
            )
        end
        println()
    end

    println("All generated recettes agree with the published operators.")
    return nothing
end

main()
