#!/usr/bin/env julia

import Pkg
const FLEXOPT_ROOT = normpath(joinpath(@__DIR__, ".."))
Pkg.activate(FLEXOPT_ROOT)

using JLD2
using KernelAbstractions: CPU
using LinearAlgebra
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

interpolation(points, offset, order) = (
    ptsSpace=points, ptsTime=1,
    offsetSpace=offset, offsetTime=0.0,
    YorderBspace=order, YorderBtime=-1,
)

const SCHEMES = [
    (
        name="FD3", points=3, order_b=-1, supplementary_order=0,
        interpolation=interpolation(1, 1.0, -1),
    ),
    (
        name="OPT3", points=3, order_b=1, supplementary_order=2,
        interpolation=interpolation(1, 1.0, -1),
    ),
    (
        name="FD5", points=5, order_b=-1, supplementary_order=0,
        interpolation=interpolation(1, 2.0, -1),
    ),
    (
        name="OPT5-ordinary-hat-supp0", points=5, order_b=1,
        supplementary_order=0,
        interpolation=interpolation(1, 2.0, -1),
    ),
]

const CASES = [
    (
        name="homogeneous",
        field_wave=(2, 0),
        material_wave=(0, 0),
        phase=0.0,
        kappa0=2.0,
        amplitude=0.0,
    ),
    (
        name="same_wave_phase0",
        field_wave=(2, 1),
        material_wave=(2, 1),
        phase=0.0,
        kappa0=2.0,
        amplitude=0.35,
    ),
    (
        name="same_wave_phase_pi2",
        field_wave=(2, 1),
        material_wave=(2, 1),
        phase=pi / 2,
        kappa0=2.0,
        amplitude=0.35,
    ),
    (
        name="short_material_phase0",
        field_wave=(2, 1),
        material_wave=(6, 3),
        phase=0.0,
        kappa0=2.0,
        amplitude=0.35,
    ),
    (
        name="short_material_phase_pi2",
        field_wave=(2, 1),
        material_wave=(6, 3),
        phase=pi / 2,
        kappa0=2.0,
        amplitude=0.35,
    ),
]

function make_recipe(dimension, scheme, h)
    equation = dimension == 1 ? "1DpoissonHetero" : "2DpoissonHetero"
    Δ = dimension == 1 ? h : ntuple(_ -> h, dimension)
    return makeOPTsemiSymbolic(Dict{String,Any}(
        "famousEquationType" => equation,
        "Δ" => Δ,
        "orderBtime" => 0,
        "orderBspace" => scheme.order_b,
        "pointsInSpace" => scheme.points,
        "pointsInTime" => 1,
        "supplementaryOrder" => scheme.supplementary_order,
        "fieldItpl" => scheme.interpolation,
        "materItpl" => scheme.interpolation,
        "nuGeometryMode" => :middle,
        "recipe_backend" => CPU(),
    ))["recette"]
end

function waves(case, dimension)
    k = collect(case.field_wave[1:dimension])
    q = collect(case.material_wave[1:dimension])
    return Float64.(k), Float64.(q)
end

function manufactured(case, dimension, coordinates)
    k, q = waves(case, dimension)
    θt = dot(k, coordinates)
    θκ = dot(q, coordinates) + case.phase
    field = cos(θt)
    kappa = case.kappa0 + case.amplitude * cos(θκ)
    force = case.amplitude * dot(q, k) * sin(θκ) * sin(θt) -
        kappa * dot(k, k) * cos(θt)
    return field, kappa, force
end

function coefficient_value(expression, mapping)
    substituted = Symbolics.substitute(expression, mapping)
    return Float64(Symbolics.value(substituted))
end

function periodic_point(index, offset, n)
    return ntuple(d -> mod1(index[d] + offset[d], n), length(index))
end

function point_coordinates(index, h)
    return Float64[(value - 1) * h for value in index]
end

function local_offsets(recipe, geometry=1)
    nodes = vec(recipe.nodes[geometry])
    centre = nodes[recipe.centresIndices[geometry]]
    return [Tuple(node - centre) for node in nodes]
end

function material_mapping(varM, offsets, index, n, h, case, dimension)
    mapping = Dict{Any,Any}()
    for local_index in eachindex(offsets)
        point = periodic_point(index, offsets[local_index], n)
        _, kappa, _ = manufactured(
            case, dimension, point_coordinates(point, h))
        for variable_index in axes(varM, 1)
            symbol = varM[variable_index, local_index]
            text = string(symbol)
            if occursin("κ", text)
                mapping[symbol] = kappa
            elseif text == "c" || occursin("c(", text)
                mapping[symbol] = 1.0
            else
                mapping[symbol] = 1.0
            end
        end
    end
    return mapping
end

function residual_error(recipe, case, dimension, n)
    h = 2pi / n
    offsets = local_offsets(recipe)
    lhs = recipe.lhs
    rhs = recipe.rhs
    geometry = 1
    residuals = Float64[]
    force_values = Float64[]

    for index in CartesianIndices(ntuple(_ -> n, dimension))
        index_tuple = Tuple(index)
        mapping = material_mapping(
            lhs.varM, offsets, index_tuple, n, h, case, dimension)
        discrete_left = 0.0
        discrete_right = 0.0
        for local_index in eachindex(offsets)
            point = periodic_point(
                index_tuple, offsets[local_index], n)
            field, _, force = manufactured(
                case, dimension, point_coordinates(point, h))
            a = coefficient_value(
                lhs.Ajiννᶜ[local_index, 1, 1, geometry], mapping)
            gamma = coefficient_value(
                rhs.Γjiννᶜ[local_index, 1, 1, geometry],
                Dict{Any,Any}())
            discrete_left += a * field
            discrete_right += gamma * force
        end
        push!(residuals, discrete_left - discrete_right)
        _, _, centre_force = manufactured(
            case, dimension, point_coordinates(index_tuple, h))
        push!(force_values, centre_force)
    end
    absolute = norm(residuals) / sqrt(length(residuals))
    relative = absolute /
        max(norm(force_values) / sqrt(length(force_values)), eps(Float64))
    return absolute, relative
end

function observed_orders(errors, sizes)
    orders = fill(NaN, size(errors))
    for scheme in axes(errors, 3), case in axes(errors, 2)
        for grid in 2:size(errors, 1)
            orders[grid, case, scheme] =
                log(errors[grid - 1, case, scheme] /
                    errors[grid, case, scheme]) /
                log(sizes[grid] / sizes[grid - 1])
        end
    end
    return orders
end

function run_dimension(dimension)
    sizes = dimension == 1 ?
        [16, 24, 32, 48, 64, 96, 128, 192] :
        [8, 12, 16, 24, 32, 48, 64]
    absolute_errors =
        fill(NaN, length(sizes), length(CASES), length(SCHEMES))
    relative_errors = similar(absolute_errors)
    recipe_seconds = fill(NaN, length(sizes), length(SCHEMES))
    residual_seconds = fill(NaN, size(absolute_errors))

    for (grid, n) in pairs(sizes), (scheme_index, scheme) in pairs(SCHEMES)
        h = 2pi / n
        @info "Poisson FD/OPT recipe" dimension n scheme=scheme.name
        recipe = nothing
        recipe_seconds[grid, scheme_index] = @elapsed begin
            recipe = make_recipe(dimension, scheme, h)
        end
        for (case_index, case) in pairs(CASES)
            residual_seconds[grid, case_index, scheme_index] = @elapsed begin
                absolute_errors[grid, case_index, scheme_index],
                relative_errors[grid, case_index, scheme_index] =
                    residual_error(recipe, case, dimension, n)
            end
        end
    end

    absolute_orders = observed_orders(absolute_errors, sizes)
    relative_orders = observed_orders(relative_errors, sizes)
    return (; dimension, sizes, absolute_errors, relative_errors,
        absolute_orders, relative_orders, recipe_seconds, residual_seconds)
end

function main()
    dimensions = get(ENV, "FLEXOPT_POISSON_DIMENSIONS", "1,2")
    selected = parse.(Int, split(dimensions, ","))
    results = Dict("dimension_$dimension" => run_dimension(dimension)
        for dimension in selected)
    output_dir = joinpath(FLEXOPT_ROOT, "scripts", "tmp",
        "famous_equation_periodic_benchmarks")
    mkpath(output_dir)
    dimension_label = join(string.(selected), "d_") * "d"
    output_file = joinpath(output_dir,
        "poisson_fd_vs_opt_$(dimension_label).jld2")
    jldsave(output_file;
        results,
        schemes=SCHEMES,
        cases=CASES,
        error_definition=(
            absolute="RMS(A*T-Gamma*f)",
            relative="RMS(A*T-Gamma*f)/RMS(f)",
        ),
        boundary_mode="periodic local wrapping; no boundary closure",
    )
    for dimension in selected
        result = results["dimension_$dimension"]
        println("\nPoisson ", dimension, "-D")
        for (scheme_index, scheme) in pairs(SCHEMES)
            finite_orders = filter(isfinite,
                vec(result.relative_orders[:, :, scheme_index]))
            @printf("%-32s median residual order %7.3f\n",
                scheme.name, median(finite_orders))
        end
    end
    println("Saved: ", output_file)
    return output_file
end

abspath(PROGRAM_FILE) == abspath(@__FILE__) && main()
