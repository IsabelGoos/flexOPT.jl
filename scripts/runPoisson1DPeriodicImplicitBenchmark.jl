#!/usr/bin/env julia
#
# Periodic Poisson benchmark for node-centred and implicit midpoint recipes.
# Periodicity gives one equation per midpoint without introducing an
# artificial Dirichlet closure for even-point stencils.

import Pkg

const FLEXOPT_ROOT = normpath(joinpath(@__DIR__, ".."))
Pkg.activate(FLEXOPT_ROOT)

using KernelAbstractions: CPU
using JLD2
using LinearAlgebra
using Printf
using SparseArrays
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
    offsetSpace=offset, offsetTime=1,
    YorderBspace=order, YorderBtime=1,
)

function configuration(
    name, points;
    family="OPT",
    mu_description="",
    order_b=points - 2,
    supplementary_order=2,
    field_points=1,
    field_offset=(points - 1) / 2,
    material_points=field_points,
    material_offset=field_offset,
    field_order=-1,
    material_order=field_order,
    taylor_inverse_mode=:scaled_svd,
)
    return (; name, family, mu_description, points, order_b, supplementary_order,
        taylor_inverse_mode,
        fieldItpl=interpolation(field_points, field_offset, field_order),
        materItpl=interpolation(material_points, material_offset, material_order))
end

function configurations()
    return [
        configuration("FD3", 3; family="FD", mu_description="central",
            order_b=-1, supplementary_order=0, field_offset=1.0),
        configuration("FD4-half", 4; family="FD", mu_description="single half-centre",
            order_b=-1, supplementary_order=0, field_offset=1.5),
        configuration("FD5", 5; family="FD", mu_description="central",
            order_b=-1, supplementary_order=0, field_offset=2.0),
        configuration("FD6-half", 6; family="FD", mu_description="single half-centre",
            order_b=-1, supplementary_order=0, field_offset=2.5),

        configuration("OPT3-central", 3; mu_description="field centre; material centre",
            field_offset=1.0),
        configuration("OPT3-central-material-all3", 3;
            mu_description="field centre; material nodes 1:3",
            field_offset=1.0, material_points=3, material_offset=0.0,
            field_order=-1, material_order=1),
        configuration("OPT3-allmu", 3; mu_description="field/material nodes 1:3",
            field_points=3, field_offset=0.0,
            material_points=3, material_offset=0.0,
            field_order=1, material_order=1),

        configuration("OPT4-field23-material23", 4;
            mu_description="implicit field nodes 2:3; material nodes 2:3",
            field_points=2, field_offset=1.0,
            material_points=2, material_offset=1.0,
            field_order=1, material_order=1),
        configuration("OPT4-field23-material-all4", 4;
            mu_description="implicit field nodes 2:3; material nodes 1:4",
            field_points=2, field_offset=1.0,
            material_points=4, material_offset=0.0,
            field_order=1, material_order=1),
        configuration("OPT4-field23-material-stagger5", 4;
            mu_description="implicit field nodes 2:3; material half-nodes 0.5:4.5",
            field_points=2, field_offset=1.0,
            material_points=5, material_offset=-0.5,
            field_order=1, material_order=1),

        configuration("OPT5-central", 5; mu_description="field centre; material centre",
            field_offset=2.0),
        configuration("OPT5-central-moore-penrose", 5;
            mu_description="field centre; material centre; direct Moore-Penrose SVD",
            field_offset=2.0, taylor_inverse_mode=:moore_penrose_svd),
        configuration("OPT5-central-material234", 5;
            mu_description="field centre; material nodes 2:4",
            field_offset=2.0, material_points=3, material_offset=1.0,
            field_order=-1, material_order=1),
        configuration("OPT5-central-material-all5", 5;
            mu_description="field centre; material nodes 1:5",
            field_offset=2.0, material_points=5, material_offset=0.0,
            field_order=-1, material_order=1),
        configuration("OPT5-allmu", 5; mu_description="field/material nodes 1:5",
            field_points=5, field_offset=0.0,
            material_points=5, material_offset=0.0,
            field_order=1, material_order=1),

        configuration("OPT6-field34-material34", 6;
            mu_description="implicit field nodes 3:4; material nodes 3:4",
            field_points=2, field_offset=2.0,
            material_points=2, material_offset=2.0,
            field_order=1, material_order=1),
        configuration("OPT6-field34-material-all6", 6;
            mu_description="implicit field nodes 3:4; material nodes 1:6",
            field_points=2, field_offset=2.0,
            material_points=6, material_offset=0.0,
            field_order=1, material_order=1),
        configuration("OPT6-field34-material-stagger7", 6;
            mu_description="implicit field nodes 3:4; material half-nodes 0.5:6.5",
            field_points=2, field_offset=2.0,
            material_points=7, material_offset=-0.5,
            field_order=1, material_order=1),
    ]
end

function make_recipe(config; delta=1.0)
    return makeOPTsemiSymbolic(Dict{String,Any}(
        "famousEquationType" => "1DpoissonHetero",
        # Build once in dimensionless local coordinates. The physical source
        # is multiplied by h^2 below.
        "Δ" => delta,
        "orderBtime" => 1,
        "orderBspace" => config.order_b,
        "pointsInSpace" => config.points,
        "pointsInTime" => 1,
        "supplementaryOrder" => config.supplementary_order,
        "fieldItpl" => config.fieldItpl,
        "materItpl" => config.materItpl,
        "nuGeometryMode" => :middle,
        "taylorInverseMode" => config.taylor_inverse_mode,
        "recipe_backend" => CPU(),
    ))["recette"]
end

function relative_scaling_error(target, reference, scale)
    denominator = max(norm(target), eps(Float64))
    return norm(target - scale .* reference) / denominator
end

function delta_scaling_check()
    selected_names = [
        "FD5",
        "OPT4-field23-material-all4",
        "OPT5-central",
        "OPT6-field34-material-all6",
    ]
    selected = filter(config -> config.name in selected_names, configurations())
    deltas = [0.5, 2pi / 52]
    println("Checking full recipe scaling relative to Delta=1")
    for config in selected
        reference = prepare_numeric_recipe(make_recipe(config; delta=1.0))
        println("\n", config.name)
        for delta in deltas
            current = prepare_numeric_recipe(make_recipe(config; delta=delta))
            scale_a = dot(vec(current.material_tensor), vec(reference.material_tensor)) /
                      dot(vec(reference.material_tensor), vec(reference.material_tensor))
            scale_gamma = dot(current.gamma, reference.gamma) /
                          dot(reference.gamma, reference.gamma)
            error_a = relative_scaling_error(
                current.material_tensor, reference.material_tensor, scale_a)
            error_gamma = relative_scaling_error(
                current.gamma, reference.gamma, scale_gamma)
            ratio_error = abs(scale_gamma / scale_a - delta^2) / delta^2
            @printf(
                "  Delta=%g  sA=% .8e  sGamma=% .8e  relA=%.3e  relGamma=%.3e  ratio-error=%.3e\n",
                delta, scale_a, scale_gamma, error_a, error_gamma, ratio_error,
            )
        end
    end
end

periodic_index(index, n) = mod1(index, n)

function numeric_value(expression, mapping)
    return Num2Float64(Symbolics.substitute(expression, mapping))
end

function prepare_numeric_recipe(recipe)
    geometry = 1
    nodes = recipe.nodes[geometry]
    centre = nodes[recipe.centresIndices[geometry]][1]
    offsets = Int[node[1] - centre for node in nodes]
    lhs = recipe.lhs
    rhs = recipe.rhs
    number_points = length(nodes)

    # The Poisson operator is linear in kappa. Precompute its response to
    # every local material basis vector so the large parameter sweep does not
    # perform Symbolics substitutions inside every grid row.
    material_tensor = zeros(Float64, number_points, number_points)
    for material_point in 1:number_points
        mapping = Dict{Any,Any}(
            lhs.varM[1, point][] => (point == material_point ? 1.0 : 0.0)
            for point in 1:number_points
        )
        for field_point in 1:number_points
            material_tensor[field_point, material_point] = numeric_value(
                lhs.Ajiννᶜ[field_point, 1, 1, geometry], mapping)
        end
    end

    force_mapping = Dict{Any,Any}(entry[] => 1.0 for entry in vec(rhs.varF))
    gamma = Float64[
        numeric_value(rhs.Γjiννᶜ[field_point, 1, 1, geometry], force_mapping)
        for field_point in 1:number_points
    ]
    return (; offsets, material_tensor, gamma)
end

function assemble_periodic(prepared, beta)
    n = length(beta)
    offsets = prepared.offsets
    number_points = length(offsets)
    A = spzeros(Float64, n, n)
    Gamma = spzeros(Float64, n, n)

    for row in 1:n
        local_beta = Float64[
            beta[periodic_index(row + offsets[point], n)]
            for point in 1:number_points
        ]
        a = prepared.material_tensor * local_beta
        for field_point in 1:number_points
            column = periodic_index(row + offsets[field_point], n)
            A[row, column] += a[field_point]
            Gamma[row, column] += prepared.gamma[field_point]
        end
    end
    return A, Gamma
end

function benchmark_cases()
    cases = NamedTuple[]
    push!(cases, (
        name="homogeneous",
        kT=2, phiT=0.0, kappa0=1.0, amplitude=0.0,
        kKappa=0, phiKappa=0.0,
        period_ratio=Inf, phase_shift=0.0,
    ))
    for (label, kKappa, phiKappa) in [
        ("same_period_phase0", 2, 0.0),
        ("same_period_phase_pi4", 2, pi / 4),
        ("same_period_phase_pi2", 2, pi / 2),
        ("material_twice_period_phase0", 1, 0.0),
        ("material_twice_period_phase_pi2", 1, pi / 2),
        ("material_half_period_phase0", 4, 0.0),
        ("material_half_period_phase_pi2", 4, pi / 2),
        ("material_quarter_period_phase_pi4", 8, pi / 4),
    ]
        push!(cases, (
            name=label,
            kT=2, phiT=0.0, kappa0=2.0, amplitude=0.35,
            kKappa=kKappa, phiKappa=phiKappa,
            period_ratio=2 / kKappa,
            phase_shift=phiKappa,
        ))
    end
    return cases
end

function manufactured_fields(case, x)
    theta_t = case.kT .* x .+ case.phiT
    exact = cos.(theta_t)
    if case.amplitude == 0
        beta = fill(case.kappa0, length(x))
        force = -case.kappa0 .* case.kT^2 .* cos.(theta_t)
    else
        theta_kappa = case.kKappa .* x .+ case.phiKappa
        beta = case.kappa0 .+ case.amplitude .* cos.(theta_kappa)
        # f = d/dx(kappa*dT/dx) = kappa'*T' + kappa*T''.
        force = case.amplitude * case.kKappa * case.kT .*
                sin.(theta_kappa) .* sin.(theta_t) .-
                beta .* case.kT^2 .* cos.(theta_t)
    end
    return exact, beta, force
end

function solve_periodic(prepared, n, case; source_scale=1.0)
    length_domain = 2pi
    h = length_domain / n
    x = collect(0:n-1) .* h
    exact, beta, force = manufactured_fields(case, x)

    A, Gamma = assemble_periodic(prepared, beta)
    b = Gamma * (source_scale .* force)

    # Periodic Poisson has a constant null vector. Even-point recipes may
    # possess an additional checkerboard null mode; retain their very large
    # errors/NaNs as a diagnostic instead of silently regularising them.
    A[end, :] .= 1 / n
    b[end] = sum(exact) / n
    solution = try
        A \ b
    catch exception
        if exception isa LinearAlgebra.SingularException
            return NaN
        end
        rethrow()
    end
    return norm(solution - exact) / sqrt(n)
end

function observed_orders(errors, dx)
    orders = fill(NaN, size(errors))
    for scheme in axes(errors, 3), case_index in axes(errors, 2)
        for grid_index in 2:size(errors, 1)
            coarse = errors[grid_index - 1, case_index, scheme]
            fine = errors[grid_index, case_index, scheme]
            if isfinite(coarse) && isfinite(fine) && coarse > 0 && fine > 0
                orders[grid_index, case_index, scheme] =
                    log(coarse / fine) /
                    log(dx[grid_index - 1] / dx[grid_index])
            end
        end
    end
    return orders
end

function main()
    sizes = [16, 24, 32, 48, 64, 96, 128, 192, 256, 384, 512]
    dx = 2pi ./ sizes
    configs = configurations()
    cases = benchmark_cases()

    # Faithful production mode: reconstruct C^l_eta, A and Gamma at each
    # actual Delta. This intentionally does not reuse a Delta=1 recipe.
    prepared = Matrix{Any}(undef, length(sizes), length(configs))
    for (grid_index, delta) in pairs(dx), (scheme_index, config) in pairs(configs)
        @info "Preparing per-Delta recipe" config=config.name delta
        prepared[grid_index, scheme_index] =
            prepare_numeric_recipe(make_recipe(config; delta=delta))
    end

    errors = fill(NaN, length(sizes), length(cases), length(configs))
    for (scheme_index, config) in pairs(configs)
        @info "Periodic convergence scheme" config=config.name
        for (case_index, case) in pairs(cases), (grid_index, n) in pairs(sizes)
            errors[grid_index, case_index, scheme_index] =
                solve_periodic(prepared[grid_index, scheme_index], n, case)
        end
    end
    orders = observed_orders(errors, dx)

    scheme_names = getproperty.(configs, :name)
    case_names = getproperty.(cases, :name)
    formula_T = "T(x) = cos(kT*x + phiT)"
    formula_kappa = "kappa(x) = kappa0 + amplitude*cos(kKappa*x + phiKappa)"
    formula_force = "f(x) = amplitude*kKappa*kT*sin(kKappa*x+phiKappa)*sin(kT*x+phiT) - kappa(x)*kT^2*cos(kT*x+phiT)"
    reference_orders = [2.0, 4.0, 6.0]

    output_dir = joinpath(FLEXOPT_ROOT, "scripts", "tmp", "poisson1d_periodic")
    mkpath(output_dir)
    output_file = joinpath(
        output_dir, "poisson1d_periodic_convergence_per_delta.jld2")
    jldsave(output_file;
        sizes, dx, errors, orders,
        scheme_names, case_names,
        scheme_metadata=configs,
        case_metadata=cases,
        formula_T, formula_kappa, formula_force,
        reference_orders,
        domain_length=2pi,
        error_definition="norm(T_numerical-T_exact)/sqrt(N)",
        array_layout="errors[grid_index, case_index, scheme_index]",
        recipe_delta_mode="A, Gamma and C^l_eta reconstructed at every stored dx",
    )

    println("\nSaved periodic convergence benchmark: ", output_file)
    println("Array layout: errors[grid, case, scheme]")
    for (scheme_index, config) in pairs(configs)
        finite_orders = filter(isfinite, vec(orders[:, :, scheme_index]))
        median_order = isempty(finite_orders) ? NaN : median(finite_orders)
        @printf("%-42s median observed order %8.3f\n",
            config.name, median_order)
    end
    return output_file
end

if get(ENV, "FLEXOPT_DELTA_SCALING_CHECK", "0") == "1"
    delta_scaling_check()
else
    main()
end
