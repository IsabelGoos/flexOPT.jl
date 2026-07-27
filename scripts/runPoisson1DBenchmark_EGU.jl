#!/usr/bin/env julia
#
# Usage:
#   julia --project=.. -t auto scripts/runPoisson1DBenchmark_EGU.jl
#
# A short smoke test can be requested with:
#   FLEXOPT_BENCHMARK_QUICK=1 julia --project=.. -t auto scripts/runPoisson1DBenchmark_EGU.jl
#
# This benchmark follows manuscript equations (18) and (54):
#   * fieldItpl/materItpl select the Y_mu interpolation in equation (18);
#   * makeOPTsemiSymbolic constructs A in equation (54), including the two
#     Taylor reconstruction coefficients and the W*Y*Y*K*K integrals.

import Pkg

const FLEXOPT_ROOT = normpath(joinpath(@__DIR__, ".."))
Pkg.activate(FLEXOPT_ROOT)

using Base.Threads
using CairoMakie
using JLD2
using KernelAbstractions: CPU
using LinearAlgebra
using SparseArrays
using Symbolics

include(joinpath(FLEXOPT_ROOT, "src", "batchFiles", "batchGPU.jl"))
include(joinpath(FLEXOPT_ROOT, "src", "commonBatchs.jl"))
include(joinpath(FLEXOPT_ROOT, "src", "planet1D.jl"))
include(joinpath(FLEXOPT_ROOT, "src", "GeoPoints.jl"))
include(joinpath(FLEXOPT_ROOT, "src", "flexOPT.jl"))

using .commonBatchs
using .flexOPT
using .GeoPoints
using .planet1D

CairoMakie.activate!()

const EQUATION = "1DpoissonHetero"
const REPRESENTATION = "matrixfree"
const RECIPE_BACKEND = CPU()
const QUICK = get(ENV, "FLEXOPT_BENCHMARK_QUICK", "0") == "1"

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
    points_in_space;
    order_b_space=1,
    supplementary_order=max(points_in_space - 1, 2),
    field_points=points_in_space,
    field_offset=0.0,
    material_points=field_points,
    material_offset=field_offset,
    interpolation_order=1,
)
    return (
        name=name,
        orderBtime=1,
        orderBspace=order_b_space,
        pointsInSpace=points_in_space,
        pointsInTime=1,
        supplementaryOrder=supplementary_order,
        fieldItpl=interpolation(field_points, field_offset, interpolation_order),
        materItpl=interpolation(material_points, material_offset, interpolation_order),
    )
end

function benchmark_configurations()
    configs = [
        # Conventional boxcar references. FD3 is always retained as baseline.
        configuration(
            "FD3", 3;
            order_b_space=-1,
            supplementary_order=0,
            field_points=1,
            field_offset=1.0,
            material_points=1,
            material_offset=1.0,
            interpolation_order=-1,
        ),
        configuration(
            "FD4", 4;
            order_b_space=-1,
            supplementary_order=0,
            field_points=1,
            field_offset=1.5,
            material_points=1,
            material_offset=1.5,
            interpolation_order=-1,
        ),
        configuration(
            "FD5", 5;
            order_b_space=-1,
            supplementary_order=0,
            field_points=1,
            field_offset=2.0,
            material_points=1,
            material_offset=2.0,
            interpolation_order=-1,
        ),

        # Proven seismic-style OPT3 baseline: hat W and one central expansion.
        configuration(
            "OPT3-normal", 3;
            supplementary_order=2,
            field_points=1,
            field_offset=1.0,
            material_points=1,
            material_offset=1.0,
            interpolation_order=-1,
        ),

        # Equation-(18) interpolation variants.
        configuration(
            "OPT3-staggered", 3;
            supplementary_order=2,
            field_points=3,
            field_offset=0.0,
            material_points=4,
            material_offset=-0.5,
            interpolation_order=1,
        ),
        # Five μ centres over the unchanged three-node trial stencil:
        # 1, 1.5, 2, 2.5, 3.
        configuration(
            "OPT3-dense-mu5", 3;
            supplementary_order=2,
            field_points=5,
            field_offset=0.0,
            material_points=5,
            material_offset=0.0,
            interpolation_order=1,
        ),
        configuration("OPT4-collocated", 4; supplementary_order=2),
        configuration(
            "OPT4-staggered", 4;
            supplementary_order=2,
            field_points=4,
            field_offset=0.0,
            material_points=5,
            material_offset=-0.5,
            interpolation_order=1,
        ),
        configuration("OPT5-collocated", 5; supplementary_order=2),
        configuration(
            "OPT5-staggered", 5;
            supplementary_order=2,
            field_points=5,
            field_offset=0.0,
            material_points=6,
            material_offset=-0.5,
            interpolation_order=1,
        ),
    ]
    return configs
end

function benchmark_cases(x)
    ∂x = Differential(x)
    raw_cases = [
        (name="homogeneous", u=cos(x), beta=1.0),
        (name="same_lambda", u=cos(x), beta=sin(x) + 2),
        (name="twice_lambda", u=cos(x), beta=sin(x / 2) + 2),
        (name="shifted_pi_3", u=cos(x), beta=sin(x + pi / 3) + 2),
        (name="lambda_2", u=cos(x), beta=cos(x)^2 + 1),
        (name="quadratic", u=cos(x), beta=x^2 + 1),
    ]
    cases = [
        merge(case, (; force=mySimplify(∂x(mySimplify(case.beta * ∂x(case.u))))))
        for case in raw_cases
    ]
    # The second quick case is required to exercise material interpolation:
    # for the homogeneous first case, collocated and half-shifted beta centres
    # are mathematically indistinguishable.
    return QUICK ? cases[1:2] : cases
end

evaluate_expression(expr, x, grid) =
    Float64[Symbolics.value(substitute(expr, Dict(x => coordinate))) for coordinate in grid]

function grid_and_model(case, x, log_h_inverse, L)
    requested_dx = exp(-log_h_inverse)
    nx = max(3, floor(Int, L / requested_dx) + 1)
    dx = L / (nx - 1)
    grid = collect(range(0.0, L; length=nx))
    beta = evaluate_expression(case.beta, x, grid)
    exact = evaluate_expression(case.u, x, grid)
    force = evaluate_expression(case.force, x, grid)
    return (; nx, dx, grid, beta, exact, force)
end

function make_recipe(config, dx)
    params = Dict{String,Any}(
        "famousEquationType" => EQUATION,
        "Δ" => dx,
        "orderBtime" => config.orderBtime,
        "orderBspace" => config.orderBspace,
        "pointsInSpace" => config.pointsInSpace,
        "pointsInTime" => config.pointsInTime,
        "supplementaryOrder" => config.supplementaryOrder,
        "fieldItpl" => config.fieldItpl,
        "materItpl" => config.materItpl,
        "recipe_backend" => RECIPE_BACKEND,
    )
    return makeOPTsemiSymbolic(params)
end

"""
Solve A*u = Gamma*f with exact Dirichlet data at the two extreme nodes.

Replacing the complete first and last equations is intentional: the OPT
interior operator is still equation (54), while the boundary equations are
u(x_left)=u_exact(x_left) and u(x_right)=u_exact(x_right).
"""
function solve_poisson(recipe, grid_data, config)
    # Static geometry appends its own singleton time coordinate. Supplying the
    # time-aware getModelPoints result here would append that dimension twice.
    model_points = (grid_data.nx,)
    model_family = (
        models=[grid_data.beta],
        modelPoints=model_points,
        Δ=grid_data.dx,
        modelName="poisson_$(config.name)_$(grid_data.nx)",
    )

    numerical_params = Dict{String,Any}(
        "optRec" => recipe,
        "modelFam" => model_family,
        "absorbingBoundaries" => nothing,
        "maskedRegionInSpace" => nothing,
        "backend" => :cpu,
        "representation" => REPRESENTATION,
        "compatibility_outputs" => false,
    )
    numerical = numericalOperatorConstruction(numerical_params)
    prepared = prepareLinearSystem(numerical["numOperators"])

    prepared.NField == 1 || error("Expected one Poisson field, got $(prepared.NField)")
    prepared.NForceField == 1 || error("Expected one force field, got $(prepared.NForceField)")
    prepared.timePointsUsedForOneStep == 1 ||
        error("Static Poisson recipe unexpectedly uses time marching")
    hasproperty(prepared, :R_force) ||
        error("The prepared system does not expose the manuscript Γ source operator")
    expected_gamma_size = (grid_data.nx, grid_data.nx)
    size(prepared.R_force) == expected_gamma_size ||
        error("Unexpected Γ dimensions $(size(prepared.R_force)); expected $expected_gamma_size")

    # OPT must redistribute the external source through Γ (manuscript eq. 55).
    # Guard against accidentally replacing Γ*f by the pointwise source vector.
    if startswith(config.name, "OPT")
        gamma_off_diagonal = copy(prepared.R_force)
        gamma_off_diagonal -= spdiagm(0 => diag(prepared.R_force))
        dropzeros!(gamma_off_diagonal)
        nnz(gamma_off_diagonal) > 0 ||
            error("OPT source operator Γ is diagonal; source redistribution was lost")
    end

    known_inputs = vcat(
        vec(prepared.known_lhs_template),
        vec(reshape(grid_data.force, :, 1, 1)),
    )
    A = sparse(prepared.A_template)
    b = copy(prepared.b_template)
    prepared.b_fun!(b, known_inputs)

    # Independent left and right boundary values (not the former hard-coded 1,1).
    left_value = first(grid_data.exact)
    right_value = last(grid_data.exact)
    A[1, :] .= 0
    A[1, 1] = 1
    b[1] = left_value
    A[end, :] .= 0
    A[end, end] = 1
    b[end] = right_value

    numerical_solution = real.(lu(A) \ b)
    isapprox(first(numerical_solution), left_value; atol=100eps(Float64)) ||
        error("Left Dirichlet condition was not imposed")
    isapprox(last(numerical_solution), right_value; atol=100eps(Float64)) ||
        error("Right Dirichlet condition was not imposed")
    return numerical_solution
end

rms_error(numerical, exact) = norm(numerical - exact) / sqrt(length(exact))

function save_results(output_dir, log_h_inverse, cases, configs, dxs, errors)
    mkpath(output_dir)
    checkpoint = joinpath(output_dir, "poisson1d_benchmark_EGU.jld2")
    config_names = getproperty.(configs, :name)
    case_names = getproperty.(cases, :name)
    jldsave(checkpoint; log_h_inverse, dxs, errors, config_names, case_names)

    figure = Figure(size=(1200, 800))
    for (case_index, case) in enumerate(cases)
        row = (case_index - 1) ÷ 3 + 1
        column = (case_index - 1) % 3 + 1
        axis = Axis(
            figure[row, column];
            title=case.name,
            xlabel="grid spacing h",
            ylabel="RMS error",
            xscale=log10,
            yscale=log10,
        )
        for (config_index, config) in enumerate(configs)
            scatterlines!(
                axis,
                dxs,
                errors[:, case_index, config_index];
                marker=:circle,
                label=config.name,
            )
        end
        axislegend(axis; position=:rb, labelsize=9)
    end
    figure_file = joinpath(output_dir, "poisson1d_benchmark_EGU.png")
    save(figure_file, figure)
    return checkpoint, figure_file
end

function main()
    @variables x
    L = 10pi
    log_h_inverse = QUICK ? [0.0, 0.5] : collect(0.0:0.5:3.0)
    configs = benchmark_configurations()
    cases = benchmark_cases(x)

    # A recipe depends on h and the discretisation, not on beta or the exact
    # solution. Build it once, serially, then distribute numerical assemblies
    # and solves across cases. This also avoids concurrent cache construction.
    recipes = Matrix{Any}(undef, length(log_h_inverse), length(configs))
    dxs = zeros(length(log_h_inverse))
    for h_index in eachindex(log_h_inverse)
        representative = grid_and_model(first(cases), x, log_h_inverse[h_index], L)
        dxs[h_index] = representative.dx
        for config_index in eachindex(configs)
            config = configs[config_index]
            @info "Constructing OPT recipe" h_index dx=representative.dx config=config.name
            recipes[h_index, config_index] = make_recipe(config, representative.dx)
        end
    end

    jobs = [
        (h_index=h_index, case_index=case_index, config_index=config_index)
        for config_index in eachindex(configs)
        for case_index in eachindex(cases)
        for h_index in eachindex(log_h_inverse)
    ]
    errors = fill(NaN, length(log_h_inverse), length(cases), length(configs))

    # numericalOperatorConstruction compiles Symbolics runtime functions whose
    # generated material symbols are not thread-local. Concurrent construction
    # has been observed to cross-contaminate jobs (for example, κ₂ leaking into
    # a one-material Poisson operator), producing non-reproducible errors.
    # Keep assembly and solves serial until that compiler path is made thread-safe.
    for job_index in eachindex(jobs)
        job = jobs[job_index]
        case = cases[job.case_index]
        config = configs[job.config_index]
        data = grid_and_model(case, x, log_h_inverse[job.h_index], L)
        solution = solve_poisson(recipes[job.h_index, job.config_index], data, config)
        errors[job.h_index, job.case_index, job.config_index] =
            rms_error(solution, data.exact)
        @info "Completed benchmark" case=case.name config=config.name nx=data.nx error=errors[job.h_index, job.case_index, job.config_index]
    end

    all(isfinite, errors) || error("At least one benchmark produced a non-finite error")
    output_dir = joinpath(FLEXOPT_ROOT, "scripts", "tmp", "poisson1d_EGU")
    checkpoint, figure_file =
        save_results(output_dir, log_h_inverse, cases, configs, dxs, errors)
    @info "Benchmark complete" checkpoint figure_file
    return errors
end

main()
