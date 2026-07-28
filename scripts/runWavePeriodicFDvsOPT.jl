#!/usr/bin/env julia
#
# Boundary-free manufactured-solution and homogeneous Fourier benchmarks for
# the wave equations in famousEquations.jl.  Every spatial resolution receives
# its own flexOPT recipe: no Δ=1 reference recipe is rescaled afterwards.

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
include(joinpath(@__DIR__, "famousEquationBenchmarkMatrix.jl"))

using .commonBatchs
using .flexOPT
using .FamousEquationBenchmarkMatrix

const WAVE_EQUATIONS = filter(e -> e.physics != :poisson, EQUATIONS)
const FIELD_WAVE = (2.0, 1.0)
const BASE = (rho=2.0, mu=3.0, lambda=4.0, velocity=2.0)
const CFL_RUN = 0.18
const SCHEMA_VERSION = "wave_periodic_fd_opt_v2"

value_for(spec::Number, ::Symbol, default) = Float64(spec)
value_for(spec::NamedTuple, key::Symbol, default) =
    hasproperty(spec, key) ? Float64(getproperty(spec, key)) : Float64(default)
value_for(spec, ::Symbol, default) = Float64(default)

function active_for(scenario, key)
    scenario.active === :all && return true
    isempty(scenario.active) && return false
    return key in scenario.active
end

function material_value_gradient(scenario, key, x, dimension)
    baseline = Float64(getproperty(BASE, key))
    active_for(scenario, key) || return baseline, zeros(dimension)
    q = collect(Float64.(scenario.material_wave[1:dimension]))
    phase = value_for(scenario.phase, key, 0.0)
    fraction = value_for(scenario.amplitude_fraction, key, 0.0)
    angle = dot(q, x) + phase
    return (
        baseline * (1 + fraction * cos(angle)),
        -baseline * fraction .* q .* sin(angle),
    )
end

function characteristic_speed(equation)
    equation.physics in (:sh_frequency, :sh_time) &&
        return sqrt(BASE.mu / BASE.rho)
    equation.physics == :acoustic && return BASE.velocity
    equation.physics == :elastic &&
        return sqrt((BASE.lambda + 2 * BASE.mu) / BASE.rho)
    error("Unsupported physics $(equation.physics)")
end

function interpolation(scheme, equation)
    return (
        ptsSpace=1,
        ptsTime=1,
        offsetSpace=(scheme.points_space - 1) / 2,
        offsetTime=equation.has_time ? (scheme.points_time - 1) / 2 : 0.0,
        YorderBspace=scheme.interpolation_order,
        YorderBtime=equation.has_time ? scheme.interpolation_order : -1,
    )
end

function make_recipe(equation, scheme, dx; dt=nothing)
    deltas = equation.has_time ?
        Tuple(vcat(fill(Float64(dx), equation.space_dimension), Float64(dt))) :
        Tuple(fill(Float64(dx), equation.space_dimension))
    params = Dict{String,Any}(
        "famousEquationType" => equation.equation,
        "Δ" => deltas,
        "orderBtime" => equation.has_time ? scheme.order_b_time : 0,
        "orderBspace" => scheme.order_b_space,
        "pointsInSpace" => scheme.points_space,
        "pointsInTime" => equation.has_time ? scheme.points_time : 1,
        "supplementaryOrder" => scheme.supplementary_order,
        "fieldItpl" => interpolation(scheme, equation),
        "materItpl" => interpolation(scheme, equation),
        "nuGeometryMode" => :middle,
        "hierarchicalTestFunctions" => false,
        "evenOrderHalfShiftMode" => :none,
        "recipe_backend" => CPU(),
    )
    return makeOPTsemiSymbolic(params)["recette"]
end

function local_offsets(recipe)
    nodes = vec(recipe.nodes[1])
    centre = nodes[recipe.centresIndices[1]]
    return [Tuple(Int.(node - centre)) for node in nodes]
end

struct NumericRecipe
    lhs
    rhs
    lhs_varm
    offsets
    nfields::Int
end

function compile_recipe(recipe)
    return NumericRecipe(
        flexOPT._compile_coefficient_recipes(recipe.lhs.Ajiννᶜ),
        flexOPT._compile_coefficient_recipes(recipe.rhs.Γjiννᶜ),
        recipe.lhs.varM,
        local_offsets(recipe),
        recipe.numbersOfTheSystem.numbersOfTheSystemL.NtypeofFields,
    )
end

# Some gallery force definitions use the literal `1` as extvars.  Symbolics
# then labels the corresponding local placeholder (for example `1₁`) even
# though it is not a physical material variable.  Such placeholders must be
# assigned one when Γ is evaluated.
function evaluate_coefficient(recipe, mapping)
    if hasproperty(recipe, :vars)
        isempty(recipe.vars) && return recipe.constant
        return recipe.f([get(mapping, variable, 1.0)
                         for variable in recipe.vars])
    end
    return complex(
        evaluate_coefficient(recipe.real_recipe, mapping),
        evaluate_coefficient(recipe.imag_recipe, mapping),
    )
end

function material_mapping(numeric, equation, scenario, x)
    mapping = Dict{Any,Any}()
    d = equation.space_dimension
    for local_index in eachindex(numeric.offsets)
        spatial_x = x[1:d] .+ collect(numeric.offsets[local_index][1:d])
        rho, _ = material_value_gradient(scenario, :rho, spatial_x, d)
        mu, _ = material_value_gradient(scenario, :mu, spatial_x, d)
        lambda, _ = material_value_gradient(scenario, :lambda, spatial_x, d)
        velocity, _ = material_value_gradient(scenario, :velocity, spatial_x, d)
        for variable_index in axes(numeric.lhs_varm, 1)
            symbol = numeric.lhs_varm[variable_index, local_index]
            text = string(symbol)
            if occursin("ρ", text)
                mapping[symbol] = rho
            elseif occursin("μ", text)
                mapping[symbol] = mu
            elseif occursin("λ", text)
                mapping[symbol] = lambda
            elseif startswith(text, "v") || occursin("v(", text)
                mapping[symbol] = velocity
            elseif occursin("ω", text)
                mapping[symbol] = 0.8 * sqrt(BASE.mu / BASE.rho) *
                    norm(FIELD_WAVE[1:d])
            else
                mapping[symbol] = 1.0
            end
        end
    end
    return mapping
end

function manufactured(equation, scenario, coordinates; branch=:P)
    d = equation.space_dimension
    x = coordinates[1:d]
    k = collect(FIELD_WAVE[1:d])
    k2 = dot(k, k)
    rho, _ = material_value_gradient(scenario, :rho, x, d)
    mu, gradmu = material_value_gradient(scenario, :mu, x, d)
    lambda, gradlambda = material_value_gradient(scenario, :lambda, x, d)
    velocity, _ = material_value_gradient(scenario, :velocity, x, d)

    if equation.physics == :sh_frequency
        theta = dot(k, x)
        omega = 0.8 * sqrt(BASE.mu / BASE.rho) * norm(k)
        u = [cos(theta)]
        f = [rho * omega^2 * cos(theta) -
             dot(gradmu, k) * sin(theta) - mu * k2 * cos(theta)]
    elseif equation.physics == :sh_time
        omega = 0.8 * sqrt(BASE.mu / BASE.rho) * norm(k)
        theta = dot(k, x) - omega * coordinates[end]
        u = [cos(theta)]
        f = [-rho * omega^2 * cos(theta) +
             dot(gradmu, k) * sin(theta) + mu * k2 * cos(theta)]
    elseif equation.physics == :acoustic
        omega = 0.8 * BASE.velocity * norm(k)
        theta = dot(k, x) - omega * coordinates[end]
        u = [cos(theta)]
        f = [(-omega^2 + velocity^2 * k2) * cos(theta)]
    elseif equation.physics == :elastic
        khat = k / norm(k)
        p = branch == :P ? khat : [-khat[2], khat[1]]
        speed = branch == :P ?
            sqrt((BASE.lambda + 2 * BASE.mu) / BASE.rho) :
            sqrt(BASE.mu / BASE.rho)
        omega = 0.8 * speed * norm(k)
        theta = dot(k, x) - omega * coordinates[end]
        stiffness = mu * k2 .* p + (lambda + mu) .* k .* dot(k, p)
        gradient = gradlambda .* dot(k, p) +
            dot(gradmu, k) .* p + k .* dot(gradmu, p)
        u = p .* cos(theta)
        f = (-rho * omega^2 .* p + stiffness) .* cos(theta) +
            gradient .* sin(theta)
    else
        error("Unsupported physics")
    end
    return u, f
end

function periodic_coordinates(index, n, dx, dt, equation)
    d = equation.space_dimension
    x = [(mod(index[j] - 1, n)) * dx for j in 1:d]
    equation.has_time || return x
    nt = max(8, n)
    push!(x, mod(index[end] - 1, nt) * dt)
    return x
end

function sample_indices(n, equation; max_samples=768)
    dimensions = equation.space_dimension + equation.has_time
    lengths = equation.has_time ?
        Tuple(vcat(fill(n, equation.space_dimension), max(8, n))) :
        Tuple(fill(n, equation.space_dimension))
    total = prod(lengths)
    stride = max(1, ceil(Int, (total / max_samples)^(1 / dimensions)))
    axes_sample = [1:stride:lengths[j] for j in 1:dimensions]
    return Iterators.product(axes_sample...)
end

function residual_error(numeric, equation, scenario, n, dx, dt; branch=:P)
    residual2 = 0.0
    force2 = 0.0
    count = 0
    d = equation.space_dimension
    for index in sample_indices(n, equation)
        centre = periodic_coordinates(index, n, dx, dt, equation)
        # material_mapping expects physical offsets, not integer offsets.
        physical = copy(centre)
        mapping = Dict{Any,Any}()
        for local_index in eachindex(numeric.offsets)
            offset = numeric.offsets[local_index]
            point = centre .+ [
                offset[j] * (j <= d ? dx : dt) for j in eachindex(offset)]
            rho, _ = material_value_gradient(scenario, :rho, point[1:d], d)
            mu, _ = material_value_gradient(scenario, :mu, point[1:d], d)
            lambda, _ = material_value_gradient(scenario, :lambda, point[1:d], d)
            velocity, _ = material_value_gradient(scenario, :velocity, point[1:d], d)
            for variable_index in axes(numeric.lhs_varm, 1)
                symbol = numeric.lhs_varm[variable_index, local_index]
                text = string(symbol)
                mapping[symbol] = occursin("ρ", text) ? rho :
                    occursin("μ", text) ? mu :
                    occursin("λ", text) ? lambda :
                    (startswith(text, "v") || occursin("v(", text)) ? velocity :
                    occursin("ω", text) ?
                        0.8 * sqrt(BASE.mu / BASE.rho) * norm(FIELD_WAVE[1:d]) :
                    1.0
            end
        end
        lhs_value = zeros(Float64, numeric.nfields)
        rhs_value = zeros(Float64, numeric.nfields)
        gamma_volume = zeros(Float64, numeric.nfields)
        for local_index in eachindex(numeric.offsets)
            offset = numeric.offsets[local_index]
            point = centre .+ [
                offset[j] * (j <= d ? dx : dt) for j in eachindex(offset)]
            field, force = manufactured(equation, scenario, point; branch)
            for equation_index in 1:numeric.nfields
                for field_index in 1:numeric.nfields
                    lhs_value[equation_index] +=
                        evaluate_coefficient(
                            numeric.lhs[local_index, field_index, equation_index, 1],
                            mapping) * field[field_index]
                end
                for force_index in 1:numeric.nfields
                    gamma = evaluate_coefficient(
                        numeric.rhs[local_index, force_index, equation_index, 1],
                        Dict{Any,Any}())
                    rhs_value[equation_index] += gamma * force[force_index]
                    gamma_volume[equation_index] += abs(gamma)
                end
            end
        end
        _, centre_force = manufactured(equation, scenario, centre; branch)
        for field in 1:numeric.nfields
            scale = max(gamma_volume[field], eps(Float64))
            residual2 += ((lhs_value[field] - rhs_value[field]) / scale)^2
            force2 += centre_force[field]^2
            count += 1
        end
    end
    absolute = sqrt(residual2 / count)
    relative = absolute / max(sqrt(force2 / count), eps(Float64))
    return absolute, relative, count
end

function numeric_symbol(numeric, equation, wave_angles, mapping)
    nf = numeric.nfields
    lhs = zeros(ComplexF64, nf, nf)
    rhs = zeros(ComplexF64, nf, nf)
    for local_index in eachindex(numeric.offsets)
        phase = cis(dot(wave_angles, numeric.offsets[local_index]))
        for i in 1:nf, j in 1:nf
            lhs[i, j] += evaluate_coefficient(
                numeric.lhs[local_index, j, i, 1], mapping) * phase
            rhs[i, j] += evaluate_coefficient(
                numeric.rhs[local_index, j, i, 1], Dict{Any,Any}()) * phase
        end
    end
    return lhs, rhs
end

function homogeneous_mapping(numeric, equation)
    scenario = first(MATERIAL_SCENARIOS)
    coordinates = zeros(equation.space_dimension + equation.has_time)
    return material_mapping(numeric, equation, scenario, coordinates)
end

function amplification_roots(numeric, equation, spatial_angles, mapping)
    nf = numeric.nfields
    blocks = Dict(-1 => zeros(ComplexF64, nf, nf),
                   0 => zeros(ComplexF64, nf, nf),
                   1 => zeros(ComplexF64, nf, nf))
    d = equation.space_dimension
    for local_index in eachindex(numeric.offsets)
        offset = numeric.offsets[local_index]
        phase = cis(dot(spatial_angles, offset[1:d]))
        time_offset = offset[end]
        for i in 1:nf, j in 1:nf
            blocks[time_offset][i, j] += evaluate_coefficient(
                numeric.lhs[local_index, j, i, 1], mapping) * phase
        end
    end
    plus, zero, minus = blocks[1], blocks[0], blocks[-1]
    left = [-zero -minus; Matrix{ComplexF64}(I, nf, nf) zeros(ComplexF64, nf, nf)]
    right = [plus zeros(ComplexF64, nf, nf); zeros(ComplexF64, nf, nf) Matrix{ComplexF64}(I, nf, nf)]
    return eigvals(left, right)
end

function max_amplification(numeric, equation; samples=13)
    mapping = homogeneous_mapping(numeric, equation)
    angles = range(0, pi; length=samples)
    maximum_modulus = 0.0
    if equation.space_dimension == 1
        for theta in angles
            maximum_modulus = max(maximum_modulus,
                maximum(abs.(amplification_roots(
                    numeric, equation, [theta], mapping))))
        end
    else
        for theta_x in angles, theta_y in angles
            maximum_modulus = max(maximum_modulus,
                maximum(abs.(amplification_roots(
                    numeric, equation, [theta_x, theta_y], mapping))))
        end
    end
    return maximum_modulus
end

function cfl_limit(equation, scheme; quick=false)
    candidates = quick ? collect(0.1:0.2:0.9) : collect(0.05:0.1:1.25)
    last_stable = 0.0
    first_unstable = Inf
    rows = NamedTuple[]
    for cfl in candidates
        dt = cfl / characteristic_speed(equation)
        elapsed = @elapsed recipe = make_recipe(equation, scheme, 1.0; dt)
        numeric = compile_recipe(recipe)
        amp = max_amplification(numeric, equation; samples=quick ? 7 : 13)
        stable = amp <= 1 + 2e-7
        push!(rows, (; cfl, amp, stable, recipe_seconds=elapsed))
        if stable
            last_stable = cfl
        else
            first_unstable = min(first_unstable, cfl)
            isfinite(first_unstable) && last_stable > 0 && break
        end
    end
    if isfinite(first_unstable) && last_stable > 0
        for _ in 1:(quick ? 2 : 5)
            cfl = (last_stable + first_unstable) / 2
            dt = cfl / characteristic_speed(equation)
            elapsed = @elapsed recipe = make_recipe(equation, scheme, 1.0; dt)
            numeric = compile_recipe(recipe)
            amp = max_amplification(numeric, equation; samples=quick ? 7 : 13)
            stable = amp <= 1 + 2e-7
            push!(rows, (; cfl, amp, stable, recipe_seconds=elapsed))
            stable ? (last_stable = cfl) : (first_unstable = cfl)
        end
    end
    return last_stable, rows
end

function physical_root(roots, exact_angle)
    candidates = filter(z -> isfinite(real(z)) && isfinite(imag(z)) &&
        abs(abs(z) - 1) < 0.08 && angle(z) >= -1e-10, roots)
    isempty(candidates) && return NaN
    return angle(candidates[argmin(abs.(angle.(candidates) .- exact_angle))])
end

function dispersion_rows_for(equation, scheme, cfl; quick=false)
    cfl_use = cfl > 0 ? 0.8cfl : CFL_RUN
    dt = cfl_use / characteristic_speed(equation)
    recipe = make_recipe(equation, scheme, 1.0; dt)
    numeric = compile_recipe(recipe)
    mapping = homogeneous_mapping(numeric, equation)
    ppws = quick ? [5.0, 8.0, 12.0, 20.0] :
        [4.0, 5.0, 6.0, 8.0, 10.0, 12.0, 16.0, 24.0, 32.0]
    directions = equation.space_dimension == 1 ?
        [(name="axis", unit=[1.0])] :
        [(name="axis", unit=[1.0, 0.0]),
         (name="diagonal", unit=[inv(sqrt(2)), inv(sqrt(2))])]
    branches = equation.physics == :elastic ? (:P, :S) : equation.branches
    rows = NamedTuple[]
    for direction in directions, ppw in ppws, branch in branches
        theta = 2pi / ppw
        spatial_angles = theta .* direction.unit
        roots = amplification_roots(numeric, equation, spatial_angles, mapping)
        speed = branch == :S ? sqrt(BASE.mu / BASE.rho) :
            branch == :P && equation.physics == :elastic ?
                sqrt((BASE.lambda + 2 * BASE.mu) / BASE.rho) :
                characteristic_speed(equation)
        exact_angle = speed * theta * dt
        numerical_angle = physical_root(roots, exact_angle)
        ratio = numerical_angle / exact_angle
        push!(rows, (
            equation=equation.label,
            equation_id=equation.equation,
            scheme=scheme.name,
            branch=String(branch),
            direction=direction.name,
            points_per_wavelength=ppw,
            cfl=cfl_use,
            numerical_phase_velocity=ratio * speed,
            exact_phase_velocity=speed,
            relative_phase_error=ratio - 1,
        ))
    end
    return rows
end

function observed_orders!(rows)
    groups = Dict{Tuple,Vector{Int}}()
    for (index, row) in pairs(rows)
        key = (row.equation_id, row.scheme, row.scenario, row.branch)
        push!(get!(groups, key, Int[]), index)
    end
    orders = fill(NaN, length(rows))
    for indices in values(groups)
        sort!(indices; by=i -> rows[i].n)
        for pair in 2:length(indices)
            old, new = indices[pair - 1], indices[pair]
            orders[new] = log(rows[old].absolute_error /
                rows[new].absolute_error) / log(rows[new].n / rows[old].n)
        end
    end
    return [merge(rows[i], (; observed_order=orders[i])) for i in eachindex(rows)]
end

function selected_equations()
    requested = split(get(ENV, "FLEXOPT_WAVE_EQUATIONS", ""), ",")
    isempty(first(requested)) && return WAVE_EQUATIONS
    return filter(e -> e.equation in requested, WAVE_EQUATIONS)
end

function selected_schemes()
    requested = split(get(ENV, "FLEXOPT_WAVE_SCHEMES", ""), ",")
    isempty(first(requested)) && return SCHEMES
    return filter(s -> s.name in requested, SCHEMES)
end

function selected_sizes(equation, quick)
    sizes = space_sizes(equation)
    quick && return equation.space_dimension == 1 ? [16, 32] : [8, 16]
    cap = parse(Int, get(ENV, "FLEXOPT_WAVE_MAX_N", string(maximum(sizes))))
    return filter(<=(cap), sizes)
end

function run_convergence(equations, schemes; quick=false)
    rows = NamedTuple[]
    timing_rows = NamedTuple[]
    for equation in equations, scheme in schemes
        scenarios = material_scenarios(equation)
        branches = equation.physics == :elastic ? (:P, :S) : equation.branches
        for n in selected_sizes(equation, quick)
            dx = 2pi / n
            dt = equation.has_time ? CFL_RUN * dx / characteristic_speed(equation) : nothing
            @info "Wave periodic recipe" equation=equation.label scheme=scheme.name n dx dt
            recipe = nothing
            recipe_seconds = @elapsed recipe = make_recipe(equation, scheme, dx; dt)
            numeric = compile_recipe(recipe)
            # Time a production-like repeated local-symbol evaluation separately
            mapping = homogeneous_mapping(numeric, equation)
            angles = fill(0.37, equation.space_dimension + equation.has_time)
            numeric_symbol(numeric, equation, angles, mapping) # warm-up
            repeats = 20
            application_seconds = @elapsed for _ in 1:repeats
                numeric_symbol(numeric, equation, angles, mapping)
            end
            push!(timing_rows, (
                equation=equation.label, equation_id=equation.equation,
                scheme=scheme.name, n, dx,
                local_space_time_points=length(numeric.offsets),
                recipe_seconds,
                symbol_application_seconds=application_seconds / repeats,
            ))
            for scenario in scenarios, branch in branches
                elapsed = @elapsed absolute, relative, samples =
                    residual_error(numeric, equation, scenario, n, dx, dt; branch)
                push!(rows, (
                    equation=equation.label, equation_id=equation.equation,
                    physics=String(equation.physics), scheme=scheme.name,
                    family=String(scheme.family), scenario=scenario.name,
                    relation=String(scenario.relation), branch=String(branch),
                    n, dx, dt=equation.has_time ? dt : NaN,
                    absolute_error=absolute, relative_error=relative,
                    samples, residual_seconds=elapsed,
                ))
            end
        end
    end
    return observed_orders!(rows), timing_rows
end

function main()
    quick = get(ENV, "FLEXOPT_WAVE_QUICK", "0") == "1"
    skip_convergence = get(ENV, "FLEXOPT_WAVE_SKIP_CONVERGENCE", "0") == "1"
    skip_fourier = get(ENV, "FLEXOPT_WAVE_SKIP_FOURIER", "0") == "1"
    equations = selected_equations()
    schemes = selected_schemes()
    convergence_rows, timing_rows = skip_convergence ?
        (NamedTuple[], NamedTuple[]) :
        run_convergence(equations, schemes; quick)

    cfl_rows = NamedTuple[]
    cfl_scan_rows = NamedTuple[]
    dispersion_rows = NamedTuple[]
    if !skip_fourier
        for equation in filter(e -> e.has_time, equations), scheme in schemes
            @info "CFL scan" equation=equation.label scheme=scheme.name
            cfl, scan = cfl_limit(equation, scheme; quick)
            for item in scan
                push!(cfl_scan_rows, merge((
                    equation=equation.label, equation_id=equation.equation,
                    scheme=scheme.name,
                ), item))
            end
            branches = equation.physics == :elastic ? (:P, :S) : equation.branches
            for branch in branches
                push!(cfl_rows, (
                    equation=equation.label, equation_id=equation.equation,
                    scheme=scheme.name, branch=String(branch), cfl_max=cfl,
                ))
            end
            append!(dispersion_rows,
                dispersion_rows_for(equation, scheme, cfl; quick))
        end
    end

    output_dir = joinpath(FLEXOPT_ROOT, "scripts", "tmp",
        "famous_equation_periodic_benchmarks")
    mkpath(output_dir)
    output_file = joinpath(output_dir, quick ?
        "wave_fd_vs_opt_quick.jld2" : "wave_fd_vs_opt.jld2")
    jldsave(output_file;
        schema_version=SCHEMA_VERSION,
        convergence_rows,
        timing_rows,
        cfl_rows,
        cfl_scan_rows,
        dispersion_rows,
        schemes,
        equations,
        base_material=BASE,
        field_wave=FIELD_WAVE,
        cfl_used_for_convergence=CFL_RUN,
        boundary_mode="periodic manufactured local residual; no boundary closure",
        error_definition=(
            absolute="RMS((A*u-Gamma*f)/sum(abs(Gamma)))",
            relative="absolute/RMS(f)",
            note="Each Δ has an independently constructed recette.",
        ),
        cpu_definition=(
            recipe="wall time for makeOPTsemiSymbolic",
            residual="wall time for sampled manufactured residual",
            application="warm local homogeneous Fourier-symbol evaluation",
        ),
    )
    println("\nSaved: ", output_file)
    println("Convergence rows: ", length(convergence_rows))
    println("Timing rows: ", length(timing_rows))
    println("CFL rows: ", length(cfl_rows))
    println("Dispersion rows: ", length(dispersion_rows))
    return output_file
end

abspath(PROGRAM_FILE) == abspath(@__FILE__) && main()
