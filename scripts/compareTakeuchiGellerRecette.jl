#!/usr/bin/env julia
#
# Operator-level validation against:
#   Geller & Takeuchi (1998), GJI 135, Eq. (18)
#   Takeuchi & Geller (2000), PEPI 119, Eqs. (23)-(25)
#
# Usage:
#   julia --project=.. scripts/compareTakeuchiGellerRecette.jl
#   FLEXOPT_RECETTE_QUICK=1 julia --project=.. scripts/compareTakeuchiGellerRecette.jl
#   FLEXOPT_RECETTE_TABLES=1 julia --project=.. scripts/compareTakeuchiGellerRecette.jl
#   FLEXOPT_RECETTE_VERBOSE=1 julia --project=.. scripts/compareTakeuchiGellerRecette.jl
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
const PRINT_TABLES = get(ENV, "FLEXOPT_RECETTE_TABLES", QUICK ? "1" : "0") == "1"
const VERBOSE_GENERATOR = get(ENV, "FLEXOPT_RECETTE_VERBOSE", "0") == "1"

interpolation(offset) = (
    ptsSpace=1,
    ptsTime=1,
    offsetSpace=offset,
    offsetTime=1,
    YorderBspace=-1,
    YorderBtime=-1,
)

function recipe_parameters(
    equation,
    spacings,
    supplementary_order;
    points_in_space=3,
    interpolation_offset=(points_in_space - 1) / 2,
)
    return Dict{String,Any}(
        "famousEquationType" => equation,
        "Δ" => spacings,
        "orderBtime" => 1,
        "orderBspace" => 1,
        "pointsInSpace" => points_in_space,
        "pointsInTime" => 3,
        "supplementaryOrder" => supplementary_order,
        "fieldItpl" => interpolation(interpolation_offset),
        "materItpl" => interpolation(interpolation_offset),
        "nuGeometryMode" => :middle,
        "recipe_backend" => CPU(),
    )
end

function generated_recipe(
    equation,
    spacings;
    supplementary_order=2,
    points_in_space=3,
    interpolation_offset=(points_in_space - 1) / 2,
)
    key = (
        equation,
        (
            Tuple(spacings)...,
            supplementary_order,
            points_in_space,
            interpolation_offset,
        ),
    )
    return get!(RECIPE_CACHE, key) do
        construct = () -> makeOPTsemiSymbolic(
            recipe_parameters(
                equation,
                spacings,
                supplementary_order;
                points_in_space,
                interpolation_offset,
            ),
        )["recette"]
        if VERBOSE_GENERATOR
            construct()
        else
            redirect_stdout(devnull) do
                construct()
            end
        end
    end
end

function material_value(symbol, values)
    name = string(symbol)
    occursin("ρ", name) && return values.rho
    (occursin("λ", name) || occursin("lambda", name)) && return values.lambda
    occursin("μ", name) && return values.mu
    occursin("v", name) && return values.velocity
    error("No material substitution was supplied for flexOPT symbol $name")
end

function numerical_coefficients(
    recipe,
    values;
    equation_index=1,
    field_index=1,
    scale=1.0,
)
    replacements = Dict{Any,Any}()
    for symbol in vec(recipe.lhs.varM)
        raw_symbol = symbol isa Num ? symbol[] : symbol
        replacements[raw_symbol] = material_value(symbol, values)
    end

    geometry = 1
    symbolic = recipe.lhs.Ajiννᶜ[
        :,
        equation_index,
        field_index,
        geometry,
    ]
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

function paper_1998_components(delta_x, delta_t)
    rho_coefficients = Dict{Tuple{Int,Int},Float64}()
    mu_coefficients = Dict{Tuple{Int,Int},Float64}()
    smoothing = (1.0, 10.0, 1.0)
    second = (1.0, -2.0, 1.0)
    for (ix, x) in enumerate(-1:1), (it, t) in enumerate(-1:1)
        rho_coefficients[(x, t)] =
            smoothing[ix] * second[it] / (12 * delta_t^2)
        mu_coefficients[(x, t)] =
            -second[ix] * smoothing[it] / (12 * delta_x^2)
    end
    return (; rho=rho_coefficients, mu=mu_coefficients)
end

function paper_2000_components(delta_x, delta_z, delta_t)
    rho_coefficients = Dict{Tuple{Int,Int,Int},Float64}()
    mu_coefficients = Dict{Tuple{Int,Int,Int},Float64}()
    smoothing = (1.0, 10.0, 1.0)
    second = (1.0, -2.0, 1.0)
    for (ix, x) in enumerate(-1:1),
        (iz, z) in enumerate(-1:1),
        (it, t) in enumerate(-1:1)

        rho_coefficients[(x, z, t)] =
            smoothing[ix] * smoothing[iz] * second[it] /
            (144 * delta_t^2)
        mu_coefficients[(x, z, t)] =
            -second[ix] * smoothing[iz] * smoothing[it] /
            (144 * delta_x^2) -
            smoothing[ix] * second[iz] * smoothing[it] /
            (144 * delta_z^2)
    end
    return (; rho=rho_coefficients, mu=mu_coefficients)
end

function paper_2000_psv_components(delta_x, delta_z, delta_t)
    offsets = collect(Iterators.product(-1:2, -1:2, -1:1))
    zero_stencil() = Dict(offset => 0.0 for offset in offsets)

    smoothing = Dict(-1 => 1.0, 0 => 10.0, 1 => 1.0)
    second = Dict(-1 => 1.0, 0 => -2.0, 1 => 1.0)
    # Eq. (52): [-5, -3, 9, -1]/12 at offsets [-1, 0, 1, 2].
    first = Dict(-1 => -5.0, 0 => -3.0, 1 => 9.0, 2 => -1.0)

    mass = zero_stencil()
    stiffness_x = zero_stencil()
    stiffness_z = zero_stencil()
    mixed = zero_stencil()
    for offset in offsets
        x, z, t = offset
        if haskey(smoothing, x) && haskey(smoothing, z)
            mass[offset] =
                smoothing[x] * smoothing[z] * second[t] /
                (144 * delta_t^2)
            stiffness_x[offset] =
                second[x] * smoothing[z] * smoothing[t] /
                (144 * delta_x^2)
            stiffness_z[offset] =
                smoothing[x] * second[z] * smoothing[t] /
                (144 * delta_z^2)
        end
        mixed[offset] =
            first[x] * first[z] * smoothing[t] /
            (1728 * delta_x * delta_z)
    end

    linear_combination(terms...) = Dict(
        offset => sum(weight * stencil[offset] for (weight, stencil) in terms)
        for offset in offsets
    )
    zeros = zero_stencil()
    return Dict(
        (1, 1) => (
            rho=mass,
            lambda=linear_combination((-1.0, stiffness_x)),
            mu=linear_combination((-2.0, stiffness_x), (-1.0, stiffness_z)),
        ),
        (1, 2) => (
            rho=zeros,
            lambda=linear_combination((-1.0, mixed)),
            mu=linear_combination((-1.0, mixed)),
        ),
        (2, 1) => (
            rho=zeros,
            lambda=linear_combination((-1.0, mixed)),
            mu=linear_combination((-1.0, mixed)),
        ),
        (2, 2) => (
            rho=mass,
            lambda=linear_combination((-1.0, stiffness_z)),
            mu=linear_combination((-1.0, stiffness_x), (-2.0, stiffness_z)),
        ),
    )
end

function combine_components(components, rho, mu)
    return Dict(
        offset => rho * components.rho[offset] + mu * components.mu[offset]
        for offset in keys(components.rho)
    )
end

function subtract_coefficients(left, right)
    Set(keys(left)) == Set(keys(right)) ||
        error("Cannot subtract different stencil supports")
    return Dict(offset => left[offset] - right[offset] for offset in keys(left))
end

function combine_elastic_components(components, rho, lambda, mu)
    return Dict(
        offset =>
            rho * components.rho[offset] +
            lambda * components.lambda[offset] +
            mu * components.mu[offset]
        for offset in keys(components.rho)
    )
end

function compare_coefficients(
    label,
    generated,
    published;
    print_table=false,
)
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
    if print_table
        println("\n", label)
        println("  offset          flexOPT              paper                difference")
        for offset in sort!(collect(keys(published)))
            @printf(
                "  %-14s  % .12e  % .12e  % .3e\n",
                string(offset),
                generated[offset],
                published[offset],
                generated[offset] - published[offset],
            )
        end
    end
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
    recipe = generated_recipe(
        "1DsismoTimeHomo",
        (delta_x, delta_t),
        supplementary_order=0,
    )
    scale = inv(delta_x^2)
    generated_components = (
        rho=numerical_coefficients(
            recipe,
            (rho=1.0, mu=0.0, velocity=0.0);
            scale,
        ),
        mu=numerical_coefficients(
            recipe,
            (rho=0.0, mu=1.0, velocity=0.0);
            scale,
        ),
    )
    published_components = paper_1998_components(delta_x, delta_t)
    compare_coefficients(
        "1998 rho coefficients",
        generated_components.rho,
        published_components.rho,
        print_table=PRINT_TABLES,
    )
    compare_coefficients(
        "1998 mu coefficients",
        generated_components.mu,
        published_components.mu,
        print_table=PRINT_TABLES,
    )
    generated = combine_components(generated_components, rho, mu)
    published = combine_components(published_components, rho, mu)
    return compare_coefficients(
        "1998 combined, Δ=($delta_x,$delta_t)",
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
    # For u_tt-v^2*laplacian(u), the v=0 evaluation isolates the mass
    # coefficient. The difference between v=1 and v=0 isolates the rigidity
    # coefficient after restoring rho*(v^2)=mu.
    # flexOPT integrates over x, z, and t. Its raw weak-form coefficients
    # therefore carry Δ^3 from the integration measure, while the two
    # derivatives contribute Δ^-2. The paper divides out that remaining Δ,
    # so converting the dimensionless recette to its convention requires
    # the complete Δ^-3 factor.
    physical_scale = inv(delta_x^3)
    generated_rho = numerical_coefficients(
        recipe,
        (rho=1.0, mu=0.0, velocity=0.0);
        scale=physical_scale,
    )
    generated_unit_velocity = numerical_coefficients(
        recipe,
        (rho=1.0, mu=1.0, velocity=1.0);
        scale=physical_scale,
    )
    generated_components = (
        rho=generated_rho,
        mu=subtract_coefficients(generated_unit_velocity, generated_rho),
    )
    published_components = paper_2000_components(
        delta_x,
        delta_z,
        delta_t,
    )
    compare_coefficients(
        "2000 rho coefficients",
        generated_components.rho,
        published_components.rho,
        print_table=PRINT_TABLES,
    )
    compare_coefficients(
        "2000 mu coefficients",
        generated_components.mu,
        published_components.mu,
        print_table=PRINT_TABLES,
    )
    generated = combine_components(generated_components, rho, mu)
    published = combine_components(published_components, rho, mu)
    return compare_coefficients(
        "2000 combined, Δ=($delta_x,$delta_z,$delta_t)",
        generated,
        published,
    )
end

function validate_2000_psv(
    delta_x,
    delta_z,
    delta_t,
    rho,
    lambda,
    mu,
)
    delta_x == delta_z == delta_t || error(
        "The current P-SV comparison uses one common grid scale; " *
        "got Δx=$delta_x, Δz=$delta_z, and Δt=$delta_t",
    )
    recipe = generated_recipe(
        "2DsismoTimeIsoHomo",
        (delta_x, delta_z, delta_t);
        points_in_space=4,
        supplementary_order=2,
        interpolation_offset=1.0,
    )
    physical_scale = inv(delta_x^3)
    published_blocks = paper_2000_psv_components(
        delta_x,
        delta_z,
        delta_t,
    )

    basis_values = (
        rho=(rho=1.0, lambda=0.0, mu=0.0, velocity=0.0),
        lambda=(rho=0.0, lambda=1.0, mu=0.0, velocity=0.0),
        mu=(rho=0.0, lambda=0.0, mu=1.0, velocity=0.0),
    )
    # The paper deliberately mixes supports: diagonal P-SV blocks use the
    # three-point Eqs. (23)-(25), already validated above, while only the
    # off-diagonal mixed derivatives use the four-point Eq. (54). Since
    # flexOPT currently selects pointsInSpace globally per recette, compare
    # Eq. (54) on its two off-diagonal blocks here.
    for (equation_index, field_index) in ((1, 2), (2, 1))
        block = (equation_index, field_index)
        generated = (
            rho=numerical_coefficients(
                recipe,
                basis_values.rho;
                equation_index,
                field_index,
                scale=physical_scale,
            ),
            lambda=numerical_coefficients(
                recipe,
                basis_values.lambda;
                equation_index,
                field_index,
                scale=physical_scale,
            ),
            mu=numerical_coefficients(
                recipe,
                basis_values.mu;
                equation_index,
                field_index,
                scale=physical_scale,
            ),
        )
        published = published_blocks[block]
        for component in (:rho, :lambda, :mu)
            compare_coefficients(
                "2000 P-SV[$equation_index,$field_index] $component",
                getproperty(generated, component),
                getproperty(published, component);
                print_table=PRINT_TABLES,
            )
        end
        compare_coefficients(
            "2000 P-SV[$equation_index,$field_index] combined",
            combine_elastic_components(generated, rho, lambda, mu),
            combine_elastic_components(published, rho, lambda, mu),
        )
    end
    return nothing
end

function validate_2000_eq54(delta_x, delta_z, delta_t, lambda)
    error(
        "Eq. (54) cannot yet be compared as a flexOPT recette: the current " *
        "dependency analyzer does not expose an isolated mixed derivative, " *
        "and the coupled elastic recipe uses one global spatial support.",
    )
end

function main()
    println("Takeuchi-Geller local-recette validation")
    println(
        "Comparing rho and mu interior coefficients separately, then combined; " *
        "no time marching is used.\n",
    )

    material_pairs = QUICK ? (
        (rho=1.0, lambda=1.0, mu=1.0),
    ) : (
        (rho=1.0, lambda=1.0, mu=1.0),
        (rho=2.7, lambda=20.0, mu=31.0),
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
        @printf(
            "ρ = %.6g, λ = %.6g, μ = %.6g\n",
            material.rho,
            material.lambda,
            material.mu,
        )
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

    println(
        "All enabled generated-recette comparisons agree with the published " *
        "operators.",
    )
    println(
        "Eq. (54) is not enabled: flexOPT currently reconstructs the coupled " *
        "elastic system globally and does not expose the paper's separately " *
        "assembled four-point mixed-derivative operator.",
    )
    return nothing
end

main()
