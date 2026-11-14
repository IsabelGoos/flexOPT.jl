using JLD2
using Symbolics
using LinearAlgebra

const ROOT = normpath(joinpath(@__DIR__, ".."))
include(joinpath(ROOT, "src", "commonBatchs.jl"))
include(joinpath(ROOT, "src", "flexOPT.jl"))
using .commonBatchs, .flexOPT

const RECIPES = (
    volume=joinpath(ROOT, "data", "semiSymbolic", "elastic2D_OPT3_38ef42fe.jld2"),
    shifted=joinpath(ROOT, "data", "semiSymbolic",
        "elastic2D_OPT3_surface_shifted_full_46fe875f.jld2"),
    clipped=joinpath(ROOT, "data", "semiSymbolic",
        "elastic2D_OPT3_surface_clipped_available_a2b2a2fb.jld2"),
)

function numeric_coefficient(expression, substitutions)
    value = Symbolics.substitute(expression, substitutions)
    Float64(Symbolics.value(value))
end

function recipe_residual(recipe, ux, uz; rho=0.0, lambda=0.0, mu=0.0)
    lhs = recipe.lhs
    substitutions = Dict{Any,Any}()
    material_values = (rho, lambda, mu)
    for ivar in axes(lhs.varM, 1), inode in axes(lhs.varM, 2)
        substitutions[lhs.varM[ivar, inode][]] = material_values[ivar]
    end
    nodes = vec(recipe.nodes[1])
    fields = (ux, uz)
    residual = zeros(Float64, 2)
    for equation in 1:2, field in 1:2, inode in eachindex(nodes)
        point = nodes[inode]
        x, z, t = Float64.(Tuple(point))
        coefficient = numeric_coefficient(
            lhs.Ajiννᶜ[inode, field, equation, 1], substitutions)
        residual[equation] += coefficient * fields[field](x, z, t)
    end
    residual
end

zero_field(x, z, t) = 0.0

function audit_recipe(path)
    recipe = load(path, "recette")
    tests = (
        lambda_same=recipe_residual(recipe, (x,z,t)->x, zero_field; lambda=1),
        lambda_cross=recipe_residual(recipe, zero_field, (x,z,t)->z; lambda=1),
        mu_same=recipe_residual(recipe, (x,z,t)->z, zero_field; mu=1),
        mu_cross=recipe_residual(recipe, zero_field, (x,z,t)->x; mu=1),
    )
    comparisons = (
        lambda=tests.lambda_same - tests.lambda_cross,
        mu=tests.mu_same - tests.mu_cross,
    )
    scale = maximum(abs, reduce(vcat, values(tests)))
    mismatch = maximum(abs, reduce(vcat, values(comparisons)))
    relative_mismatch = mismatch / max(scale, eps())
    rigid_rotation = recipe_residual(recipe,
        (x,z,t)->z, (x,z,t)->-x; mu=1)
    lambda, mu = 2.0, 3.0
    traction_free_uniform = recipe_residual(recipe,
        (x,z,t)->x,
        (x,z,t)->-(lambda / (lambda + 2mu)) * z;
        lambda, mu)
    (; comparisons, scale, mismatch, relative_mismatch,
       rigid_rotation, traction_free_uniform,
       test_derivative_orders=recipe.testDerivativeOrders,
       nodes=size(recipe.nodes[1]))
end

results = map(audit_recipe, RECIPES)
for name in keys(results)
    println("\n=== ", name, " ===")
    display(getproperty(results, name))
end

@assert results.volume.mismatch < 5e-6
for result in (results.shifted, results.clipped)
    @assert result.relative_mismatch < 5e-6
    @assert norm(result.rigid_rotation) < 5e-6
    @assert norm(result.traction_free_uniform) < 2e-5
end

println("\nOPT3 weak cross-partial patch audit passed")
