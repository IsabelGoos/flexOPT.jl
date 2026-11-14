#!/usr/bin/env julia
#
# helpers/testTestFunctionMoments.jl -- checks src/motorsOPT/testFunctionMoments.jl
# against the results derived with sympy and confirmed in flexOPT (Nov 2026).
#
#   julia helpers/testTestFunctionMoments.jl      (no package needed)

include(joinpath(@__DIR__, "..", "src", "motorsOPT", "testFunctionMoments.jl"))
using .TestFunctionMoments
using Test

const Q = Rational{BigInt}

@testset "test-function moments" begin
    @testset "B-spline moments" begin
        hat1 = hatTestFunction(1)
        M = testFunctionMoments(hat1; mmax=6)
        @test M[1] == 1 && M[3] == 1 // 6 && M[5] == 1 // 15 && M[7] == 1 // 28
        @test M[2] == 0 && M[4] == 0                       # symmetric
        # cubic, d = 1: M₂ = κ₂ = 1/3, M₄ = κ₄ + 3κ₂² = -1/30 + 1/3
        Mc = testFunctionMoments(centredBspline(3, 1); mmax=4)
        @test Mc[3] == 1 // 3 && Mc[5] == 3 // 10
        @test testFunctionMoments(hatTestFunction(2); mmax=4)[3] == 4 // 6
    end

    @testset "predicted orders" begin
        n3, n5, n7 = defaultStencilNodes.((3, 5, 7))
        @test predictedOrder(n3, hatTestFunction(1)).order == 2
        @test predictedOrder(n5, hatTestFunction(1)).order == 4
        @test predictedOrder(n7, hatTestFunction(1)).order == 6
        # constant κ only: (0, b) rows -- the hat is O(h^(P+1))
        lead3 = predictedOrder(n3, hatTestFunction(1)).leading
        @test all(ab -> ab[1] > 0, first.(lead3))        # no (0, b) term: κ-constant part is O(h⁴)
        mix5 = mixTestFunctions(Q[4 // 15, 11 // 15], [hatTestFunction(1), hatTestFunction(2)])
        @test predictedOrder(n5, mix5).order == 6
        mix7 = mixTestFunctions(Q[40 // 231, 191 // 231], [hatTestFunction(1), hatTestFunction(2)])
        @test predictedOrder(n7, mix7).order == 8
    end

    @testset "optimal weights" begin
        n3, n5, n7 = defaultStencilNodes.((3, 5, 7))
        r5 = optimalTestFunctionWeights(n5, [hatTestFunction(1), hatTestFunction(2)])
        @test r5.order == 6 && r5.weights == Q[4 // 15, 11 // 15] && isempty(r5.freeDirections)
        r7 = optimalTestFunctionWeights(n7, [hatTestFunction(1), hatTestFunction(2)])
        @test r7.order == 8 && r7.weights == Q[40 // 231, 191 // 231]
        # cubic weight forced to zero when added to hat₁ + hat₂ (P = 5)
        r5c = optimalTestFunctionWeights(n5, [hatTestFunction(1), hatTestFunction(2), centredBspline(3, 1)])
        @test r5c.order == 6 && r5c.weights == Q[4 // 15, 11 // 15, 0]
        # P = 3: nothing beats order 2 in heterogeneous media
        r3 = optimalTestFunctionWeights(n3, [hatTestFunction(1), centredBspline(2, 1 // 2), centredBspline(3, 1 // 2)])
        @test r3.order == 2
        # P = 3 with exact κ and f: still order 2
        @test optimalTestFunctionWeights(n3, [hatTestFunction(1), centredBspline(3, 1 // 2)];
            kappa=:exact, source=:exact).order == 2
    end

    @testset "mass term (SH frequency)" begin
        n3, n5, n7 = defaultStencilNodes.((3, 5, 7))
        @test predictedOrder(n3, hatTestFunction(1); mass=true).order == 2
        @test predictedOrder(n5, hatTestFunction(1); mass=true).order == 4
        r5 = optimalTestFunctionWeights(n5, [hatTestFunction(1), hatTestFunction(2)]; mass=true)
        @test r5.order == 6 && r5.weights == Q[4 // 15, 11 // 15]
        r7 = optimalTestFunctionWeights(n7, [hatTestFunction(1), hatTestFunction(2)]; mass=true)
        @test r7.order == 8 && r7.weights == Q[40 // 231, 191 // 231]
        # the mass rows alone first appear at h^(P+1)
        mix5 = mixTestFunctions(Q[4 // 15, 11 // 15], [hatTestFunction(1), hatTestFunction(2)])
        lead = predictedOrder(n5, mix5; mass=true).leading
        @test any(p -> first(p) isa Tuple{Symbol,Int,Int}, lead)
    end

    @testset "flexOPT knots" begin
        @test flexOPTRefPoints([-2, 0, 2], 5) == [1, 3, 5]
        @test flexOPTRefPoints([-2, 0, 2], 7) == [2, 4, 6]
        @test_throws ArgumentError flexOPTRefPoints([-3, 0, 3], 5)
    end
end

# A readable summary
for P in (3, 5, 7)
    nodes = defaultStencilNodes(P)
    r = optimalTestFunctionWeights(nodes, [hatTestFunction(1), hatTestFunction(2)])
    println("P = $P: hat₁ order ", predictedOrder(nodes, hatTestFunction(1)).order,
        "; best hat₁/hat₂ mixture: order ≥ ", r.order, ", weights ", r.weights)
end
