#!/usr/bin/env julia
#
# helpers/diagnoseMuBoth.jl   (Claude-authored; does not modify src/)
#
# Why do multi-centre μ sets break the Poisson recipe only when they are
# applied to the field AND the material ("both")?
#
# In 1-D the hierarchical Taylor inverse reproduces degree ≤ P-1 exactly and
# exhausts the P nodal values, so the reconstruction at any μ is the same
# interpolating polynomial. With Σ_μ Y_μ = 1 the operator should therefore
# be IDENTICAL for every μ set and every role. This prints the local A row
# (and Γ row) at Δ = 1 for
#   κ ≡ 1          (then the material μ set cannot matter at all), and
#   κ = 1 + 0.1 s  (linear in the local coordinate s: also reproduced exactly),
# for central / field / material / both. Any difference from "central"
# beyond round-off is a bug in the recipe assembly, not a property of μ.
#
# Usage:  julia --project=. helpers/diagnoseMuBoth.jl
#         FLEXOPT_DIAG_POINTS=3,4 FLEXOPT_DIAG_KINDS=nodes,midpoints julia --project=. helpers/diagnoseMuBoth.jl

include(joinpath(@__DIR__, "runPoissonMuBenchmark.jl"))   # main() is guarded

const DIAG_POINTS = env_ints("FLEXOPT_DIAG_POINTS", "3,4")
const DIAG_KINDS = env_syms("FLEXOPT_DIAG_KINDS", "nodes,midpoints")

function local_row(numeric, kappa_of_s)
    mapping = Dict{Any,Any}()
    for (j, offset) in enumerate(numeric.offsets)
        for v in axes(numeric.lhs_varm, 1)
            symbol = numeric.lhs_varm[v, j]
            mapping[symbol] = occursin("κ", string(symbol)) ? kappa_of_s(offset[1]) : 1.0
        end
    end
    A = [evaluate_coefficient(numeric.lhs[j, 1, 1, 1], mapping) for j in eachindex(numeric.offsets)]
    G = [evaluate_coefficient(numeric.rhs[j, 1, 1, 1], Dict{Any,Any}()) for j in eachindex(numeric.offsets)]
    return A, G
end

fmt(v) = join((@sprintf("%+.12f", x) for x in v), "  ")

for P in DIAG_POINTS
    central = mu_geometry(P, :central)
    for kind in DIAG_KINDS
        g = mu_geometry(P, kind)
        cases = [
            ("central", central, central),
            ("field-$(kind)", g, central),
            ("material-$(kind)", central, g),
            ("both-$(kind)", g, g),
        ]
        println("\n================ P = $P, μ set '$(kind)' = $(format_positions(mu_positions(P, g))) (relative to centre)")
        reference = nothing
        for (label, field, material) in cases
            scheme = mu_scheme("OPT$P-$label", :OPT, P,
                label == "central" ? "reference" : split(label, "-")[1],
                Symbol(kind), field, material)
            numeric = compile_recipe(mu_recipe("1DpoissonHetero", 1, scheme, 1.0)["recette"])
            offsets = [o[1] for o in numeric.offsets]
            for (kname, kfun) in (("κ≡1", s -> 1.0), ("κ=1+0.1s", s -> 1.0 + 0.1 * s))
                A, G = local_row(numeric, kfun)
                println(rpad(label, 22), rpad(kname, 10), "offsets ", offsets)
                println("    A  = ", fmt(A), "    ΣA = ", @sprintf("%+.3e", sum(A)))
                kname == "κ≡1" && println("    Γ  = ", fmt(G))
                if label == "central"
                    reference = kname == "κ≡1" ? (A1=A, G1=G) : merge(reference, (A2=A,))
                else
                    Aref = kname == "κ≡1" ? reference.A1 : reference.A2
                    println("    max|A - A_central| = ", @sprintf("%.3e", maximum(abs.(A - Aref))),
                        kname == "κ≡1" ? @sprintf("   max|Γ - Γ_central| = %.3e", maximum(abs.(G - reference.G1))) : "")
                end
            end
        end
    end
end
