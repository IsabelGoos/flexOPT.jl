# dev/spherical1Delastic/stages/p1_symbolic.jl
# P-SV Step P1 (notebooks/SphericalTI_PSV_oct2026.ipynb, cells psv-p1-03 … psv-p1-11).
# Symbolic checks: cheap, rerun every time; only the brute-force Cartesian
# check is cached.

@show checkTIisotropicLimitAndLoveNotation()
symPSV = deriveSpheroidalSymbolic()
println("E_φ = ", symPSV.E[3])
println("\nR_U =\n", symPSV.RU)
println("\nR_V =\n", symPSV.RV)
for l in (1, 2, 5, 20)
    @show l verifySpheroidalProjection(; l=l)
end

# ---- independent m ≠ 0 check, full Cartesian TI operator: cached
psvCartesianParams = Dict{String,Any}(
    "stage" => "p1_cartesian", "codeVersion" => CODE_VERSION["p1_cartesian"],
    "harmonics" => [:Y21, :Y32],
)
devRedo!("cartesian", psvCartesianParams, SPH_CACHE_DIR, "p1_cartesianCheck")
psvCartesianChecks = myProduceOrLoad(psvCartesianParams, SPH_CACHE_DIR, "p1_cartesianCheck") do cfg
    Dict{String,Any}(string(h) => (@time cartesianSpheroidalCheck(; harmonic=h)) for h in cfg["harmonics"])
end
for h in psvCartesianParams["harmonics"]
    @show psvCartesianChecks[string(h)]
end

# ---- fluid (Q = P/ω, Geller & Ohminato 1994)
for l in (1, 3, 20)
    @show l verifyFluidProjection(; l=l)
end

# ---- famousEquation entries and their natural weak forms
eqPSV, forcePSV = famousEquations(PSV_EQUATION)
@show eqPSV.fields eqPSV.vars forcePSV.exprs
wfPSV = naturalWeakForm(eqPSV)
for (ie, name) in enumerate(("U-equation", "V-equation"))
    println("\n$name — volume terms (test derivative order => expression):")
    for t in wfPSV.volume_terms[ie]
        println("   ", t.test_derivative_orders, "  =>  ", t.expression)
    end
    println("$name — boundary fluxes:")
    for b in wfPSV.boundary_fluxes[ie]
        println("   ", b.normal_coordinate, "  =>  ", b.expression)
    end
end
@show verifyPSVEntry(; l=3) verifyPSVEntry(; l=10)
eqQ, forceQ = famousEquations(FLUID_EQUATION)
@show eqQ.exprs eqQ.vars
wfQ = naturalWeakForm(eqQ)
for t in wfQ.volume_terms[1]; println("   ", t.test_derivative_orders, "  =>  ", t.expression); end
for b in wfQ.boundary_fluxes[1]; println("   flux ", b.normal_coordinate, "  =>  ", b.expression); end
