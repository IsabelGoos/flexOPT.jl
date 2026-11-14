# dev/spherical1Delastic/stages/step1_symbolic.jl
# Ported from notebooks/SphericalTI_SH_oct2026.ipynb, Step 1 (cells sh-step1-03 … sh-step1-11).
# Symbolic checks: cheap, rerun every time; only the brute-force Cartesian
# check (1c) is cached.

# ---- 1a: Cartesian TI equation
@show checkTIisotropicLimitAndLoveNotation()
eq3D, force3D = famousEquations("3DsismoFreqTIradial")
@show eq3D.fields eq3D.vars eq3D.coordinates
for i in 1:3
    s = string(eq3D.exprs[i])
    println("component $i : ", length(s), " characters\n   ", first(s, 400), " …\n")
end

# ---- 1b–1c: spherical components, normalised toroidal harmonic, projection
sym = deriveToroidalSymbolic()
println("E_r  = ", sym.E[1])
println("E_θ  = ", sym.E[2])
println("\nR(r) = ∫ C*·(ρω²u + ∇·σ) dΩ =\n", sym.R)
println("\nexpected =\n", sym.expected)
println("\nsymbolic R - expected (expanded, ×r³) = ",
    Symbolics.expand(sym.R * x^3 - sym.expected * x^3))
for l in (1, 2, 5, 20)
    @show l verifyToroidalProjection(; l=l)
end

# ---- independent check with m ≠ 0 (brute force in Cartesian coordinates): cached
cartesianParams = Dict{String,Any}(
    "stage" => "step1_cartesian",
    "codeVersion" => CODE_VERSION["step1_cartesian"],
    "harmonics" => [:Y21, :Y32],
)
devRedo!("cartesian", cartesianParams, SPH_CACHE_DIR, "step1_cartesianCheck")
cartesianChecks = myProduceOrLoad(cartesianParams, SPH_CACHE_DIR, "step1_cartesianCheck") do cfg
    Dict{String,Any}(string(h) => (@time cartesianToroidalCheck(; harmonic=h)) for h in cfg["harmonics"])
end
for h in cartesianParams["harmonics"]
    @show cartesianChecks[string(h)]
end

# ---- 1d: the radial equation as a famousEquation entry
eqSH, forceSH = famousEquations("1DsphericalSH_TI_Freq")
@show eqSH.exprs eqSH.fields eqSH.vars eqSH.coordinates forceSH.exprs
wf = naturalWeakForm(eqSH)
println("\nvolume terms (expression, derivative order on the test function):")
for t in wf.volume_terms[1]
    println("   ", t.test_derivative_orders, "  =>  ", t.expression)
end
println("\nboundary fluxes (normal coordinate, expression):")
for b in wf.boundary_fluxes[1]
    println("   ", b.normal_coordinate, "  =>  ", b.expression)
end
println("\nexpected: (1,) => -Lr²∂W + Lr W ;  (0,) => Lr ∂W, ω²ρr² W, -Nₗ W")

# Check that the geometric fields reproduce the radial operator
let s = sym.symbols
    @variables ρr²(x) Lr²(x) Lr(x) Nₗ(x) ω
    e = substitute(eqSH.exprs, Dict(ρr² => s.ρ * x^2, Lr² => s.L * x^2, Lr => s.L * x,
        Nₗ => s.L + (s.Λ - 2) * s.N))
    e = substitute(e, Dict(ω => s.ω))
    diff = Symbolics.expand(expand_derivatives(e) - x^2 * sym.expected)
    println("x·(entry(geometric fields) - r²·expected) = ", Symbolics.expand(x * diff))
end
