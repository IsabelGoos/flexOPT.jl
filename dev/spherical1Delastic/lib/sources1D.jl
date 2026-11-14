# dev/spherical1Delastic/lib/sources1D.jl
#
# Point sources at the north pole for the 1-D spherical machinery of
# lib/elastic1D.jl: forces AND moment tensors (first moment of the force
# density), the radial source term ("Γ") of every scheme, attenuation
# (complex moduli as in classic DSM, Q varying like any other parameter) and
# readers for the classic DSM spectra (tish / tipsv .spc).
#
# Weak form.  A source density s(x) = f δ(x - x_s) - ∇·(M δ(x - x_s)) acts on
# a test function w through
#     ⟨s, w̄⟩ = f · w̄(x_s) + M : ∇w̄(x_s).
# With w = X(r) conj(P|B|C)_lm and the source on the polar axis, every field
# (P → U, B → V, C → W) receives the radial functional
#     F · X(r_s) + D · ∂_r X(r_s) + Rr · X(r_s)/r_s ,
# with the exact (F, D, Rr) of `poleSourceCoefficients` for |m| ≤ 2
# (sympy, every l: proto/mt_source.py, pole_check.py; matches tish `calg2`
# and tipsv to 0.1-0.2 % once the l = 1 toroidal term is kept).
#
# Units (as `sourceLoads`): K_ξ x = -(1/R) [F·w0 + (D/R)·w1 + Rr/(R ξs)·w0],
# with w0, w1 the discrete X(ξs), ∂ξX(ξs) functionals.  M in (g/cm³)(km/s)²km³
# = 1e25 dyn·cm, f in (g/cm³)(km/s)²km² -- the classic DSM output units.

using LinearAlgebra, SparseArrays, Printf

# ---------------------------------------------------------------------------
# 1. pole coefficients
# ---------------------------------------------------------------------------

"Point source: force f = (f_r, f_θ, f_φ), moment tensor M = (Mrr, Mrθ, Mrφ, Mθθ, Mθφ, Mφφ)."
pointSource(; f=(0.0, 0.0, 0.0), M=(0.0, 0.0, 0.0, 0.0, 0.0, 0.0)) =
    (f=Tuple(Float64.(collect(f))), M=Tuple(Float64.(collect(M))))

"""
    poleSourceCoefficients(l, m, src) -> (P=(F, D, Rr), B=(F, D, Rr), C=(F, D, Rr))

Exact radial functionals of a point source at the north pole (θ = 0, the
local θ̂ along φ = 0) for the fully normalised harmonics of lib/elastic1D.jl
(Condon–Shortley, P = r̂Y, B = k⁻¹∇₁Y, C = -k⁻¹r̂×∇₁Y): field `P|B|C`
receives F·X(r_s) + D·∂_rX(r_s) + Rr·X(r_s)/r_s.  Zero for |m| > 2.
With c = √((2l+1)/4π), k = √(l(l+1)), q = √((l-1)(l+2)):
  m = 0   P: F = c f_r, D = c Mrr, Rr = c (Mθθ+Mφφ);  B: Rr = -c k/2 (Mθθ+Mφφ)
  |m| = 1 P: Rr = c k/2 (-m Mrθ + i Mrφ)
          B: F = c/2 (-m f_θ + i f_φ),  D = -Rr = c/2 (-m Mrθ + i Mrφ)
          C: F = c/2 ( m f_φ + i f_θ),  D = -Rr = c/2 ( m Mrφ + i Mrθ)
  |m| = 2 B: Rr = c q/4 (Mθθ - Mφφ - 2i sgn(m) Mθφ)
          C: Rr = c q/4 (-2Mθφ - i sgn(m)(Mθθ - Mφφ))
"""
function poleSourceCoefficients(l::Integer, m::Integer, src)
    z = zero(ComplexF64)
    PF = PD = PR = BF = BD = BR = CF = CD = CR = z
    abs(m) > min(l, 2) && return (P=(z, z, z), B=(z, z, z), C=(z, z, z))
    c = sqrt((2l + 1) / (4π)); k = sqrt(l * (l + 1)); q = sqrt(max((l - 1) * (l + 2), 0))
    fr, fθ, fφ = src.f
    Mrr, Mrθ, Mrφ, Mθθ, Mθφ, Mφφ = src.M
    if m == 0
        PF = c * fr; PD = c * Mrr; PR = c * (Mθθ + Mφφ)
        BR = -c * k / 2 * (Mθθ + Mφφ)
    elseif abs(m) == 1
        PR = c * k / 2 * (-m * Mrθ + im * Mrφ)
        BF = c / 2 * (-m * fθ + im * fφ); BD = c / 2 * (-m * Mrθ + im * Mrφ); BR = -BD
        CF = c / 2 * (m * fφ + im * fθ);  CD = c / 2 * (m * Mrφ + im * Mrθ);  CR = -CD
    else
        s = sign(m)
        BR = c * q / 4 * (Mθθ - Mφφ - 2im * s * Mθφ)
        CR = c * q / 4 * (-2Mθφ - im * s * (Mθθ - Mφφ))
    end
    if l == 0                       # B and C do not exist for l = 0
        BF = BD = BR = CF = CD = CR = z
    end
    return (P=(PF, PD, PR), B=(BF, BD, BR), C=(CF, CD, CR))
end

# ---------------------------------------------------------------------------
# 2. radial source weights ("Γ")
# ---------------------------------------------------------------------------

"Default source stencil order p (p+1 nodes) of a scheme: its formal order of accuracy."
sourceOrderOf(scheme) = occursin("OPT5", scheme) ? 4 : occursin("OPT4", scheme) ? 3 : 2

"""
    sourceRadialWeights(sys, rs_km; method=:lagrange, p=2, gamma=nothing)

(B0, B1, ξs): nd × nf matrices; column f carries, on the dofs of field f of
the solid zone containing ξs, the weights of the functionals
X ↦ X(ξs) (B0) and X ↦ ∂ξX(ξs) (B1).
  :hat      linear-spline values and slopes (Galerkin DSM: exact for the
            linear trial space, but ∂ξX is piecewise constant, O(h) and
            erratic in the position of ξs inside its element);
  :lagrange Fornberg weights on the p+1 nodes nearest ξs inside the zone
            (one-sided near the zone ends): the dual "derivative with
            respect to the source position" form, O(h^p);
  :gamma    the native flexOPT RHS recipe (`sourceGammaStencils`): its node
            stencils, Lagrange-interpolated (order p) to an off-node ξs.
"""
function sourceRadialWeights(sys::ElasticGlobalSystem, rs_km; method=:lagrange, p=2, gamma=nothing)
    ξs = rs_km / sys.R
    nd = size(sys.T, 1); nf = sys.kind === :PSV ? 2 : 1
    cands = [iz for (iz, seg) in enumerate(sys.segments)
             if !seg.fluid && seg.ξ[1] <= ξs <= seg.ξ[end]]
    isempty(cands) && error("source radius $(rs_km) km is not inside a solid zone")
    iz = last(cands)                       # on a discontinuity: the upper zone
    ξ = sys.segments[iz].ξ; n = length(ξ)
    nodes = Int[]; w0 = Float64[]; w1 = Float64[]
    if method === :hat
        j = clamp(searchsortedlast(ξ, ξs), 1, n - 1); h = ξ[j+1] - ξ[j]
        nodes = [j, j + 1]; w0 = [(ξ[j+1] - ξs) / h, (ξs - ξ[j]) / h]; w1 = [-1 / h, 1 / h]
    elseif method === :lagrange
        idx = _nearestStencil(ξ, ξs, p)
        nodes = collect(idx)
        w0 = fornbergWeights(ξs, ξ[idx], 0); w1 = fornbergWeights(ξs, ξ[idx], 1)
    elseif method === :gamma
        isnothing(gamma) && error("method = :gamma needs gamma = sourceGammaStencils(config)")
        h = ξ[2] - ξ[1]
        idx = _nearestStencil(ξ, ξs, p)
        ℓ = fornbergWeights(ξs, ξ[idx], 0)          # δ(ξ - ξs) ≈ Σ ℓ_k δ(ξ - ξ_k)
        acc0 = Dict{Int,Float64}(); acc1 = Dict{Int,Float64}()
        for (a, kk) in enumerate(idx), (o, g0, g1) in zip(gamma.offsets, gamma.g0, gamma.g1)
            j = kk + o
            1 <= j <= n || error("Γ stencil leaves the source zone: move the source or use :lagrange")
            acc0[j] = get(acc0, j, 0.0) + ℓ[a] * g0
            acc1[j] = get(acc1, j, 0.0) + ℓ[a] * g1 / h
        end
        nodes = sort(collect(union(keys(acc0), keys(acc1))))
        w0 = [get(acc0, j, 0.0) for j in nodes]; w1 = [get(acc1, j, 0.0) for j in nodes]
    else
        error("method must be :hat, :lagrange or :gamma")
    end
    B0 = zeros(nd, nf); B1 = zeros(nd, nf)
    for f in 1:nf, (a, i) in enumerate(nodes)
        B0[sys.dofmap[iz][f, i], f] += w0[a]
        B1[sys.dofmap[iz][f, i], f] += w1[a]
    end
    return (B0=B0, B1=B1, ξs=ξs, zone=iz, nodes=nodes, w0=w0, w1=w1)
end

"p+1 consecutive node indices centred on ξs (uniform zone spacing), clipped to 1:n."
function _nearestStencil(ξ, ξs, p)
    n = length(ξ)
    n >= p + 1 || error("zone has $(n) nodes, stencil needs $(p + 1)")
    s = (ξs - ξ[1]) / (ξ[2] - ξ[1]) + 1          # fractional node index
    lo = clamp(round(Int, s - p / 2), 1, n - p)
    return lo:lo+p
end

# --- native flexOPT Γ: an SH entry whose RHS carries a moment-like flux -------

const SH_SOURCE_EQUATION = "1DsphericalSH_TI_Freq_src"

"""
famousEquation `"1DsphericalSH_TI_Freq_src"`: the SH operator of
`"1DsphericalSH_TI_Freq"` with the r²-weighted source term of a point force
AND a moment-like flux,

    ∂x(Lr² ∂x W) - ∂x(Lr W) + Lr ∂x W + ω² ρr² W - Nₗ W  =  fr² - ∂x(gr²)

(fr² = r²F + rG₀, gr² = r²G₁ for the functional F X + G₁ X' + G₀ X/r).  Only
used to read flexOPT's own right-hand-side recipe Γ (`prepared.R_force`):
the stencils that the scheme applies to δ and to -∂δ.
"""
function flexOPT.famousEquation(::Val{:eq_1DsphericalSH_TI_Freq_src})
    @variables ω
    @variables W(x) ρr²(x) Lr²(x) Lr(x) Nₗ(x) fr²(x) gr²(x)
    expr = ∂x(Lr² * ∂x(W)) - ∂x(Lr * W) + Lr * ∂x(W) + ω^2 * ρr² * W - Nₗ * W
    exprs = expr
    fields = W
    vars = (ρr², Lr², Lr, Nₗ, ω)
    extexprs = fr² - ∂x(gr²)
    extfields = (fr², gr²)
    # one RHS variable only: the coefficients of fr² and ∂x(gr²) are pure numbers,
    # and src/fullyNumericalOPT/makeCostFunctions.jl (frequency branch, ~l. 857)
    # pads rhs.CartesianDependencies with NtypeofFields (= 1) columns
    extvars = fr²
    coordinates = (x,)
    ∂, ∂² = usefulPartials(coordinates)
    return exprs, fields, vars, extexprs, extfields, extvars, coordinates, ∂, ∂²
end

"""
    sourceGammaStencils(config; n=25, verbose=true)

Native flexOPT Γ of `config` (Δ = 1, scaled coordinate): the RHS recipe of
`"1DsphericalSH_TI_Freq_src"` applied to a unit nodal value of fr² (→ δ) and
of gr² (→ -∂δ) at the middle node of an n-node segment, read off
`prepared.R_force` through evaluateLinearSystem.  Calibrated by the moment
conditions of the exact functionals (Σ g0 = 1, Σ o·g1 = 1 in node units), so
that only the SHAPE of the stencils is compared; the raw sums are returned
too.  `offsets` are relative to the source node.
"""
function sourceGammaStencils(config; n=25, verbose=true)
    rec = shRecipe(withRecipeExtras(config, Dict("equationRevisionSource" => "v1"));
                   equation=SH_SOURCE_EQUATION)
    models = Any[ones(n), ones(n), zeros(n), ones(n), 1.0]   # any LHS values: Γ does not depend on them
    modelFam = (models=models, modelName="sphericalSHsrc", modelPoints=(n), Δ=(1.0))
    params = Dict{String,Any}("optRec" => rec, "modelFam" => modelFam,
                              "absorbingBoundaries" => nothing, "maskedRegionInSpace" => nothing,
                              "representation" => "matrixfree")
    numOps = numericalOperatorConstruction(params)["numOperators"]
    prepared = prepareLinearSystem(numOps)
    kf = copy(prepared.known_lhs_template)
    F = copy(prepared.known_rhs_template)
    size(F, 2) == 2 || error("expected 2 force fields (fr², gr²), got $(size(F, 2))")
    size(F, 3) == 1 || @warn "force template has $(size(F, 3)) time slots; using the last"
    k = (n + 1) ÷ 2
    col(field) = begin
        F .= 0.0; F[k, field, end] = 1.0
        _, b = evaluateLinearSystem(prepared, kf, F; sparse_output=false)
        collect(b)
    end
    v0 = col(1); v1 = col(2)
    o = collect(1:length(v0)) .- k
    keep = findall(i -> abs(v0[i]) > 1e-13 * maximum(abs.(v0)) || abs(v1[i]) > 1e-13 * maximum(abs.(v1)), eachindex(v0))
    s0 = sum(v0); s1 = sum(o .* v1)
    g0 = v0[keep] ./ s0; g1 = v1[keep] ./ s1
    out = (name=config.name, offsets=o[keep], g0=g0, g1=g1, rawSum0=s0, rawMoment1=s1,
           moments0=[sum(o[keep] .^ q .* g0) for q in 0:5],
           moments1=[sum(o[keep] .^ q .* g1) for q in 0:5])
    if verbose
        @printf("%-8s Γ[δ]   offsets %s\n", config.name, string(out.offsets))
        @printf("         g0 = %s   (raw Σ = %.6g)\n", join((@sprintf("%+.6f", x) for x in g0), " "), s0)
        @printf("         g1 = %s   (raw Σo·g = %.6g)\n", join((@sprintf("%+.6f", x) for x in g1), " "), s1)
        @printf("         moments Σ o^q g0, q=0..5: %s\n", join((@sprintf("%+.2e", x) for x in out.moments0), " "))
        @printf("         moments Σ o^q g1, q=0..5: %s\n", join((@sprintf("%+.2e", x) for x in out.moments1), " "))
    end
    return out
end

"Fornberg stencils of the same functionals on p+1 nodes centred on the source node."
function fornbergSourceStencil(p)
    o = collect(-(p ÷ 2):(p - p ÷ 2))
    return (name="Fornberg p=$(p)", offsets=o, g0=fornbergWeights(0.0, Float64.(o), 0),
            g1=fornbergWeights(0.0, Float64.(o), 1))
end

# ---------------------------------------------------------------------------
# 3. exact check: δ and δ' sources in a homogeneous traction-free SH shell
# ---------------------------------------------------------------------------

_dsbj(l, z) = l / z * sphericalBesselJ(l, z) - sphericalBesselJ(l + 1, z)
_dsby(l, z) = l / z * sphericalBesselY(l, z) - sphericalBesselY(l + 1, z)

"""
    shellSourceGreen(Ω, l, ξs; a, ρ, β)

Exact surface toroidal response W(1) of the traction-free shell a ≤ ξ ≤ 1
(isotropic, ρ, β; weak operator Ω²T - H of the K_ξ systems) to the
functionals ℓ(X) = X(ξs) ("force") and ℓ(X) = X'(ξs) ("dprime"), i.e.
K w = -ℓ.  W = -W_lo(ξs or ')·W_hi(1) / (ξs² μ Wr), Wr = W_lo W_hi' - W_lo' W_hi.
"""
function shellSourceGreen(Ω, l, ξs; a, ρ, β)
    μ = ρ * β^2; k = Ω / β
    J(x) = sphericalBesselJ(l, k * x);  Jd(x) = k * _dsbj(l, k * x)
    Y(x) = sphericalBesselY(l, k * x);  Yd(x) = k * _dsby(l, k * x)
    trJ(x) = Jd(x) - J(x) / x;  trY(x) = Yd(x) - Y(x) / x
    cJ, cY = trY(a), -trJ(a);     Wlo(x) = cJ * J(x) + cY * Y(x); Wlod(x) = cJ * Jd(x) + cY * Yd(x)
    dJ, dY = trY(1.0), -trJ(1.0); Whi(x) = dJ * J(x) + dY * Y(x); Whid(x) = dJ * Jd(x) + dY * Yd(x)
    C = ξs^2 * μ * (Wlo(ξs) * Whid(ξs) - Wlod(ξs) * Whi(ξs))
    return (force=-Wlo(ξs) * Whi(1.0) / C, dprime=-Wlod(ξs) * Whi(1.0) / C)
end

"""
    shellSourceConvergence(schemes; configs, methods, l=10, kΩ=20.3, ρ=3.3, β=4.5,
        rcmb=3480.0, depth=590.0, nIntervals=[50,100,200,400], recipes=Dict())

Relative error of the surface response to δ (force) and δ' (moment) sources
in the homogeneous SH shell, per scheme × radial source method × grid.
`methods`: vector of (label, method, p) with method :hat | :lagrange, or
(label, :gamma, p, stencil) for a fixed Γ stencil, or (label, :gammaOwn, nothing)
for each scheme's own Γ (`gammas[scheme]`, p = sourceOrderOf(scheme)); p = nothing
→ sourceOrderOf(scheme).
"""
function shellSourceConvergence(schemes; configs, methods, l=10, kΩ=20.3, ρ=3.3, β=4.5,
                                 R=6371.0, rcmb=3480.0, depth=590.0,
                                 nIntervals=[50, 100, 200, 400], recipes=Dict{Any,Any}(),
                                 gammas=Dict{String,Any}(), nearRows=:fallback)
    model = layeredPSVModel(; R=R, radii=[rcmb, R], ρ=[ρ], Vp=[sqrt(3) * β], Vs=[β])
    a = rcmb / R; ξs = (R - depth) / R; Ω = kΩ * β
    ex = shellSourceGreen(Ω, l, ξs; a=a, ρ=ρ, β=β)
    @printf("exact (l=%d, Ω=%.3f, ξs=%.5f): force %.8e  dprime %.8e\n", l, Ω, ξs, ex.force, ex.dprime)
    rows = NamedTuple[]
    for n in nIntervals
        grid = elasticRadialGrid(model; h=(R - rcmb) / n, region=:outerSolid, minIntervals=4)
        for sc in schemes
            sys = elasticGlobalSystem(:SH, model, grid; scheme=sc, configs=configs, recipes=recipes,
                                      nearRows=nearRows)
            Kf = lu(elasticOperator(sys, Ω, l)); Ws = surfaceDofs(sys)[1]
            for mt in methods
                label, meth, p = mt[1], mt[2], mt[3]
                isnothing(p) && (p = sourceOrderOf(sc))
                γ = length(mt) >= 4 ? mt[4] : nothing
                if meth === :gammaOwn
                    haskey(gammas, sc) || continue
                    meth = :gamma; γ = gammas[sc]
                end
                (meth === :gamma && isnothing(γ)) && continue
                w = try
                    sourceRadialWeights(sys, R - depth; method=meth, p=p, gamma=γ)
                catch err
                    @warn "skipped" sc label err; continue
                end
                eF = abs((Kf \ (-w.B0[:, 1]))[Ws] - ex.force) / abs(ex.force)
                eD = abs((Kf \ (-w.B1[:, 1]))[Ws] - ex.dprime) / abs(ex.dprime)
                push!(rows, (scheme=sc, method=label, n=n, h_km=(R - rcmb) / n, errForce=eF, errDprime=eD))
                @printf("%-8s %-14s n=%4d  force %.3e  dprime %.3e\n", sc, label, n, eF, eD)
            end
        end
    end
    return rows
end

"Observed orders log2(e(n)/e(2n)) for each (scheme, method)."
function sourceConvergenceOrders(rows)
    for key in unique([(r.scheme, r.method) for r in rows])
        rr = sort(filter(r -> (r.scheme, r.method) == key, rows); by=r -> r.n)
        length(rr) < 2 && continue
        oF = [log2(rr[i].errForce / rr[i+1].errForce) for i in 1:length(rr)-1]
        oD = [log2(rr[i].errDprime / rr[i+1].errDprime) for i in 1:length(rr)-1]
        @printf("%-8s %-14s orders  force %s | dprime %s\n", key[1], key[2],
                join((@sprintf("%5.2f", o) for o in oF), " "), join((@sprintf("%5.2f", o) for o in oD), " "))
    end
end

# ---------------------------------------------------------------------------
# 4. attenuation: complex moduli as in classic DSM (tish/tipsv calcoef)
# ---------------------------------------------------------------------------
#   c(ω, Q) = (1 + ln(ω/2π)/(πQ) + i/(2Q))² = 1 + a(ω)/Q + b(ω)/Q²,
#   a = 2x + i, b = (x + i/2)², x = ln(ω/2π)/π   (reference 1 Hz; ω = 0 → x = 0)
# tipsv split of the TI moduli: κ-parts A-4N/3, C-4N/3, F+2N/3 with Qκ,
# μ-parts 4N/3, 4N/3, -2N/3, L, N with Qμ; fluids: 1/λ → 1/(c_κ λ)
# = (1/λ)(1 - a/Qκ + (a²-b)/Qκ² + O(Qκ⁻³)).  The operators are linear in the
# material fields, so two extra basis sets (moduli weighted by 1/Q and 1/Q²)
# give the exact complex-modulus operator at every ω with Q(r) arbitrary.

"(a(ω), b(ω)) of the classic DSM attenuation factor; ω real (rad/s)."
function dsmAttenuationAB(ω)
    x = ω > 0 ? log(ω / (2π)) / π : 0.0
    return (a=2x + im, b=(x + im / 2)^2)
end

"Material function for the 1/Q^order basis sets (see the comment above)."
function anelasticMaterials(order::Integer)
    return function (model, iz, ξ)
        m = elasticMaterials(model, iz, ξ)
        qμ = _poly(model.C_Qμ[iz, :], ξ); qκ = _poly(model.C_Qκ[iz, :], ξ)
        iμ = qμ > 0 ? qμ^(-order) : 0.0
        iκ = qκ > 0 ? qκ^(-order) : 0.0
        N43 = 4m.N / 3
        return (ρ=0.0, A=(m.A - N43) * iκ + N43 * iμ, C=(m.C - N43) * iκ + N43 * iμ,
                F=(m.F + m.N * 2 / 3) * iκ - (m.N * 2 / 3) * iμ, L=m.L * iμ, N=m.N * iμ,
                λ=m.λ, λinv=iκ / m.λ, ρinv=0.0)
    end
end

"""
    anelasticSystems(kind, model, grid; kw...)  -> (q1, q2)

The 1/Q and 1/Q² basis systems (same assembly as `elasticGlobalSystem`,
with `anelasticMaterials(1|2)`); use with `anelasticOperator`.
"""
anelasticSystems(kind, model, grid; kw...) =
    (q1=elasticGlobalSystem(kind, model, grid; materials=anelasticMaterials(1), kw...),
     q2=elasticGlobalSystem(kind, model, grid; materials=anelasticMaterials(2), kw...))

"""
    anelasticOperator(sys, q, Ω, ω, l)

K(Ω, l) with complex moduli at the real frequency ω (classic DSM):
  Ω²(T - a T₁ + (a² - b) T₂) + Σ_key (H + a H₁ + b H₂)  + Ω Cpl
(T₁, T₂ carry only the fluid 1/λ terms, H₁, H₂ only the solid moduli).
`q === nothing` → elastic.
"""
function anelasticOperator(sys::ElasticGlobalSystem, q, Ω, ω, l)
    isnothing(q) && return elasticOperator(sys, Ω, l)
    a, b = dsmAttenuationAB(ω)
    Λ = l * (l + 1); k = sqrt(Λ)
    Hs(s) = s.H0 .+ k .* s.Hk .+ Λ .* s.HΛ
    return Ω^2 .* (sys.T .- a .* q.q1.T .+ (a^2 - b) .* q.q2.T) .+
           Hs(sys) .+ a .* Hs(q.q1) .+ b .* Hs(q.q2) .+ Ω .* sys.Cpl
end

# ---------------------------------------------------------------------------
# 5. synthetics for point sources (forces and moment tensors) at the pole
# ---------------------------------------------------------------------------

"""
    pointSourceSynthetics(psv, sh, R; rs_km, sources, radial=(method=:lagrange, p=2),
        anelastic=nothing, receiversθdeg=10:10:170, φr=π/2, tlen=3276.8, fmax=0.05,
        omegaiTlen=1e-4, lmaxOf=ω -> ceil(Int, ω*R/3.0) + 30, lcap=1500, verbose=true)

Surface (u_r, u_θ, u_φ) on the meridian φ = φr for every source in
`sources` (Dict name => pointSource(...)), source exactly on the polar axis
(no 1e-5 clamp: the analytic pole coefficients are used, |m| ≤ 2).  One LU per
(ω, l) serves all sources and all radial functionals (columns U·X, U·X',
V·X, V·X' for P-SV; W·X, W·X' for SH), frequencies in parallel over Julia
threads (`julia -t auto`).  `parts` selects :psv and/or :sh; `split = true`
also returns "name|psv" and "name|sh".  `anelastic = (psv=(q1=..,q2=..), sh=...)`
switches on the complex moduli.  `omegaiTlen = 1e-4`: ω_I·tlen = ln 10⁴
(classic DSM uses 1e-2; that leaves up to ~10 % of wrapped antipodal surface
waves in the first 2000 s, see p6).
"""
function pointSourceSynthetics(psv::ElasticGlobalSystem, sh, R; rs_km, sources,
                               radial=(method=:lagrange, p=2), anelastic=nothing,
                               receiversθdeg=10:10:170, φr=π / 2,
                               tlen=3276.8, fmax=0.05, omegaiTlen=1e-4,
                               lmaxOf=ω -> ceil(Int, ω * R / 3.0) + 30, lcap=1500,
                               verbose=true, parts=(:psv, :sh), split=false)
    ωI = -log(omegaiTlen) / tlen
    kmax = floor(Int, fmax * tlen)
    ωs = [2π * k / tlen for k in 0:kmax]
    θr = deg2rad.(collect(receiversθdeg)); nr = length(θr)
    lmaxAll = min(lcap, maximum(lmaxOf(ω) for ω in ωs))
    recTabs = [normalisedLegendre(lmaxAll, 2, θ) for θ in θr]
    names = collect(keys(sources))
    rad = (; method=:lagrange, p=2, gamma=nothing, radial...)
    wP = sourceRadialWeights(psv, rs_km; method=rad.method, p=rad.p, gamma=rad.gamma)
    ξs = wP.ξs
    LP = complex.(-(1 / R) .* hcat(wP.B0[:, 1], wP.B1[:, 1], wP.B0[:, 2], wP.B1[:, 2]))
    LS = nothing
    if !isnothing(sh)
        wS = sourceRadialWeights(sh, rs_km; method=rad.method, p=rad.p, gamma=rad.gamma)
        LS = complex.(-(1 / R) .* hcat(wS.B0[:, 1], wS.B1[:, 1]))
    end
    # receiver sums per source, per load column j, per l, per receiver
    Z() = [zeros(ComplexF64, lmaxAll + 1, nr) for _ in 1:4]
    S = Dict(nm => (Y=Z(), Bθ=Z(), Bφ=Z(), Cθ=Z(), Cφ=Z()) for nm in names)
    for nm in names, l in 0:lmaxAll, m in -min(l, 2):min(l, 2)
        c = poleSourceCoefficients(l, m, sources[nm])
        sP = (c.P[1] + c.P[3] / (R * ξs), c.P[2] / R, c.B[1] + c.B[3] / (R * ξs), c.B[2] / R)
        sC = (c.C[1] + c.C[3] / (R * ξs), c.C[2] / R)
        (all(iszero, sP) && all(iszero, sC)) && continue
        d = S[nm]
        for (ir, tab) in enumerate(recTabs)
            Y, Bθ, Bφ, Cθ, Cφ = vectorHarmonics(tab, l, m, φr)
            for j in 1:4
                d.Y[j][l+1, ir] += sP[j] * Y; d.Bθ[j][l+1, ir] += sP[j] * Bθ; d.Bφ[j][l+1, ir] += sP[j] * Bφ
            end
            for j in 1:2
                d.Cθ[j][l+1, ir] += sC[j] * Cθ; d.Cφ[j][l+1, ir] += sC[j] * Cφ
            end
        end
    end
    Us, Vs = surfaceDofs(psv)
    Ws = isnothing(sh) ? 0 : surfaceDofs(sh)[1]
    qP = isnothing(anelastic) ? nothing : anelastic.psv
    qS = (isnothing(anelastic) || isnothing(sh)) ? nothing : anelastic.sh
    doP = :psv in parts; doS = (:sh in parts) && !isnothing(sh)
    nω = length(ωs)
    mk() = Dict(nm => (ur=zeros(ComplexF64, nω, nr), uθ=zeros(ComplexF64, nω, nr),
                       uφ=zeros(ComplexF64, nω, nr)) for nm in names)
    outP = mk(); outS = mk()
    # frequencies are independent: one thread per ω (start Julia with -t auto);
    # each thread writes only its own row iω.  BLAS threads would oversubscribe.
    Threads.nthreads() > 1 && LinearAlgebra.BLAS.set_num_threads(1)
    ndone = Threads.Atomic{Int}(0); t0 = time()
    Threads.@threads :dynamic for iω in 1:nω
        ω = ωs[iω]
        Ωc = (ω - im * ωI) * R
        lmax = min(lmaxAll, lmaxOf(ω))
        for l in 0:lmax
            if doP
                G = lu(anelasticOperator(psv, qP, Ωc, ω, l)) \ LP
                gU = G[Us, :]; gV = G[Vs, :]
                for nm in names
                    d = S[nm]; o = outP[nm]
                    for j in 1:4
                        @views o.ur[iω, :] .+= gU[j] .* d.Y[j][l+1, :]
                        @views o.uθ[iω, :] .+= gV[j] .* d.Bθ[j][l+1, :]
                        @views o.uφ[iω, :] .+= gV[j] .* d.Bφ[j][l+1, :]
                    end
                end
            end
            if doS && l >= 1
                gW = (lu(anelasticOperator(sh, qS, Ωc, ω, l)) \ LS)[Ws, :]
                for nm in names
                    d = S[nm]; o = outS[nm]
                    for j in 1:2
                        @views o.uθ[iω, :] .+= gW[j] .* d.Cθ[j][l+1, :]
                        @views o.uφ[iω, :] .+= gW[j] .* d.Cφ[j][l+1, :]
                    end
                end
            end
        end
        k = Threads.atomic_add!(ndone, 1) + 1
        verbose && (k % 20 == 0 || k == nω) &&
            @printf("  %s: %d/%d frequencies, %.1f s (%d threads)\n", psv.label, k, nω,
                    time() - t0, Threads.nthreads())
    end
    spec(o) = (ω=ωs, ωI=ωI, tlen=tlen, θdeg=collect(receiversθdeg), φr=φr, ur=o.ur, uθ=o.uθ, uφ=o.uφ)
    res = Dict{String,Any}()
    for nm in names
        res[nm] = spec((ur=outP[nm].ur .+ outS[nm].ur, uθ=outP[nm].uθ .+ outS[nm].uθ,
                        uφ=outP[nm].uφ .+ outS[nm].uφ))
        if split
            res[nm * "|psv"] = spec(outP[nm]); res[nm * "|sh"] = spec(outS[nm])
        end
    end
    return res
end

# ---------------------------------------------------------------------------
# 6. classic DSM spectra (tish / tipsv .spc, VERSION172 format)
# ---------------------------------------------------------------------------

"(tlen, np, ωI, idx, u::Matrix{ComplexF64} (nω × 3: r, θ, φ)) of one .spc file."
function readClassicSpc(path)
    L = readlines(path)
    tlen = parse(Float64, strip(L[2])); np = parse(Int, strip(L[3]))
    ωI = parse(Float64, split(L[4])[1])
    idx = Int[]; rows = Vector{ComplexF64}[]
    k = 6
    while k + 2 <= length(L) && !isempty(strip(L[k]))
        a = split(L[k]); b = split(L[k+1]); c = split(L[k+2])
        push!(idx, parse(Int, a[1]))
        push!(rows, [complex(parse(Float64, a[2]), parse(Float64, a[3])),
                     complex(parse(Float64, b[1]), parse(Float64, b[2])),
                     complex(parse(Float64, c[1]), parse(Float64, c[2]))])
        k += 3
    end
    return (tlen=tlen, np=np, ωI=ωI, idx=idx, u=permutedims(reduce(hcat, rows)))
end

"""
    classicSpectra(dir, θdeg; parts=("sh", "psv"))

Sum of the tish (sh_DDD.spc) and tipsv (psv_DDD.spc) spectra in `dir`, in the
spec format of `pointSourceSynthetics` (usable by `elasticTimeSeries`).
"""
function classicSpectra(dir, θdeg; parts=("sh", "psv"))
    θdeg = collect(θdeg); nr = length(θdeg)
    s1 = readClassicSpc(joinpath(dir, @sprintf("%s_%03d.spc", parts[1], round(Int, θdeg[1]))))
    nω = length(s1.idx)
    ur = zeros(ComplexF64, nω, nr); uθ = zeros(ComplexF64, nω, nr); uφ = zeros(ComplexF64, nω, nr)
    for (ir, θ) in enumerate(θdeg), part in parts
        s = readClassicSpc(joinpath(dir, @sprintf("%s_%03d.spc", part, round(Int, θ))))
        ur[:, ir] .+= s.u[:, 1]; uθ[:, ir] .+= s.u[:, 2]; uφ[:, ir] .+= s.u[:, 3]
    end
    ωs = 2π .* s1.idx ./ s1.tlen
    return (ω=ωs, ωI=s1.ωI, tlen=s1.tlen, θdeg=θdeg, ur=ur, uθ=uθ, uφ=uφ)
end

"Relative L2 spectral misfit per component, ‖A - B‖/‖B‖ over the common frequencies."
function spectralMisfit(A, B)
    n = min(length(A.ω), length(B.ω))
    maximum(abs.(A.ω[1:n] .- B.ω[1:n])) < 1e-9 || error("frequency grids differ")
    abs(A.ωI - B.ωI) < 1e-9 * max(A.ωI, 1e-12) || @warn "ω_I differs" A.ωI B.ωI
    f(c) = (a = getfield(A, c)[1:n, :]; b = getfield(B, c)[1:n, :];
            norm(b) == 0 ? NaN : norm(a .- b) / norm(b))
    return (Z=f(:ur), R=f(:uθ), T=f(:uφ))
end
