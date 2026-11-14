# dev/spherical1Delastic/lib/elastic1D.jl
#
# Generic 1-D spherical elastic machinery (SH, P-SV, fluid in Q = P/ω) used
# by the P-SV steps (P2, P5) and reusable for SH.  Everything is written in
# ξ = r/R with Ω = ωR, on the K_ξ scale of conventional DSM:
#
#     K(Ω, l) = Ω² T + H0 + k Hk + Λ HΛ + Ω Cpl ,   k = √Λ, Λ = l(l+1)
#
# (linearity of every recipe in its material fields: each basis matrix is
# built ONCE per zone and scheme; Cpl = Geller & Ohminato 1994 solid–fluid
# coupling).  Validated in a Python prototype (conventional DSM version)
# against exact eigenfrequencies of a fluid core + solid mantle (5×5
# determinant) and inner core + fluid + mantle (8×8): real spectrum, h²
# convergence, grid starting at r = 0 without special treatment.
#
# Radial equations (r²-weighted, see lib/sphericalTI_PSV.jl and
# lib/sphericalTI_SH.jl):
#   SH   ∂(r²L W') - ∂(rL W) + rL W' + (ω²ρr² - L - (Λ-2)N) W
#   P-SV ∂(r²C U' + 2rF U - krF V) - 2rF U' - krL V' + (ω²ρr² - 4A + 4N - ΛL) U + k(2A+L-2N) V
#        ∂(r²L V' - rL V + krL U) + krF U' + rL V' + (ω²ρr² - ΛA - L + 2N) V + k(2A+L-2N) U
#   fluid ∂(r²/ρ q') + (ω²r²/λ - Λ/ρ) q
# Coupling (G&O eq. 35, divided by R):  solid below / fluid above (ICB):
#   +Ω ξ_b² in (U-row, q-col) and (q-row, U-col);  fluid below / solid above
#   (CMB): -Ω ξ_a² in both.  The source: K x = -(1/R)·X(ξs)·c (point force).
#
# Requires lib/sphericalTI_SH*.jl and lib/sphericalTI_PSV.jl.

using LinearAlgebra, SparseArrays, Printf

const SH_EQUATION = "1DsphericalSH_TI_Freq"
const PSV_EQUATION = "1DsphericalPSV_TI_Freq"
const FLUID_EQUATION = "1DsphericalFluidQ_Freq"
# bump when the corresponding famousEquation entry changes (enters the recipe cache key)
const EQUATION_REVISIONS = Dict(SH_EQUATION => "v1", PSV_EQUATION => "v1", FLUID_EQUATION => "v1")

# ---------------------------------------------------------------------------
# materials, zones, grids
# ---------------------------------------------------------------------------

"""
    elasticMaterials(model, iz, ξ)

TI parameters from the planet1D zone polynomials at ξ = r/R:
A = ρVph², C = ρVpv², L = ρVsv², N = ρVsh², F = η(A - 2L); in a fluid
λ = ρVpv² (incompressibility).
"""
function elasticMaterials(model, iz, ξ)
    ρ = _poly(model.C_ρ[iz, :], ξ)
    vpv = _poly(model.C_Vpv[iz, :], ξ); vph = _poly(model.C_Vph[iz, :], ξ)
    vsv = _poly(model.C_Vsv[iz, :], ξ); vsh = _poly(model.C_Vsh[iz, :], ξ)
    η = _poly(model.C_η[iz, :], ξ)
    A = ρ * vph^2; C = ρ * vpv^2; L = ρ * vsv^2; N = ρ * vsh^2
    return (ρ=ρ, A=A, C=C, F=η * (A - 2L), L=L, N=N, λ=ρ * vpv^2)
end

"""
    elasticRadialGrid(model; h=10.0, rmin=0.0, minIntervals=12, region=:all)

Segments bottom → top with one uniform h per zone (h_zone = width/ceil(width/h),
at least `minIntervals`), each tagged `fluid`.  `region = :all` (P-SV: inner
core from rmin, fluid outer core, mantle, crust) or `:outerSolid` (SH).
A fluid top layer (ocean) is not supported yet.
"""
function elasticRadialGrid(model; h=10.0, rmin=0.0, minIntervals=12, region=:all)
    R = model.averagedPlanetRadiusInKilometer
    zones = region === :all ? (1:model.nzone) : outerSolidZones(model)
    model.solid_or_fluid[last(zones)] == "S" ||
        error("a fluid top layer (ocean) is not supported yet")
    segs = NamedTuple[]
    for iz in zones
        rb = max(model.bottomRadius[iz], rmin); rt = model.topRadius[iz]
        rt > rb || continue
        hz = h isa AbstractDict ? h[iz] : h isa AbstractVector ? h[iz] : h
        n = max(ceil(Int, (rt - rb) / hz), minIntervals)
        seg = _segment(iz, collect(range(rb, rt; length=n + 1)), R)
        push!(segs, merge(seg, (fluid=model.solid_or_fluid[iz] == "F",)))
    end
    return (R=R, segments=segs, region=region)
end

# ---------------------------------------------------------------------------
# recipes and multi-field flexOPT operators
# ---------------------------------------------------------------------------

"Recipe of `config` for `equation` (CPU, cached; the equation revision enters the key)."
elasticRecipe(config, equation) =
    shRecipe(withRecipeExtras(config, Dict("equationRevisionElastic" => EQUATION_REVISIONS[equation]));
             equation=equation)

"Recipe cache in memory: Dict((name, equation) => recipe)."
function getRecipe!(recipes::AbstractDict, config, equation)
    return get!(recipes, (config.name, equation)) do
        elasticRecipe(config, equation)
    end
end

"""
    flexOPTOperator(optRec, models, nfields)

Dense LHS operator with rows AND columns in field-major order
((field-1)·n + node).  flexOPT returns residual rows interleaved
(equation fastest) and unknown columns field-major (prepareNumericalLinearSystem).
"""
function flexOPTOperator(optRec, models, nfields::Integer)
    A = flexOPTSHMatrix(optRec, models)
    nfields == 1 && return A
    n = size(A, 1) ÷ nfields
    perm = [(ip - 1) * nfields + ie for ie in 1:nfields for ip in 1:n]
    return A[perm, :]
end

# ---------------------------------------------------------------------------
# zone model sets (scaled local coordinate s, unit spacing) for the bases
# ---------------------------------------------------------------------------

"""
    zoneModelSets(kind, model, seg)

flexOPT `models` (in each equation's `vars` order) for the four basis
matrices T (Ω²), H0, Hk (× k), HΛ (× Λ), in the scaled coordinate
s = (ξ - ξ_bottom)/h: second-order fields unchanged, first-order × h,
zeroth-order × h².  Also returns the indices of the first-order fields
(negated when the segment is mirrored).
"""
function zoneModelSets(kind::Symbol, model, seg; materials=elasticMaterials)
    ξ = seg.ξ; h = seg.hξ; n = length(ξ); z = zeros(n)
    m = [materials(model, seg.zone, x) for x in ξ]
    g(f) = getfield.(m, f)
    ρ, A, C, F, L, N = g(:ρ), g(:A), g(:C), g(:F), g(:L), g(:N)
    λinv, ρinv = _λinv.(m), _ρinv.(m)
    if kind === :SH          # (ρr², Lr², Lr, Nₗ, ω)
        sets = (T=Any[h^2 .* ρ .* ξ .^ 2, z, z, z, 1.0],
                H0=Any[z, ξ .^ 2 .* L, h .* ξ .* L, h^2 .* (L .- 2N), 0.0],
                Hk=nothing,
                HΛ=Any[z, z, z, h^2 .* N, 0.0])
        return sets, [3]
    elseif kind === :PSV     # (ρr², Cr², Fr, Lr², Lr, kFr, kLr, mUU, mVV, mUV, ω)
        sets = (T=Any[h^2 .* ρ .* ξ .^ 2, z, z, z, z, z, z, z, z, z, 1.0],
                H0=Any[z, ξ .^ 2 .* C, h .* ξ .* F, ξ .^ 2 .* L, h .* ξ .* L, z, z,
                       h^2 .* (-4A .+ 4N), h^2 .* (-L .+ 2N), z, 0.0],
                Hk=Any[z, z, z, z, z, h .* ξ .* F, h .* ξ .* L, z, z,
                       h^2 .* (2A .+ L .- 2N), 0.0],
                HΛ=Any[z, z, z, z, z, z, z, h^2 .* (-L), h^2 .* (-A), z, 0.0])
        return sets, [3, 5, 6, 7]
    elseif kind === :fluid   # (Pr², κr², Mₗ, ω)
        sets = (T=Any[z, h^2 .* ξ .^ 2 .* λinv, z, 1.0],
                H0=Any[ξ .^ 2 .* ρinv, z, z, 0.0],
                Hk=nothing,
                HΛ=Any[z, z, h^2 .* ρinv, 0.0])
        return sets, Int[]
    end
    error("kind must be :SH, :PSV or :fluid")
end

# fluid coefficients 1/λ, 1/ρ; a material function may override them
# (attenuation bases, lib/sources1D.jl); elasticMaterials does not
_λinv(m) = hasproperty(m, :λinv) ? m.λinv : 1 / m.λ
_ρinv(m) = hasproperty(m, :ρinv) ? m.ρinv : 1 / m.ρ

_mirror(m, firstIdx) = Any[(i in firstIdx ? -1 : 1) .* (x isa AbstractArray ? reverse(x) : x)
                            for (i, x) in enumerate(m)]
_reverseBlocks(v, nf, n) = vcat([reverse(v[(f-1)*n+1:f*n]) for f in 1:nf]...)
_equationOf(kind) = kind === :SH ? SH_EQUATION : kind === :PSV ? PSV_EQUATION : FLUID_EQUATION
_nfields(kind) = kind === :PSV ? 2 : 1

"""
    nearBoundaryConfig(base; j, P=base.pointsInSpace, degree=P-1)

Recipe for the row j nodes inside a segment end (j = 0: the end row), with
the interior scheme's reference size P and Taylor exactness (degree P-1, 4
for OPT5) instead of the OPT3 ones:
  j = 0   the clipped half-element recipe (`shClippedBoundaryConfig`, p = P:
          test centre ν = P, reference nodes 1..P+1, test integral clipped
          at the end);
  j ≥ 1   test centre ν = P - j, reference nodes 1..P (all inside the
          segment, one-sided), full test support.
Written for the top end; the bottom end uses the mirrored segment.
"""
function nearBoundaryConfig(base; j::Integer, P::Integer=base.pointsInSpace, degree::Integer=P - 1)
    base.orderBspace == -1 && error("near-boundary recipes need a weak-form base (OPT3/4/5)")
    # end row: OPT3 keeps its p = 2 half element (p = 3, degree 2 converges worse,
    # p7 2026-10-28); P ≥ 4 uses the scheme's own size and degree
    j == 0 && return P <= 3 ? shClippedBoundaryConfig(base; p=2) :
                              shClippedBoundaryConfig(base; p=P, degree=degree)
    0 < j < P - 1 || error("j = $(j) must be in 1:$(P - 2)")
    ν = P - j
    return (name="$(base.name)_nb$(P)_j$(j)_d$(degree)", orderBtime=base.orderBtime,
            orderBspace=base.orderBspace, pointsInSpace=P,
            pointsInTime=base.pointsInTime, supplementaryOrder=base.supplementaryOrder,
            fieldItpl=(ptsSpace=P, ptsTime=1, offsetSpace=0.0, offsetTime=1,
                       YorderBspace=1, YorderBtime=base.fieldItpl.YorderBtime),
            materItpl=merge(base.materItpl, (offsetSpace=Float64(ν - 1),)),
            extra=Dict{String,Any}(
                "nuCentre" => (ν,),
                "trialFunctionRefPoints" => (Tuple(1:P),),
                "testIntegrationBounds" => ((1.0, Float64(P)),),
                "exactTaylorTotalDegree" => degree))
end

"""
    flexOPTZoneBasis(kind, config, recipes, model, seg; configs, nearRows=:fallback,
                     clipP=2, fallbackName="OPT3")

Basis matrices (Dict :T, :H0, :Hk, :HΛ; field-major, K_ξ scale) of one zone:
the scheme's recipe in the interior, OPT3 rows where its stencil cannot
reach (`nearRows = :fallback`) or moving ν (`:movingNu`), native clipped
half-element rows (p = clipP, from OPT3 parameters) at both ends of every
field, mirrored recipe for the bottom; or (`:clippedFamily`) the scheme's own
near-boundary recipes for every row its stencil cannot serve, end row included
(`nearBoundaryConfig`).  Rows are divided by the recipe
scale (σh; σ = 2 for strong-form FD, 1 otherwise).
"""
function flexOPTZoneBasis(kind::Symbol, config, recipes, model, seg; configs,
                          nearRows=:fallback, clipP=2, fallbackName="OPT3",
                          materials=elasticMaterials)
    nearRows in (:fallback, :movingNu, :clippedFamily) ||
        error("nearRows must be :fallback, :movingNu or :clippedFamily")
    eq = _equationOf(kind); nf = _nfields(kind); n = length(seg.ξ)
    sets, firstIdx = zoneModelSets(kind, model, seg; materials=materials)
    base = getRecipe!(recipes, config, eq)
    rec = nearRows === :movingNu ?
        elasticRecipe(withRecipeExtras(config, Dict("nuGeometryMode" => :all)), eq) : base
    fbConfig = only(filter(c -> c.name == fallbackName, configs))
    fbRec = getRecipe!(recipes, fbConfig, eq)
    σ = config.orderBspace == -1 ? 2.0 : 1.0
    offs = recipeOffsets(base)
    Jtop = max(maximum(offs), 1); Jbot = max(-minimum(offs), 1)   # rows to replace at each end
    bad = union(collect(1:(-minimum(offs))), collect((n - maximum(offs) + 1):n))
    inner = nearRows === :fallback ? setdiff(bad, [1, n]) : Int[]
    # :clippedFamily -- row j (0 = end row) from the end: the scheme's own
    # near-boundary recipe (nearBoundaryConfig), bottom rows from the mirrored segment
    famRecs = nearRows === :clippedFamily ?
        [getRecipe!(recipes, nearBoundaryConfig(config; j=j), eq) for j in 0:max(Jtop, Jbot)-1] :
        [getRecipe!(recipes, shClippedBoundaryConfig(fbConfig; p=clipP), eq)]
    nTop = nearRows === :clippedFamily ? Jtop : 1
    nBot = nearRows === :clippedFamily ? Jbot : 1
    n >= 2 * max(nTop, nBot) + 1 || error("segment too short ($(n) nodes) for the near-boundary rows")
    rowsOf(nodes) = [(f - 1) * n + i for f in 1:nf for i in nodes]
    scale = fill(σ * seg.hξ, nf * n)
    scale[rowsOf(inner)] .= seg.hξ
    scale[rowsOf(vcat(collect(1:nBot), collect(n-nTop+1:n)))] .= seg.hξ
    out = Dict{Symbol,Matrix{Float64}}()
    for key in (:T, :H0, :Hk, :HΛ)
        m = getfield(sets, key)
        if isnothing(m)
            out[key] = zeros(nf * n, nf * n); continue
        end
        A = flexOPTOperator(rec, m, nf)
        if !isempty(inner)
            Fb = flexOPTOperator(fbRec, m, nf)
            r = rowsOf(inner); A[r, :] .= Fb[r, :]
        end
        for (jj, fr) in enumerate(famRecs)
            j = jj - 1
            top = flexOPTOperator(fr, m, nf)
            bot = flexOPTOperator(fr, _mirror(m, firstIdx), nf)
            for f in 1:nf
                j < nTop && (A[(f-1)*n+n-j, :] .= top[(f-1)*n+n-j, :])
                j < nBot && (A[(f-1)*n+1+j, :] .= _reverseBlocks(bot[(f-1)*n+n-j, :], nf, n))
            end
        end
        out[key] = A ./ scale
    end
    return out
end

# ---------------------------------------------------------------------------
# conventional DSM zone bases (linear splines, exact Gauss integration)
# ---------------------------------------------------------------------------

"""
    dsmZoneBasis(kind, model, seg; nq=8)

Conventional DSM basis matrices (same keys and ordering as `flexOPTZoneBasis`):
element-by-element Galerkin of the natural weak forms above (natural
half-element rows at both ends).
"""
function dsmZoneBasis(kind::Symbol, model, seg; nq=8, materials=elasticMaterials)
    nf = _nfields(kind); ξn = seg.ξ; n = length(ξn)
    B = Dict(key => zeros(nf * n, nf * n) for key in (:T, :H0, :Hk, :HΛ))
    q, w = gaussLegendre(nq)
    iU(i) = i; iV(i) = n + i
    for e in 1:n-1
        a, b = ξn[e], ξn[e+1]; he = b - a
        for (qk, wk) in zip(q, w)
            x = (a + b) / 2 + he / 2 * qk; wt = wk * he / 2
            M = materials(model, seg.zone, x)
            X = ((b - x) / he, (x - a) / he); dX = (-1 / he, 1 / he)
            for i in 1:2, j in 1:2
                I, J = e + i - 1, e + j - 1
                Xi, Xj, dXi, dXj = X[i], X[j], dX[i], dX[j]
                if kind === :SH
                    B[:T][I, J] += wt * M.ρ * x^2 * Xi * Xj
                    B[:H0][I, J] += wt * (-dXi * (x^2 * M.L * dXj - x * M.L * Xj) +
                                          Xi * (x * M.L * dXj + (-M.L + 2M.N) * Xj))
                    B[:HΛ][I, J] += wt * (-M.N * Xi * Xj)
                elseif kind === :fluid
                    B[:T][I, J] += wt * x^2 * _λinv(M) * Xi * Xj
                    B[:H0][I, J] += wt * (-(x^2 * _ρinv(M)) * dXi * dXj)
                    B[:HΛ][I, J] += wt * (-_ρinv(M) * Xi * Xj)
                else
                    A_, C_, F_, L_, N_ = M.A, M.C, M.F, M.L, M.N
                    B[:T][iU(I), iU(J)] += wt * M.ρ * x^2 * Xi * Xj
                    B[:T][iV(I), iV(J)] += wt * M.ρ * x^2 * Xi * Xj
                    # U-test, U-trial
                    B[:H0][iU(I), iU(J)] += wt * (-dXi * (x^2 * C_ * dXj + 2x * F_ * Xj) +
                                                  Xi * (-2x * F_ * dXj + (-4A_ + 4N_) * Xj))
                    B[:HΛ][iU(I), iU(J)] += wt * (-L_ * Xi * Xj)
                    # U-test, V-trial (× k)
                    B[:Hk][iU(I), iV(J)] += wt * (dXi * x * F_ * Xj +
                                                  Xi * (-x * L_ * dXj + (2A_ + L_ - 2N_) * Xj))
                    # V-test, V-trial
                    B[:H0][iV(I), iV(J)] += wt * (-dXi * (x^2 * L_ * dXj - x * L_ * Xj) +
                                                  Xi * (x * L_ * dXj + (-L_ + 2N_) * Xj))
                    B[:HΛ][iV(I), iV(J)] += wt * (-A_ * Xi * Xj)
                    # V-test, U-trial (× k)
                    B[:Hk][iV(I), iU(J)] += wt * (-dXi * x * L_ * Xj +
                                                  Xi * (x * F_ * dXj + (2A_ + L_ - 2N_) * Xj))
                end
            end
        end
    end
    return B
end

# ---------------------------------------------------------------------------
# global systems
# ---------------------------------------------------------------------------

struct ElasticGlobalSystem
    kind::Symbol                       # :SH or :PSV
    T::SparseMatrixCSC{Float64,Int}
    H0::SparseMatrixCSC{Float64,Int}
    Hk::SparseMatrixCSC{Float64,Int}
    HΛ::SparseMatrixCSC{Float64,Int}
    Cpl::SparseMatrixCSC{Float64,Int}
    dofmap::Vector{Matrix{Int}}        # per zone: dofmap[z][field, node] (1 = U or W or q, 2 = V)
    segments::Vector{Any}
    R::Float64
    label::String
end

"K(Ω, l) = Ω²T + H0 + k Hk + Λ HΛ + Ω Cpl (sparse, complex Ω allowed)."
function elasticOperator(sys::ElasticGlobalSystem, Ω, l)
    Λ = l * (l + 1); k = sqrt(Λ)
    return Ω^2 .* sys.T .+ sys.H0 .+ k .* sys.Hk .+ Λ .* sys.HΛ .+ Ω .* sys.Cpl
end

"""
    elasticGlobalSystem(kind, model, grid; scheme="DSMconv", configs=nothing,
                        recipes=Dict(), nearRows=:fallback, label=scheme)

Assemble the zones of `grid` (kind = :SH: outer solid region; :PSV: all
zones, solid (U, V) and fluid (q)).  Same-type neighbours share the interface
node (half-element rows added: continuity of the fields and natural
continuity of the fluxes); solid–fluid interfaces carry separate unknowns
and the G&O coupling terms.
"""
function elasticGlobalSystem(kind::Symbol, model, grid; scheme="DSMconv", configs=nothing,
                             recipes=Dict{Any,Any}(), nearRows=:fallback, label=scheme,
                             materials=elasticMaterials)
    kind in (:SH, :PSV) || error("kind must be :SH or :PSV")
    segs = grid.segments
    dofmap = Matrix{Int}[]
    nd = 0
    prevFluid = nothing; prevLast = Int[]
    zoneKinds = Symbol[]
    for seg in segs
        zk = seg.fluid ? :fluid : kind
        kind === :SH && seg.fluid && error("SH grids must contain solid zones only")
        nf = _nfields(zk); n = length(seg.ξ)
        D = zeros(Int, nf, n)
        for i in 1:n, f in 1:nf
            if i == 1 && prevFluid === seg.fluid && !isempty(prevLast)
                D[f, i] = prevLast[f]
            else
                nd += 1; D[f, i] = nd
            end
        end
        push!(dofmap, D); push!(zoneKinds, zk)
        prevFluid = seg.fluid; prevLast = D[:, n]
    end
    trip = Dict(key => (Int[], Int[], Float64[]) for key in (:T, :H0, :Hk, :HΛ, :Cpl))
    if scheme != "DSMconv"
        isnothing(configs) && error("flexOPT schemes need `configs`")
        config = only(filter(c -> c.name == scheme, configs))
    end
    for (iz, seg) in enumerate(segs)
        zk = zoneKinds[iz]
        Bz = scheme == "DSMconv" ? dsmZoneBasis(zk, model, seg; materials=materials) :
             flexOPTZoneBasis(zk, config, recipes, model, seg; configs=configs, nearRows=nearRows,
                              materials=materials)
        g = vec(permutedims(dofmap[iz]))        # local field-major index → global dof
        for key in (:T, :H0, :Hk, :HΛ)
            I, J, V = findnz(sparse(Bz[key]))
            append!(trip[key][1], g[I]); append!(trip[key][2], g[J]); append!(trip[key][3], V)
        end
    end
    # solid–fluid coupling (G&O 1994): one interface per change of type
    for iz in 2:length(segs)
        lo, up = segs[iz-1], segs[iz]
        lo.fluid == up.fluid && continue
        kind === :PSV || continue
        ξa = up.ξ[1]
        if lo.fluid && !up.fluid          # CMB-like: fluid below, solid above
            q = dofmap[iz-1][1, end]; u = dofmap[iz][1, 1]; s = -ξa^2
        else                              # ICB-like: solid below, fluid above
            u = dofmap[iz-1][1, end]; q = dofmap[iz][1, 1]; s = +ξa^2
        end
        push!(trip[:Cpl][1], u, q); push!(trip[:Cpl][2], q, u); push!(trip[:Cpl][3], s, s)
    end
    mk(key) = sparse(trip[key][1], trip[key][2], trip[key][3], nd, nd)   # duplicates are summed
    return ElasticGlobalSystem(kind, mk(:T), mk(:H0), mk(:Hk), mk(:HΛ), mk(:Cpl),
                               dofmap, collect(segs), grid.R, label)
end

"Load vectors -(1/R)·X(ξs) on each solid field at the source radius (one column per field)."
function sourceLoads(sys::ElasticGlobalSystem, rs_km)
    ξs = rs_km / sys.R
    nd = size(sys.T, 1); nf = sys.kind === :PSV ? 2 : 1
    B = zeros(nd, nf)
    for (iz, seg) in enumerate(sys.segments)
        seg.fluid && continue
        (seg.ξ[1] <= ξs <= seg.ξ[end]) || continue
        j = clamp(searchsortedlast(seg.ξ, ξs), 1, length(seg.ξ) - 1)
        h = seg.ξ[j+1] - seg.ξ[j]
        wts = ((seg.ξ[j+1] - ξs) / h, (ξs - seg.ξ[j]) / h)
        for f in 1:nf
            B[sys.dofmap[iz][f, j], f] = -wts[1] / sys.R
            B[sys.dofmap[iz][f, j+1], f] = -wts[2] / sys.R
        end
        return B
    end
    error("source radius $(rs_km) km is not inside a solid zone")
end

"Global dofs of the surface node (U, V) or (W,)."
surfaceDofs(sys::ElasticGlobalSystem) = Tuple(sys.dofmap[end][:, end])

# ---------------------------------------------------------------------------
# harmonics: P, B, C at (θ, φ) from a normalised Legendre table
# ---------------------------------------------------------------------------

"(Y, B_θ, B_φ, C_θ, C_φ) of Y_lm, B_lm = k⁻¹∇₁Y_lm, C_lm = -k⁻¹ r̂×∇₁Y_lm at (θ, φ)."
function vectorHarmonics(tab, l::Integer, m::Integer, φ::Real)
    am = abs(m)
    sgn = (m < 0 && isodd(am)) ? -1.0 : 1.0
    X = sgn * tab.X[l+1, am+1]; D = sgn * tab.D[l+1, am+1]
    e = cis(m * φ)
    Y = X * e
    l == 0 && return (Y, zero(ComplexF64), zero(ComplexF64), zero(ComplexF64), zero(ComplexF64))
    k = sqrt(l * (l + 1))
    dθY = D * e; ims = (im * m) * (X / tab.sinθ) * e
    return (Y, dθY / k, ims / k, ims / k, -dθY / k)
end

# ---------------------------------------------------------------------------
# 3-component synthetics
# ---------------------------------------------------------------------------

"""
    elasticSynthetics(psv, sh, R; rs_km, forces, θs=0.0, φs=0.0,
        receiversθdeg=10:10:180, φr=0.0, tlen=3276.8, fmax=0.05, omegaiTlen=1e-2,
        lmaxOf=ω -> ceil(Int, ω*R/3.0) + 30, lcap=1500, mmax=2, verbose=true)

Surface displacement (u_r = vertical, u_θ = radial, u_φ = transverse on
φr = 0) for point forces f = (f_r, f_θ, f_φ) at (rs, θs, φs), flat in ω.
`forces` is a Dict name => (f_r, f_θ, f_φ); one LU per (ω, l) serves all of
them.  `sh` may be `nothing` (P-SV only).  |m| ≤ mmax (exact for sources on
the polar axis).
"""
function elasticSynthetics(psv::ElasticGlobalSystem, sh, R; rs_km, forces,
                           θs=0.0, φs=0.0, receiversθdeg=10:10:180, φr=0.0,
                           tlen=3276.8, fmax=0.05, omegaiTlen=1e-2,
                           lmaxOf=ω -> ceil(Int, ω * R / 3.0) + 30, lcap=1500, mmax=2,
                           verbose=true)
    ωI = -log(omegaiTlen) / tlen
    kmax = floor(Int, fmax * tlen)
    ωs = [2π * k / tlen for k in 0:kmax]
    θr = deg2rad.(collect(receiversθdeg)); nr = length(θr)
    lmaxAll = min(lcap, maximum(lmaxOf(ω) for ω in ωs))
    recTabs = [normalisedLegendre(lmaxAll, mmax, θ) for θ in θr]
    srcTab = normalisedLegendre(lmaxAll, mmax, θs)
    names = collect(keys(forces))
    # receiver sums Σ_m c_lm × harmonic(receiver), per force, l, receiver
    S = Dict(nm => Dict(s => zeros(ComplexF64, lmaxAll + 1, nr)
                        for s in (:rP, :rB, :θP, :θB, :φP, :φB, :θC, :φC)) for nm in names)
    for nm in names
        f = forces[nm]
        for l in 0:lmaxAll, m in -min(l, mmax):min(l, mmax)
            Ys, Bθs, Bφs, Cθs, Cφs = vectorHarmonics(srcTab, l, m, φs)
            cP = f[1] * conj(Ys)
            cB = f[2] * conj(Bθs) + f[3] * conj(Bφs)
            cC = f[2] * conj(Cθs) + f[3] * conj(Cφs)
            (iszero(cP) && iszero(cB) && iszero(cC)) && continue
            for (ir, tab) in enumerate(recTabs)
                Y, Bθ, Bφ, Cθ, Cφ = vectorHarmonics(tab, l, m, φr)
                d = S[nm]
                d[:rP][l+1, ir] += cP * Y;  d[:rB][l+1, ir] += cB * Y
                d[:θP][l+1, ir] += cP * Bθ; d[:θB][l+1, ir] += cB * Bθ
                d[:φP][l+1, ir] += cP * Bφ; d[:φB][l+1, ir] += cB * Bφ
                d[:θC][l+1, ir] += cC * Cθ; d[:φC][l+1, ir] += cC * Cφ
            end
        end
    end
    bP = complex.(sourceLoads(psv, rs_km)); Us, Vs = surfaceDofs(psv)
    bS = isnothing(sh) ? nothing : complex.(sourceLoads(sh, rs_km))
    Ws = isnothing(sh) ? 0 : surfaceDofs(sh)[1]
    out = Dict(nm => (ur=zeros(ComplexF64, length(ωs), nr), uθ=zeros(ComplexF64, length(ωs), nr),
                      uφ=zeros(ComplexF64, length(ωs), nr)) for nm in names)
    t0 = time()
    for (iω, ω) in enumerate(ωs)
        Ωc = (ω - im * ωI) * R
        lmax = min(lmaxAll, lmaxOf(ω))
        for l in 0:lmax
            G = lu(elasticOperator(psv, Ωc, l)) \ bP          # columns: U-load, V-load
            gUU, gUV, gVU, gVV = G[Us, 1], G[Us, 2], G[Vs, 1], G[Vs, 2]
            gW = zero(ComplexF64)
            if !isnothing(sh) && l >= 1
                gW = (lu(elasticOperator(sh, Ωc, l)) \ bS)[Ws, 1]
            end
            for nm in names
                d = S[nm]; o = out[nm]
                @views o.ur[iω, :] .+= gUU .* d[:rP][l+1, :] .+ gUV .* d[:rB][l+1, :]
                @views o.uθ[iω, :] .+= gVU .* d[:θP][l+1, :] .+ gVV .* d[:θB][l+1, :] .+ gW .* d[:θC][l+1, :]
                @views o.uφ[iω, :] .+= gVU .* d[:φP][l+1, :] .+ gVV .* d[:φB][l+1, :] .+ gW .* d[:φC][l+1, :]
            end
        end
        verbose && (iω % 20 == 0 || iω == length(ωs)) &&
            @printf("  %s: f = %.4f Hz (%d/%d), lmax = %d, %.1f s\n", psv.label,
                    ω / 2π, iω, length(ωs), lmax, time() - t0)
    end
    return Dict(nm => (ω=ωs, ωI=ωI, tlen=tlen, θdeg=collect(receiversθdeg),
                       ur=out[nm].ur, uθ=out[nm].uθ, uφ=out[nm].uφ) for nm in names)
end

"Time series of the three components (Z = u_r, R = u_θ, T = u_φ); see shTimeSeries."
elasticTimeSeries(spec; np=4096, taper=0.2) =
    (Z=shTimeSeries(spec; np=np, taper=taper, component=:ur),
     R=shTimeSeries(spec; np=np, taper=taper, component=:uθ),
     T=shTimeSeries(spec; np=np, taper=taper, component=:uφ))

"""
    operatorDistance(Bflex, Bdsm, n, nf; margin=4)

Per basis key: least-squares scale s (Bflex ≈ s Bdsm; 1 expected on the K_ξ
scale) and relative difference on interior rows (`margin` nodes excluded at
each end of every field block).
"""
function operatorDistance(Bflex, Bdsm, n, nf; margin=4)
    rows = [(f - 1) * n + i for f in 1:nf for i in (1+margin):(n-margin)]
    out = Dict{Symbol,Any}()
    for key in (:T, :H0, :Hk, :HΛ)
        a = vec(Bflex[key][rows, :]); b = vec(Bdsm[key][rows, :])
        nb = norm(b)
        if nb == 0
            out[key] = (scale=NaN, relDiff=norm(a)); continue
        end
        s = dot(a, b) / dot(b, b)
        out[key] = (scale=s, relDiff=norm(a .- s .* b) / norm(s .* b))
    end
    return out
end
