# helpers/sphericalTI_SH_step2.jl
#
# Step 2 of the 1-D spherical TI SH workflow (notebooks/SphericalTI_SH_oct2026.ipynb):
# radial models and grids from src/planet1D, flexOPT operators for the
# famousEquation entry "1DsphericalSH_TI_Freq" (helpers/sphericalTI_SH.jl),
# and the conventional DSM SH matrices as reference.
#
# Conventions
# -----------
# * ξ = r / R (R = planet1D averaged planet radius, km), Ω = ω R (km/s units).
#   In ξ the radial operator keeps exactly the same form as in r with ω → Ω.
# * Zones: only the outermost continuous solid region (what SH sees).  Each
#   zone is a separate smooth segment; the interface radius appears once at
#   the top of the lower zone and once at the bottom of the upper zone
#   (duplicated nodes, as in planet1D.VerticalGridStructure) -- Step 4 will
#   overlap the operators there.
# * Grids: `mode=:uniform` (default) uses h_zone = width / ceil(width / h),
#   i.e. one uniform h per zone, generally different from zone to zone;
#   `mode=:dsm` uses planet1D.VerticalGridStructure (Geller & Takeuchi 1995,
#   eq. 6.2, controlled by `re` and ωmax).
# * Scaled local coordinate per zone: s = (ξ - ξ_bottom)/h_zone, so EVERY zone
#   and every h uses the SAME recipe built with Δ = 1.  The factors of h go
#   into the geometric material fields:
#        ∂s(ξ²L ∂sW) - ∂s(hξL W) + hξL ∂sW + (Ω² h²ρξ² - h²Nₗ) W  = h² × (ξ-form)
#   i.e. models (ρr², Lr², Lr, Nₗ, ω) = (h²ρξ², ξ²L, hξL, h²Nₗ, Ω), and the
#   flexOPT matrix in s is expected to be ≈ h × K_ξ, K_ξ = Ω²T - H.
#
# Requires helpers/sphericalTI_SH.jl to be included first.

using LinearAlgebra, Printf

# ---------------------------------------------------------------------------
# models (planet1D)
# ---------------------------------------------------------------------------

_coef4(v::Number) = [Float64(v), 0.0, 0.0, 0.0]
_coef4(v::AbstractVector) = (length(v) == 4 || error("need 4 polynomial coefficients"); Float64.(v))
_poly(c, x) = c[1] + c[2] * x + c[3] * x^2 + c[4] * x^3

"""
    layeredSHModel(; R=6371.0, radii=[3480.0, 6371.0], ρ=[3.3], Vsv=[4.5], Vsh=[4.5], Qμ=[600.0])

In-memory planet1D model (DSM-SH format) with `length(radii)-1` solid zones
between consecutive `radii` (km).  Each of ρ (g/cm³), Vsv, Vsh (km/s) per
zone is either a number or 4 polynomial coefficients in r/R (planet1D
convention).  Built with planet1D's own `"DSM-SH"` constructor, so it goes
through the same code path as PREM.
"""
function layeredSHModel(; R=6371.0, radii=[3480.0, 6371.0], ρ=[3.3], Vsv=[4.5],
                        Vsh=Vsv, Qμ=fill(600.0, length(radii) - 1))
    nz = length(radii) - 1
    all(length(v) == nz for v in (ρ, Vsv, Vsh, Qμ)) ||
        error("ρ, Vsv, Vsh, Qμ need one entry per zone ($nz)")
    arr = Float64[]
    for i in 1:nz
        append!(arr, [radii[i], radii[i+1]])
        append!(arr, _coef4(ρ[i])); append!(arr, _coef4(Vsv[i]))
        append!(arr, _coef4(Vsh[i])); push!(arr, Float64(Qμ[i]))
    end
    return planet1D.planet1DPSVmodel(nz, "DSM-SH", Float64(R), arr)
end

"Indices (bottom → top) of the outermost continuous solid zones."
function outerSolidZones(model)
    sf = model.solid_or_fluid
    top = findlast(==("S"), sf)
    isnothing(top) && error("model has no solid zone")
    bottom = top
    while bottom > 1 && sf[bottom-1] == "S"
        bottom -= 1
    end
    return bottom:top
end

# ---------------------------------------------------------------------------
# grids
# ---------------------------------------------------------------------------

"""
    shRadialGrid(model; mode=:uniform, h=50.0, re=1e-2, ωmax=nothing,
                 rmin=0.0, minIntervals=12)

Segments (one per outer-solid zone, bottom → top) with fields
`zone, r (km), ξ, h_km, hξ`.

* `mode=:uniform`: `h` in km — a number, or a per-zone vector/Dict indexed by
  the planet1D zone number; h_zone = width/ceil(width/h), at least
  `minIntervals` intervals per zone (thin zones get a smaller h).
* `mode=:dsm`: planet1D.VerticalGridStructure(model, re, ωmax) (ωmax in rad/s).
* `rmin` (km): DSM-style cut-off; nodes below rmin are dropped.
"""
function shRadialGrid(model; mode=:uniform, h=50.0, re=1e-2, ωmax=nothing,
                      rmin=0.0, minIntervals=12)
    R = model.averagedPlanetRadiusInKilometer
    zones = outerSolidZones(model)
    segments = NamedTuple[]
    if mode === :uniform
        for iz in zones
            rb = max(model.bottomRadius[iz], rmin)
            rt = model.topRadius[iz]
            rt > rb || continue
            hz = h isa AbstractDict ? h[iz] : h isa AbstractVector ? h[iz] : h
            n = max(ceil(Int, (rt - rb) / hz), minIntervals)
            r = collect(range(rb, rt; length=n + 1))
            push!(segments, _segment(iz, r, R))
        end
    elseif mode === :dsm
        isnothing(ωmax) && error("mode=:dsm needs ωmax (rad/s)")
        g = planet1D.VerticalGridStructure(model, Float64(re), Float64(ωmax))
        rall = g.RadiusInKm
        cuts = [0; findall(k -> rall[k+1] == rall[k], 1:length(rall)-1); length(rall)]
        pieces = [rall[cuts[i]+1:cuts[i+1]] for i in 1:length(cuts)-1]
        length(pieces) == model.nzone ||
            error("could not split the planet1D grid into $(model.nzone) zones " *
                  "(found $(length(pieces)))")
        for iz in zones
            r = filter(>=(rmin), pieces[iz])
            length(r) >= 2 || continue
            push!(segments, _segment(iz, r, R))
        end
    else
        error("mode must be :uniform or :dsm")
    end
    return (R=R, segments=segments, mode=mode)
end

function _segment(iz, r, R)
    hk = r[2] - r[1]
    maximum(abs.(diff(r) .- hk)) <= 1e-7 * max(1.0, abs(hk)) ||
        error("zone $iz: grid is not uniform")
    return (zone=iz, r=r, ξ=r ./ R, h_km=hk, hξ=hk / R)
end

# ---------------------------------------------------------------------------
# materials and flexOPT models on a segment
# ---------------------------------------------------------------------------

"ρ, L = ρVsv², N = ρVsh² from the zone's planet1D polynomials at ξ."
function shZoneMaterials(model, iz, ξ)
    ρ = _poly(model.C_ρ[iz, :], ξ)
    vsv = _poly(model.C_Vsv[iz, :], ξ)
    vsh = _poly(model.C_Vsh[iz, :], ξ)
    return (ρ=ρ, L=ρ * vsv^2, N=ρ * vsh^2)
end

"""
    shScaledModels(model, seg, l, Ω)

flexOPT `models` (order of `vars` of "1DsphericalSH_TI_Freq") in the scaled
local coordinate s of `seg`:  [h²ρξ², ξ²L, hξL, h²(L+(Λ-2)N), Ω].
"""
function shScaledModels(model, seg, l, Ω)
    h = seg.hξ
    Λ = l * (l + 1)
    m = [shZoneMaterials(model, seg.zone, ξ) for ξ in seg.ξ]
    ρ = getfield.(m, :ρ); L = getfield.(m, :L); N = getfield.(m, :N)
    ξ = seg.ξ
    return Any[h^2 .* ρ .* ξ .^ 2, L .* ξ .^ 2, h .* L .* ξ,
               h^2 .* (L .+ (Λ - 2) .* N), Float64(Ω)]
end

# ---------------------------------------------------------------------------
# conventional DSM reference (linear splines, exact Gauss integration)
# ---------------------------------------------------------------------------

function gaussLegendre(n::Integer)
    β = [k / sqrt(4k^2 - 1) for k in 1:n-1]
    E = eigen(SymTridiagonal(zeros(n), β))
    return E.values, 2 .* E.vectors[1, :] .^ 2
end

"""
    dsmSHConventional(model, seg, l, Ω; nq=8)

Element-by-element linear-spline Galerkin matrices on `seg` in ξ:
    T_ij = ∫ ρ ξ² X_i X_j dξ
    H_ij = ∫ L (ξX_i' - X_i)(ξX_j' - X_j) dξ + (Λ-2) ∫ N X_i X_j dξ
Gauss–Legendre with `nq` points per element (exact for the cubic planet1D
profiles).  Returns (K = Ω²T - H, T, H).
"""
function dsmSHConventional(model, seg, l, Ω; nq=8)
    ξn = seg.ξ
    n = length(ξn)
    T = zeros(n, n); H = zeros(n, n)
    Λ = l * (l + 1)
    q, w = gaussLegendre(nq)
    for e in 1:n-1
        a, b = ξn[e], ξn[e+1]
        he = b - a
        for (qk, wk) in zip(q, w)
            ξ = (a + b) / 2 + he / 2 * qk
            wt = wk * he / 2
            mat = shZoneMaterials(model, seg.zone, ξ)
            X = ((b - ξ) / he, (ξ - a) / he)
            dX = (-1 / he, 1 / he)
            for i in 1:2, j in 1:2
                I, J = e + i - 1, e + j - 1
                T[I, J] += wt * mat.ρ * ξ^2 * X[i] * X[j]
                H[I, J] += wt * (mat.L * (ξ * dX[i] - X[i]) * (ξ * dX[j] - X[j]) +
                                 (Λ - 2) * mat.N * X[i] * X[j])
            end
        end
    end
    return (K=Ω^2 .* T .- H, T=T, H=H)
end

# ---------------------------------------------------------------------------
# flexOPT recipes and matrices
# ---------------------------------------------------------------------------

# Part of the recipe cache key: bump it whenever the famousEquation entry in
# helpers/sphericalTI_SH.jl changes, so that stale cached recipes are not reused
# (myProduceOrLoad hashes the parameters, not the equation definition).
const SH_EQUATION_REVISION = "1DsphericalSH_TI_Freq-v1"  # (recipe_backend=:cpu also enters the cache key)

"""
    shRecipe(config; equation="1DsphericalSH_TI_Freq", form=nothing)

Semi-symbolic recipe at Δ = 1 (scaled local coordinate), cached with
myProduceOrLoad.  `config` is a NamedTuple as in seismo1Dbenchmark.ipynb.
Strong form is used automatically when orderBspace == -1 (FD).
"""
function shRecipe(config; equation="1DsphericalSH_TI_Freq", form=nothing)
    form = isnothing(form) ? (config.orderBspace == -1 ? :strong : :weak) : form
    p = Dict{String,Any}(
        "famousEquationType" => equation,
        "Δ" => 1.0,
        "orderBtime" => config.orderBtime,
        "orderBspace" => config.orderBspace,
        "pointsInSpace" => config.pointsInSpace,
        "pointsInTime" => config.pointsInTime,
        "supplementaryOrder" => config.supplementaryOrder,
        "fieldItpl" => config.fieldItpl,
        "materItpl" => config.materItpl,
        "variationalForm" => form,
        "equationRevision" => SH_EQUATION_REVISION,
        # Recipes MUST be built on the CPU: on a GPU backend (Metal here)
        # makeOPTrecette.jl computes the coefficients in Float32
        # (coefficientFloat, l. ~500), which leaves a ~1e-7 inexactness in the
        # stencils, amplified as ~eps/(kh)^2 in the residual (Step 2b-ii').
        "recipe_backend" => :cpu,
    )
    # optional extra recipe parameters (e.g. the clipped boundary recipes of
    # Step 3a': nuCentre, trialFunctionRefPoints, testIntegrationBounds, ...)
    hasproperty(config, :extra) && merge!(p, Dict{String,Any}(config.extra))
    return myProduceOrLoad(makeOPTsemiSymbolic, p, "semiSymbolic")
end

"""
    flexOPTSHMatrix(optRec, models; modelName)

Dense LHS matrix of the radial SH operator on one segment, through the
standard path numericalOperatorConstruction → prepareLinearSystem →
evaluateLinearSystem (as in seismo1Dbenchmark.ipynb).
"""
function flexOPTSHMatrix(optRec, models; modelName="sphericalSH")
    Nx = length(models[1])
    modelFam = (models=models, modelName=modelName, modelPoints=(Nx), Δ=(1.0))
    params = Dict{String,Any}("optRec" => optRec, "modelFam" => modelFam,
                              "absorbingBoundaries" => nothing,
                              "maskedRegionInSpace" => nothing,
                              # as in all the working notebooks: the default
                              # (:auto -> sparse) is not accepted by
                              # _residual_numerical_operator in src/
                              "representation" => "matrixfree")
    numOps = numericalOperatorConstruction(params)["numOperators"]
    prepared = prepareLinearSystem(numOps)
    knownField = copy(prepared.known_lhs_template)
    knownForce = copy(prepared.known_rhs_template)
    knownForce .= 0.0
    A, _ = evaluateLinearSystem(prepared, knownField, knownForce; sparse_output=false)
    return Matrix(A)
end

_interior(n, config) = (config.pointsInSpace):(n - config.pointsInSpace + 1)

"""
    compareWithDSM(A, K, rows)

Least-squares scale s with A ≈ s K on the interior `rows`, and the relative
matrix difference ‖A - sK‖/‖sK‖ there.
"""
function compareWithDSM(A, K, rows)
    a = vec(A[rows, :]); k = vec(K[rows, :])
    s = dot(a, k) / dot(k, k)
    return (scale=s, relDiff=norm(a .- s .* k) / norm(s .* k))
end

"""
    residualOfExact(A, M, W, rows)

Geller–Takeuchi-type accuracy measure: max|A W| / max|M W| on interior rows,
with W the exact solution sampled at the nodes and M the Ω²-part of the same
operator (M = A(Ω) - A(0)).  Scale-free.  O(h²) for conventional DSM and FD3;
OPT3 is also O(h²), but with a much smaller error constant than FD3.
"""
residualOfExact(A, M, W, rows) =
    maximum(abs.(A[rows, :] * W)) / maximum(abs.(M[rows, :] * W))

# ---------------------------------------------------------------------------
# exact solution: homogeneous isotropic solid, W = j_l(Ω ξ / β)
# ---------------------------------------------------------------------------

"Spherical Bessel j_l(x): upward recurrence for x > l, Miller downward otherwise."
function sphericalBesselJ(l::Integer, x::Real)
    x == 0 && return l == 0 ? 1.0 : 0.0
    j0 = sin(x) / x
    l == 0 && return j0
    if x > l
        jm, j = j0, sin(x) / x^2 - cos(x) / x
        for n in 1:l-1
            jm, j = j, (2n + 1) / x * j - jm
        end
        return j
    end
    nstart = l + 20 + ceil(Int, sqrt(40l))
    jp, j, jl = 0.0, 1e-30, 0.0
    for n in nstart:-1:1
        jm = (2n + 1) / x * j - jp
        jp, j = j, jm
        n - 1 == l && (jl = j)
        if abs(j) > 1e250
            jp /= 1e250; j /= 1e250; jl /= 1e250
        end
    end
    return jl * (j0 / j)
end

# ---------------------------------------------------------------------------
# drivers
# ---------------------------------------------------------------------------

"""
    shConvergenceHomogeneous(configs; l=10, ρ=3.3, β=4.5, R=6371.0, rcmb=3480.0,
                             wavelengths=30.0, nIntervals=[50,100,200,400],
                             grid=:uniform, reList=[1e-1,3e-2,1e-2,3e-3])

Homogeneous isotropic shell rcmb ≤ r ≤ R (one planet1D zone), Ω chosen so
that the shell holds `wavelengths`/(2π) radial wavelengths-ish
(Ω = wavelengths·β/(1-ξ_cmb)).  For each grid (uniform: `nIntervals`;
dsm: `reList`) and each config: scale s/h of A ≈ s K_DSM, relative matrix
difference, and the residual of the exact solution j_l for flexOPT and for
conventional DSM.  Returns a vector of NamedTuples (also printed).
"""
function shConvergenceHomogeneous(configs; l=10, ρ=3.3, β=4.5, R=6371.0,
                                  rcmb=3480.0, wavelengths=30.0,
                                  nIntervals=[50, 100, 200, 400],
                                  grid=:uniform, reList=[1e-1, 3e-2, 1e-2, 3e-3])
    model = layeredSHModel(; R=R, radii=[rcmb, R], ρ=[ρ], Vsv=[β], Vsh=[β])
    ξcmb = rcmb / R
    Ω = wavelengths * β / (1 - ξcmb)
    recipes = Dict(c.name => shRecipe(c) for c in configs)
    gridParams = grid === :uniform ? [(h=(R - rcmb) / n,) for n in nIntervals] :
                 [(re=re,) for re in reList]
    out = NamedTuple[]
    for gp in gridParams
        g = grid === :uniform ?
            shRadialGrid(model; mode=:uniform, h=gp.h) :
            shRadialGrid(model; mode=:dsm, re=gp.re, ωmax=Ω / R)
        seg = only(g.segments)
        n = length(seg.ξ)
        W = [sphericalBesselJ(l, Ω * ξ / β) for ξ in seg.ξ]
        dsm = dsmSHConventional(model, seg, l, Ω)
        dsm0 = dsmSHConventional(model, seg, l, 0.0)
        for c in configs
            rows = _interior(n, c)
            A = flexOPTSHMatrix(recipes[c.name], shScaledModels(model, seg, l, Ω))
            A0 = flexOPTSHMatrix(recipes[c.name], shScaledModels(model, seg, l, 0.0))
            cmp = compareWithDSM(A, dsm.K, rows)
            row = (config=c.name, grid=grid, gridParam=first(values(gp)),
                   n=n - 1, h_km=seg.h_km, scale_over_h=cmp.scale / seg.hξ,
                   relDiffDSM=cmp.relDiff,
                   resFlexOPT=residualOfExact(A, A .- A0, W, rows),
                   resDSM=residualOfExact(dsm.K, dsm.K .- dsm0.K, W, rows))
            push!(out, row)
            @printf("%-8s %-8s %10.4g  n=%5d  h=%8.3f km  s/h=% .4f  ‖A-sK‖/‖sK‖=%.3e  res(flexOPT)=%.3e  res(DSM)=%.3e\n",
                    row.config, string(grid), row.gridParam, row.n, row.h_km,
                    row.scale_over_h, row.relDiffDSM, row.resFlexOPT, row.resDSM)
        end
    end
    return out
end

"""
    shZoneComparison(model, configs; l=10, Ω=..., h=50.0, mode=:uniform, re=1e-2, rmin=0.0)

Zone-by-zone (outer solid region) comparison of flexOPT vs conventional DSM
matrices for a real planet1D model (e.g. PREM), at one (l, Ω).  No exact
solution here: reports s/h and the relative matrix difference per zone.
"""
function shZoneComparison(model, configs; l=10, Ω=200.0, h=50.0, mode=:uniform,
                          re=1e-2, rmin=0.0)
    R = model.averagedPlanetRadiusInKilometer
    recipes = Dict(c.name => shRecipe(c) for c in configs)
    g = shRadialGrid(model; mode=mode, h=h, re=re, ωmax=Ω / R, rmin=rmin)
    out = NamedTuple[]
    for seg in g.segments
        dsm = dsmSHConventional(model, seg, l, Ω)
        for c in configs
            n = length(seg.ξ)
            rows = _interior(n, c)
            isempty(rows) && continue
            A = flexOPTSHMatrix(recipes[c.name], shScaledModels(model, seg, l, Ω))
            cmp = compareWithDSM(A, dsm.K, rows)
            row = (zone=seg.zone, rbot=seg.r[1], rtop=seg.r[end], n=n - 1,
                   h_km=seg.h_km, config=c.name, scale_over_h=cmp.scale / seg.hξ,
                   relDiffDSM=cmp.relDiff)
            push!(out, row)
            @printf("zone %2d [%7.1f,%7.1f] km  n=%4d  h=%7.3f km  %-8s s/h=% .4f  ‖A-sK‖/‖sK‖=%.3e\n",
                    row.zone, row.rbot, row.rtop, row.n, row.h_km, row.config,
                    row.scale_over_h, row.relDiffDSM)
        end
    end
    return out
end

# ---------------------------------------------------------------------------
# recipe exactness diagnostic (constant unit material, polynomial fields)
# ---------------------------------------------------------------------------

"""
    shStencilExactness(configs; n=41)

Pure-recipe check in the scaled coordinate s (Δ = 1), with ONE unit material
field at a time (all others zero) on n nodes, applied to W = s^k:

* `:Lr²`  (operator ∂s(∂sW) = W''):  A·1 and A·s must vanish exactly;
  A·s², A·s³ must be proportional to W'' (fitted scale σ)
* `:Lr`   (operator -∂s(W) + ∂sW ≡ 0): A·s^k must vanish for every k

Returned errors are relative to max|A|·max|W| on interior rows.  Anything
much above ~1e-13 is an inexactness of the recipe itself; in the spherical
residual it is amplified by 1/(kh)², which would explain a floor that grows
as h decreases.
"""
function shStencilExactness(configs; n=41)
    s = collect(0.0:n-1)
    out = NamedTuple[]
    for c in configs
        rec = shRecipe(c)
        rows = _interior(n, c)
        unit(i) = Any[(j == i ? ones(n) : zeros(n)) for j in 1:4]
        A2 = flexOPTSHMatrix(rec, push!(unit(2), 0.0))   # Lr² = 1
        A1 = flexOPTSHMatrix(rec, push!(unit(3), 0.0))   # Lr  = 1
        rel(A, W, target) = maximum(abs.(A[rows, :] * W .- target[rows])) /
                            (maximum(abs.(A[rows, :])) * maximum(abs.(W)))
        # The Lr part must vanish identically, so it cannot be normalised by
        # itself: compare it with the size of ONE of its two terms (|∫X'W| ~ 1/2)
        # i.e. report max|A_Lr W| / (max|W| · 1/2) and max|A_Lr| directly.
        absLr(W) = maximum(abs.(A1[rows, :] * W)) / (0.5 * maximum(abs.(W)))
        W2 = s .^ 2; W3 = s .^ 3
        d2 = 2.0 .* ones(n); d3 = 6.0 .* s
        σ = dot(A2[rows, :] * W3, d3[rows]) / dot(d3[rows], d3[rows])
        row = (config=c.name, σ=σ,
               Lr²_const=rel(A2, ones(n), zeros(n)),
               Lr²_linear=rel(A2, s, zeros(n)),
               Lr²_quadratic=rel(A2, W2, σ .* d2),
               Lr²_cubic=rel(A2, W3, σ .* d3),
               Lr_maxabsA=maximum(abs.(A1[rows, :])),
               Lr_const=absLr(ones(n)),
               Lr_linear=absLr(s),
               Lr_cubic=absLr(W3))
        push!(out, row)
        @printf("%-8s σ=% .9f | Lr²: 1→%.1e s→%.1e s²→%.1e s³→%.1e | Lr: max|A|=%.1e  1→%.1e s→%.1e s³→%.1e\n",
                row.config, row.σ, row.Lr²_const, row.Lr²_linear, row.Lr²_quadratic,
                row.Lr²_cubic, row.Lr_maxabsA, row.Lr_const, row.Lr_linear, row.Lr_cubic)
        println("          recipe backend: ", rec["recette"].recipe_backend)
    end
    return out
end
