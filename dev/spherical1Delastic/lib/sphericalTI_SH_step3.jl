# helpers/sphericalTI_SH_step3.jl
#
# Step 3 of the 1-D spherical TI SH workflow (notebooks/SphericalTI_SH_oct2026.ipynb):
# boundaries.  Benchmark = homogeneous isotropic mantle shell a ≤ r ≤ b,
# traction-free at BOTH ends (SH sees the CMB as a free surface too).
#
#   exact:  W = A j_l(kr) + B y_l(kr),  traction ∝ W' - W/r = 0 at r = a, b
#           ⇔  D(k) = g_j(ka) g_y(kb) - g_j(kb) g_y(ka) = 0,
#              g_f(x) = x f_{l-1}(x) - (l+2) f_l(x),        Ω = β k  (ξ = r/R)
#   discrete: A(Ω) = A₀ + Ω² M  (ω enters only through ω²ρr²), so the
#           eigenfrequencies are the generalized eigenvalues of A₀ v = -Ω² M v.
#
# Boundary rows (each scheme × each closure):
#   * rows whose stencil would reach outside the segment (OPT4: 1 on one side,
#     2 on the other; OPT5: 2 + 2) cannot use the scheme's own recipe;
#     the non-end ones are taken from OPT3 (same s ≈ h·K scale) -- the
#     simplest form of "overlapping" two recipes;
#   * the END rows get one of two closures:
#       :dsmHalf        Galerkin natural condition: the conventional DSM row
#                       (exact half element, zero traction flux), rescaled to
#                       the scheme's convention;
#       (:traction, q)  traction collocation  Lr²·∂sW - Lr·W = 0  with
#                       one-sided Fornberg weights on q+1 nodes; the row has no
#                       Ω² part (M row = 0).
#
# Requires helpers/sphericalTI_SH.jl and helpers/sphericalTI_SH_step2.jl.

using LinearAlgebra, Printf

# ---------------------------------------------------------------------------
# exact toroidal eigenfrequencies of a homogeneous shell
# ---------------------------------------------------------------------------

"Spherical Bessel y_l(x) by upward recurrence (stable for y)."
function sphericalBesselY(l::Integer, x::Real)
    y0 = -cos(x) / x
    l == 0 && return y0
    y1 = -cos(x) / x^2 - sin(x) / x
    for n in 1:l-1
        y0, y1 = y1, (2n + 1) / x * y1 - y0
    end
    return y1
end

_gTraction(f, l, x) = x * f(l - 1, x) - (l + 2) * f(l, x)

"D(k) for a traction-free shell a ≤ ξ ≤ b (ξ = r/R)."
shellToroidalDeterminant(k, l, a, b) =
    _gTraction(sphericalBesselJ, l, k * a) * _gTraction(sphericalBesselY, l, k * b) -
    _gTraction(sphericalBesselJ, l, k * b) * _gTraction(sphericalBesselY, l, k * a)

"""
    shellToroidalRoots(l, a, b; kmax, dk=0.01, kmin=0.2)

Roots k (in 1/ξ units) of D(k) in (kmin, kmax): sign changes on a dk grid,
then bisection.  Ω_n = β k_n.  (l ≥ 2; l = 1 has the rigid-rotation mode.)
"""
function shellToroidalRoots(l, a, b; kmax, dk=0.01, kmin=0.2)
    f(k) = shellToroidalDeterminant(k, l, a, b)
    ks = collect(kmin:dk:kmax)
    vals = f.(ks)
    roots = Float64[]
    for i in 1:length(ks)-1
        sign(vals[i]) == sign(vals[i+1]) && continue
        lo, hi = ks[i], ks[i+1]
        flo = vals[i]
        for _ in 1:100
            mid = (lo + hi) / 2
            fm = f(mid)
            if sign(fm) == sign(flo)
                lo, flo = mid, fm
            else
                hi = mid
            end
        end
        push!(roots, (lo + hi) / 2)
    end
    return roots
end

# ---------------------------------------------------------------------------
# closures
# ---------------------------------------------------------------------------

"Fornberg (1988) weights for the m-th derivative at x0 on nodes x."
function fornbergWeights(x0::Real, x::AbstractVector, m::Integer=1)
    n = length(x)
    c = zeros(n, m + 1)
    c1 = 1.0; c4 = x[1] - x0; c[1, 1] = 1.0
    for i in 2:n
        mn = min(i - 1, m); c2 = 1.0; c5 = c4; c4 = x[i] - x0
        for j in 1:i-1
            c3 = x[i] - x[j]; c2 *= c3
            if j == i - 1
                for k in mn:-1:1
                    c[i, k+1] = c1 * (k * c[i-1, k] - c5 * c[i-1, k+1]) / c2
                end
                c[i, 1] = -c1 * c5 * c[i-1, 1] / c2
            end
            for k in mn:-1:1
                c[j, k+1] = (c4 * c[j, k+1] - k * c[j, k]) / c3
            end
            c[j, 1] = c4 * c[j, 1] / c3
        end
        c1 = c2
    end
    return c[:, m+1]
end

"Stencil offsets (relative to ν) of a recipe, from its nodes and centre."
function recipeOffsets(rec)
    r = rec["recette"]
    nodes = r.nodes[1]
    c = nodes[r.centresIndices[1]]
    return [round(Int, Tuple(n .- c)[1]) for n in nodes]
end

"""
    shTractionRow(models, iend, q)

Row (length n) of the traction collocation  Lr²·∂sW - Lr·W  at node `iend`
(1 or n) in the scaled coordinate s (unit spacing), with one-sided Fornberg
weights on q+1 nodes.  `models` as from `shScaledModels`.
"""
function shTractionRow(models, iend, q)
    Lr², Lr = models[2], models[3]
    n = length(Lr²)
    idx = iend == 1 ? (1:q+1) : (n-q:n)
    w = fornbergWeights(Float64(iend), Float64.(collect(idx)))
    row = zeros(n)
    row[idx] .+= Lr²[iend] .* w
    row[iend] -= Lr[iend]
    return row
end

# ---------------------------------------------------------------------------
# native flexOPT boundary rows: clipped (half-element) recipe   [Step 3a']
# ---------------------------------------------------------------------------

"""
    shClippedBoundaryConfig(base; p, degree=p-1)

1-D analogue of the `clipped_available` surface recipe of
scripts/earth2D_OPT3_freesurface.jl: p one-sided trial nodes 1..p, the
test-function centre on the last node (`nuCentre = p`, the free surface),
linear one-sided field interpolation, material taken at the surface node, and
the test function integrated only over the solid side
(`testIntegrationBounds = (1, p)`) -- i.e. a half element with the natural
(zero-traction) flux.  `base` is a weak-form volume config (OPT3).
p = 2 should reproduce the conventional DSM half row.
"""
function shClippedBoundaryConfig(base; p::Integer, degree::Integer=p - 1, refPoints=:extended)
    base.orderBspace == -1 && error("the clipped boundary recipe needs a weak-form base (e.g. OPT3)")
    return (name="$(base.name)_clip$(p)_$(refPoints)", orderBtime=base.orderBtime,
            orderBspace=base.orderBspace, pointsInSpace=p,
            pointsInTime=base.pointsInTime, supplementaryOrder=base.supplementaryOrder,
            fieldItpl=(ptsSpace=p, ptsTime=1, offsetSpace=0.0, offsetTime=1,
                       YorderBspace=1, YorderBtime=base.fieldItpl.YorderBtime),
            materItpl=merge(base.materItpl, (offsetSpace=Float64(p - 1),)),
            extra=Dict{String,Any}(
                "nuCentre" => (p,),
                # :extended (default, Step 3a''): reference nodes 1..p+1 so that the
                # test centre ν = p is interior to the reference family and only
                # testIntegrationBounds clips it (a true half hat).  :as2D reproduces
                # the 2D surface choice (1..max(p,3)), for which ν = p is the LAST
                # reference node when p ≥ 3 (truncated end family -- Step 3a′ showed
                # O(h) errors for p = 3, 5 with it).
                "trialFunctionRefPoints" => (refPoints === :as2D ? Tuple(1:max(p, 3)) :
                                             Tuple(1:p+1),),
                "testIntegrationBounds" => ((1.0, Float64(p)),),
                "exactTaylorTotalDegree" => degree))
end

"Models of `seg` mirrored (s → -s): arrays reversed and the Lr field negated."
function _mirroredModels(m)
    return Any[reverse(m[1]), reverse(m[2]), -reverse(m[3]), reverse(m[4]), m[5]]
end

"""
    shClippedBoundaryRows(clipRec, model, seg, l)

Top (free surface, node n) and bottom (CMB, node 1) rows (A₀ and Ω²-part)
from the clipped recipe.  The bottom row uses the same recipe on the mirrored
segment, mapped back.
"""
function shClippedBoundaryRows(clipRec, model, seg, l)
    m0 = shScaledModels(model, seg, l, 0.0)
    m1 = shScaledModels(model, seg, l, 1.0)
    n = length(m0[1])
    A0 = flexOPTSHMatrix(clipRec, m0); A1 = flexOPTSHMatrix(clipRec, m1)
    R0 = flexOPTSHMatrix(clipRec, _mirroredModels(m0))
    R1 = flexOPTSHMatrix(clipRec, _mirroredModels(m1))
    top = (A0=A0[n, :], M=A1[n, :] .- A0[n, :])
    bottom = (A0=reverse(R0[n, :]), M=reverse(R1[n, :] .- R0[n, :]))
    return (top=top, bottom=bottom)
end

"""
    compareClippedWithDSMHalf(clipRec, model, seg, l, Ω)

Relative difference (after the best scale) between the clipped boundary rows
at Ω and the conventional DSM half rows: ≈ 0 expected for p = 2.
"""
function compareClippedWithDSMHalf(clipRec, model, seg, l, Ω)
    rows = shClippedBoundaryRows(clipRec, model, seg, l)
    K = dsmSHConventional(model, seg, l, Ω).K
    n = size(K, 1)
    res = NamedTuple[]
    for (name, r, kr) in (("top", rows.top, K[n, :]), ("bottom", rows.bottom, K[1, :]))
        a = r.A0 .+ Ω^2 .* r.M
        s = dot(a, kr) / dot(kr, kr)
        push!(res, (side=name, scale_over_h=s / seg.hξ, relDiff=norm(a .- s .* kr) / norm(s .* kr),
                    support=findall(!iszero, a)))
    end
    return res
end

# ---------------------------------------------------------------------------
# assembly with closures, eigenfrequencies
# ---------------------------------------------------------------------------

"config with extra recipe parameters added (kept in the cache key)"
function withRecipeExtras(config, extras::AbstractDict)
    old = hasproperty(config, :extra) ? Dict{String,Any}(config.extra) : Dict{String,Any}()
    return merge(config, (extra=merge(old, Dict{String,Any}(extras)),))
end

"""
    shShellSystem(config, recipes, model, seg, l; closure=:dsmHalf, fallback="OPT3",
                  nearRows=:fallback, clipRefPoints=:extended)

(A₀, M) for one segment with both ends closed.  `recipes` is a Dict name →
recipe (must contain `fallback` when the scheme's stencil is wider than 3 and
nearRows = :fallback).
"""
function shShellSystem(config, recipes, model, seg, l; closure=:dsmHalf, fallback="OPT3",
                       clipBase="OPT3", clipConfig=nothing, nearRows=:fallback,
                       clipRefPoints=:extended)
    # nearRows = :fallback  -> rows the stencil cannot reach are taken from OPT3
    #          = :movingNu  -> recipe built with nuGeometryMode = :all, so that
    #            prepareNumericalOperatorGeometry moves ν inward near the ends and
    #            every row except the two end rows keeps the scheme's own order
    rec = nearRows === :movingNu ?
        shRecipe(withRecipeExtras(config, Dict("nuGeometryMode" => :all))) :
        recipes[config.name]
    m0 = shScaledModels(model, seg, l, 0.0)
    A0 = flexOPTSHMatrix(rec, m0)
    M = flexOPTSHMatrix(rec, shScaledModels(model, seg, l, 1.0)) .- A0
    n = size(A0, 1)

    # rows the scheme's stencil cannot serve (outside nodes dropped)
    offs = recipeOffsets(rec)
    lowBad = 1:(-minimum(offs))
    highBad = (n - maximum(offs) + 1):n
    bad = union(collect(lowBad), collect(highBad))
    inner = nearRows === :movingNu ? Int[] : setdiff(bad, [1, n])
    if !isempty(inner)
        fr = recipes[fallback]
        F0 = flexOPTSHMatrix(fr, m0)
        FM = flexOPTSHMatrix(fr, shScaledModels(model, seg, l, 1.0)) .- F0
        A0[inner, :] .= F0[inner, :]
        M[inner, :] .= FM[inner, :]
    end

    # end rows
    if closure === :dsmHalf
        K0 = dsmSHConventional(model, seg, l, 0.0).K
        KM = dsmSHConventional(model, seg, l, 1.0).K .- K0
        rows = _interior(n, config)
        a = vec(A0[rows, :] .+ M[rows, :]); k = vec(K0[rows, :] .+ KM[rows, :])
        s = dot(a, k) / dot(k, k)               # the scheme's scale (≈ h or 2h)
        for i in (1, n)
            A0[i, :] .= s .* K0[i, :]
            M[i, :] .= s .* KM[i, :]
        end
    elseif closure isa Tuple && first(closure) === :traction
        q = closure[2]
        for (i, nb) in ((1, 2), (n, n - 1))
            row = shTractionRow(m0, i, q)
            wgt = norm(A0[nb, :]) / norm(row)    # conditioning only
            A0[i, :] .= wgt .* row
            M[i, :] .= 0.0
        end
    elseif closure isa Tuple && first(closure) === :clipped
        # native flexOPT half-element rows; whole-row scaling does not change
        # the generalized eigenvalues, so no rescaling is needed
        haskey(recipes, clipBase) || error("closure (:clipped, p) needs the base recipe $clipBase")
        clipRec = shRecipe(shClippedBoundaryConfig(clipConfig; p=closure[2], refPoints=clipRefPoints))
        rows = shClippedBoundaryRows(clipRec, model, seg, l)
        A0[n, :] .= rows.top.A0;    M[n, :] .= rows.top.M
        A0[1, :] .= rows.bottom.A0; M[1, :] .= rows.bottom.M
    else
        error("closure must be :dsmHalf, (:traction, q) or (:clipped, p)")
    end
    return A0, M
end

"""
    shEigenfrequencies(A0, M)

Real positive Ω from A₀ v = -Ω² M v (generalized eigenvalues; constraint rows
with M = 0 give infinite eigenvalues and are dropped; complex pairs from
non-symmetric stencils are dropped and counted).
"""
function shEigenfrequencies(A0, M)
    λ = eigvals(A0, -M)
    λ = λ[isfinite.(λ)]
    good = [abs(imag(z)) <= 1e-8 * abs(real(z)) && real(z) > 0 for z in λ]
    Ω = sort(sqrt.(real.(λ[good])))
    return (Ω=Ω, nComplex=count(z -> abs(imag(z)) > 1e-8 * abs(real(z)), λ))
end

"""
    shShellEigenBenchmark(configs; l=10, ρ=3.3, β=4.5, R=6371.0, rcmb=3480.0,
        nIntervals=[50,100,200,400], closures=[:dsmHalf, (:traction,2), (:traction,4)],
        nModes=8)

Relative eigenfrequency errors (Ω_h - Ω)/Ω of the first `nModes` toroidal
overtones of the shell, for every scheme × closure × grid, plus conventional
DSM (natural BC) as reference.  Returns a vector of NamedTuples.
"""
function shShellEigenBenchmark(configs; l=10, ρ=3.3, β=4.5, R=6371.0, rcmb=3480.0,
                               nIntervals=[50, 100, 200, 400],
                               closures=[:dsmHalf, (:traction, 2), (:traction, 4)],
                               nModes=8, nearRows=:fallback, clipRefPoints=:extended)
    model = layeredSHModel(; R=R, radii=[rcmb, R], ρ=[ρ], Vsv=[β], Vsh=[β])
    a, b = rcmb / R, 1.0
    kexact = shellToroidalRoots(l, a, b; kmax=80.0)
    length(kexact) >= nModes || error("not enough analytic roots below kmax")
    Ωexact = β .* kexact[1:nModes]
    @printf("analytic Ω (ξ units, l=%d): %s\n", l, join((@sprintf("%.6f", w) for w in Ωexact), "  "))
    recipes = Dict(c.name => shRecipe(c) for c in configs)
    clipConfig = only(filter(c -> c.name == "OPT3", configs))
    out = NamedTuple[]
    report(name, cl, n, h, Ω, nC) = begin
        k = min(nModes, length(Ω))
        err = (Ω[1:k] .- Ωexact[1:k]) ./ Ωexact[1:k]
        push!(out, (scheme=name, closure=string(cl), n=n, h_km=h, relErr=err, nComplex=nC))
        @printf("%-8s %-14s n=%4d h=%7.3f km  nC=%d | %s\n", name, string(cl), n, h, nC,
                join((@sprintf("%+.2e", e) for e in err), " "))
        any(abs.(err) .> 0.1) && println("   ⚠ large error: spurious or missing mode?")
    end
    for n in nIntervals
        seg = only(shRadialGrid(model; mode=:uniform, h=(R - rcmb) / n).segments)
        K0 = dsmSHConventional(model, seg, l, 0.0).K
        KM = dsmSHConventional(model, seg, l, 1.0).K .- K0
        ev = shEigenfrequencies(K0, KM)
        report("DSMconv", :natural, n, seg.h_km, ev.Ω, ev.nComplex)
        for c in configs, cl in closures
            try
                A0, M = shShellSystem(c, recipes, model, seg, l; closure=cl,
                                      clipConfig=clipConfig, nearRows=nearRows,
                                      clipRefPoints=clipRefPoints)
                ev = shEigenfrequencies(A0, M)
                report(nearRows === :movingNu ? c.name * "+mν" : c.name, cl, n, seg.h_km, ev.Ω, ev.nComplex)
            catch err
                @printf("%-8s %-14s n=%4d  FAILED: %s\n", c.name, string(cl), n,
                        first(sprint(showerror, err), 300))
            end
        end
    end
    return out
end
