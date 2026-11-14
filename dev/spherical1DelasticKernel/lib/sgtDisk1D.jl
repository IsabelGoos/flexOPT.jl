# dev/spherical1DelasticKernel/lib/sgtDisk1D.jl
#
# Strain Green's tensor "disks": for a point source on the north pole at
# radius rs, the displacement and strain fields on the meridian half-plane
# (r, θ), φ = 0, for every elementary source component
#   family = :moment  (Mrr, Mrθ, Mrφ, Mθθ, Mθφ, Mφφ)   → TSGT side (source S)
#   family = :force   (fr, fθ, fφ)                      → RSGT side (receiver R,
#                                                          by reciprocity)
# Only the non-zero (source, component) pairs are stored (`sgtPairs`).
# Any other azimuth, source position or mechanism is reconstructed by
# rotation (lib/sgtGeometry1D.jl) -- this replaces the hand-written
# h1..h10 / H1..H20 tables of Zhao & Chevrot (2011, App. B) and of the old
# DSM-Kernel SGTsh/SGTpsv codes.
#
# Uses (dev/spherical1Delastic/lib): ElasticGlobalSystem, elasticOperator,
# anelasticOperator, sourceRadialWeights, poleSourceCoefficients,
# elasticMaterials, fornbergWeights, _nearestStencil, _poly.
#
# Units: the g/cm³-km-s system of lib/sources1D.jl (M in (g/cm³)(km/s)²km³
# = 1e25 dyn·cm, f in (g/cm³)(km/s)²km²); displacement in km, strain
# dimensionless, r in km.

using LinearAlgebra, SparseArrays, Printf

# ---------------------------------------------------------------------------
# radial evaluation of a global solution vector at an arbitrary radius
# ---------------------------------------------------------------------------

"""
    radialEvaluator(sys, r_km; p=2, side=:upper)

Functionals that give, from a solution vector of `sys`, the radial values
and r-derivatives at r_km: Lagrange/Fornberg weights on the p+1 nodes of
the zone containing r_km nearest to it (the same stencils as
`sourceRadialWeights(...; method=:lagrange)`, so that source and receiver
functionals are identical → exact discrete reciprocity).  On a zone
boundary `side` picks the zone above (:upper) or below (:lower).  In a
fluid zone a 2nd derivative is added (p+2 nodes, needed for ∂r u_r).
Returns `nothing` if r_km is outside the region of `sys`.
"""
function radialEvaluator(sys, r_km; p::Integer=2, side::Symbol=:upper)
    ξr = r_km / sys.R
    cands = [iz for (iz, seg) in enumerate(sys.segments)
             if seg.ξ[1] - 1e-12 <= ξr <= seg.ξ[end] + 1e-12]
    isempty(cands) && return nothing
    iz = side === :upper ? last(cands) : first(cands)
    seg = sys.segments[iz]
    ξ = seg.ξ
    fluid = seg.fluid
    q = fluid ? max(p, 3) : p                      # fluid: needs q''
    idx = _nearestStencil(ξ, clamp(ξr, ξ[1], ξ[end]), q)
    w0 = fornbergWeights(ξr, ξ[idx], 0)
    w1 = fornbergWeights(ξr, ξ[idx], 1)
    w2 = fluid ? fornbergWeights(ξr, ξ[idx], 2) : zeros(length(idx))
    nf = size(sys.dofmap[iz], 1)
    dofs = [sys.dofmap[iz][f, i] for f in 1:nf, i in idx]     # nf × npts
    return (zone=iz, planetZone=seg.zone, fluid=fluid, r=r_km, ξ=ξr, R=sys.R,
            dofs=dofs, w0=w0, w1=w1, w2=w2)
end

"Material values (ρ, dρ/dr, λ) at an evaluator's radius (fluid relations)."
function _fluidMaterial(model, ev)
    m = elasticMaterials(model, ev.planetZone, ev.ξ)
    c = model.C_ρ[ev.planetZone, :]
    dρdξ = c[2] + 2c[3] * ev.ξ + 3c[4] * ev.ξ^2
    return (ρ=m.ρ, dρ=dρdξ / ev.R, λ=m.λ)
end

"""
    psvRadial(ev, G, l, ωc, fm) -> (U, dU, V, dV)  (each a row over the columns of G)

Radial functions of the P-SV solution(s) G (nd × ncol) at the evaluator.
Solid: U, V and their r-derivatives.  Fluid (unknown q = P/ω of Geller &
Ohminato 1994): u = -∇q/(ρω), i.e.
  U = -q_r/(ρω),  V = -k q/(ρωr),
  U' = -q_rr/(ρω) + ρ' q_r/(ρ²ω),  V' = -k (q_r/(ρωr) - q/(ρωr²) - ρ' q/(ρ²ωr)),
with ω = ωc (complex), ρ' = dρ/dr; div u = ω q/λ holds to discretisation
accuracy (checked in k1).
"""
function psvRadial(ev, G, l, ωc, fm)
    R = ev.R
    if !ev.fluid
        U = vec(ev.w0' * G[ev.dofs[1, :], :]); dU = vec(ev.w1' * G[ev.dofs[1, :], :]) ./ R
        V = vec(ev.w0' * G[ev.dofs[2, :], :]); dV = vec(ev.w1' * G[ev.dofs[2, :], :]) ./ R
        return (U, dU, V, dV)
    end
    qv = vec(ev.w0' * G[ev.dofs[1, :], :])
    qr = vec(ev.w1' * G[ev.dofs[1, :], :]) ./ R
    qrr = vec(ev.w2' * G[ev.dofs[1, :], :]) ./ R^2
    k = sqrt(l * (l + 1.0)); r = ev.r; ρ = fm.ρ; ρr = fm.dρ; w = ωc
    U = -qr ./ (ρ * w)
    dU = -qrr ./ (ρ * w) .+ ρr .* qr ./ (ρ^2 * w)
    V = -k .* qv ./ (ρ * w * r)
    dV = -k .* (qr ./ (ρ * w * r) .- qv ./ (ρ * w * r^2) .- ρr .* qv ./ (ρ^2 * w * r))
    return (U, dU, V, dV)
end

"Toroidal radial functions (W, dW) of SH solution(s) G at the evaluator (zeros if `ev === nothing`)."
function shRadial(ev, G, ncol)
    isnothing(ev) && return (zeros(ComplexF64, ncol), zeros(ComplexF64, ncol))
    W = vec(ev.w0' * G[ev.dofs[1, :], :]); dW = vec(ev.w1' * G[ev.dofs[1, :], :]) ./ ev.R
    return (W, dW)
end

# ---------------------------------------------------------------------------
# angular matrices of one (l, m): the strain is linear in
#   g = (∂rU, U/r, V/r, W/r, ∂rV, ∂rW)   and the displacement in (U, V, W),
# so ε_c(r, θ) = Σ_j g_j(r) A[c, j, θ],  u_c = Σ_j (U, V, W)_j B[c, j, θ]
# (A from modeFields at r = 1 with unit radial vectors: every undifferentiated
# U, V, W enters the strain divided by r exactly once).
# ---------------------------------------------------------------------------

const _GSLOT = (2, 1, 3, 5, 4, 6)   # position in rad = (U, dU, V, dV, W, dW) of g_j

function angularMatrices!(A, B, tabs, l, m)
    for (iθ, tab) in enumerate(tabs)
        for j in 1:6
            rad = ntuple(k -> k == _GSLOT[j] ? one(ComplexF64) : zero(ComplexF64), 6)
            _, e = modeFields(rad, 1.0, tab, l, m)
            for c in 1:6
                A[c, j, iθ] = e[c]
            end
        end
        for (j, k) in enumerate((1, 3, 5))
            rad = ntuple(q -> q == k ? one(ComplexF64) : zero(ComplexF64), 6)
            u, _ = modeFields(rad, 1.0, tab, l, m)
            for c in 1:3
                B[c, j, iθ] = u[c]
            end
        end
    end
    return A, B
end

# ---------------------------------------------------------------------------
# disk radii: the radial grid of the dev/spherical1Delastic system
# ---------------------------------------------------------------------------

"""
    diskRadii(sys; stride=1, rmin=0.0, rmax=Inf)

Radii at which a disk is stored: in every zone (segment) of `sys` a uniform
subset of its nodes, every `stride`-th node (spacing ≈ stride·h; when
stride does not divide the number of intervals the subset is a uniform grid
of the same zone, evaluated exactly by `radialEvaluator`).  Both ends of
every zone are kept, so each discontinuity appears twice (the lower zone
with side = :lower, the upper one with side = :upper) and the
interpolation never crosses it.  The centre r = 0 is skipped (1/r in the
formulas; the interpolants extrapolate to it).  Returns (r_km, side, seg).

Known limitation: the strain at a node uses one-sided Fornberg derivatives at
the ends of each zone, so it is only O(h) accurate right at the free surface
and at discontinuities (a traction-continuity closure would fix it; TODO).
"""
function diskRadii(sys; stride::Integer=1, rmin=0.0, rmax=Inf)
    r = Float64[]; side = Symbol[]; segi = Int[]
    for (iz, seg) in enumerate(sys.segments)
        a, b = seg.r[1], seg.r[end]
        (b < rmin || a > rmax) && continue
        n = length(seg.r) - 1
        ns = max(1, round(Int, n / stride))
        rr = n % ns == 0 ? seg.r[1:(n ÷ ns):end] : collect(range(a, b; length=ns + 1))
        for (k, x) in enumerate(rr)
            (x < rmin - 1e-9 || x > rmax + 1e-9) && continue
            x <= 1e-9 && continue          # r = 0: the formulas contain 1/r (fields are regular; extrapolated)
            push!(r, x); push!(side, k == length(rr) ? :lower : :upper); push!(segi, iz)
        end
    end
    return (r_km=r, side=side, seg=segi)
end

# ---------------------------------------------------------------------------
# the disk
# ---------------------------------------------------------------------------

"""
    computeSGTDisk(psv, sh, model; family, rs_km, radii=diskRadii(psv), θdeg, tlen, fmax,
        omegaiTlen=1e-4, fmin=0.0, ωlist=nothing, p=2, anelastic=nothing,
        lmaxOf=ω -> ceil(Int, ω*R/3.0) + 30, lcap=2000, verbose=true)

Frequency-domain disk (NamedTuple) of a pole source at radius `rs_km`
(`family = :moment` or `:force`) on `radii` (default: every node of the P-SV
grid, `diskRadii`; a plain vector of radii in km is also accepted, side
:upper) and colatitudes
`θdeg` (degrees, φ = 0).  Frequencies ω_k = 2πk/tlen for fmin ≤ f ≤ fmax
(or `ωlist`), complex frequency ω - iω_I with ω_I = -ln(omegaiTlen)/tlen
(lib/sources1D.jl convention).  One sparse LU per (ω, l) and system serves
all source components and all radii.  `anelastic = (psv=(q1,q2), sh=(q1,q2))`
switches on the complex moduli of classic DSM.

`adjoint` (default: true for :force): solve with Kᵀ instead of K.  The
receiver-side Green's function is the solution of the ADJOINT discrete
problem; for a symmetric operator (DSMconv) nothing changes, for a
non-symmetric one (OPT3, OPT5: Petrov–Galerkin rows, clipped / fallback
boundary rows) it makes reciprocity exact (u_n(R) = w_Rᵀ K⁻¹ w_S both ways)
and the kernels the exact discrete gradient δu = (K⁻ᵀ w_R)ᵀ δK (K⁻¹ w_S).
(PREM, OPT3, h = 10 km: reciprocity 3.6e-2 with K, vs 4e-12 for DSMconv.)

Fields: strain[iω, ipair, ir, iθ] for `pairs.strain`, disp[...] for
`pairs.disp`, lmax[iω], tail[iω, ir] = max over the last 10 degrees and
over θ of |strain contribution of degree l| / max_θ |strain| at that radius
(truncation diagnostic; at r = rs it stays O(1) because the strain is
singular at θ = 0 -- exclude the source neighbourhood there).
"""
function computeSGTDisk(psv, sh, model; family::Symbol, rs_km, radii=diskRadii(psv), θdeg, tlen, fmax,
                        omegaiTlen=1e-4, fmin=0.0, ωlist=nothing, p::Integer=2,
                        anelastic=nothing, lmaxOf=nothing, lcap=2000, adjoint=(family === :force),
                        verbose=true)
    if radii isa AbstractVector
        r_km = collect(Float64, radii); sides = fill(:upper, length(r_km))
    else
        r_km = collect(Float64, radii.r_km); sides = collect(radii.side)
    end
    R = psv.R
    lmaxOf = isnothing(lmaxOf) ? (ω -> ceil(Int, ω * R / 3.0) + 30) : lmaxOf
    ωI = -log(omegaiTlen) / tlen
    ωs = isnothing(ωlist) ? [2π * k / tlen for k in ceil(Int, fmin * tlen):floor(Int, fmax * tlen)] :
         collect(Float64, ωlist)
    θdeg = collect(Float64, θdeg)
    nr, nθ, nω = length(r_km), length(θdeg), length(ωs)
    pairs = sgtPairs(family)
    nsrc = family === :moment ? 6 : 3
    srcs = [unitSource(family, s) for s in 1:nsrc]
    # source loads (columns: U·X, U·X', V·X, V·X' | W·X, W·X')
    wP = sourceRadialWeights(psv, rs_km; method=:lagrange, p=p)
    ξs = wP.ξs
    LP = complex.(-(1 / R) .* hcat(wP.B0[:, 1], wP.B1[:, 1], wP.B0[:, 2], wP.B1[:, 2]))
    LS = nothing
    hasSH = !isnothing(sh)
    if hasSH
        try
            wS = sourceRadialWeights(sh, rs_km; method=:lagrange, p=p)
            LS = complex.(-(1 / R) .* hcat(wS.B0[:, 1], wS.B1[:, 1]))
        catch
            hasSH = false            # source outside the SH region: no toroidal field
        end
    end
    evP = [radialEvaluator(psv, r; p=p, side=sd) for (r, sd) in zip(r_km, sides)]
    any(isnothing, evP) && error("a disk radius is outside the P-SV system")
    # SH: same side; a radius on the bottom of the SH region (CMB) seen from
    # below (fluid side) has no toroidal field
    evS = [(hasSH && !evP[i].fluid) ? radialEvaluator(sh, r_km[i]; p=p, side=sides[i]) : nothing
           for i in eachindex(r_km)]
    # an SH evaluator in a zone that is fluid in P-SV cannot exist (SH = outer solid)
    fms = [ev.fluid ? _fluidMaterial(model, ev) : nothing for ev in evP]
    lmaxAll = min(lcap, maximum(lmaxOf(ω) for ω in ωs))
    tabs = [legendreTable2(lmaxAll, 2, deg2rad(θ)) for θ in θdeg]
    # pole coefficients per (l, m, source): sP (4), sC (2)
    coefs = Dict{Tuple{Int,Int,Int},Tuple{NTuple{4,ComplexF64},NTuple{2,ComplexF64}}}()
    for l in 0:lmaxAll, m in -min(l, 2):min(l, 2), s in 1:nsrc
        c = poleSourceCoefficients(l, m, srcs[s])
        sP = (c.P[1] + c.P[3] / (R * ξs), c.P[2] / R, c.B[1] + c.B[3] / (R * ξs), c.B[2] / R)
        sC = (c.C[1] + c.C[3] / (R * ξs), c.C[2] / R)
        (all(iszero, sP) && all(iszero, sC)) && continue
        coefs[(l, m, s)] = (ComplexF64.(sP), ComplexF64.(sC))
    end
    # pair lookup: for source s, list of (pair index, component)
    sPairs = [[(ip, c) for (ip, (ss, c)) in enumerate(pairs.strain) if ss == s] for s in 1:nsrc]
    dPairs = [[(ip, c) for (ip, (ss, c)) in enumerate(pairs.disp) if ss == s] for s in 1:nsrc]
    strain = zeros(ComplexF64, nω, length(pairs.strain), nr, nθ)
    disp = zeros(ComplexF64, nω, length(pairs.disp), nr, nθ)
    lmaxs = zeros(Int, nω); tails = zeros(nω, nr)
    qP = isnothing(anelastic) ? nothing : anelastic.psv
    qS = (isnothing(anelastic) || !hasSH) ? nothing : anelastic.sh
    radU = zeros(ComplexF64, 4, nr); raddU = similar(radU); radV = similar(radU); raddV = similar(radU)
    radW = zeros(ComplexF64, 2, nr); raddW = similar(radW)
    Gstr = zeros(ComplexF64, nr, 6); Gdis = zeros(ComplexF64, nr, 3)
    Aθ = zeros(ComplexF64, 6, 6, nθ); Bθ = zeros(ComplexF64, 3, 3, nθ)
    t0 = time()
    for (iω, ω) in enumerate(ωs)
        ωc = ω - im * ωI; Ωc = ωc * R
        lmax = min(lmaxAll, lmaxOf(ω)); lmaxs[iω] = lmax
        Sacc = zeros(ComplexF64, length(pairs.strain), nr, nθ)
        Lacc = similar(Sacc)                     # contribution of the current degree l
        lastIncr = zeros(nr)                      # max over the last 10 degrees of max_θ |Δε|
        Dacc = zeros(ComplexF64, length(pairs.disp), nr, nθ)
        for l in 0:lmax
            Kp = anelasticOperator(psv, qP, Ωc, ω, l)
            Gp = lu(adjoint ? sparse(transpose(Kp)) : Kp) \ LP
            Gs = if hasSH && l >= 1
                Ks = anelasticOperator(sh, qS, Ωc, ω, l)
                lu(adjoint ? sparse(transpose(Ks)) : Ks) \ LS
            else
                nothing
            end
            fill!(Lacc, 0)
            # radial factors per load column, then per source:
            #   strain:  g = (∂rU, U/r, V/r, W/r, ∂rV, ∂rW),  displacement: (U, V, W)
            for ir in 1:nr
                U, dU, V, dV = psvRadial(evP[ir], Gp, l, ωc, fms[ir])
                W, dW = isnothing(Gs) ? (zeros(ComplexF64, 2), zeros(ComplexF64, 2)) :
                        shRadial(evS[ir], Gs, 2)
                radU[:, ir] .= U; raddU[:, ir] .= dU; radV[:, ir] .= V; raddV[:, ir] .= dV
                radW[:, ir] .= W; raddW[:, ir] .= dW
            end
            for m in -min(l, 2):min(l, 2)
                any(s -> haskey(coefs, (l, m, s)), 1:nsrc) || continue
                angularMatrices!(Aθ, Bθ, tabs, l, m)          # Aθ[c, j, iθ], Bθ[c, j, iθ]
                for s in 1:nsrc
                    haskey(coefs, (l, m, s)) || continue
                    sP, sC = coefs[(l, m, s)]
                    for ir in 1:nr
                        r = r_km[ir]
                        U = transpose(view(radU, :, ir)) * collect(sP); dU = transpose(view(raddU, :, ir)) * collect(sP)
                        V = transpose(view(radV, :, ir)) * collect(sP); dV = transpose(view(raddV, :, ir)) * collect(sP)
                        W = transpose(view(radW, :, ir)) * collect(sC); dW = transpose(view(raddW, :, ir)) * collect(sC)
                        Gstr[ir, 1] = dU; Gstr[ir, 2] = U / r; Gstr[ir, 3] = V / r
                        Gstr[ir, 4] = W / r; Gstr[ir, 5] = dV; Gstr[ir, 6] = dW
                        Gdis[ir, 1] = U; Gdis[ir, 2] = V; Gdis[ir, 3] = W
                    end
                    for (ip, c) in sPairs[s]
                        mul!(view(Lacc, ip, :, :), Gstr, view(Aθ, c, :, :), true, true)
                    end
                    for (ip, c) in dPairs[s]
                        mul!(view(Dacc, ip, :, :), Gdis, view(Bθ, c, :, :), true, true)
                    end
                end
            end
            Sacc .+= Lacc
            if l > lmax - 10
                for ir in 1:nr
                    lastIncr[ir] = max(lastIncr[ir], maximum(abs, view(Lacc, :, ir, :)))
                end
            end
        end
        strain[iω, :, :, :] .= Sacc; disp[iω, :, :, :] .= Dacc
        for ir in 1:nr
            tot = maximum(abs, view(Sacc, :, ir, :))
            tails[iω, ir] = tot > 0 ? lastIncr[ir] / tot : 0.0
        end
        verbose && (iω % 10 == 0 || iω == nω) &&
            @printf("  disk %s rs=%.1f km: f = %.4f Hz (%d/%d), lmax = %d, median tail = %.1e, %.1f s\n",
                    family, rs_km, ω / 2π, iω, nω, lmax, sort(tails[iω, :])[(nr + 1) ÷ 2], time() - t0)
    end
    return (family=family, rs_km=Float64(rs_km), R=R, r_km=r_km, side=sides,
            seg=[ev.zone for ev in evP], θdeg=θdeg, ω=ωs, ωI=ωI,
            tlen=Float64(tlen), pairs=pairs, strain=strain, disp=disp, lmax=lmaxs, tail=tails,
            p=p, anelastic=!isnothing(anelastic), adjoint=adjoint, label=psv.label)
end

# ---------------------------------------------------------------------------
# evaluation of a disk on the pole frame
# ---------------------------------------------------------------------------

"Index of radius r_km in the disk (must be one of its radii)."
function diskRadiusIndex(disk, r_km; tol=1e-6)
    i = findfirst(x -> abs(x - r_km) <= tol * max(1.0, abs(r_km)), disk.r_km)
    isnothing(i) && error("radius $(r_km) km is not a disk radius")
    return i
end

"Lagrange weights (order `order`) for θ (degrees) on the disk's θ grid: (indices, weights)."
function diskθWeights(disk, θdeg; order::Integer=5)
    g = disk.θdeg; n = length(g)
    j = searchsortedlast(g, θdeg)
    if 1 <= j <= n && abs(g[j] - θdeg) <= 1e-12 * max(1.0, abs(θdeg))
        return (j:j, [1.0])
    end
    n >= order + 1 || error("θ grid too short for order $(order)")
    lo = clamp(j - (order - 1) ÷ 2, 1, n - order)
    idx = lo:lo+order
    return (idx, fornbergWeights(θdeg, g[idx], 0))
end

"""
    rotateSourceAboutPole(family, a, φ)

Components (pole-local basis r̂ = z, θ̂ = x, φ̂ = y at φ = 0) of the source
`a` rotated by -φ about the polar axis: the field of `a` at azimuth φ equals
(in the local basis) the field of the rotated source at φ = 0.
"""
function rotateSourceAboutPole(family::Symbol, a, φ)
    c, s = cos(φ), sin(φ)
    # rotation by -φ about z in the (r, θ, φ) ~ (z, x, y) ordering
    Q = [1.0 0 0; 0 c s; 0 -s c]
    if family === :moment
        return mtVector(Q * mtMatrix(a) * Q')
    else
        return Tuple(Q * collect(a))
    end
end

"""
    poleField(disk, iω, a, ir, θdeg, φ; which=:strain, order=5)

Field of the source `a` (6 MT or 3 force components, pole-local basis) at
(r = disk.r_km[ir], θ, φ) in the local basis (r̂, θ̂, φ̂) there:
strain as a 3×3 matrix, or displacement as a 3-vector.
"""
function poleField(disk, iω::Integer, a, ir::Integer, θdeg, φ; which::Symbol=:strain, order::Integer=5)
    b = rotateSourceAboutPole(disk.family, a, φ)
    idx, w = diskθWeights(disk, θdeg; order=order)
    if which === :strain
        e = zeros(ComplexF64, 6)
        for (ip, (s, c)) in enumerate(disk.pairs.strain)
            iszero(b[s]) && continue
            v = zero(ComplexF64)
            for (a_, i) in enumerate(idx)
                v += w[a_] * disk.strain[iω, ip, ir, i]
            end
            e[c] += b[s] * v
        end
        return strainMatrix(e)
    else
        u = zeros(ComplexF64, 3)
        for (ip, (s, c)) in enumerate(disk.pairs.disp)
            iszero(b[s]) && continue
            v = zero(ComplexF64)
            for (a_, i) in enumerate(idx)
                v += w[a_] * disk.disp[iω, ip, ir, i]
            end
            u[c] += b[s] * v
        end
        return u
    end
end

# ---------------------------------------------------------------------------
# time domain (same transform as shTimeSeries: one-sided spectrum, cosine
# taper over the last `taper` fraction of the band, exp(ω_I t) undamping)
# ---------------------------------------------------------------------------

function _taperWeights(nω, taper)
    w = ones(nω)
    if taper > 0
        k0 = floor(Int, (1 - taper) * (nω - 1))
        for k in k0:nω-1
            w[k+1] = 0.5 * (1 + cos(π * (k - k0) / max(1, nω - 1 - k0)))
        end
    end
    return w
end

"""
    diskSnapshot(disk, which, ip, t; taper=0.2)

Real time-domain field of the stored pair `ip` (`which = :strain | :disp`)
on the whole (r, θ) disk at time t (s), source time function = the flat
spectrum of the disk (Dirac moment / force), band-limited by `taper`.
The disk must hold the frequencies 2πk/tlen from k = 0 (fmin = 0).
"""
function diskSnapshot(disk, which::Symbol, ip::Integer, t::Real; taper=0.2)
    A = which === :strain ? disk.strain : disk.disp
    nω = length(disk.ω)
    abs(disk.ω[1]) < 1e-12 || error("diskSnapshot needs the zero frequency (fmin = 0)")
    w = _taperWeights(nω, taper)
    Δω = 2π / disk.tlen
    out = 0.5 .* real.(w[1] .* A[1, ip, :, :])
    for k in 2:nω
        out .+= real.(w[k] * cis(disk.ω[k] * t) .* A[k, ip, :, :])
    end
    return out .* (Δω / π * exp(disk.ωI * t))
end
