# dev/raytheory/lib/rayTauP.jl
#
# Fast path for laterally homogeneous models (VerticalModel): for a given ray
# parameter p the ray is known in closed form up to 1-D integrals,
#
#   Spherical1D (η = r/v, p in s/rad):
#       Δ = ∫ p / (r √(η² - p²)) dr ,  T = ∫ η² / (r √(η² - p²)) dr
#   Cartesian   (η = 1/v, p in s/km):
#       X = ∫ p / √(η² - p²) dz ,      T = ∫ η² / √(η² - p²) dz
#
# walked leg by leg with the SAME event logic as the tracer
# (boundaryDecision), so the phase vocabulary is shared. Every segment is
# integrated with Gauss-Legendre after ζ = a + (b-a)(1 - cos πu)/2, which
# removes the inverse-square-root singularity at turning points.
# Arrivals T(Δ) come from root finding on Δ(p) (all branches, triplications
# included), like TauP.

"Gauss-Legendre nodes/weights on [0, 1] (Golub-Welsch, stdlib only)"
function gaussLegendre01(n::Int)
    β = [i / sqrt(4i^2 - 1) for i in 1:n-1]
    E = eigen(SymTridiagonal(zeros(n), β))
    return (E.values .+ 1) ./ 2, E.vectors[1, :] .^ 2
end
const GL_NODES = Ref{Tuple{Vector{Float64},Vector{Float64}}}(gaussLegendre01(64))

etaOf(m::VerticalModel{Spherical1D}, k, ζ, mode) = ζ / zoneVelocity(m, k, ζ, mode)
etaOf(m::VerticalModel, k, ζ, mode) = 1 / zoneVelocity(m, k, ζ, mode)
weightOf(::VerticalModel{Spherical1D}, ζ) = 1 / ζ
weightOf(::VerticalModel, ζ) = 1.0

"(ΔX, T) of one segment [a, b] of zone k (η ≥ p on it)"
function segmentIntegrals(m::VerticalModel, k, a, b, p, mode)
    (b - a) <= 0 && return 0.0, 0.0
    u, w = GL_NODES[]
    D = 0.0; Tt = 0.0
    for i in eachindex(u)
        ζ = a + (b - a) * (1 - cos(π * u[i])) / 2
        dζ = (b - a) * (π / 2) * sin(π * u[i])
        η = etaOf(m, k, ζ, mode)
        s2 = η^2 - p^2
        s2 <= 0 && continue
        c = weightOf(m, ζ) * dζ * w[i] / sqrt(s2)
        D += p * c
        Tt += η^2 * c
    end
    return D, Tt
end

"first ζ between `from` and `to` (exclusive of `from`) where η(ζ) ≤ p, or nothing"
function findTurning(m::VerticalModel, k, from, to, p, mode; nscan=256)
    f(ζ) = etaOf(m, k, ζ, mode) - p
    prev = from
    for i in 1:nscan
        ζ = from + (to - from) * i / nscan
        if f(ζ) <= 0
            lo, hi = prev, ζ           # f(lo) > 0, f(hi) <= 0
            for _ in 1:80
                mid = (lo + hi) / 2
                f(mid) > 0 ? (lo = mid) : (hi = mid)
            end
            return (lo + hi) / 2
        end
        prev = ζ
    end
    return nothing
end

"""
    taupWalk(m, phase, srcDepth, p, startDir; totalReflection=true) -> (ok, Δ, T, legs)

Δ in radians (Spherical1D) or km (Cartesian). `legs` lists the (ζa, ζb,
zone, mode) segments (for plotting / debugging). With `totalReflection`, a
ray that cannot be transmitted through a discontinuity which is not one of
the phase's events (post-critical incidence) is totally reflected in the
same mode — this is the back-branch of the upper-mantle triplications, as
in TauP. Explicit events never fall back.
"""
function taupWalk(m::VerticalModel, phase::Phase, srcDepth, p, startDir::Symbol; maxSeg=2000, totalReflection=true)
    ζ = zetaOfDepth(m, srcDepth)
    dir = startDir
    k = zoneAt(m, ζ; side=dir === :down ? :below : :above)
    mode = phase.mode
    j = 1
    Δ = 0.0; T = 0.0
    legs = Tuple{Float64,Float64,Int,Symbol}[]
    etaOf(m, k, ζ, mode) > p || return (false, NaN, NaN, legs)
    for _ in 1:maxSeg
        if dir === :down
            bot = m.zones[k].bottom
            ζt = findTurning(m, k, ζ, bot, p, mode)
            if ζt !== nothing
                d, t = segmentIntegrals(m, k, ζt, ζ, p, mode)
                Δ += d; T += t; push!(legs, (ζ, ζt, k, mode))
                ζ = ζt; dir = :up
                continue
            end
            d, t = segmentIntegrals(m, k, bot, ζ, p, mode)
            Δ += d; T += t; push!(legs, (ζ, bot, k, mode))
            ζ = bot; b = k
        else
            top = m.zones[k].top
            ζt = findTurning(m, k, ζ, top, p, mode)
            if ζt !== nothing
                d, t = segmentIntegrals(m, k, ζ, ζt, p, mode)
                Δ += d; T += t; push!(legs, (ζ, ζt, k, mode))
                ζ = ζt; dir = :down
                continue
            end
            d, t = segmentIntegrals(m, k, ζ, top, p, mode)
            Δ += d; T += t; push!(legs, (ζ, top, k, mode))
            ζ = top; b = k + 1
        end
        jBefore = j
        action, newMode, j = boundaryDecision(m, phase, j, b, dir, mode)
        if action === :arrive
            return (true, Δ, T, legs)
        elseif action === :dead
            return (false, Δ, T, legs)
        end
        kOut = (action === :reflect) ? k : (dir === :down ? k - 1 : k + 1)
        (kOut < 1 || kOut > nzones(m)) && return (false, Δ, T, legs)
        if !(zoneVelocity(m, kOut, ζ, newMode) > 0 && etaOf(m, kOut, ζ, newMode) > p)
            # evanescent: an explicit event fails; a default transmission
            # (same mode, not an event) becomes a total reflection
            (action === :pass || (action === :transmit && totalReflection && newMode === mode && j == jBefore)) || return (false, Δ, T, legs)
            action = :reflect; kOut = k
        end
        mode = newMode; k = kOut
        action === :reflect && (dir = dir === :down ? :up : :down)
    end
    return (false, Δ, T, legs)
end

"""
Δ(p) on (0, pmax): `np` uniform samples, then midpoints are inserted where
two neighbours differ by more than `dmax` in distance (default 0.25° or
1 km) or where the ray changes between alive and dead (branch ends),
until the step in p is ~1e-9 pmax. Returns (ps, Δs, ok) sorted in p.
"""
function sampleDeltaOfP(m::VerticalModel, ph::Phase, srcDepth, dir, pmax; np=400, dmax=nothing, maxPass=40)
    dmax === nothing && (dmax = m.geom isa Spherical1D ? 0.25 : 1.0)
    f(p) = (r = taupWalk(m, ph, srcDepth, p, dir); r[1] ? toDistance(m, r[2]) : NaN)
    ps = collect(range(1e-6, 1 - 1e-9; length=np)) .* pmax
    Ds = f.(ps)
    for _ in 1:maxPass
        newp = Float64[]
        for i in 1:length(ps)-1
            (ps[i+1] - ps[i]) < 1e-9 * pmax && continue
            a, b = Ds[i], Ds[i+1]
            if isnan(a) != isnan(b) || (!isnan(a) && abs(a - b) > dmax)
                push!(newp, (ps[i] + ps[i+1]) / 2)
            end
        end
        isempty(newp) && break
        append!(ps, newp); append!(Ds, f.(newp))
        o = sortperm(ps); ps = ps[o]; Ds = Ds[o]
    end
    return ps, Ds, .!isnan.(Ds)
end

struct Arrival
    phase::String
    distance::Float64      # degrees (Spherical1D) or km (Cartesian)
    time::Float64          # s
    p::Float64             # s/rad (Spherical1D) or s/km
    startDir::Symbol
    takeoff::Float64       # degrees from the downward vertical
end
Base.show(io::IO, a::Arrival) = @printf(io, "%-8s Δ=%8.3f  T=%9.3f s  p=%10.5f  takeoff=%6.2f° (%s)",
    a.phase, a.distance, a.time, a.p, a.takeoff, a.startDir)

toDistance(::VerticalModel{Spherical1D}, Δ) = rad2deg(Δ)
toDistance(::VerticalModel, Δ) = Δ

"""
    taupArrivals(m, phase, srcDepth, dist; np=1500) -> Vector{Arrival}

All arrivals of `phase` at epicentral distance `dist` (degrees in
Spherical1D, km in Cartesian), by scanning p and refining every sign change
of Δ(p) - dist. Spherical: also Δ' = 360k ± dist (PKKP, major arcs).
"""
function taupArrivals(m::VerticalModel, phase, srcDepth, dist; np=400, dmax=nothing)
    ph = asPhase(phase)
    ζs = zetaOfDepth(m, srcDepth)
    dirs = ph.startDirection === :any ? (:down, :up) : (ph.startDirection,)
    out = Arrival[]
    for dir in dirs
        k = zoneAt(m, ζs; side=dir === :down ? :below : :above)
        pmax = etaOf(m, k, ζs, ph.mode)
        ps, Ds, ok = sampleDeltaOfP(m, ph, srcDepth, dir, pmax; np=np, dmax=dmax)
        targets = m.geom isa Spherical1D ? [dist, 360 - dist, 360 + dist, 720 - dist] : [dist]
        for D in unique(targets), i in 1:length(ps)-1
            (ok[i] && ok[i+1]) || continue
            f1, f2 = Ds[i] - D, Ds[i+1] - D
            f1 * f2 <= 0 || continue
            lo, hi, flo = ps[i], ps[i+1], f1
            for _ in 1:80
                mid = (lo + hi) / 2
                o, Δ, _, _ = taupWalk(m, ph, srcDepth, mid, dir)
                o || break
                fm = toDistance(m, Δ) - D
                (fm * flo > 0) ? (lo = mid; flo = fm) : (hi = mid)
            end
            pa = (lo + hi) / 2
            o, Δ, T, _ = taupWalk(m, ph, srcDepth, pa, dir)
            tol = m.geom isa Spherical1D ? 1e-6 : 1e-6 * max(1.0, abs(D))
            (o && abs(toDistance(m, Δ) - D) < tol) || continue          # jump, not a root
            v = zoneVelocity(m, k, ζs, ph.mode)
            s = m.geom isa Spherical1D ? pa * v / ζs : pa * v
            i0 = asin(clamp(s, -1, 1))
            push!(out, Arrival(ph.name, dist, T, pa, dir, rad2deg(dir === :down ? i0 : π - i0)))
        end
    end
    sort!(out; by=a -> a.time)
    # remove duplicates from neighbouring brackets
    keep = Arrival[]
    for a in out
        any(b -> abs(b.time - a.time) < 1e-6 && abs(b.p - a.p) < 1e-8, keep) || push!(keep, a)
    end
    return keep
end

"""
    travelTimeCurve(m, phase, srcDepth; np=600) -> Vector of branches (Δs, Ts, ps)

T(Δ) curves for plotting, triplications included: Δ(p), T(p) on the adaptive
p sampling, split wherever the ray dies or Δ jumps. Distances are folded
into [0, 180]° for Spherical1D.
"""
function travelTimeCurve(m::VerticalModel, phase, srcDepth; np=600, dmax=nothing)
    ph = asPhase(phase)
    ζs = zetaOfDepth(m, srcDepth)
    branchesOut = Tuple{Vector{Float64},Vector{Float64},Vector{Float64}}[]
    dirs = ph.startDirection === :any ? (:down, :up) : (ph.startDirection,)
    for dir in dirs
        k = zoneAt(m, ζs; side=dir === :down ? :below : :above)
        pmax = etaOf(m, k, ζs, ph.mode)
        ps, Ds, ok = sampleDeltaOfP(m, ph, srcDepth, dir, pmax; np=np, dmax=dmax)
        cur = (Float64[], Float64[], Float64[])
        jump = m.geom isa Spherical1D ? 2.0 : 10.0
        for i in eachindex(ps)
            if !ok[i] || (!isempty(cur[1]) && abs(Ds[i] - cur[1][end]) > jump)
                length(cur[1]) > 1 && push!(branchesOut, cur)
                cur = (Float64[], Float64[], Float64[])
                ok[i] || continue
            end
            _, Δ, T, _ = taupWalk(m, ph, srcDepth, ps[i], dir)
            D = toDistance(m, Δ)
            if m.geom isa Spherical1D
                D = mod(D, 360.0); D > 180 && (D = 360 - D)
            end
            push!(cur[1], D); push!(cur[2], T); push!(cur[3], ps[i])
        end
        length(cur[1]) > 1 && push!(branchesOut, cur)
    end
    return branchesOut
end
