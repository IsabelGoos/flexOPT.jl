# dev/raytheory/lib/rayTracer.jl
#
# Characteristics of H integrated with RK4 in T, with exact (bisection)
# location of every boundary crossing, where the phase's event list decides
# what happens (layer 3: interfaceSlowness). The same code runs for
# Spherical1D, Cartesian2D and Cartesian3D; only rayRHS / the medium API differ.

"""
What to do when a ray travelling `dir` in `mode` reaches boundary b with the
phase's next event index j. Returns (action, newMode, newj) with action in
:transmit, :reflect, :pass (not a discontinuity), :arrive, :dead.
"""
function boundaryDecision(med, phase::Phase, j::Int, b::Int, dir::Symbol, mode::Symbol)
    isModelBottom(med, b) && return (:dead, mode, j)
    isDiscontinuous(med, b) || return (:pass, mode, j)
    if j <= length(phase.events)
        e = phase.events[j]
        if (e.direction === :any || e.direction === dir) && matchesBoundary(med, b, e.boundary)
            return (e.action, e.mode, j + 1)
        end
    end
    if isFreeSurface(med, b)
        return (j > length(phase.events) && dir === :up) ? (:arrive, mode, j) : (:dead, mode, j)
    end
    isFluidSolid(med, b) && return (:dead, mode, j)
    return (:transmit, mode, j)
end

"""
Largest safe RK4 step (s) in region k: a quarter of the time to cross the
zone at the local velocity, and 0.2 / |∇v| (the turning time scale), so that
thin, strongly graded zones (e.g. the 10-m regolith zones of the Mars models)
are never stepped over with an extrapolated polynomial.
"""
function stepLimit(med::VerticalModel, k::Int, x, mode)
    md = mode in (:P, :S) ? mode : :P
    v = abs(velocity(med, k, x, md)); v > 0 || (md = :P; v = abs(velocity(med, k, x, :P)))
    z = med.zones[k]
    g = abs(zoneVelocityDeriv(med, k, zetaOf(med, x), md))
    return min(0.25 * (z.top - z.bottom) / v, g > 0 ? 0.2 / g : Inf)
end
stepLimit(med, k, x, mode) = Inf

"One hit of a boundary along a ray"
struct RayHit
    boundary::Symbol
    action::Symbol
    modeIn::Symbol
    modeOut::Symbol
    x::Vector{Float64}
    T::Float64
end

struct RayResult
    phase::Phase
    ok::Bool
    reason::Symbol              # :arrived, or why it died
    T::Float64
    x0::Vector{Float64}
    p0::Vector{Float64}
    xEnd::Vector{Float64}
    pEnd::Vector{Float64}
    path::Vector{Vector{Float64}}
    pathT::Vector{Float64}
    pathMode::Vector{Symbol}
    hits::Vector{RayHit}
    maxResidual::Float64        # max |H| along the ray (accuracy check)
end

function rk4(H, geom, med, k, x, p, mode, h)
    k1x, k1p = rayRHS(H, geom, med, k, x, p, mode)
    k2x, k2p = rayRHS(H, geom, med, k, x .+ (h / 2) .* k1x, p .+ (h / 2) .* k1p, mode)
    k3x, k3p = rayRHS(H, geom, med, k, x .+ (h / 2) .* k2x, p .+ (h / 2) .* k2p, mode)
    k4x, k4p = rayRHS(H, geom, med, k, x .+ h .* k3x, p .+ h .* k3p, mode)
    return x .+ (h / 6) .* (k1x .+ 2k2x .+ 2k3x .+ k4x), p .+ (h / 6) .* (k1p .+ 2k2p .+ 2k3p .+ k4p)
end

"""
    traceRay(med, phase, x0, p0; H=IsotropicHamiltonian(), dT=0.5, Tmax=5000.0,
             keepPath=true, bbox=nothing)

Trace one ray of `phase` from x0 with initial canonical momentum p0 (on the
H = 0 shell of the initial mode). Coordinates: Spherical1D (r, θ) with
p = (p_r, p_θ); Cartesian (x, z) / (x, y, z) with z depth. `bbox` (Cartesian):
a vector of (min, max) per coordinate; the ray dies when it leaves it.
`totalReflection`: see taupWalk (post-critical default transmissions reflect).
"""
function traceRay(med, phase::Phase, x0, p0; H::AbstractHamiltonian=IsotropicHamiltonian(),
        dT=0.5, Tmax=5000.0, keepPath=true, bbox=nothing, region=nothing, totalReflection=true)
    geom = geometry(med)
    x = collect(Float64, x0); p = collect(Float64, p0)
    mode = phase.mode
    q = physicalSlowness(geom, x, p)
    goingDown = geom isa Spherical1D ? q[1] < 0 : q[end] > 0
    k = region === nothing ? initialRegion(med, x, goingDown) : region
    T = 0.0; j = 1
    path = [copy(x)]; pathT = [0.0]; pathMode = [mode]
    hits = RayHit[]
    maxres = abs(hamiltonianResidual(H, geom, med, k, x, p, mode))
    finish(ok, reason) = RayResult(phase, ok, reason, T, collect(Float64, x0), collect(Float64, p0),
        x, p, path, pathT, pathMode, hits, maxres)
    nstep = 0
    while T < Tmax
        nstep += 1
        nstep > 10^6 && return finish(false, :maxSteps)
        h = min(dT, Tmax - T, stepLimit(med, k, x, mode))
        x1, p1 = rk4(H, geom, med, k, x, p, mode, h)
        # earliest boundary crossed during this step
        sBest = 2.0; bBest = 0
        for b in boundaryIds(med, k)
            f0 = insideDistance(med, k, b, x)
            f1 = insideDistance(med, k, b, x1)
            (f1 < 0 && f1 < f0) || continue
            lo, hi = 0.0, 1.0
            for _ in 1:60
                mid = (lo + hi) / 2
                xm, _ = rk4(H, geom, med, k, x, p, mode, mid * h)
                insideDistance(med, k, b, xm) > 0 ? (lo = mid) : (hi = mid)
                (hi - lo) * h < 1e-12 && break
            end
            s = (lo + hi) / 2
            s < sBest && (sBest = s; bBest = b)
        end
        if bBest == 0
            x, p = x1, p1
            T += h
            if keepPath
                push!(path, copy(x)); push!(pathT, T); push!(pathMode, mode)
            end
            maxres = max(maxres, abs(hamiltonianResidual(H, geom, med, k, x, p, mode)))
            if bbox !== nothing && any(i -> !(bbox[i][1] <= x[i] <= bbox[i][2]), eachindex(bbox))
                return finish(false, :leftBox)
            end
            if geom isa Spherical1D && x[1] < 1.0
                return finish(false, :center)
            end
            continue
        end
        # move exactly onto the boundary
        x, p = rk4(H, geom, med, k, x, p, mode, sBest * h)
        T += sBest * h
        keepPath && (push!(path, copy(x)); push!(pathT, T); push!(pathMode, mode))
        dir = crossDirection(med, k, bBest)
        jBefore = j
        action, newMode, j = boundaryDecision(med, phase, j, bBest, dir, mode)
        isPass = action === :pass        # zone boundary without a jump: Snell keeps H = 0 exactly
        isPass && (action = :transmit)
        if action === :arrive
            push!(hits, RayHit(boundaryName(med, bBest), :arrive, mode, mode, copy(x), T))
            return finish(true, :arrived)
        elseif action === :dead
            push!(hits, RayHit(boundaryName(med, bBest), :dead, mode, mode, copy(x), T))
            return finish(false, Symbol("stopped_at_", boundaryName(med, bBest)))
        end
        kOut = action === :transmit ? across(med, k, bBest) : k
        n = outwardNormal(med, k, bBest, x)
        qin = physicalSlowness(geom, x, p)
        qout = interfaceSlowness(H, geom, med, kOut, x, qin, n, newMode, action)
        if qout === nothing && (isPass || (totalReflection && action === :transmit && newMode === mode && j == jBefore))
            # post-critical at a discontinuity that is not an event: total reflection
            action = :reflect; kOut = k
            qout = interfaceSlowness(H, geom, med, kOut, x, qin, n, newMode, action)
        end
        if qout === nothing
            push!(hits, RayHit(boundaryName(med, bBest), :evanescent, mode, newMode, copy(x), T))
            return finish(false, Symbol("evanescent_at_", boundaryName(med, bBest)))
        end
        (isPass && action === :transmit) || push!(hits, RayHit(boundaryName(med, bBest), action, mode, newMode, copy(x), T))
        p = canonicalMomentum(geom, x, qout)
        mode = newMode
        k = kOut
        keepPath && (push!(path, copy(x)); push!(pathT, T); push!(pathMode, mode))
    end
    return finish(false, :Tmax)
end

# ---------------------------------------------------------------- initial conditions
"""
    sphericalStart(med, depth, takeoff, mode) -> (x0, p0)

Source at `depth` km, takeoff angle (rad) measured from the downward
vertical (0 = straight down, π/2 = horizontal, π = straight up).
p_θ = r sin(i)/v is the ray parameter in s/rad.
"""
function sphericalStart(med::VerticalModel{Spherical1D}, depth, takeoff, mode)
    r = med.surface - depth
    k = zoneAt(med, r; side=cos(takeoff) > 0 ? :below : :above)
    v = zoneVelocity(med, k, r, mode)
    return [r, 0.0], [-cos(takeoff) / v, r * sin(takeoff) / v]
end
"start with a given ray parameter p (s/rad) and direction (:down / :up)"
function sphericalStartP(med::VerticalModel{Spherical1D}, depth, p, mode, dir::Symbol)
    r = med.surface - depth
    k = zoneAt(med, r; side=dir === :down ? :below : :above)
    v = zoneVelocity(med, k, r, mode)
    s = 1 / v^2 - p^2 / r^2
    s < 0 && return nothing
    return [r, 0.0], [dir === :down ? -sqrt(s) : sqrt(s), p]
end

"Cartesian 2-D: source (x, z), takeoff from the downward vertical toward +x"
function cartesianStart(med, xs::AbstractVector, takeoff, mode; azimuth=0.0)
    goingDown = cos(takeoff) > 0
    k = initialRegion(med, xs, goingDown)
    v = velocity(med, k, xs, mode)
    if length(xs) == 2
        return collect(Float64, xs), [sin(takeoff), cos(takeoff)] ./ v
    else
        return collect(Float64, xs), [sin(takeoff) * cos(azimuth), sin(takeoff) * sin(azimuth), cos(takeoff)] ./ v
    end
end
rayStart(::IsotropicHamiltonian, med, xs, takeoff, mode; azimuth=0.0) = cartesianStart(med, xs, takeoff, mode; azimuth=azimuth)
