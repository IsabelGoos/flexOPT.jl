# dev/raytheory/lib/rayShooting.jl
#
# Two-point ray tracing by shooting with the general tracer (any medium).
# A fan of takeoff angles is traced; every sign change of
# (endpoint - receiver) between two neighbouring successful rays is refined
# by bisection. All branches (triplications, multipathing) are returned.

struct ShotArrival
    phase::String
    time::Float64
    takeoff::Float64          # degrees from the downward vertical
    ray::RayResult
    misfit::Float64           # endpoint error (deg or km)
end
Base.show(io::IO, a::ShotArrival) = @printf(io, "%-8s T=%9.3f s  takeoff=%7.3f°  misfit=%.1e  maxRes=%.1e",
    a.phase, a.time, a.takeoff, a.misfit, a.ray.maxResidual)

endpointDistance(::Spherical1D, ray) = rad2deg(ray.xEnd[2] - ray.x0[2])
endpointDistance(::Cartesian2D, ray) = ray.xEnd[1]

"""
    shootSpherical(med, phase, depth, Δdeg; nfan=721, dT=0.5, kw...)

Rays from a source at `depth` km, receiver at the surface at Δdeg (minor arc).
"""
function shootSpherical(med::VerticalModel{Spherical1D}, phase, depth, Δdeg; nfan=721, dT=0.5, Tmax=4000.0, tol=1e-7)
    ph = asPhase(phase)
    i1, i2 = ph.startDirection === :down ? (0.0, π / 2) : ph.startDirection === :up ? (π / 2, π) : (0.0, π)
    is = range(i1 + 1e-4, i2 - 1e-4; length=nfan)
    shoot(i) = traceRay(med, ph, sphericalStart(med, depth, i, ph.mode)...; dT=dT, Tmax=Tmax, keepPath=false)
    return fanAndBisect(shoot, is, r -> endpointDistance(Spherical1D(), r), Δdeg, ph, tol; keep=i -> traceRay(med, ph, sphericalStart(med, depth, i, ph.mode)...; dT=dT, Tmax=Tmax))
end

"""
    shootCartesian2D(med, phase, xs, X; nfan=721, ...)

Source at xs = (x, z), receiver at the surface (z = 0) at x = X.
Takeoffs span both sides (rays to -x too). For `H = ChristoffelHamiltonian()`
the takeoff is the PHASE-normal angle (the ray leaves along the group velocity).
"""
function shootCartesian2D(med, phase, xs, X; nfan=721, dT=0.01, Tmax=200.0, tol=1e-9, bbox=nothing, H::AbstractHamiltonian=IsotropicHamiltonian())
    ph = asPhase(phase)
    i1, i2 = ph.startDirection === :down ? (-π / 2, π / 2) : ph.startDirection === :up ? (π / 2, 3π / 2) : (-π, π)
    is = range(i1 + 1e-4, i2 - 1e-4; length=nfan)
    shoot(i; path=false) = traceRay(med, ph, rayStart(H, med, xs, i, ph.mode)...; H=H, dT=dT, Tmax=Tmax, keepPath=path, bbox=bbox)
    return fanAndBisect(shoot, is, r -> r.xEnd[1], X, ph, tol; keep=i -> shoot(i; path=true))
end

function fanAndBisect(shoot, is, endpoint, target, ph, tol; keep)
    rays = [shoot(i) for i in is]
    f = [r.ok ? endpoint(r) - target : NaN for r in rays]
    out = ShotArrival[]
    for n in 1:length(is)-1
        (rays[n].ok && rays[n+1].ok) || continue
        f[n] * f[n+1] <= 0 || continue
        # same event history on both sides (otherwise a branch jump, not a root)
        sig(r) = [(h.boundary, h.action, h.modeOut) for h in r.hits]
        sig(rays[n]) == sig(rays[n+1]) || continue
        lo, hi, flo = is[n], is[n+1], f[n]
        for _ in 1:100
            mid = (lo + hi) / 2
            r = shoot(mid)
            r.ok || break
            fm = endpoint(r) - target
            (fm * flo > 0) ? (lo = mid; flo = fm) : (hi = mid)
            abs(hi - lo) < 1e-14 && break
        end
        ia = (lo + hi) / 2
        r = keep(ia)
        (r.ok && abs(endpoint(r) - target) < sqrt(tol)) || continue
        push!(out, ShotArrival(ph.name, r.T, rad2deg(ia), r, abs(endpoint(r) - target)))
    end
    sort!(out; by=a -> a.time)
    return out
end

"""
    shootCartesian3D(med, phase, xs, xr; guess=nothing, maxit=30)

3-D two-point problem: Newton on (takeoff, azimuth) with a finite-difference
Jacobian of the surface endpoint (x, y). Initial guess: the 2-D fan in the
vertical plane through source and receiver (exact for laterally homogeneous
media), unless `guess = (takeoff, azimuth)` in radians.
"""
function shootCartesian3D(med, phase, xs, xr; guess=nothing, dT=0.01, Tmax=200.0, maxit=30, tol=1e-8, bbox=nothing, H::AbstractHamiltonian=IsotropicHamiltonian())
    ph = asPhase(phase)
    az0 = atan(xr[2] - xs[2], xr[1] - xs[1])
    shoot(i, a; path=false) = traceRay(med, ph, rayStart(H, med, xs, i, ph.mode; azimuth=a)...; H=H, dT=dT, Tmax=Tmax, keepPath=path, bbox=bbox)
    if guess === nothing
        # fan in the vertical plane
        best = nothing; bestErr = Inf
        for i in range(1e-3, π / 2 - 1e-3; length=181)
            r = shoot(i, az0)
            r.ok || continue
            e = hypot(r.xEnd[1] - xr[1], r.xEnd[2] - xr[2])
            e < bestErr && (bestErr = e; best = i)
        end
        best === nothing && return nothing
        u = [best, az0]
    else
        u = collect(Float64, guess)
    end
    r = shoot(u...)
    for it in 1:maxit
        r.ok || return nothing
        F = [r.xEnd[1] - xr[1], r.xEnd[2] - xr[2]]
        norm(F) < tol && break
        J = zeros(2, 2); δ = 1e-6
        for c in 1:2
            du = copy(u); du[c] += δ
            rr = shoot(du...)
            rr.ok || return nothing
            J[:, c] = ([rr.xEnd[1] - xr[1], rr.xEnd[2] - xr[2]] .- F) ./ δ
        end
        step = J \ F
        λ = 1.0
        while λ > 1e-4
            un = u .- λ .* step
            rn = shoot(un...)
            if rn.ok && hypot(rn.xEnd[1] - xr[1], rn.xEnd[2] - xr[2]) < norm(F)
                u = un; r = rn; break
            end
            λ /= 2
        end
        λ <= 1e-4 && break
    end
    r = shoot(u...; path=true)
    return ShotArrival(ph.name, r.T, rad2deg(u[1]), r, hypot(r.xEnd[1] - xr[1], r.xEnd[2] - xr[2]))
end
