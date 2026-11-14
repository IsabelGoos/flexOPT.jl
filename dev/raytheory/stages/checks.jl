# dev/raytheory/stages/checks.jl
#
# All numerical checks of the ray-theory core, without any plotting package
# (step `check`; s1/s2/s3 call the relevant part and then plot). Every check
# compares two INDEPENDENT computations:
#   tracer (RK4 on Hamilton's equations)  vs  τ-p integrals  vs  closed form / TauP / Fermat.

struct RTCheck
    name::String
    value::Float64
    reference::Float64
    tol::Float64
end
passed(c::RTCheck) = abs(c.value - c.reference) <= c.tol
function printChecks(cs::Vector{RTCheck}; title="")
    isempty(title) || println("\n--- ", title, " ---")
    for c in cs
        @printf("  %-58s %16.9f  ref %16.9f  |Δ| %.1e  %s\n", c.name, c.value, c.reference,
            abs(c.value - c.reference), passed(c) ? "ok" : "FAIL")
    end
    @printf("  %d / %d passed\n", count(passed, cs), length(cs))
    return cs
end

# ------------------------------------------------------------------ Cartesian (isotropic)
function cartesianModels()
    v0, g = 2.0, 0.5
    grad = layeredCartesian(Cartesian2D(); tops=[0.0], vp=[v0], vs=[v0 / sqrt(3)], gp=[g], gs=[g / sqrt(3)], zmax=1e3, name="gradient")
    H, vp1, vs1, vp2, vs2 = 10.0, 5.0, 3.0, 7.0, 4.0
    layer = layeredCartesian(Cartesian2D(); tops=[0.0, H], vp=[vp1, vp2], vs=[vs1, vs2], zmax=1e3, names=[:moho], name="layer")
    dip = deg2rad(15)
    nd = [-sin(dip), cos(dip)]
    dipping = PiecewiseMedium(Cartesian2D(),
        [PlaneInterface([0.0, 0.0], [0.0, 1.0], :surface; freeSurface=true), PlaneInterface([0.0, 8.0], nd, :refl)],
        [linearRegion([0, 0], 4.0, [0, 0], 2.3, [0, 0]), linearRegion([0, 0], 6.0, [0, 0], 3.5, [0, 0])])
    n3 = normalize([0.2, -0.25, 1.0])
    dipping3 = PiecewiseMedium(Cartesian3D(),
        [PlaneInterface([0.0, 0, 0], [0.0, 0, 1], :surface; freeSurface=true), PlaneInterface([0.0, 0, 9.0], n3, :refl)],
        [linearRegion([0, 0, 0], 4.0, [0, 0, 0], 2.3, [0, 0, 0]), linearRegion([0, 0, 0], 6.0, [0, 0, 0], 3.5, [0, 0, 0])])
    grad3 = layeredCartesian(Cartesian3D(); tops=[0.0], vp=[v0], vs=[v0 / sqrt(3)], gp=[g], gs=[g / sqrt(3)], zmax=1e3)
    return (; v0, g, grad, H, vp1, vs1, layer, nd, dipping, n3, dipping3, grad3)
end

function checkCartesian()
    M = cartesianModels()
    cs = RTCheck[]
    # 1. constant gradient: T = (2/g) asinh(g X / 2 v0)
    for X in (5.0, 20.0, 60.0)
        Tan = 2 / M.g * asinh(M.g * X / (2M.v0))
        push!(cs, RTCheck("gradient v0+gz, X=$X km: τ-p", taupArrivals(M.grad, "P", 0.0, X)[1].time, Tan, 1e-7))
        push!(cs, RTCheck("gradient v0+gz, X=$X km: tracer (shooting)", shootCartesian2D(M.grad, "P", [0.0, 0.0], X; nfan=361, dT=0.02, Tmax=100.0)[1].time, Tan, 1e-7))
    end
    # 2. layer over half-space: PmP (image source) and PmS (Snell solved by bisection)
    H, vp1, vs1 = M.H, M.vp1, M.vs1
    for X in (5.0, 30.0)
        f(p) = H * p * vp1 / sqrt(1 - (p * vp1)^2) + H * p * vs1 / sqrt(1 - (p * vs1)^2) - X
        lo, hi = 0.0, (1 - 1e-12) / vp1
        for _ in 1:200
            mid = (lo + hi) / 2; f(mid) > 0 ? (hi = mid) : (lo = mid)
        end
        p = (lo + hi) / 2
        Tps = H / (vp1 * sqrt(1 - (p * vp1)^2)) + H / (vs1 * sqrt(1 - (p * vs1)^2))
        for (ph, Tan) in (("PmP", sqrt(X^2 + 4H^2) / vp1), ("PmS", Tps))
            push!(cs, RTCheck("layer $ph, X=$X km: τ-p", taupArrivals(M.layer, ph, 0.0, X)[1].time, Tan, 1e-7))
            push!(cs, RTCheck("layer $ph, X=$X km: tracer", shootCartesian2D(M.layer, ph, [0.0, 0.0], X; nfan=361, dT=0.05, Tmax=100.0)[1].time, Tan, 1e-7))
        end
    end
    # 3. dipping reflector 2-D (image source)
    xs = [0.0, 0.0]
    img = xs .- 2 * dot(M.nd, xs .- [0.0, 8.0]) .* M.nd
    for X in (-10.0, 5.0, 20.0)
        s = shootCartesian2D(M.dipping, "Pv{refl}P", xs, X; nfan=361, dT=0.05, Tmax=60.0, bbox=[(-200.0, 200.0), (-1.0, 200.0)])
        push!(cs, RTCheck("dipping (15°) reflector 2-D, X=$X: tracer", s[1].time, norm([X, 0.0] .- img) / 4.0, 1e-7))
    end
    # 4. dipping reflector 3-D with oblique strike (Newton shooting on takeoff+azimuth)
    xs3 = [0.0, 0, 0]
    img3 = xs3 .- 2 * dot(M.n3, xs3 .- [0.0, 0, 9.0]) .* M.n3
    for xr in ([10.0, 5.0, 0.0], [-6.0, 12.0, 0.0])
        s = shootCartesian3D(M.dipping3, "Pv{refl}P", xs3, xr; dT=0.02, Tmax=60.0, bbox=[(-200.0, 200.0), (-200.0, 200.0), (-1.0, 200.0)])
        push!(cs, RTCheck("oblique dipping reflector 3-D, xr=$(xr[1:2]): tracer", s.time, norm(xr .- img3) / 4.0, 1e-7))
    end
    # 5. 3-D gradient, receiver off-axis
    X = hypot(12, 9)
    s = shootCartesian3D(M.grad3, "P", [0.0, 0, 0], [12.0, -9.0, 0.0]; dT=0.02, Tmax=60.0)
    push!(cs, RTCheck("gradient 3-D, xr=(12,-9): tracer", s.time, 2 / M.g * asinh(M.g * X / (2M.v0)), 1e-7))
    return cs
end

# ------------------------------------------------------------------ Spherical 1-D
"tracer vs τ-p at the same ray parameter (two independent formulations)"
function checkSphericalTracer(m; cases=[("P", 0.0, 60.0), ("PP", 0.0, 140.0), ("PcP", 0.0, 50.0), ("ScS", 0.0, 50.0),
        ("PKIKP", 0.0, 140.0), ("PKiKP", 0.0, 140.0), ("SKS", 0.0, 140.0), ("PKP", 0.0, 150.0), ("P410s", 0.0, 60.0),
        ("Pms", 0.0, 60.0), ("pP", 100.0, 40.0), ("sP", 100.0, 40.0), ("S", 500.0, 70.0)], dT=0.25)
    cs = RTCheck[]
    for (ph, dep, d) in cases
        for a in taupArrivals(m, ph, dep, d)
            x0, p0 = sphericalStartP(m, dep, a.p, asPhase(ph).mode, a.startDir)
            r = traceRay(m, asPhase(ph), x0, p0; dT=dT, keepPath=false)
            push!(cs, RTCheck(@sprintf("%s h=%g Δ=%g p=%.2f: tracer T vs τ-p", ph, dep, d, a.p), r.T, a.time, 1e-5))
            fold(D) = (D = mod(D, 360.0); D > 180 ? 360 - D : D)      # major-arc roots (e.g. PKP at 210°)
            push!(cs, RTCheck(@sprintf("%s h=%g Δ=%g p=%.2f: tracer Δ vs τ-p", ph, dep, d, a.p), fold(rad2deg(r.xEnd[2])), a.distance, 1e-6))
        end
    end
    return cs
end

"""
τ-p arrivals vs the TauP table (ObsPy TauP on the SAME model, reference/).
Each TauP arrival is matched by ray parameter. TauP also lists 'arrivals'
whose ray grazes a zone boundary (p = η there, the same p at many
distances: e.g. PS with p = 568.13 s/rad, the P leg grazing PREM's 600-km
boundary, where the published coefficients leave a 4e-3 km/s step that TauP
treats as a reflector). They are branch ends, not geometric rays of our
model, and are counted separately ("grazing"), not as failures.
"""
function checkTauP(m, csvfile; tol=5e-3)
    rows = [split(l, ',') for l in readlines(csvfile)[2:end]]
    groups = Dict{Tuple{String,Float64,Float64},Vector{Tuple{Float64,Float64}}}()
    for r in rows
        push!(get!(groups, (String(r[1]), parse(Float64, r[3]), parse(Float64, r[4])), Tuple{Float64,Float64}[]),
            (parse(Float64, r[5]), parse(Float64, r[6])))
    end
    # slownesses of rays grazing a zone boundary (η at every boundary, P and S)
    crit = Float64[]
    for k in 2:nzones(m), kk in (k - 1, k), md in (:P, :S)
        b = m.boundaries[k]
        v = zoneVelocity(m, kk, b.zeta, md); v > 0 && push!(crit, b.zeta / v)
    end
    stats = Dict{String,Vector{Float64}}()   # phase => [n, matched, branchEnd, unmatched, maxdT]
    unmatched = String[]
    for (key, tp) in sort(collect(groups))
        ph, z, d = key
        ours = taupArrivals(m, ph, z, d)
        st = get!(stats, ph, zeros(5))
        for (T, p) in tp
            st[1] += 1
            i = isempty(ours) ? 0 : argmin([abs(a.p - p) for a in ours])
            if i > 0 && abs(ours[i].p - p) < 2.0 && abs(ours[i].time - T) < 0.05
                st[2] += 1; st[5] = max(st[5], abs(ours[i].time - T))
            elseif any(c -> abs(c - p) < 0.05, crit)
                st[3] += 1
            else
                st[4] += 1
                push!(unmatched, @sprintf("%s h=%g Δ=%g T=%.3f p=%.3f", ph, z, d, T, p))
            end
        end
    end
    @printf("  %-7s %7s %8s %11s %10s %10s\n", "phase", "TauP", "matched", "grazing", "unmatched", "max|ΔT| s")
    cs = RTCheck[]
    for ph in sort(collect(keys(stats)))
        st = stats[ph]
        @printf("  %-7s %7d %8d %11d %10d %10.4f\n", ph, st[1], st[2], st[3], st[4], st[5])
        push!(cs, RTCheck("TauP $ph: max |T - T_TauP| over matched arrivals", st[5], 0.0, tol))
        push!(cs, RTCheck("TauP $ph: unmatched geometric arrivals", st[4], 0.0, 0.0))
    end
    isempty(unmatched) || println("  unmatched: ", join(unmatched, "; "))
    return cs
end

# ------------------------------------------------------------------ anisotropy
"exact VTI qP phase velocity (Thomsen 1986), θ from the symmetry axis — independent of rayChristoffel.jl"
function vtiPhaseVelocity(A, θ)
    C11, C33, C13, C44 = A[1, 1], A[3, 3], A[1, 3], A[4, 4]
    s2 = sin(θ)^2
    D = sqrt((C33 - C44)^2 + 2 * (2 * (C13 + C44)^2 - (C33 - C44) * (C11 + C33 - 2C44)) * s2 +
             ((C11 + C33 - 2C44)^2 - 4 * (C13 + C44)^2) * s2^2)
    return sqrt((C33 + C44 + (C11 - C33) * s2 + D) / 2)
end
function vtiGroup(A, θ)
    v = vtiPhaseVelocity(A, θ); h = 1e-6
    dv = (vtiPhaseVelocity(A, θ + h) - vtiPhaseVelocity(A, θ - h)) / 2h
    return hypot(v, dv), θ + atan(dv / v)         # |V_g|, group angle
end
"straight-ray qP time between two points of a homogeneous VTI medium (closed form + 1-D root)"
function vtiTime(A, Δx, Δz)
    φ = atan(abs(Δx), abs(Δz))
    lo, hi = 0.0, π / 2
    for _ in 1:100
        mid = (lo + hi) / 2; vtiGroup(A, mid)[2] < φ ? (lo = mid) : (hi = mid)
    end
    return hypot(Δx, Δz) / vtiGroup(A, (lo + hi) / 2)[1]
end

anisoSurface2D() = PlaneInterface([0.0, 0], [0.0, 1], :surface; freeSurface=true)
const QP_UP = Phase("qP", :qP, :up, PhaseEvent[])

function checkAnisotropy()
    H = ChristoffelHamiltonian()
    cs = RTCheck[]
    Ai = isotropicA(6.0, 3.5)
    med = PiecewiseMedium(Cartesian2D(), [anisoSurface2D()], [ARegion(Ai)])
    s = shootCartesian2D(med, QP_UP, [0.0, 10.0], 7.0; H=H, nfan=181, dT=0.2, Tmax=20.0)
    push!(cs, RTCheck("isotropic tensor through Christoffel, qP", s[1].time, hypot(7, 10) / 6, 1e-9))
    Av = vtiA(3.0, 1.5, 0.25, 0.10, 0.2)
    medv = PiecewiseMedium(Cartesian2D(), [anisoSurface2D()], [ARegion(Av)])
    for X in (5.0, 15.0, 40.0)
        s = shootCartesian2D(medv, QP_UP, [0.0, 10.0], X; H=H, nfan=181, dT=0.2, Tmax=60.0)
        push!(cs, RTCheck("VTI (ε=.25, δ=.1) qP, X=$X: tracer vs exact group velocity", s[1].time, vtiTime(Av, X, 10.0), 1e-8))
    end
    med2 = PiecewiseMedium(Cartesian2D(), [anisoSurface2D(), PlaneInterface([0.0, 8.0], [0.0, 1.0], :i1)],
        [ARegion(isotropicA(2.5, 1.4)), ARegion(Av)])
    for X in (6.0, 25.0)
        s = shootCartesian2D(med2, QP_UP, [0.0, 20.0], X; H=H, nfan=181, dT=0.2, Tmax=60.0)
        f(xc) = vtiTime(Av, xc, 12.0) + hypot(X - xc, 8.0) / 2.5
        a, b = 0.0, X
        for _ in 1:200
            m1 = a + 0.382(b - a); m2 = a + 0.618(b - a); f(m1) < f(m2) ? (b = m2) : (a = m1)
        end
        push!(cs, RTCheck("VTI → isotropic transmission, X=$X: tracer vs Fermat", s[1].time, f((a + b) / 2), 1e-8))
    end
    med4 = PiecewiseMedium(Cartesian2D(), [anisoSurface2D(), PlaneInterface([0.0, 10.0], [0.0, 1.0], :moho)],
        [ARegion(isotropicA(5.0, 3.0)), ARegion(isotropicA(7.0, 4.0))])
    ph = Phase("qPv{moho}qS", :qP, :down, [PhaseEvent(:moho, :reflect, :qS1, :down)])
    s = shootCartesian2D(med4, ph, [0.0, 0.0], 5.0; H=H, nfan=181, dT=0.2, Tmax=60.0)
    tiso = shootCartesian2D(cartesianModels().layer, "PmS", [0.0, 0.0], 5.0; nfan=361, dT=0.05, Tmax=100.0)[1].time
    push!(cs, RTCheck("P→S reflection: Christoffel vs isotropic Hamiltonian", s[1].time, tiso, 1e-8))
    med3 = PiecewiseMedium(Cartesian3D(), [PlaneInterface([0.0, 0, 0], [0.0, 0, 1], :surface; freeSurface=true)], [ARegion(Av)])
    s = shootCartesian3D(med3, QP_UP, [0.0, 0, 10.0], [9.0, 12.0, 0.0]; H=H, dT=0.2, Tmax=60.0, guess=(π - 0.8, atan(12, 9)))
    push!(cs, RTCheck("VTI qP 3-D, receiver (9,12): tracer vs exact", s.time, vtiTime(Av, 15.0, 10.0), 1e-8))
    return cs
end

# ------------------------------------------------------------------ gridded media (Marmousi)
"""
Gridded-model checks: (1) tracer on a GRIDDED constant-gradient model vs the
closed form (tests the B-spline + its gradient, edges included), (2) eikonal
convergence on the same model, (3) on the smoothed model `gp` (Marmousi):
earliest traced ray vs Richardson-extrapolated eikonal first arrival at
surface receivers `Xs` (only where a ray IS the first arrival).
"""
function checkGrid(gp=nothing; xs=3.0, Xs=6.5:0.5:8.5, tol=3e-3)
    cs = RTCheck[]
    v0, g, h = 2.0, 0.5, 0.025
    G = GridField2D([v0 + g * (j - 1) * h for i in 1:401, j in 1:201], 0.0, 0.0, h, h)
    med, bbox = gridMedium2D(G)
    for X in (2.0, 5.0, 9.0)
        s = shootCartesian2D(med, "P", [0.5, 0.0], 0.5 + X; nfan=721, dT=0.01, Tmax=20.0, bbox=bbox)
        push!(cs, RTCheck("gridded gradient, X=$X km: tracer vs closed form", s[1].time, 2 / g * asinh(g * X / (2v0)), 1e-7))
    end
    e = [(T = eikonalOnGrid(G, 0.5, 0.0; h=hh); T[1][argmin(abs.(T[2] .- 9.5)), 1]) for hh in (0.0125, 0.00625)]
    push!(cs, RTCheck("gridded gradient, X=9 km: eikonal (Richardson 12.5/6.25 m)", 2 * e[2] - e[1], 2 / g * asinh(g * 9 / (2v0)), 5e-3))
    gp === nothing && return cs
    med, bbox = gridMedium2D(gp)
    E = [eikonalOnGrid(gp, xs, 0.0; h=hh) for hh in (0.005, 0.0025)]
    for X in Xs
        s = shootCartesian2D(med, "P", [xs, 0.0], X; nfan=2881, dT=0.005, Tmax=10.0, bbox=bbox)
        isempty(s) && continue
        e = [T[argmin(abs.(xg .- X)), 1] for (T, xg, _) in E]
        push!(cs, RTCheck("Marmousi, offset $(X - xs) km: first ray vs eikonal (Richardson)", s[1].time, 2 * e[2] - e[1], tol))
    end
    return cs
end
