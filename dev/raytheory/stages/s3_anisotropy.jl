# dev/raytheory/stages/s3_anisotropy.jl
# Layer 2 (Christoffel dispersion relation) and anisotropic interface matching:
#   checks: isotropic tensor, homogeneous VTI vs exact group velocity, VTI→iso
#   transmission vs Fermat, P→S reflection, 3-D VTI;
#   figures: slowness and wave (group-velocity) curves, qP rays + wavefronts in
#   VTI and in a 30°-tilted TI (rays ∦ phase normals), anisotropic Snell.
printChecks(checkAnisotropy(); title="anisotropy (Christoffel Hamiltonian)")

Hc = ChristoffelHamiltonian()
Av = vtiA(3.0, 1.5, 0.25, 0.10, 0.2)            # Thomsen: vp0, vs0, ε, δ, γ
At = rotateA(Av, rotY(deg2rad(30)))             # same, symmetry axis tilted 30° toward +x
fig = Figure(size=(1600, 1000))

# (a) slowness curves and (b) wave surfaces of the three branches
axS = Axis(fig[1, 1]; title="VTI slowness curves (s/km)", xlabel="p_x", ylabel="p_z", aspect=DataAspect(), yreversed=true)
axW = Axis(fig[2, 1]; title="VTI wave surfaces: group velocity (km/s)", xlabel="V_x", ylabel="V_z", aspect=DataAspect(), yreversed=true)
θs = range(0, 2π; length=1441)
for (mi, md) in enumerate((:qP, :qS1, :qS2))
    P = Vector{Point2f}(); W = Vector{Point2f}()
    for θ in θs
        n = [sin(θ), 0.0, cos(θ)]
        λ, V = branches(christoffel(Av, n))
        p = n ./ sqrt(λ[mi])
        vg = groupVelocity(Av, V[:, mi], p)
        push!(P, Point2f(p[1], p[3])); push!(W, Point2f(vg[1], vg[3]))
    end
    lines!(axS, P; color=RT_MODE_COLOR[md], label=String(md))
    scatter!(axW, W; color=RT_MODE_COLOR[md], markersize=1.5)
end
axislegend(axS; position=:rt, framevisible=false)

"rays of a fan of phase-normal angles from a point source, wavefronts at times Ts"
function anisoFan!(ax, med, xs, Ts; nray=72, Tmax=maximum(Ts), mode=:qP, normals=true)
    rays = RayResult[]
    for ψ in range(0, 2π; length=nray + 1)[1:end-1]
        x0, p0 = christoffelStart(med, xs, [sin(ψ), cos(ψ)], mode)
        r = traceRay(med, Phase(String(mode), mode, :any, PhaseEvent[]), x0, p0; H=Hc, dT=0.05, Tmax=Tmax)
        push!(rays, r)
        plotRayPath!(ax, r; linewidth=0.6, alpha=0.6)
    end
    for Tk in Ts
        pts = Point2f[]
        for r in rays
            i = findfirst(t -> t >= Tk - 1e-9, r.pathT)
            i === nothing || push!(pts, Point2f(r.path[i][1], r.path[i][end]))
        end
        length(pts) > 2 && lines!(ax, pts; color=:black, linewidth=1)
    end
    if normals                       # phase normal (red) vs ray (blue) on a few rays
        for r in rays[1:6:end]
            i = findfirst(t -> t >= Ts[end] - 1e-9, r.pathT)
            i === nothing && continue
            x = r.path[i]
            d = r.p0 ./ norm(r.p0) .* 2.5
            lines!(ax, [x[1], x[1] + d[1]], [x[end], x[end] + d[end]]; color=:red, linewidth=2)
        end
    end
    return rays
end

surf2 = PlaneInterface([0.0, 0], [0.0, 1], :surface; freeSurface=true)
ax = Axis(fig[1, 2]; title="VTI qP: rays (blue) and wavefronts; red = phase normal", xlabel="x (km)", ylabel="z (km)",
    aspect=DataAspect(), yreversed=true)
anisoFan!(ax, PiecewiseMedium(Cartesian2D(), [surf2], [ARegion(Av)]), [0.0, 25.0], [2.0, 4.0, 6.0])
ax = Axis(fig[1, 3]; title="TTI (axis tilted 30°) qP: asymmetric wavefronts", xlabel="x (km)", ylabel="z (km)",
    aspect=DataAspect(), yreversed=true)
anisoFan!(ax, PiecewiseMedium(Cartesian2D(), [surf2], [ARegion(At)]), [0.0, 25.0], [2.0, 4.0, 6.0])

# (d) anisotropic Snell: qP from a VTI half-space into an isotropic layer
med2 = PiecewiseMedium(Cartesian2D(), [surf2, PlaneInterface([0.0, 8.0], [0.0, 1.0], :i1)],
    [ARegion(isotropicA(2.5, 1.4)), ARegion(Av)])
ax = Axis(fig[2, 2]; title="qP: VTI half-space → isotropic layer (8 km), receivers every 5 km", xlabel="x (km)", ylabel="z (km)",
    aspect=DataAspect(), yreversed=true)
hlines!(ax, [8.0]; color=:gray40)
for X in 0:5:40
    s = shootCartesian2D(med2, QP_UP, [0.0, 20.0], Float64(X); H=Hc, nfan=181, dT=0.1, Tmax=60.0)
    isempty(s) && continue
    plotRayPath!(ax, s[1].ray; linewidth=1.5)
    @printf("  VTI→iso X = %4.1f km: T = %.6f s\n", X, s[1].time)
end
# (e) TTI: travel times to the left and right
ax = Axis(fig[2, 3]; title="qP travel time to the surface, source at 25 km", xlabel="x (km)", ylabel="T (s)")
for (lab, A, c) in (("VTI", Av, 1), ("TTI 30°", At, 2), ("isotropic (vp0)", isotropicA(3.0, 1.5), 3))
    med = PiecewiseMedium(Cartesian2D(), [surf2], [ARegion(A)])
    Xs = Float64[]; Ts = Float64[]
    for X in -40.0:4.0:40.0
        s = shootCartesian2D(med, QP_UP, [0.0, 25.0], X; H=Hc, nfan=181, dT=0.2, Tmax=60.0)
        isempty(s) || (push!(Xs, X); push!(Ts, s[1].time))
    end
    lines!(ax, Xs, Ts; color=rtColor(c), label=lab)
end
axislegend(ax; position=:ct, framevisible=false)
devSave(fig, "anisotropy")
