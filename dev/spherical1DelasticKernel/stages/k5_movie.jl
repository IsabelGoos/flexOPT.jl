# dev/spherical1DelasticKernel/stages/k5_movie.jl
# K5: time-domain waveform-kernel movie on the great-circle section, WHOLE Earth (mantle,
# fluid outer core, inner core), short period.
#
#   δu_Z(R, t) = ∫ K(Q, t) δln α(Q) dV,
#   K(Q, t) = (Δω/π) e^{ω_I t} Re Σ_k w'_k K_k(Q) e^{iω_k t}      (same synthesis as diskSnapshot)
#
# with K_k the frequency-domain δln α (and δln β) kernel densities of lib/sgtKernels1D.jl
# (TI parameters, η = F/(A-2L) kept; in the fluid core L = N = 0 → δln α = δln κ kernel).
# Source: moment tensor k5M at k5Depth km below (lat 0, lon 0); receiver: Z at (0, k5Δ°).
# Section: the equator, longitudes k5λ (S at 0°, R at k5Δ), radii = the disk radii
# (every k5Stride-th node of the P-SV grid, whole Earth, both sides of each discontinuity).
#
# Memory/speed: the disks are NOT stored for all frequencies.  For each frequency (in
# parallel over threads) the :moment and :force disks are computed on the exact θ grid of
# the section (no θ interpolation), the kernel K_k is evaluated on the section and cached
# in data/dev/spherical1DelasticKernel/k5/ (one small file per frequency; an interrupted run
# resumes), then the disks are dropped.  At the end the snapshots K(Q, t) for t in k5Times
# are synthesised and saved (Float32) together with the reference trace u_Z(t).
#
# Run with threads (one frequency per thread; BLAS is set to 1 thread):
#   julia --project=. -t auto dev/spherical1DelasticKernel/run.jl k5          # quick test, T ≥ 40 s
#   julia --project=. -t auto dev/spherical1DelasticKernel/run.jl k5 --full   # the movie, T ≥ 10 s
# Rough cost (extrapolated from k4): a few hours with ~10 threads, ≈ 1-1.5 GB RAM per thread.
# Lower k5fmax (e.g. 0.05) / coarser k5θstep for a quick test.

using LinearAlgebra
BLAS.set_num_threads(1)

k5Scheme = "OPT3"; k5h = 2.5                 # radial grid of the systems (km); OPT3 needs ~h ≤ 2.5 km at 10 s
k5Stride = 4                                 # section radii every k5Stride·h = 10 km
k5θstep = 0.25                               # section / disk θ step (°)
k5λ = -15.0:k5θstep:165.0                    # section longitudes (°)
k5Depth = 150.0; k5Δ = 150.0
k5M = (1.0, 0.6, -0.8, -0.3, 0.4, -0.7)      # (Mrr, Mrθ, Mrφ, Mθθ, Mθφ, Mφφ) × 1e25 dyn·cm (same MT as p6/p8)
k5tlen = 2048.0; k5fmax = 0.1; k5Oi = 1e-3   # periods ≥ 10 s (flat to 0.08 Hz, cosine taper to 0.1 Hz)
k5Taper = 0.2
k5Times = 1120.0:0.5:1240.0                  # movie frames (s)
k5Tag = "PKP150"
if !DEV.full    # quick test (minutes): T ≥ 40 s, coarse grids; the real movie needs --full
    k5h = 10.0; k5Stride = 2; k5θstep = 0.5; k5λ = -15.0:k5θstep:165.0; k5fmax = 0.025
    k5Tag = "PKP150_quick"
    println("k5 QUICK TEST (T ≥ 40 s, h = 10 km); add --full for the T ≥ 10 s movie")
end

prem = planet1D.getSet1Dmodel!(:PREM)
Rp = prem.averagedPlanetRadiusInKilometer
k5S = PointGeo(Rp - k5Depth, π / 2, 0.0)
k5R = PointGeo(Rp, π / 2, deg2rad(k5Δ))
k5n = receiverDirection(:Z, k5S, k5R)
sys = kerSystems("PREM", k5Scheme, k5h)
psv, sh = sys["psv"], sys["sh"]
k5Radii = diskRadii(psv; stride=k5Stride)
k5θdisk = collect(0.0:k5θstep:maximum(abs.(collect(k5λ) .- [0.0 k5Δ])) + 1e-9)
nr, nλ = length(k5Radii.r_km), length(k5λ)
k5ωs = [2π * k / k5tlen for k in 0:floor(Int, k5fmax * k5tlen)]
nω = length(k5ωs)
@printf("k5: %s h = %.1f km, %d section radii × %d longitudes, %d frequencies (f ≤ %.3f Hz), θ step %.2f°, %d threads\n",
        k5Scheme, k5h, nr, nλ, nω, k5fmax, k5θstep, Threads.nthreads())

# ---- geometry of the section (independent of r and ω): for every longitude, the θ index in
# each disk, the source rotated about its pole (b) and the basis change T to Q's local frame
function k5Geometry(family, P, a)
    Rs, ap = sourceToPole(family, a, P)
    out = Vector{Any}(undef, nλ)
    for (j, λ) in enumerate(k5λ)
        Q = PointGeo(Rp, π / 2, deg2rad(λ))
        θp, φp, T = frameOf(Rs, Q)
        iθ = round(Int, θp / k5θstep) + 1
        abs(k5θdisk[iθ] - θp) < 1e-6 || error("section point not on the disk θ grid ($θp°)")
        out[j] = (iθ=iθ, b=rotateSourceAboutPole(family, ap, φp), T=T)
    end
    return out
end
k5GS = k5Geometry(:moment, k5S, k5M)
k5GR = k5Geometry(:force, k5R, k5n)

# materials per section radius (zone side as stored in the disk)
k5Mat = [begin
    iz = _zoneOfRadius(prem, r; side=sd)
    m = elasticMaterials(prem, iz, r / Rp)
    η = abs(m.A - 2m.L) > 0 ? m.F / (m.A - 2m.L) : 1.0
    (A=m.A, C=m.C, F=m.F, L=m.L, N=m.N, ρ=m.ρ, η=η)
end for (r, sd) in zip(k5Radii.r_km, k5Radii.side)]

"Field (local basis at Q) from one disk at one frequency: strain 3×3 and displacement 3-vector."
function k5Field(disk, g, ir)
    e = zeros(ComplexF64, 6); u = zeros(ComplexF64, 3)
    for (ip, (s, c)) in enumerate(disk.pairs.strain)
        iszero(g.b[s]) || (e[c] += g.b[s] * disk.strain[1, ip, ir, g.iθ])
    end
    for (ip, (s, c)) in enumerate(disk.pairs.disp)
        iszero(g.b[s]) || (u[c] += g.b[s] * disk.disp[1, ip, ir, g.iθ])
    end
    return g.T * strainMatrix(e) * g.T', g.T * u
end

"Kernel densities δln α, δln β of one frequency on the section (and the trace sample u_Z(R))."
function k5Frequency(k)
    ω = k5ωs[k]
    kw = (θdeg=k5θdisk, radii=k5Radii, tlen=k5tlen, fmax=k5fmax, omegaiTlen=k5Oi, ωlist=[ω],
          p=sourceOrderOf(k5Scheme), verbose=false)
    dS = computeSGTDisk(psv, sh, prem; family=:moment, rs_km=k5S.r, kw...)
    dR = computeSGTDisk(psv, sh, prem; family=:force, rs_km=k5R.r, kw...)
    ωc = ω - im * dS.ωI
    Kα = zeros(ComplexF32, nr, nλ); Kβ = zeros(ComplexF32, nr, nλ)
    for ir in 1:nr, j in 1:nλ
        εS, uS = k5Field(dS, k5GS[j], ir)
        εR, uR = k5Field(dR, k5GR[j], ir)
        K = kernelDensities(εR, εS, uR, uS, ωc)
        m = k5Mat[ir]
        Kα[ir, j] = 2m.A * K.A + 2m.C * K.C + 2m.η * m.A * K.F
        Kβ[ir, j] = 2m.L * K.L + 2m.N * K.N - 4m.η * m.L * K.F
    end
    # reference trace sample: displacement at R (Z) of the source, from the :moment disk
    uZ = ComplexF64(dot(collect(k5n), sgtAt(dS, 1, k5S, k5M, k5R; which=:disp)))
    return (Kα=Kα, Kβ=Kβ, uZ=uZ, ωI=dS.ωI, lmax=dS.lmax[1], tailS=maximum(dS.tail), tailR=maximum(dR.tail))
end

k5Params = Dict{String,Any}("stage" => "k5_freq", "codeVersion" => CODE_VERSION_K["k5_freq"],
    "scheme" => k5Scheme, "h" => k5h, "stride" => k5Stride, "θstep" => k5θstep,
    "λ" => (first(k5λ), last(k5λ)), "depth" => k5Depth, "Δ" => k5Δ, "M" => k5M,
    "tlen" => k5tlen, "fmax" => k5fmax, "omegaiTlen" => k5Oi)
k5Dir = datadir(KER_CACHE_DIR, "k5_" * commonBatchs.hash_parameters(k5Params))
mkpath(k5Dir)
"redo" in DEV.redo || "all" in DEV.redo ? foreach(f -> rm(joinpath(k5Dir, f)), readdir(k5Dir)) : nothing
k5File(k) = joinpath(k5Dir, @sprintf("freq_%04d.jld2", k))
println("per-frequency cache: ", k5Dir)

k5Todo = [k for k in 1:nω if !isfile(k5File(k))]
println(length(k5Todo), " of ", nω, " frequencies to compute")
k5Lock = ReentrantLock(); k5Done = Threads.Atomic{Int}(0); k5t0 = time()
Threads.@threads :greedy for k in reverse(k5Todo)         # one frequency at a time, high (slowest) first
    res = k5Frequency(k)
    lock(k5Lock) do
        jldsave(k5File(k); Kα=res.Kα, Kβ=res.Kβ, uZ=res.uZ, ωI=res.ωI, lmax=res.lmax,
                tailS=res.tailS, tailR=res.tailR, f=k5ωs[k] / 2π)
        n = Threads.atomic_add!(k5Done, 1) + 1
        @printf("  f = %.4f Hz done (%d/%d), lmax %d, tail S %.1e R %.1e, %.0f s elapsed\n",
                k5ωs[k] / 2π, n, length(k5Todo), res.lmax, res.tailS, res.tailR, time() - k5t0)
    end
end

# ---- synthesis of the snapshots and of the trace
println("\n---- synthesis of ", length(k5Times), " snapshots")
Δω = 2π / k5tlen
w = _taperWeights(nω, k5Taper); w′ = copy(w); w′[1] /= 2
ωI = load(k5File(1), "ωI")
Kα_t = zeros(Float32, length(k5Times), nr, nλ); Kβ_t = similar(Kα_t)
uZω = zeros(ComplexF64, nω)
for k in 1:nω
    d = load(k5File(k))
    uZω[k] = d["uZ"]
    for (it, t) in enumerate(k5Times)
        c = ComplexF32((Δω / π) * w′[k] * cis(k5ωs[k] * t) * exp(ωI * t))
        @views Kα_t[it, :, :] .+= real.(c .* d["Kα"])
        @views Kβ_t[it, :, :] .+= real.(c .* d["Kβ"])
    end
end
tTrace = collect(0.0:0.25:k5tlen-0.25)
uZt = [(Δω / π) * real(sum(w′[k] * uZω[k] * cis(k5ωs[k] * t) for k in 1:nω)) * exp(ωI * t) for t in tTrace]
# saved as Float16 / scale (scale = 99.9 % quantile of |K| over all frames, values clipped at
# ±6e4·scale near the source and receiver) so that each file stays ≈ 230 MB; the full
# Float32/complex per-frequency fields stay in the cache directory above.
k5Save(name, K) = begin
    sc = Float32(quantile_abs(K, 0.999))
    file = datadir(KER_CACHE_DIR, "k5_movie_$(k5Tag)_$(name).jld2")
    jldsave(file; K=Float16.(clamp.(K ./ sc, -6f4, 6f4)), scale=sc, times=collect(k5Times),
            r_km=k5Radii.r_km, side=String.(k5Radii.side), lon_deg=collect(k5λ), tTrace=tTrace, uZ=uZt,
            source_depth_km=k5Depth, delta_deg=k5Δ, M=collect(k5M), scheme=k5Scheme, h_km=k5h,
            fmax=k5fmax, taper=k5Taper, tlen=k5tlen,
            units="K × scale: δu_Z (km; M in 1e25 dyn·cm) per unit δln $(name) per km³")
    println("movie data -> ", file, "  (", round(filesize(file) / 1e6; digits=1), " MB)")
end
k5Save("alpha", Kα_t)
k5Save("beta", Kβ_t)

# ---- quick-look figure: trace + one snapshot (the movie itself is rendered from k5Out)
let it = argmin(abs.(collect(k5Times) .- 1190.0))
    fig = Figure(size=(1000, 900))
    ax1 = Axis(fig[1, 1]; xlabel="time (s)", title="u_Z at Δ = $(k5Δ)°, source $(k5Depth) km, T ≥ $(round(Int, 1 / k5fmax)) s")
    sel = (tTrace .>= first(k5Times) - 40) .& (tTrace .<= last(k5Times) + 40)
    lines!(ax1, tTrace[sel], uZt[sel]; color=:black); vlines!(ax1, [k5Times[it]]; color=:red)
    ax2 = Axis(fig[2:3, 1]; aspect=DataAspect(), title="δln α waveform kernel, t = $(k5Times[it]) s")
    θ = deg2rad.(90 .- collect(k5λ))
    x = [r * cos(t) for r in k5Radii.r_km, t in θ]; y = [r * sin(t) for r in k5Radii.r_km, t in θ]
    v = Kα_t[it, :, :]; lim = quantile_abs(v, 0.98)
    scatter!(ax2, vec(x), vec(y); color=vec(v), colormap=:balance, colorrange=(-lim, lim), markersize=2)
    devSave(fig, "k5_quicklook_$(k5Tag)")
end
