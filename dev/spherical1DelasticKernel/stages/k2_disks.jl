# dev/spherical1DelasticKernel/stages/k2_disks.jl
# K2: production-like PREM disks on the dev/spherical1Delastic grid, cached:
#   :moment disk for a source at depth k2Depth  (TSGT side, 20 pairs + 10 displacement pairs)
#   :force  disk for a receiver at the surface  (RSGT side, 10 + 5 pairs)
# Storage = nω × npairs × nr × nθ complex numbers (printed); the accuracy of the
# on-the-fly fields is set by k2Stride (radial subset of the grid), k2θstep and the
# interpolation method (see k1 (d)).  Figures: time snapshots of ε_θφ for f_φ at the
# receiver (≈ ZC2011 Fig. 9, h_Tθφ) and ε_rr for f_r (≈ Fig. 10, h_Zrr).

k2Scheme = "OPT3"; k2h = 10.0
k2Depth = 100.0
k2tlen = 1638.4; k2fmax = 0.0125; k2Oi = 1e-4
k2Stride = 4; k2θstep = 1.0
k2Times = [400.0, 800.0, 1200.0]
prem = planet1D.getSet1Dmodel!(:PREM)
Rp = prem.averagedPlanetRadiusInKilometer
sys = kerSystems("PREM", k2Scheme, k2h)
psv, sh = sys["psv"], sys["sh"]
k2Params(family, rs) = Dict{String,Any}("stage" => "k2_disks", "codeVersion" => CODE_VERSION_K["k2_disks"],
    "family" => String(family), "rs" => rs, "scheme" => k2Scheme, "h" => k2h, "tlen" => k2tlen,
    "fmax" => k2fmax, "omegaiTlen" => k2Oi, "stride" => k2Stride, "θstep" => k2θstep)
k2Disks = Dict{Symbol,Any}()
for (family, rs) in ((:moment, Rp - k2Depth), (:force, Rp))
    params = k2Params(family, rs)
    devRedo!("disks", params, KER_CACHE_DIR, "k2_disks")
    # DrWatson (tagsave) only saves Dicts: wrap the disk NamedTuple
    k2Disks[family] = myProduceOrLoad(params, KER_CACHE_DIR, "k2_disks") do cfg
        Dict("disk" => computeSGTDisk(psv, sh, prem; family=family, rs_km=rs, radii=diskRadii(psv; stride=k2Stride),
                       θdeg=0:k2θstep:180, tlen=k2tlen, fmax=k2fmax, omegaiTlen=k2Oi,
                       p=sourceOrderOf(k2Scheme)))
    end["disk"]
    d = k2Disks[family]
    nb = 16 * (length(d.strain) + length(d.disp))
    @printf("%s disk (rs = %.1f km): %d ω × %d+%d pairs × %d r × %d θ = %.1f MB; lmax %d..%d\n",
            family, rs, length(d.ω), length(d.pairs.strain), length(d.pairs.disp), length(d.r_km),
            length(d.θdeg), nb / 1e6, minimum(d.lmax), maximum(d.lmax))
    t = d.tail[end, :]
    @printf("   truncation tail at the highest frequency: median %.1e, max %.1e (at r = %.0f km)\n",
            sort(t)[(length(t) + 1) ÷ 2], maximum(t), d.r_km[argmax(t)])
end

"Polar snapshot figure of one stored pair (source on the pole at the top)."
function k2SnapshotFigure(d, which, ip, times, title)
    fig = Figure(size=(420 * length(times), 460))
    θ = deg2rad.(d.θdeg)
    x = [r * sin(t) for r in d.r_km, t in θ]; y = [r * cos(t) for r in d.r_km, t in θ]
    for (k, tt) in enumerate(times)
        v = diskSnapshot(d, which, ip, tt)
        lim = maximum(abs, v) / 3
        ax = Axis(fig[1, k]; aspect=DataAspect(), title="t = $(Int(tt)) s", xlabel="km", ylabel="km")
        scatter!(ax, vec(x), vec(y); color=vec(v), colormap=:balance, colorrange=(-lim, lim), markersize=3)
    end
    Label(fig[0, :], title)
    return fig
end
dR = k2Disks[:force]
ipT = findfirst(==((3, 6)), dR.pairs.strain)          # f_φ → ε_θφ
ipZ = findfirst(==((1, 1)), dR.pairs.strain)          # f_r → ε_rr
devSave(k2SnapshotFigure(dR, :strain, ipT, k2Times, "RSGT ε_θφ for f_φ at the receiver (≈ ZC2011 Fig. 9, h_Tθφ)"), "rsgt_T_thetaphi_snapshots")
devSave(k2SnapshotFigure(dR, :strain, ipZ, k2Times, "RSGT ε_rr for f_r at the receiver (≈ ZC2011 Fig. 10, h_Zrr)"), "rsgt_Z_rr_snapshots")
