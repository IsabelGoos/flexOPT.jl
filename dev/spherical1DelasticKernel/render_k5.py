"""Render the k5 waveform-kernel movie (dev/spherical1DelasticKernel/stages/k5_movie.jl output)
in the style of the old DSM-Kernel GMT movies (PKP130.mov): trace with a time cursor on top,
kernel snapshot on the great-circle section below, fixed colour scale, one frame per 0.5 s
shown 1 s.

usage: python3 render_k5.py k5_movie_PKP150_alpha.jld2 out.mp4 [--sat VALUE] [--fps-model 1]
       python3 render_k5.py --synthetic out.mp4      (layout test with a fake field)
"""
import sys, os, subprocess, tempfile, argparse
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.collections import QuadMesh

# colour scale sampled from the old movie's colour bar (red - orange - yellow - white - green - blue)
_CB = [(-5.0, (255, 17, 19)), (-4.2, (255, 75, 19)), (-3.4, (255, 115, 14)), (-2.5, (255, 177, 12)),
       (-1.9, (255, 220, 81)), (-1.2, (255, 241, 197)), (-0.4, (252, 252, 252)), (0.4, (252, 252, 252)),
       (0.8, (165, 241, 206)), (1.4, (56, 213, 149)), (2.0, (17, 179, 136)), (2.6, (18, 153, 203)),
       (3.3, (16, 132, 255)), (4.0, (20, 101, 255)), (5.0, (16, 60, 255))]
CMAP = LinearSegmentedColormap.from_list(
    "dsmk", [((v + 5) / 10, tuple(c / 255 for c in rgb)) for v, rgb in _CB])
CMAP.set_under((166 / 255, 11 / 255, 11 / 255)); CMAP.set_over((13 / 255, 40 / 255, 200 / 255))

R_EARTH = 6371.0
R_CMB, R_ICB = 3480.0, 1221.5


def read_jld2(path):
    import h5py
    f = h5py.File(path, "r")
    g = lambda k: f[k][()]
    K = np.asarray(g("K")).astype(np.float32)            # h5py: reversed (Julia column-major) dims
    K = np.transpose(K, (2, 1, 0)) if K.shape[0] != len(g("times")) else K   # -> (nt, nr, nλ)
    return dict(K=K * float(g("scale")), times=np.asarray(g("times")), r=np.asarray(g("r_km")),
                lon=np.asarray(g("lon_deg")), tt=np.asarray(g("tTrace")), u=np.asarray(g("uZ")),
                delta=float(g("delta_deg")), depth=float(g("source_depth_km")),
                fmax=float(g("fmax")))


def synthetic():
    lon = np.arange(-15, 165.01, 0.5); r = np.concatenate([np.arange(10, 6371, 20.0), [6371.0]])
    times = np.arange(1120, 1240.01, 0.5)
    L, Rr = np.meshgrid(lon, r)
    K = np.array([np.sin(0.2 * L + 0.01 * Rr - 0.3 * t) * np.exp(-((Rr - 2500 - 10 * (t - 1120)) / 400) ** 2)
                  for t in times], np.float32) * 3e-28
    tt = np.arange(0, 2048, 0.25)
    u = 1e-7 * np.exp(-((tt - 1190) / 4) ** 2) * np.sin(tt) - 0.6e-7 * np.exp(-((tt - 1215) / 5) ** 2) * np.cos(tt)
    return dict(K=K, times=times, r=r, lon=lon, tt=tt, u=u, delta=150.0, depth=150.0, fmax=0.1)


def edges(v):
    v = np.asarray(v, float)
    e = np.empty(len(v) + 1)
    e[1:-1] = 0.5 * (v[1:] + v[:-1]); e[0] = v[0] - 0.5 * (v[1] - v[0]); e[-1] = v[-1] + 0.5 * (v[-1] - v[-2])
    return e


def nice_sat(D, user=None):
    if user:
        s = user
    else:
        lon, r, K = D["lon"], D["r"], D["K"]
        # robust: exclude the source and receiver neighbourhoods (singular)
        L, Rr = np.meshgrid(lon, r)
        far = (np.hypot((L - 0) * np.pi / 180 * Rr, R_EARTH - D["depth"] - Rr) > 600) & \
              (np.hypot((L - D["delta"]) * np.pi / 180 * Rr, R_EARTH - Rr) > 600)
        s = np.quantile(np.abs(K[:, far]), 0.995)
    n = int(np.floor(np.log10(s)))
    m = max(1, int(round(s / 10 ** n)))
    return m, n


def frame(D, it, sat, out, lonc=None):
    m, n = sat
    vmax = m * 10.0 ** n
    lon, r = D["lon"], D["r"]
    lonc = 0.5 * (lon[0] + lon[-1]) if lonc is None else lonc
    fig = plt.figure(figsize=(5.95, 8.42), dpi=200)
    fig.patch.set_facecolor("white")
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 6.2})

    # ---- trace with cursor (axes in page points, as in the old 595x842 movies)
    t = D["times"][it]
    ax = fig.add_axes([60 / 595, 1 - 185 / 842, 355 / 595, 145 / 842])
    t1, t2 = D["times"][0] - 20, D["times"][-1] + 20
    sel = (D["tt"] >= t1) & (D["tt"] <= t2)
    A = np.max(np.abs(D["u"][sel])) * 1.15
    ax.plot(D["tt"][sel], D["u"][sel], color="k", lw=0.55)
    ax.axhline(0, color="0.6", lw=0.4)
    ax.axvline(t, color="0.55", lw=0.5)
    ax.set_xlim(t1, t2); ax.set_ylim(-A, A)
    ex = int(np.floor(np.log10(A))); a1 = float(f"{A / 10 ** ex:.0f}") if A / 10 ** ex >= 1 else 1
    yt = a1 * 10 ** ex if a1 * 10 ** ex <= A else 10 ** ex
    ax.set_yticks([-yt, 0, yt]); ax.set_yticklabels([f"−{yt:.0e}".replace("e-0", "e-"), "0", f"{yt:.0e}".replace("e-0", "e-")])
    ax.set_xticks(np.arange(np.ceil(t1 / 10) * 10, t2 + 0.1, 10 if t2 - t1 <= 100 else 20))
    ax.tick_params(direction="in", length=2.5, width=0.5, pad=2)
    for s_ in ax.spines.values(): s_.set_linewidth(0.6)
    ax.set_xlabel(f"Time = {t:.1f} s", labelpad=2)

    # ---- section (fan), whole Earth, centred on lonc at the top
    axf = fig.add_axes([0.02, 0.215, 0.96, 0.52]); axf.set_aspect("equal"); axf.axis("off")
    psi_e = np.deg2rad(edges(lon) - lonc)
    r_s = np.sort(r)
    order = np.argsort(r, kind="stable")
    re = edges(r_s); re[0] = max(re[0], 0.0); re[-1] = R_EARTH
    PSI, RR = np.meshgrid(psi_e, re)
    X, Y = RR * np.sin(PSI), RR * np.cos(PSI)
    Kt = D["K"][it][order, :]
    axf.pcolormesh(X, Y, Kt, cmap=CMAP, vmin=-vmax, vmax=vmax, shading="flat", rasterized=True)
    # outline, core boundaries, side edges
    ps = np.linspace(psi_e[0], psi_e[-1], 600)
    for rr, lw in ((R_EARTH, 0.9), (R_CMB, 0.6), (R_ICB, 0.5)):
        axf.plot(rr * np.sin(ps), rr * np.cos(ps), color="k", lw=lw)
    for p in (psi_e[0], psi_e[-1]):
        axf.plot([0, R_EARTH * np.sin(p)], [0, R_EARTH * np.cos(p)], color="k", lw=0.9)
    # degree ticks on the outer arc
    for L in np.arange(np.ceil(lon[0]), lon[-1] + 1e-6, 1.0):
        p = np.deg2rad(L - lonc); major = (round(L) % 5 == 0)
        d = 90 if major else 45
        axf.plot([R_EARTH * np.sin(p), (R_EARTH + d) * np.sin(p)], [R_EARTH * np.cos(p), (R_EARTH + d) * np.cos(p)],
                 color="k", lw=0.5)
        if major:
            lab = f"{int(round(L)) % 360}°"
            axf.text((R_EARTH + 230) * np.sin(p), (R_EARTH + 230) * np.cos(p), lab, rotation=-np.rad2deg(p),
                     ha="center", va="center", fontsize=5.6)
    # radius ticks along the right edge
    pe = psi_e[-1]
    for rr in np.arange(0, R_EARTH + 1, 250.0):
        major = rr % 1000 == 0
        d = 110 if major else 55
        x0, y0 = rr * np.sin(pe), rr * np.cos(pe)
        tx, ty = np.sin(pe + np.pi / 2), np.cos(pe + np.pi / 2)   # tick direction: outward, below the edge
        axf.plot([x0, x0 + d * tx], [y0, y0 + d * ty], color="k", lw=0.5)
        if major and rr > 0:
            axf.text(x0 + 150 * tx, y0 + 150 * ty, f"{int(rr)}", ha="center", va="top", fontsize=5.6)
    lim = R_EARTH + 450
    axf.set_xlim(-lim, lim); axf.set_ylim(-700, lim)

    # ---- colour bar
    axc = fig.add_axes([150 / 595, 1 - 735 / 842, 255 / 595, 10 / 842])
    sm = plt.cm.ScalarMappable(cmap=CMAP, norm=plt.Normalize(-m, m))
    cb = fig.colorbar(sm, cax=axc, orientation="horizontal", extend="both", extendfrac=0.06)
    cb.set_ticks(np.arange(-m, m + 1, 1 if m <= 6 else 2))
    cb.ax.tick_params(length=2, width=0.5, pad=1.5, labelsize=5.6)
    cb.outline.set_linewidth(0.5)
    cb.set_label("Sensitivity", labelpad=2)
    fig.text(412 / 595, 1 - 728 / 842, rf"$10^{{{n}}}$ /km$^3$", fontsize=5.4, va="center")
    fig.savefig(out, dpi=200, facecolor="white")
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("inp"); ap.add_argument("out", nargs="?")
    ap.add_argument("--sat", type=float, default=None)
    ap.add_argument("--seconds-per-frame", type=float, default=1.0)
    ap.add_argument("--only", type=int, default=None, help="render one frame index to PNG")
    a = ap.parse_args()
    if a.inp == "--synthetic" or a.inp == "synthetic":
        D = synthetic(); out = a.out
    else:
        D = read_jld2(a.inp); out = a.out
    sat = nice_sat(D, a.sat)
    print("colour scale ±", sat[0], "x 1e", sat[1])
    if a.only is not None:
        frame(D, a.only, sat, out); print(out); return
    tmp = tempfile.mkdtemp(dir=os.path.dirname(os.path.abspath(out)))
    for it in range(len(D["times"])):
        frame(D, it, sat, os.path.join(tmp, f"f{it:04d}.png"))
    subprocess.run(["ffmpeg", "-v", "error", "-y", "-framerate", f"{1 / a.seconds_per_frame}", "-i",
                    os.path.join(tmp, "f%04d.png"), "-vf", "scale=trunc(iw/2)*2:trunc(ih/2)*2,format=yuv420p",
                    "-r", "25", "-c:v", "libx264", "-preset", "veryfast", "-crf", "18", out], check=True)
    print(out)


if __name__ == "__main__":
    main()
