"""Seafloor depth along ECr10 picked from the published depth section (pseudo-traces).
Pick: first sample whose |amplitude| > 0.25 on 3 consecutive samples, then outlier rejection
(> 150 m from a 15-trace running median) and a 3-trace median. Depth resolution ~9 m (pixel),
lateral ~2.5 shot points. Distance assumes 50 m per shot point (JCG survey, Minami et al. 2022).
"""
import numpy as np
from scipy.ndimage import median_filter
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from PIL import Image

d = np.load("ecr10_pseudotraces.npz"); a, z, sp = d["data"], d["depth_m"], d["shotpoint"]
m = np.abs(a) > 0.25
ok = m[:, :-2] & m[:, 1:-1] & m[:, 2:]
raw = z[ok.argmax(1)]
sm = median_filter(raw, 15)
good = np.abs(raw - sm) < 150
idx = np.arange(len(raw))
sf = median_filter(np.interp(idx, idx[good], raw[good]), 3)
# remove narrow downward spikes (picks falling through a weak seafloor for 1-3 traces)
m7 = median_filter(sf, 7)
spike = (sf - m7) > 40
sf[spike] = m7[spike]; good &= ~spike
dist = (sp.max() - sp) * 50.0 / 1e3        # km from the NW end
np.savetxt("ecr10_seafloor.csv", np.c_[sp, dist, sf, good],
           delimiter=",", fmt=["%.1f", "%.3f", "%.1f", "%d"],
           header="shot_point,distance_km_from_NW(50m/SP assumed),seafloor_depth_m,direct_pick(1)/interpolated(0)")
print(f"{good.mean()*100:.1f}% direct picks; depth range {sf.min():.0f}-{sf.max():.0f} m; length {dist.max():.0f} km")

img = np.array(Image.open("/mnt/user-data/uploads/1794203516413_image.png").convert("RGB"))
crop = img[170:725, 170:2312]
fig, ax = plt.subplots(2, 1, figsize=(16, 9), gridspec_kw=dict(height_ratios=[1.3, 1]))
ax[0].imshow(crop, extent=[sp[0], sp[-1], z[-1], z[0]], aspect="auto")
ax[0].plot(sp, sf, "lime", lw=1)
ax[0].plot(sp[~good], raw[~good], "y.", ms=2)
ax[0].set_ylim(2000, 0); ax[0].set_xlim(sp[0], sp[-1])
ax[0].set_title("Seafloor pick (green) on the original figure; yellow = rejected raw picks")
ax[0].set_ylabel("depth (m)")
ax[1].plot(dist, -sf, "k", lw=1)
ax[1].fill_between(dist, -sf, -2000, color="0.85")
ax[1].set_xlim(dist[0], dist[-1]); ax[1].set_ylim(-1500, 0)
ax[1].set_xlabel("distance from NW end (km, 50 m per shot point)"); ax[1].set_ylabel("elevation (m)")
ax[1].set_title("ECr10 bathymetric profile")
sec = ax[1].secondary_xaxis("top", functions=(lambda k: sp.max() - k * 1e3 / 50, lambda s: (sp.max() - s) * 50 / 1e3))
sec.set_xlabel("shot point")
fig.tight_layout(); fig.savefig("ecr10_seafloor.png", dpi=110)
