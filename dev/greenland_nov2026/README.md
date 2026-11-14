# dev/greenland_nov2026 — 3-D elastic waves through Greenland (OPT3)

Waveforms at the 16 GLISN stations of Toyokuni et al. (2018, *PEPI* 277,
81–98, Fig. 1 / Table 1) for a single force (ice slide) on the coast near
Nuuk. Claude-authored; `src/` is not modified (new code: this folder and
`helpers/greenlandIceWater.jl`).

```sh
julia --project=. -t auto dev/greenland_nov2026/run.jl step1
julia --project=. -t auto dev/greenland_nov2026/run.jl step2
julia --project=. -t auto dev/greenland_nov2026/run.jl step3
```
step1 = model + review figure, step2 = recipe + operators + propagation,
step3 = traces, record section, snapshots, CSV. Grid case as 2nd argument: `quick` (default, 50 km), `standard` (30 km), `full` (20 km).

## Model (step1)

* **Box**: `GeoPoints.constructLocalBox`, centre 72.23°N 32.42°W, x axis at
  azimuth 170° (≈ south), y ≈ east; half-width set at run time = farthest
  station (or source) + 150 km, at least 1200 km (printed by step1). The box is a true Cartesian frame, so the Earth's curvature
  is in it: the surface is ~230 km below the tangent plane at the corners;
  400 km of Earth is kept below the surface there.
* **Rock**: `getParamsAndTopo` → planet1D PREM (no ocean layer, Moho at 24 km).
* **Ice and ocean**: ETOPO 2022 60″ *ice surface* and *bedrock* (NOAA NCEI),
  downloaded once (NCSS, else OPeNDAP ASCII, else your own files in
  `dataInput/greenland/` with "surface" / "bed" in their names). Ice =
  surface − bed (assumed grounded), ice properties from the paper
  (0.920 g/cm³, 3.81, 1.94 km/s); ocean = sea water (1.03, 1.45, 0).
* **Sub-cell mixing**: a 3 km ice sheet or a 3 km ocean is much thinner
  than a cell, so every cell gets its volume fractions of air / ice /
  water / rock and the Voigt average (ρ, ρVp², ρVs²) of the non-air part.
  A cell with ≥ 50 % air is air.

Figures of step1 (`outputs/step1_<case>/`): `model_*` (maps + section),
`model3d_icesheet_*` (ETOPO ice surface over bed, exaggerated), `model3d_curved_*`
(true curved surface in the box, stations, source, slide arrow), `model3d_fence_*`
(Vs on two vertical slices). step3 adds `snapshots3d_*` (u_z draped on the surface); step4 writes
`outputs/step4_<case>/propagation_<case>.mp4` (surface u_z map + traces with a time cursor).

## Numerics (step2)

* Recipe: `makeOPTsemiSymbolic` of `famousEquations("3DsismoTimeIsoHeteroForce")`,
  OPT3 (3 points in space and time, B-spline order 1, supplementary order 2,
  hierarchical Taylor inverse), scaled coordinates Δ = (1,1,1,1) exactly like
  the Kirishima 2-D benchmark; CPU Float64; cached in `data/semiSymbolic`.
* Operators: `numericalOperatorConstruction` + `prepareLinearSystem`;
  Cerjan layers on the four sides and the bottom.
* **Free surface**: `prepareLinearSystem` applies its pinned-void / Dietrich
  closure only in 2-D. `apply_free_surface_3d` (helpers) does the same in
  3-D on the prepared matrices: air rows u = 0, and for material points
  next to air the air coefficients are moved to the mirror point (sign +1
  normal component, −1 tangential, as `_apply_dietrich_surface_2d`). It is a
  staircase closure of a curved surface; **not yet benchmarked in 3-D**.
* Time marching: the lumped-future predictor–corrector of
  `propagateLinearSystem`, written out so that only the station nodes and
  surface snapshots are stored and the source enters through its 9 force
  columns only.
* Source: force on the Earth = −m a(t), a(τ) ∝ sin(2πτ) sin²(πτ) over
  100 s (zero net impulse), slide azimuth 245°, peak 10¹¹ N (illustrative;
  everything is linear in it).

## Limits (please read)

* **Short periods**: OPT3 needs ~4 points per shear wavelength: ≳ 60 s at
  50 km, ≳ 40 s at 30 km, ≳ 25 s at 20 km. 20 s needs dx ≈ 16 km
  (~1.6 M points, > 60 GB of operators).
* **Long periods**: a ~±1250 km box with Cerjan layers, no gravity and no
  self-gravitation is reasonable up to ~200–300 s. 300–1500 s waves are
  longer than the box (Rayleigh wavelength ~6000 km at 1500 s): use the 1-D
  spherical DSM (`application/Elastic_Spherical1D_Freq`) for that end of the band.
* Memory (rough, A + L + R sparse + one copy during the closure): quick ≈ 3 GB,
  standard ≈ 14 GB, full ≈ 40 GB; step1 prints the estimate for the actual grid.
* Pure-water cells only appear when dx < ocean depth; μ = 0 cells have not
  been tested in the elastic OPT operator (a warning is printed).
