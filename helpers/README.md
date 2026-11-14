# helpers/

Claude-authored utility code for the notebooks, kept separate from `src/`
(which Nobu writes and maintains himself). Nothing here replaces or
duplicates `src/` logic -- it only wraps/combines existing `src/` functions
for a specific notebook step. If something here turns out to be broadly
useful, it should graduate into `src/` deliberately, not by accretion.

## Contents

- `sphericalTI_SH.jl` -- 1-D spherical, frequency-domain TI workflow, SH
  (toroidal) Step 1 (notebooks/SphericalTI_SH_oct2026.ipynb): Cartesian TI
  tensor with radial axis, spherical-component divergence, projection onto
  normalised toroidal harmonics, checks, and two `famousEquation` methods
  added from outside src/: `"3DsismoFreqTIradial"` and
  `"1DsphericalSH_TI_Freq"`.

The Kirishima rebuild (notebooks/Elastic2DtopoBenchmark_oct2026.ipynb) still
only calls existing src/ functions directly.
- `sphericalTI_SH_step2.jl` -- Step 2 of the same notebook: planet1D models
  (in-memory DSM-SH zones, PREM), per-zone grids (uniform h per zone by
  default, or planet1D's GT95 `re` grid), scaled local coordinate so one
  Δ=1 recipe serves all zones/h, conventional DSM SH reference matrices,
  and the convergence / zone-by-zone comparison drivers.
- `sphericalTI_SH_step3.jl` -- Step 3: traction-free homogeneous shell,
  exact toroidal eigenfrequencies (j_l, y_l), end-row closures
  (:dsmHalf Galerkin natural row / one-sided traction collocation of order q),
  OPT3 rows for the near-boundary rows wider stencils cannot reach, and the
  generalized eigenproblem A₀v = -Ω²Mv.
- `sphericalTI_SH_step3b.jl` -- Step 3b: time-domain SH synthetics (DSM
  style) for a point force with |m| ≤ 2 normalised toroidal harmonics,
  PREM outer solid region joined at interfaces by adding half-element rows,
  per-scheme global T/H_L/H_N basis matrices, complex-frequency DSM sum over
  l and ω, inverse transform with e^{ω_I t}.
- `sphericalTI_PSV.jl` -- P-SV Step P1 (notebooks/SphericalTI_PSV_oct2026.ipynb):
  spheroidal projection of the TI operator (zonal jet), expected r²-weighted
  radial forms, Cartesian m ≠ 0 check, fluid equation in Q = P/ω
  (Geller & Ohminato 1994), famousEquation entries "1DsphericalPSV_TI_Freq"
  and "1DsphericalFluidQ_Freq".
- `runPoissonMuBenchmark.jl` -- Poisson 1-D/2-D benchmark of the interpolated
  Taylor expansion (manuscript §2.1.5): μ / μᶜ sets (nodes, midpoints, dense,
  ±1/4, ±1/2, outer, Y orders 0/1/2) on field, material or both, against FD_P
  and the one-centre OPT_P; periodic truncation residual (1-D, 2-D, reusing
  scripts/runWavePeriodicFDvsOPT.jl) and 1-D Dirichlet solution error (EGU
  cases). CPU/Float64 recipes. Outputs in scripts/tmp/poisson_mu_benchmark/.
- `optTestFunctionPrototype.jl` -- stand-alone (no package) closed-form 1-D
  natural-weak OPT_P for (κT')' = f with an arbitrary test function W
  (Lagrange interpolants on the stencil; reproduces flexOPT's 1-D residuals).
  Moment conditions on W: P=3 impossible; P=5 W = 4/15 hat₁ + 11/15 hat₂ (O6);
  P=7 W = 40/231 hat₁ + 191/231 hat₂ (O8). Derivation:
  `derivations/optTestFunctionMoments_sympy.py` (sympy).
- `runPoissonMixedTestFunction.jl` -- the same mixtures inside flexOPT
  (hat₂ via trialFunctionRefPoints, mixture assembled from component recipes
  normalised by ΣΓ), periodic residual 1-D (P=3,5,7) and 2-D tensor (P=5).
  Outputs in scripts/tmp/poisson_mixed_test_function/.
- `testTestFunctionMoments.jl` -- checks src/motorsOPT/testFunctionMoments.jl
  (exact moment theory of the OPT test function; `julia helpers/testTestFunctionMoments.jl`).
- `greenlandIceWater.jl` -- dev/greenland_nov2026: ETOPO 2022 surface/bed reading, sub-cell air/ice/water/rock mixing, 3-D free-surface closure on a prepared OPT system (pinned air rows + Dietrich mirror moves, 3-D version of src's 2-D `_pinned_void_rows_2d` / `_apply_dietrich_surface_2d`), landslide force history, station ZNE bases.
