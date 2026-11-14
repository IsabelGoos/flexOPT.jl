# dev/spherical1DelasticKernel — strain Green's tensors and Fréchet kernels (1-D spherical)

Re-implementation of the logic of Zhao & Chevrot (2011, GJI 185, 922; "ZC") and of the
old DSM-Kernel `SGTsh` / `SGTpsv` + `KernelMaker`, on top of `dev/spherical1Delastic/lib`
(global SH / P-SV + fluid systems, any scheme: DSMconv, OPT3, OPT5, …; pole sources
`poleSourceCoefficients`; attenuation). Nothing in `src/` is touched.

```
julia --project=. dev/spherical1DelasticKernel/run.jl k1      # checks (PREM, DSMconv + OPT3)
julia --project=. dev/spherical1DelasticKernel/run.jl k2      # PREM disks on the grid + snapshots
julia --project=. dev/spherical1DelasticKernel/run.jl k3      # Born kernels vs brute force
# REPL: include("dev/spherical1DelasticKernel/run.jl"); kernel1D("k3", "--redo=born")
```

## Design

1. **Disk** (`lib/sgtDisk1D.jl`, `computeSGTDisk`). Point source on the north pole at
   radius `rs`; for every frequency, one sparse LU per (ω, l) and system serves the 6
   radial load columns (U·X, U·X′, V·X, V·X′ | W·X, W·X′), hence every elementary source.
   Stored on the meridian φ = 0 over (r, θ):
   * `family = :moment` (Mrr, Mrθ, Mrφ, Mθθ, Mθφ, Mφφ) — source side (TSGT),
   * `family = :force` (fr, fθ, fφ) — receiver side (RSGT, by reciprocity).

   Only the non-zero (source component, field component) **pairs** are stored. φ = 0 is a
   mirror plane, so even sources (no or two φ indices) only give even components there:
   **20 strain pairs** for moments and **10** for forces, i.e. exactly ZC's 20 H_i and
   10 h_i, but in a basis that is generated rather than hand-written (`sgtPairs`).
   Displacements are stored too (10 / 5 pairs, for density kernels).
2. **Radii = the dev/spherical1Delastic grid** (`diskRadii(psv; stride)`): every
   `stride`-th node of each zone. Discontinuities are stored on both sides, so the
   interpolation never crosses them. The node values use the same Fornberg stencils as
   the source loads, which makes discrete reciprocity exact.
3. **On-the-fly fields** (`lib/sgtInterp1D.jl`): `diskInterpolants(disk, iω; method)`
   builds the interpolation in (r, θ), one zone at a time.
   * `:lagrange` (default): local, order 3 in r and 5 in θ.
   * `:cubic` / `:linear`: Interpolations.jl. A cubic spline is global over the zone and
     rings near the singular source field.

   Accuracy is the user's choice (stride, θ step, method); `checkInterpolation` measures it.
4. **Geometry** (`lib/sgtGeometry1D.jl`): explicit 3×3 rotations, with no azimuth tables:
   * Rs rotates S onto the pole;
   * the source tensor or vector is rotated into the pole-local basis;
   * the azimuth φ′ of Q about S is handled by rotating the source by −φ′ about the axis
     (`rotateSourceAboutPole`);
   * the field is brought back to Q's local geographic basis with
     T = E_Q′ Rs′ E(θ′, φ′).

   `sgtAt` (exact radius) and `sgtAtInterp` (any radius) work for moment tensors and for
   receiver forces.
5. **Kernels** (`lib/sgtKernels1D.jl`): frequency-domain waveform kernels
   δu_n(ω) = ∫ [ω_c² δρ u^R·u^S − ε^R : δC : ε^S] dV, written for Love's
   (ρ, A, C, F, L, N) and isotropic (λ, μ), (κ, μ). Traveltime and amplitude kernels
   (ZC eqs 9–12, Paper II) are the next step.
6. **ZC route** (`lib/sgtZC1D.jl`): ZC's own h1–h10 and H1–H20 with the *correct*
   extraction, plus (B10), (B11) and (B14). It is checked against route 4 to 1e-14, so it
   can be used to compare old SGT files function by function.

## Validation (standalone run when the lib was written; DSMconv, layered model with a fluid core)

| check | result |
|---|---|
| per-(l,m) strain vs Cartesian finite differences | 3e-9 |
| reciprocity: TSGT-side displacement vs RSGT-side strain : M, generic S, R | 1e-13 |
| :moment disk vs `pointSourceSynthetics` (independent) | 1e-7 (θ interpolation) |
| ZC route vs rotation route | 6e-15 (receiver), 4e-14 (source) |
| Born ∫K δm dV vs brute force, δvs mantle shell | ≤ 8e-4 |
| same, δρ mantle shell | 1.4e-2 → 4e-3 for h 20 → 10 km (O(h²) of DSMconv) |
| same, δvp whole fluid core (fields u = −∇q/(ρω)) | ≤ 8e-4 |
| same through the stored grid + on-the-fly (lagrange) | ≤ 1e-3 |

Known limitations / TODO:
* At the ends of each zone (free surface, discontinuities) the strain uses one-sided O(h)
  radial derivatives of the discrete solution. Measured: 2e-2 at 3 km depth and 3e-2 at
  13 km above the CMB, against ≤ 1e-3 elsewhere. A traction-continuity closure
  (σ_rr = σ_rθ = 0 at the surface) would fix it.
* Near the source, the strain converges slowly in l. `disk.tail[iω, ir]` reports the
  truncation per radius; `lmaxOf` from the surface criterion is not enough at r ≈ rs,
  small θ.
* Storage: nω × npairs × nr × nθ complex numbers; k2 prints it.

## Audit of the old DSM-Kernel code (2026-10-28)

The old pipeline was emulated line by line and fed with *exact* strains from this lib.
ZC (B10), (B11) and (B14) are **correct as printed**; the bugs are in the code.

**Receiver side (SGTpsv/SGTsh `udertorsgt` + KernelMaker `rsgt2h3freq`)**
1. At the extraction point (R on the pole, Q on φ = 0), ZC's φ_RQ = π (the direction from
   Q to R is north). `udertorsgt` assumed φ_RQ = 0 for h2 and h10, so `rsgt(2)` and
   `rsgt(10)` have the wrong sign. `rsgt(9)` was already sign-"fixed", which is
   consistent with φ_RQ = π.
2. The `! modified` factors 2 and 4 in `rsgt2h3freq` are wrong. Plain (B11) is exact
   (1e-15) once item 1 is fixed.
3. Typo for T: `h3(3) = -h3(3)` must be `-h3(2)`.
4. (B10) defines L and T relative to the **Q–R** great circle (angle φ_QR at R). Kernels
   need R and T relative to the source, i.e. a rotation by ψ = φ_SR − φ_QR at the
   receiver. KernelMaker never computes φ_QR, so its R and T kernels are wrong; with the
   rotation the result is exact (1e-14).

**Source side (`udertotsgt` + `tsgt2h4freq`)**
5. At the extraction point (S on the pole), ZC's angles are φ_QS = 0 and φ_SQ = π. The
   code assumed the opposite, so H5, H6, H7, H15, H16 (single cos or sin of one angle)
   have the wrong sign.
6. H8, H17, H18, H19, H20 are extracted with wrong combinations, not just wrong signs.
   The correct formulas are in `zcTSGT`.
7. `tsgt2h4freq` contracts H_ijpq M_pq with the off-diagonal moments counted once. For
   M_θφ this cannot be absorbed into H, because H9, H12 and H17–H20 also serve the
   diagonal sources. With 2M_rθ, 2M_rφ, 2M_θφ and the corrected H the result is exact
   (4e-14).

**DSM side (SGTsh / SGTpsv)**
8. `caldvecphi0_withplm`: `bvecdt(1,-m)` is never set (the line assigns `bvecdt(1,m)`
   twice). ∂θu_r of the m < 0 terms is missing or garbage, which breaks ε_rθ in
   TSGT/RSGT.
9. P-SV receiver: `caldveczero` has the wrong sign on the φ component for m = 1. The
   back-propagation loads use opposite sign conventions for icomp = 1 (−conj) and 2, 3
   (+conj, "NF touches"), and SH uses −conj with a sign-flipped `calbveczero`. As a
   result, the P-SV contribution to the transverse RSGT (h8, h10) has the wrong sign
   relative to the SH contribution. This is deduced from the code and not run: the new
   code takes forces from `poleSourceCoefficients` and passes the reciprocity check.
10. P-SV, l ≠ 0, field points in the **fluid core**: `ksta = kkdr + 2·iista − 1` and
    `interpolate(2, …)` assume 2 interleaved unknowns per node, but fluid zones have 1
    (`calspdr`). The outer-core TSGT/RSGT read the wrong unknowns.
11. Fluid strain (`calupfluid`) is isotropic only. That is fine for κ, but the full strain
    is −∇∇q/(ρω) (done here).
12. Pole (θ = 0, π): every horizontal harmonic and every derivative is set to zero, and
    `locallyCartesianDerivatives` only fills the diagonal there. The strain directly
    below the source or receiver is wrong.
13. The l-truncation (`calamp`) and `calmdr` depth cut use surface-amplitude criteria,
    which truncate the near-source field of the TSGT at depth.
14. Minor: `ir0` is hard-coded to 1 while the file names use `r0(r0_n)`;
    `if(rsgtswitch)` treats an integer as a logical.
