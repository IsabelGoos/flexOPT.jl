# dev/raytheory — ray theory in flexOPT

One Hamiltonian framework for every ray-theoretical arrival, in 1-D spherical
planets and in 2-D/3-D Cartesian media, isotropic or anisotropic. No phase gets
its own equation: P, S, qP, qS1, qS2, reflections and P↔S conversions all come
out of three layers.

| layer | what | where |
|---|---|---|
| 1. Hamilton–Jacobi | characteristics of H(x, p) = 0, integrated in the travel time T | `lib/rayHamiltonian.jl`, `lib/rayTracer.jl` |
| 2. dispersion relation | which branches exist: isotropic v_P, v_S, or Christoffel det[A_ijkl p_j p_l − δ_ik] = 0 (qP, qS1, qS2; group ∦ phase) | `lib/rayHamiltonian.jl`, `lib/rayChristoffel.jl` |
| 3. interface matching | tangential slowness conserved; the outgoing normal slowness solves the outgoing branch's dispersion relation (closed form if isotropic, 6×6 Stroh-type eigenproblem if anisotropic, outgoing = by GROUP velocity) | `interfaceSlowness` in both files |

A **phase** is only an initial branch plus a list of interface events
(`lib/rayPhases.jl`): `P410s` = start P, at the 410 crossed upward transmit
as S. TauP-like names are parsed: `P S PP SS PS SP pP sP sS PcP ScS PcS PKP
SKS PKKP PKiKP PKIKP PmP Pms P410s S410p Pv410P P^410P`, plus named
Cartesian interfaces: `Pv{refl}P`, `P{i1}s`, ... Anisotropic phases are built
explicitly, e.g. `Phase("qPv{moho}qS", :qP, :down, [PhaseEvent(:moho, :reflect, :qS1, :down)])`.

At a discontinuity that is not an event the ray keeps its mode, and is
**totally reflected** when transmission is post-critical (the back-branch of the
upper-mantle triplications, as in TauP). Fluid–solid boundaries are crossed
only by explicit events, so `P` never enters the core and `PKP` does.

## Geometries and the ray parameter

| geometry | coordinates | conserved | units |
|---|---|---|---|
| `Spherical1D` | (r, θ), p = (p_r, p_θ) | p_θ = r sin i / v | s/rad |
| `Cartesian2D` | (x, z), z depth ↓ | p_x = sin i / v (laterally homogeneous only) | s/km |
| `Cartesian3D` | (x, y, z) | (p_x, p_y) idem | s/km |

Media:
* `VerticalModel` (velocity depends on depth only): `readDSMModel(path)`
  reads `dataInput/model/*/*.psv.model` directly (stand-alone reader, TI
  reduced to isotropic Voigt by default), `layeredCartesian(...)` for layers
  with gradients. Discontinuities are named automatically (`:surface`,
  `:moho`, `:cmb`, `:icb`, `:d400`, `:d670`, ...; `:d410` finds the closest,
  i.e. PREM's 400).
* `PiecewiseMedium`: planar, possibly dipping interfaces in 2-D/3-D, regions
  `linearRegion` (v0 + g⋅x), `functionRegion` (any v(x)), or `ARegion`
  (anisotropic 6×6 A = C/ρ: `isotropicA`, `vtiA` from Thomsen parameters,
  `rotateA` for tilted TI).
* gridded 2-D models (`lib/rayGrid.jl`): `smoothGrid` (Gaussian σ) +
  C² cubic B-spline with linear ghost nodes at the edges (`gridEval` gives v
  and ∇v), `gridMedium2D` → a `PiecewiseMedium` under a free surface. Ray
  theory needs the smoothing: σ ≫ wavelength is the validity condition, and
  caustics/multipathing depend on σ. `eikonalOnGrid` (first-order fast
  sweeping) gives independent first-arrival times on the same model.

Two independent solvers:
* **tracer** (`traceRay`, any medium): RK4 on Hamilton's equations, boundary
  crossings located by bisection; `shootSpherical`, `shootCartesian2D`
  (fan + bisection, all branches), `shootCartesian3D` (Newton on takeoff and azimuth).
* **τ–p** (`taupArrivals`, `travelTimeCurve`, laterally homogeneous only):
  Δ(p), T(p) integrals with Gauss–Legendre after a substitution that removes
  the turning-point singularity; adaptive p sampling; same event logic.

## Running

From the repository root (the core needs no package; CairoMakie for the figures):

```sh
julia --project=. dev/raytheory/run.jl check                 # all checks, no plots (~1-2 min)
julia --project=. dev/raytheory/run.jl s1                    # Cartesian figures
julia --project=. dev/raytheory/run.jl s2 PREM ak135 DWThot  # spherical: checks + figures per model
julia --project=. dev/raytheory/run.jl s3                    # anisotropy
julia --project=. dev/raytheory/run.jl s4 s100 s200          # acoustic Marmousi, σ = 100 / 200 m
```

REPL: `include("dev/raytheory/run.jl")`, then `raytheory("s2", "vprem")`.
Figures in `dev/raytheory/outputs/<step>/`.

Using it directly:

```julia
include("dev/raytheory/lib/rayTheory.jl")
m = readDSMModel("dataInput/model/earth/PREM.psv.model")
taupArrivals(m, "P410s", 0.0, 60.0)          # fast 1-D path
a = taupArrivals(m, "SKS", 0.0, 110.0)[1]
r = traceRay(m, parsePhase("SKS"), sphericalStartP(m, 0.0, a.p, :S, :down)...; dT=0.5)   # full ray path
```

## Validation (`check`, 105 checks with the Marmousi file present)

| check | reference | agreement |
|---|---|---|
| gradient medium v0 + gz, 2-D and 3-D | T = (2/g) asinh(gX / 2v0) | 1e-8 s |
| layer over half-space, PmP and PmS (P→S at the Moho) | image source / Snell | 1e-8 s |
| 15° dipping reflector (2-D), oblique dipping plane (3-D, Newton shooting) | image source | 1e-8 s |
| PREM: tracer vs τ–p at the same p (P, PP, PcP, ScS, PKP, PKIKP, PKiKP, SKS, P410s, Pms, pP, sP, S) | two independent formulations | < 1e-5 s, < 1e-6° |
| PREM: τ–p vs TauP (ObsPy) on the same model, 17 phases × 3 depths × 13 distances, 593 TauP arrivals | `reference/taup_premDSM.csv` | all geometric arrivals found, max 1.7 ms (TauP's linear sampling) |
| isotropic tensor through Christoffel | isotropic | 1e-12 s |
| homogeneous VTI qP (ε = 0.25, δ = 0.1), 2-D and 3-D | exact Thomsen phase velocity → group velocity | 1e-9 s |
| qP: VTI half-space → isotropic layer | Fermat (minimisation over the crossing point with group velocities) | 1e-9 s |
| P→S reflection, Christoffel vs isotropic Hamiltonian | isotropic layer 3 | 1e-9 s |
| GRIDDED gradient model through the B-spline (2, 5, 9 km) | closed form | 1e-10 s |
| eikonal (fast sweeping, Richardson 12.5/6.25 m) on the same | closed form | 1.5 ms (first-order scheme) |
| Marmousi σ = 100 m, surface source, offsets 3.5–5.5 km: earliest ray | eikonal first arrival (Richardson 5/2.5 m) | 0.7–1.0 ms (eikonal's own error) |

Marmousi (`s4`): the model is `notebooks/tmp/seismicModelMarmousi.jld2`
(802 × 676, 9.2 km × 3.0 km, 450 m of water at 1.5 km/s). Where the
eikonal is earlier than every ray it is not a ray arrival: at short offsets the
direct wave grazing the free surface in the water, and at 6 km offset with
σ = 100 m energy running along the bottom edge of the model (56 ms earlier;
it disappears with σ = 200 m, where rays and eikonal agree to 0.8 ms out to
6 km). `s4` draws ray fans with eikonal wavefronts for a surface and a buried
source, ray wavefronts (folds = caustics/triplications) and T(x) of every
branch at the surface.

The only TauP arrivals not reproduced are 29 `PS` "arrivals" whose P leg grazes
PREM's 600-km zone boundary: the published PREM coefficients leave a 4e-3 km/s
step there, which TauP treats as a reflector. Steps < 1e-3 km/s are not
discontinuities here (`jumpTol` in `rayModels.jl`).

`reference/`: `writeND.jl` (DSM model → TauP .nd), `make_taup_reference.py`
(ObsPy; only needed to regenerate the CSV), `premDSM.nd`, `taup_premDSM.csv`.

## Not done yet / next

* speed: plain `Vector`s everywhere (a PREM shooting fan of 721 rays takes ~1 s);
  StaticArrays would give ×10–50.
* amplitudes: dynamic (paraxial) ray tracing for geometrical spreading, and
  Zoeppritz / anisotropic transmission coefficients at the events.
* anisotropy in `Spherical1D` (radial TI from the DSM models without Voigt
  reduction) and curved interfaces (GeoPoints manifolds with normals).
* qS1/qS2 singularities: branches are labelled by eigenvalue order, so a ray
  may jump branch where the two shear velocities cross (needs coupling ray theory).
* diffracted and head waves (Pdiff, Pn), 3-D spherical; a higher-order eikonal
  (the first-order one here is only a check).
* the eikonal is the Hamilton–Jacobi PDE of layer 1 on a grid: a natural place
  for flexOPT's own operators (OPT stencils for |∇T|² = s²) later.
