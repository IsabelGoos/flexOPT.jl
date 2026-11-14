# application/ — try flexOPT one equation at a time

This folder is the **entry point for people who want to use or test flexOPT
without reading `src/`**. There is one sub-folder per equation (or per
"output" such as sound). Each one runs from a single command, writes its
figures and numbers to its own `outputs/` folder, and compiles only what it
needs: if you never touch `Neutrino_Earth/`, nothing of it is loaded.

Nothing in `src/`, `dev/` or `notebooks/` was changed to make this folder.
It only *calls* them.

```
flexOPT/
├── src/            the engine: OPT operators, recipes, solvers, planet1D, GeoPoints …   (developers)
├── dev/            script versions of sub-projects under benchmark (stages, caches)     (developers)
├── notebooks/      exploratory and benchmark notebooks                                   (developers)
└── application/    ← you are here: one folder per equation, one command each           (users, students)
```

---

## 1. Install (once, ~15–30 min the first time)

1. Install **Julia 1.12** with [juliaup](https://github.com/JuliaLang/juliaup):
   ```sh
   curl -fsSL https://install.julialang.org | sh      # macOS / Linux
   juliaup add 1.12 && juliaup default 1.12
   ```
2. Clone the repository and instantiate the environment (versions are fixed by `Manifest.toml`):
   ```sh
   git clone <flexOPT repository URL> flexOPT
   cd flexOPT
   julia --project=. -e 'using Pkg; Pkg.instantiate(); Pkg.precompile()'
   ```
3. (Optional, for notebooks) install the Jupyter kernel as explained in the top-level `README.md`:
   `julia scripts/install_notebook_kernel.jl`.

**Every command below is run from the repository root** (the folder that
contains `Project.toml`). If you run from elsewhere, set
`FLEXOPT_ROOT=/absolute/path/to/flexOPT` in your environment.

---

## 2. Quick start: three commands to check that everything works

| command | what it checks | time* |
|---|---|---|
| `julia --project=. application/RayTheory/run.jl check` | Julia + environment (ray theory uses only the standard library) | < 1 min |
| `julia --project=. application/Poisson_1D/run.jl` | the whole OPT pipeline: recipe → operators → linear solve | 1–3 min |
| `julia --project=. application/Poisson_1D/run.jl convergence` | FD vs OPT accuracy against an exact solution | 3–10 min |

\* first run on a laptop, mostly compilation. Recipes are cached in
`data/semiSymbolic/`, so a second run is much faster.

The Poisson run prints something like
`Nx=… Δx=… error=…` and writes
`application/Poisson_1D/outputs/forward/forward_sameλ_OPT3_Nx….png` and `.csv`.

---

## 3. What is available

Status: **ready** = a self-contained application with its own code here;
**via dev/** = this folder is a thin entry point to a sub-project still being
benchmarked in `dev/` (it works, but its options and outputs may change);
**notebook** = only available as a notebook for now.

| folder | equation / physics | tasks | status |
|---|---|---|---|
| [`Poisson_1D/`](Poisson_1D/) | d/dx(κ dT/dx) = f, heterogeneous κ, exact solution | `forward`, `convergence`, `list` | **ready** (template) |
| [`Elastic_Spherical1D_Freq/`](Elastic_Spherical1D_Freq/) | 1-D spherical Earth, TI + attenuation, frequency domain (DSM-like SH, P-SV, fluid core) | `forward` (seismograms, eigenfrequencies), `partials` (strain Green's tensors, Fréchet kernels) | via `dev/spherical1Delastic`, `dev/spherical1DelasticKernel` |
| [`Elastic_Cartesian2D_Time/`](Elastic_Cartesian2D_Time/) | 2-D/3-D Cartesian elastic waves with topography (Kirishima, NIED model) vs SPECFEM2D | `step1`, `step2`, `step3` | via `dev/kirishima_oct2026` |
| [`RayTheory/`](RayTheory/) | rays (Hamiltonian), 1-D spherical and 2-D/3-D Cartesian, isotropic and anisotropic, TauP-like phase names | `check`, `s1`, `s2`, `s3` | via `dev/raytheory` |
| [`Neutrino_Earth/`](Neutrino_Earth/) | 3-flavour neutrino oscillations through a planet (flux, cross sections, detector) | — | notebook (see its README) |
| [`Sonification/`](Sonification/) | turns any trace into a WAV file (an *output*, not an equation) | default demo, or your own trace | **ready** |

Not yet wrapped here (open the notebook directly):

| topic | notebook |
|---|---|
| 2-D acoustic, OPT vs FD, Gaussian initial fields | `notebooks/homogeneous_OPT_FD_benchmark.ipynb`, `notebooks/gaussianDrift_OPT_debug.ipynb` |
| 2-D elastic benchmarks (homogeneous, heterogeneous, double couple, free surface) | `notebooks/HomogeneousElastic2DBenchmark*.ipynb`, `notebooks/HeterogeneousElastic2DBenchmark.ipynb` |
| Marmousi elastic | `notebooks/marmousi_elastic_OPT_debug.ipynb` |
| nonlinear groundwater (Richards-type) | `notebooks/groundWater.ipynb` |
| regional 2-D/3-D boxes from topography (GeoPoints) | `notebooks/getRegionalBox.ipynb`, `notebooks/propagateInside_run.ipynb` |
| how an OPT recipe is built | `notebooks/optMakerExample.ipynb` |
| small course on inverse theory | `notebooks/smallCourseInverseTheory.ipynb` |

---

## 4. How every application is organised

```
application/<Name>/
├── run.jl        the command line and REPL entry point
├── app.jl        the physics: model → recipe → operators → solve → observe   (ready applications)
├── schemes.jl    the numerical schemes offered for this equation, with their status
├── README.md     what the equation is, options, expected results
└── outputs/      figures, CSV, WAV … (git-ignored, safe to delete)
application/common/appBootstrap.jl   shared loader used by every run.jl
```

### The pipeline (the same for every PDE)

```
 model on a grid        OPT recipe               numerical operators          solve                observe
 (planet1D, GeoPoints,  makeOPTsemiSymbolic  →   numericalOperatorConstruction → prepareLinearSystem → field, error,
  images, analytic)     (cached in data/)                                        timeMarching…         seismograms,
                                                                                                        kernels, WAV
```

`Poisson_1D/app.jl` is the shortest complete example of this chain (about
60 lines of real code) and is the **template to copy** for a new equation.

### Tasks and outputs

The same equation can be used for different jobs, with different loops:

* `forward` — one solution (field, seismograms);
* `convergence` — the same solution for several spacings and schemes, with the observed order;
* `partials` — derivatives of the data with respect to the model (kernels);
* later: `inversion`.

The task is the first word on the command line. Outputs are plain files
(PNG, CSV, JLD2, WAV) so that any other tool can read them.

### Schemes and their status

Each `schemes.jl` lists, for that equation only, the schemes (`convFD3`,
`OPT3`, …) with a `status`:

* `:validated` — benchmarked against an exact or reference solution; used by default;
* `:experimental` — being benchmarked; refused unless you add `experimental=true`.

Each scheme also has a `version`. The version is part of the cache key of
the recipe, so when an operator is improved during a benchmark, bumping the
version is enough: old recipes are not reused by mistake.

`OPT3` is the default everywhere.

---

## 5. Running

### From a terminal

```sh
julia --project=. application/<Name>/run.jl <task> key=value key=value …
```

* Values are numbers, strings, or comma lists: `schemes=convFD3,OPT3`.
* Add `-t auto` after `julia` for multithreading where the application uses it.
* `FLEXOPT_GPU=metal julia --project=. …` loads the Apple GPU backend the way
  the notebooks do. The default is the CPU, which is what recipes are built on.

### From a REPL (fastest when you try many things)

Compilation is the slow part in Julia, so keep one session open per application:

```julia
julia --project=.
julia> include("application/Poisson_1D/run.jl")            # first time: compiles and runs the default
julia> poisson1D("convergence", "case=quadratic")           # same as the command line, now fast
julia> r = poisson1DForward(scheme="convFD5", case="λ_2", h=0.05)
julia> r.error, r.Nx
```

For the applications that go through `dev/`, the REPL functions are those of
the dev sub-project (`spherical(...)`, `kernel1D(...)`, `kirishima(...)`,
see `dev/README.md`).

**After editing `src/`, restart Julia.** Re-including `src/` modules in a live
session makes their exported names ambiguous.

### In a notebook

```julia
ENV["FLEXOPT_ROOT"] = "/absolute/path/to/flexOPT"     # only if the notebook is outside the repository
FLEXOPT_MODULES = [:flexOPT]
include(joinpath(ENV["FLEXOPT_ROOT"], "application", "common", "appBootstrap.jl"))
include(joinpath(flexopt_root, "application", "Poisson_1D", "app.jl"))
r = poisson1DForward(scheme="OPT3", case="sameλ")
```

---

## 6. Loading flexOPT modules yourself

`common/appBootstrap.jl` replaces the usual block of `include`s at the top of
the notebooks. Declare what you need, then include it:

```julia
FLEXOPT_MODULES = [:GeoPoints, :flexOPT]     # dependencies (commonBatchs, planet1D) are added automatically
include("application/common/appBootstrap.jl")
```

| name | file in `src/` | what it gives |
|---|---|---|
| `:flexOPT` | `flexOPT.jl` | OPT recipes, numerical operators, solvers, famous equations |
| `:commonBatchs` | `commonBatchs.jl` | caching (`myProduceOrLoad`), symbolic helpers, interpolation |
| `:planet1D` | `planet1D.jl` | 1-D planet models (PREM, ak135, …), DSM-style grids |
| `:GeoPoints` | `GeoPoints.jl` | geographic points, local boxes, topography (GMT), NIED models, images → models |
| `:fieldOutput` | `fieldOutput.jl` | derivatives and expressions of a computed field |
| `:sonification` | `sonification.jl` | traces → audio (WAV), stereo, waveform videos |
| `:neutrinoOscillation`, `:neutrinoFlux`, `:neutrinoCrossSections`, `:neutrinoDetector` | `neutrino*.jl` | oscillation probabilities, fluxes, cross sections, detector |
| `:seismicStations` | `seismicStations.jl` | NIED / FDSN station lists in a local box |
| `:elasticWave2D`, `:elasticWave3D`, `:elasticGreens2D` | `elastic*.jl` | reference FD propagators and analytical 2-D Green's functions |
| `:specfemBenchmark` | `specfemBenchmark.jl` | prepare / run / read SPECFEM2D comparisons |

Each module is loaded at most once per Julia session.

---

## 7. Where things are written

| what | where | safe to delete? |
|---|---|---|
| figures, CSV, WAV of an application | `application/<Name>/outputs/` | yes |
| OPT recipes (expensive symbolic step) | `data/semiSymbolic/` | yes, they are rebuilt (slowly) |
| dev sub-project caches | `data/dev/<sub_project>/` | yes |

All of these are git-ignored.

---

## 8. Troubleshooting

| symptom | fix |
|---|---|
| `Cannot locate flexOPT` | run from the repository root, or set `FLEXOPT_ROOT` |
| `Package X not found` | `julia --project=. -e 'using Pkg; Pkg.instantiate()'` |
| `WAV` not found (Sonification) | `julia --project=. -e 'using Pkg; Pkg.add("WAV")'` |
| `scheme '…' is experimental` | add `experimental=true`, knowing that it is still being benchmarked |
| names "ambiguous" or "already defined" after editing `src/` | restart Julia |
| Metal errors on a Mac | run without `FLEXOPT_GPU=metal` (CPU default) |
| first run very slow | normal: compilation + recipe construction; the second run reuses both |
| a figure is not displayed | figures are saved as PNG in `outputs/`; open them from there |

---

## 9. Adding a new application

1. Copy `Poisson_1D/` to `application/<Physics>_<Geometry><Dim>_<Domain>/`
   (e.g. `Acoustic_Cartesian2D_Time`, `Elastic_Spherical1D_Freq`).
2. In `schemes.jl`, list the schemes you have benchmarked for this equation, with `status` and `from`.
3. In `app.jl`, keep the five steps: model → recipe → operators → solve → observe.
   The equation must exist as a `famousEquation` in `src/motorsOPT/famousEquations.jl`
   (or be added from outside `src/`, as `helpers/sphericalTI_*.jl` do).
4. In `run.jl`, set `FLEXOPT_MODULES`, include `../common/appBootstrap.jl`, and dispatch on the task word.
5. Add a row to the table in section 3 and a short `README.md` with the expected result,
   so that anybody can tell whether their run worked.

A sub-project that is still moving can stay in `dev/` and be exposed here
with a three-line `run.jl` (see `RayTheory/run.jl`) until it is stable.
