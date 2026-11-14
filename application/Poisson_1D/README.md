# Poisson_1D — heterogeneous 1-D Poisson equation (template application)

$$\frac{d}{dx}\left(\kappa(x)\,\frac{dT}{dx}\right) = f(x), \qquad x\in[0,10\pi],\quad T \text{ fixed at both ends.}$$

The source $f$ is computed symbolically so that the exact solution is
$T(x)=\cos x$ for every conductivity model $\kappa(x)$. Every run can therefore
report its own error, which makes this the simplest place to see what the OPT
operators buy compared with conventional finite differences.

In `src/` the equation is the famous equation `"1DpoissonHetero"`
(`src/motorsOPT/famousEquations.jl`). The configurations come from
`notebooks/poisson1Dbenchmark_small.ipynb`.

## Run

```sh
julia --project=. application/Poisson_1D/run.jl list
julia --project=. application/Poisson_1D/run.jl forward                         # OPT3, κ = sin x + 2, h = e⁻²
julia --project=. application/Poisson_1D/run.jl forward scheme=convFD3 case=quadratic h=0.1
julia --project=. application/Poisson_1D/run.jl convergence                     # convFD3 vs OPT3, 7 spacings
julia --project=. application/Poisson_1D/run.jl convergence schemes=convFD3,convFD5,OPT3 case=twiceλ logInvH=0,1,2,3
julia --project=. application/Poisson_1D/run.jl forward precision=Float32
```

| option | default | meaning |
|---|---|---|
| `scheme` | `OPT3` | one of `list` (forward) |
| `schemes` | `convFD3,OPT3` | several schemes (convergence) |
| `case` | `sameλ` | κ(x): `homogeneous`, `sameλ` (sin x + 2), `twiceλ` (sin x/2 + 2), `sameλ_shifted_π_3`, `λ_2` (cos²x + 1), `quadratic` (x² + 1) |
| `h` | `e⁻² ≈ 0.135` | requested grid spacing (the real Δx divides 10π) |
| `logInvH` | `0,0.5,…,3` | spacings h = e^(−logInvH) for `convergence` |
| `precision` | `Float64` | `Float32` to try single precision in the linear system |
| `experimental` | `false` | allow schemes whose status is `:experimental` |

## Outputs

`outputs/forward/forward_<case>_<scheme>_Nx<N>.png|csv` — κ, numerical and exact T, and the pointwise error.
`outputs/convergence/convergence_<case>_<schemes>.png|csv` — error vs Δx on log axes; the terminal also prints the observed order between successive spacings.

The error is ‖T − exact‖₂ / Nx, the metric used in the notebook.

## What to expect

`OPT3` should be clearly more accurate than `convFD3` at the same spacing and
converge at a higher observed order. If `OPT3` and `convFD3` give the same
error, the recipe was probably built from a wrong configuration: check
`schemes.jl` against the notebook.

## Use from Julia

```julia
include("application/Poisson_1D/run.jl")       # or, without running anything: FLEXOPT_MODULES=[:flexOPT]; include(".../common/appBootstrap.jl"); include(".../Poisson_1D/app.jl")
r = poisson1DForward(scheme="OPT3", case="λ_2", h=0.05)
r.X, r.T, r.exact, r.error
rows = poisson1DConvergence(schemes=["convFD3", "OPT3"], case="quadratic")
poisson1DOrders(rows, "OPT3")
```

## Files

* `schemes.jl` — schemes for this equation, with `status` (`:validated` / `:experimental`), `version` (part of the recipe cache key) and `from` (where the configuration was benchmarked).
* `app.jl` — model → recipe (`makeOPTsemiSymbolic`, cached) → operators (`numericalOperatorConstruction`) → solve (`prepareLinearSystem`, `timeMarchingSchemeLinear`) → error. Copy this file to start a new application.
* `run.jl` — command line and REPL entry (`poisson1D("convergence", "case=quadratic")`).
