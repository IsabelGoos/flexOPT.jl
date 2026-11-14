# Neutrino_Earth — neutrino oscillations through a planet

Three-flavour oscillations (ν and ν̄) along paths through a 1-D or 2-D planet
whose electron density comes from a model (PREM-like layers or a StagYY
mantle-convection snapshot), folded with atmospheric fluxes (Daemonflux),
charged-current cross sections and a detector model.

Not yet wrapped as a one-command application. Use:

| what | how |
|---|---|
| step-by-step notebook | `notebooks/neutrinoPropagation.ipynb` (copy `notebooks/neutrinoPropagation.config.example.jl` to `neutrinoPropagation_local.jl` and set the paths to your StagYY files) |
| joint seismic + neutrino 2-D Earth (one figure) | `julia --project=. --threads=auto scripts/earth2D_StagYY_neutrino_seismo.jl` (needs `FLEXOPT_STAGYY_DIR`; see the header of the script) |

Modules (load them with `application/common/appBootstrap.jl`):
`:neutrinoOscillation`, `:neutrinoFlux`, `:neutrinoCrossSections`, `:neutrinoDetector`.

Planned tasks for this folder: `forward` (oscillation probabilities and event
counts for a given planet model) and `partials` (sensitivity of event counts to
electron density), with the same command-line pattern as `Poisson_1D/`.
