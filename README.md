# flexOPT.jl

flexOPT.jl is implemented in Julia. For a long time, it has been a hand-made code but from July 2026, I allowed myself to use Codex-chatGPT to prepare notebooks, to debug and to harmonise the codes.
Actually flexOPT.jl can read various types of Earth and planet 1D models, pictures and propose optimally accurate operators (Fuji & Duretz 2025) for any (linear) partial differential equations with boundary conditions, and compute neutrino oscillation using realistic flux and detector effects.

The notebooks in `notebooks/` provide demos and exploratory workflows.

## Reproducible setup

Install Julia 1.12 with [juliaup](https://github.com/JuliaLang/juliaup), clone this repository, and run this command from the clone:

```bash
julia scripts/install_notebook_kernel.jl
```

The installer:

- activates the cloned repository rather than a user-global Julia environment;
- instantiates the versions recorded in `Manifest.toml`;
- precompiles with one task to avoid IJulia precompile-task failures;
- installs a Jupyter kernel named **flexOPT 8 threads 1.12**;
- embeds the absolute path of that user’s clone in their local kernelspec.

After installation, reload VS Code, open a notebook, and choose **Select Kernel → Jupyter Kernel → flexOPT 8 threads 1.12**.

To choose another thread count, set `FLEXOPT_THREADS` while installing. For example:

```bash
FLEXOPT_THREADS=4 julia scripts/install_notebook_kernel.jl
```

The committed notebooks expect the default 8-thread kernel. A different thread count can still be selected manually.

## Notebooks outside the repository

A notebook may be stored anywhere. Select the installed flexOPT kernel. If it is launched with another kernel, set the clone location before running the notebook bootstrap:

```julia
ENV["FLEXOPT_ROOT"] = "/absolute/path/to/flexOPT"
```

The project files control package versions; thread count and clone path belong to each user’s local kernelspec and therefore are intentionally not stored in `Project.toml`.
