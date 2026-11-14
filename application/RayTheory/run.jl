# application/RayTheory/run.jl -- ray theory (Hamiltonian; isotropic and anisotropic;
# 1-D spherical planets and 2-D/3-D Cartesian media; TauP-like phase names).
#
# Thin entry point to dev/raytheory/run.jl; its core (dev/raytheory/lib) uses only the
# Julia standard library, so this is the quickest module to try.
#   julia --project=. application/RayTheory/run.jl check            # all checks, no plots
#   julia --project=. application/RayTheory/run.jl s2 PREM ak135     # 1-D spherical, several models
#   julia --project=. application/RayTheory/run.jl s3                # anisotropy

global DEV_ARGS = copy(@isdefined(APP_ARGS) ? APP_ARGS : ARGS)
Base.include(Main, joinpath(@__DIR__, "..", "..", "dev", "raytheory", "run.jl"))
