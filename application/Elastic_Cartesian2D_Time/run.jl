# application/Elastic_Cartesian2D_Time/run.jl -- 2-D (and 3-D) Cartesian elastic waves, time domain,
# with real topography and NIED velocity models (Kirishima), compared with SPECFEM2D.
#
# Thin entry point to dev/kirishima_oct2026/run.jl (the code is still being benchmarked there):
#   julia --project=. application/Elastic_Cartesian2D_Time/run.jl step1          # box, model, geometry
#   julia --project=. application/Elastic_Cartesian2D_Time/run.jl step2 heter    # FD3 / OPT3 / OPT5 vs SPECFEM2D
#   julia --project=. application/Elastic_Cartesian2D_Time/run.jl step3 homog    # flat free surface
# Options as in dev/README.md (--redo=opt, --full, --reload).
# SPECFEM2D must be installed for the reference traces (see src/specfemBenchmark.jl).

global DEV_ARGS = copy(@isdefined(APP_ARGS) ? APP_ARGS : ARGS)
Base.include(Main, joinpath(@__DIR__, "..", "..", "dev", "kirishima_oct2026", "run.jl"))
