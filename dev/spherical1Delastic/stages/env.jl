# dev/spherical1Delastic/stages/env.jl
# Ported from notebooks/SphericalTI_SH_oct2026.ipynb, cell sh-step1-01 (+ the
# helper includes of cells sh-step2-01, sh-step3-01, sh-step3b-01).
# src/ modules: loaded once per Julia session (after editing src/, restart
# Julia). lib/*.jl (the single source of the code; helpers/ only forwards
# to them): re-read with --reload. One top-level
# statement per line on purpose (bindings created by include must be
# visible to the next statement).
sphLoadEnv = !@isdefined(SPHERICAL_ENV)
sphLoadHelpers = sphLoadEnv || DEV.reload
@show VERSION Threads.nthreads() Base.active_project()
# Metal must be loaded before batchGPU.jl selects the backend.
using Metal
Metal.functional() || error("Metal.jl is loaded, but cannot access the Apple GPU")
sphLoadEnv && include(joinpath(flexopt_root, "src", "batchFiles", "batchGPU.jl"))
sphLoadEnv && include(joinpath(flexopt_root, "src", "commonBatchs.jl"))
sphLoadEnv && include(joinpath(flexopt_root, "src", "planet1D.jl"))
sphLoadEnv && planet1D.configure_input!()
sphLoadEnv && include(joinpath(flexopt_root, "src", "flexOPT.jl"))
using .commonBatchs, .planet1D, .flexOPT
using Symbolics, LinearAlgebra, Random, SparseArrays, Printf
SPH_LIB = joinpath(flexopt_root, "dev", "spherical1Delastic", "lib")
sphLoadHelpers && include(joinpath(SPH_LIB, "sphericalTI_SH.jl"))
sphLoadHelpers && include(joinpath(SPH_LIB, "sphericalTI_SH_step2.jl"))
sphLoadHelpers && include(joinpath(SPH_LIB, "sphericalTI_SH_step3.jl"))
sphLoadHelpers && include(joinpath(SPH_LIB, "sphericalTI_SH_step3b.jl"))
sphLoadHelpers && include(joinpath(SPH_LIB, "sphericalTI_PSV.jl"))
sphLoadHelpers && include(joinpath(SPH_LIB, "elastic1D.jl"))
sphLoadHelpers && include(joinpath(SPH_LIB, "elastic1D_checks.jl"))
sphLoadHelpers && include(joinpath(SPH_LIB, "sources1D.jl"))
sphLoadHelpers && include(joinpath(SPH_LIB, "nearBoundary1D.jl"))
using CairoMakie
CairoMakie.activate!(type="png")
SPHERICAL_ENV = true
