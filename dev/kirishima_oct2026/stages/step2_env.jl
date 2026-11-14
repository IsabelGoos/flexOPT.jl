# dev/kirishima_oct2026/stages/step2_env.jl
# Ported from notebooks/Elastic2DtopoBenchmark_oct2026.ipynb, cell step2-bootstrap.

kirishimaLoadEnv2 = !@isdefined(KIRISHIMA_ENV2)  # once per session (see env_step1.jl)
using Metal
kirishimaLoadEnv2 && include(joinpath(flexopt_root, "src", "batchFiles", "batchGPU.jl"))
backend isa KernelAbstractions.CPU && @warn(
    "No GPU backend is available; OPT recipes will be built on CPU")
@show backend
kirishimaLoadEnv2 && include(joinpath(flexopt_root, "src", "elasticWave2D.jl"))
kirishimaLoadEnv2 && include(joinpath(flexopt_root, "src", "flexOPT.jl"))
kirishimaLoadEnv2 && include(joinpath(flexopt_root, "src", "specfemBenchmark.jl"))
using .elasticWave2D, .flexOPT, .specfemBenchmark
using SparseArrays, Base64, JLD2
KIRISHIMA_ENV2 = true
