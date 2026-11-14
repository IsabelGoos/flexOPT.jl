# dev/greenland_nov2026/stages/env.jl
# Load src/ modules once per Julia session (after editing src/: restart Julia;
# re-including src/ modules makes their exported names ambiguous).
# Every line is its own top-level statement on purpose.
greenlandLoadEnv = !@isdefined(GREENLAND_ENV)
greenlandLoadEnv && !@isdefined(commonBatchs) && include(joinpath(flexopt_root, "src", "commonBatchs.jl"))
greenlandLoadEnv && !@isdefined(planet1D) && include(joinpath(flexopt_root, "src", "planet1D.jl"))
greenlandLoadEnv && planet1D.configure_input!()
greenlandLoadEnv && !@isdefined(GeoPoints) && include(joinpath(flexopt_root, "src", "GeoPoints.jl"))
using .commonBatchs, .planet1D, .GeoPoints
using Metal
greenlandLoadEnv && !@isdefined(backend) && include(joinpath(flexopt_root, "src", "batchFiles", "batchGPU.jl"))
greenlandLoadEnv && !@isdefined(flexOPT) && include(joinpath(flexopt_root, "src", "flexOPT.jl"))
using .flexOPT
(greenlandLoadEnv || "--reload" in (@isdefined(DEV_ARGS) ? DEV_ARGS : ARGS)) &&
    include(joinpath(flexopt_root, "helpers", "greenlandIceWater.jl"))
using .GreenlandIceWater
using CairoMakie, Statistics, LinearAlgebra, SparseArrays, StaticArrays, Printf
using DSP, Interpolations, DelimitedFiles, JLD2   # KernelAbstractions comes with batchGPU.jl
CairoMakie.activate!(type="png")
GREENLAND_ENV = true
@show VERSION Threads.nthreads() backend
