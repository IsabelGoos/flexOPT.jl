# dev/kirishima_oct2026/stages/env_step1.jl
# Ported from notebooks/Elastic2DtopoBenchmark_oct2026.ipynb, cell step1-bootstrap.

# Load the Step-1 environment once per Julia session (src/ modules cannot be
# re-included safely in the same session: after editing src/, restart Julia).
# Every line is its own top-level statement on purpose (module bindings
# created by include must be visible to the next statement).
kirishimaLoadEnv1 = !@isdefined(KIRISHIMA_ENV1)
kirishimaLoadEnv1 && include(joinpath(flexopt_root, "src", "commonBatchs.jl"))
kirishimaLoadEnv1 && include(joinpath(flexopt_root, "src", "planet1D.jl"))
kirishimaLoadEnv1 && planet1D.configure_input!()
kirishimaLoadEnv1 && include(joinpath(flexopt_root, "src", "GeoPoints.jl"))
using .commonBatchs, .planet1D, .GeoPoints
using CairoMakie, Statistics, LinearAlgebra
CairoMakie.activate!(type="png")
KIRISHIMA_ENV1 = true
@show VERSION Threads.nthreads() Base.active_project()
