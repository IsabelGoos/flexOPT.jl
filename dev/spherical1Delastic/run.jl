# dev/spherical1Delastic/run.jl -- 1-D spherical TI SH (DSM-style), script version
#
# Same code as notebooks/SphericalTI_SH_oct2026.ipynb, split into stages/*.jl;
# the expensive parts are cached with myProduceOrLoad (explicitly, in each
# stage file). Steps:
#   step1    symbolic derivation and checks           (stages/step1_symbolic.jl)
#   step2    recipes, convergence, zone comparisons   (stages/step2_operators.jl)
#   step3a   shell eigenfrequencies                   (stages/step3a_eigen.jl)
#   step3a1  3a′: clipped closures                    (stages/step3a1_clipped.jl)
#   step3a2  3a″: fallback vs moving ν                (stages/step3a2_movingnu.jl)
#   step3b   PREM SH seismograms                      (stages/step3b_seismograms.jl)
#   p1       P-SV + fluid (Q = P/ω) symbolic checks     (stages/p1_symbolic.jl)
#   p2       P-SV/fluid operators vs DSM on PREM zones, and exact
#            fluid-core + mantle eigenfrequencies       (stages/p2_checks.jl)
#   p5       PREM 3-component seismograms, f_z and f_y  (stages/p5_seismograms.jl)
#   p6       moment tensors: Γ stencils, exact δ/δ' check, PREM vs classic DSM
#            (tish + tipsv), elastic and anelastic, wrap-around  (stages/p6_sources.jl)
#   p7       near-boundary rows for OPT > 3 (own off-centred recipe family vs
#            OPT3 fallback vs moving ν): row residuals, shell eigenfrequencies,
#            δ/δ' sources, core-mantle P-SV  (stages/p7_nearboundary.jl)
#   p8       convergence curves P-SV and SH on PREM vs very fine classic DSM
#            (zero frequency included)  (stages/p8_convergence.jl)
#            -> run with threads: julia --project=. -t auto dev/spherical1Delastic/run.jl p8
#   p8t      CPU time of our schemes, not cached (run with -t 1 for the one-core
#            comparison with tmp/classic_timing)  (stages/p8t_timing.jl)
#
# From a terminal, in the repository root:
#   julia --project=. dev/spherical1Delastic/run.jl step3b
#   julia --project=. dev/spherical1Delastic/run.jl step3b --redo=spectra
#   julia --project=. dev/spherical1Delastic/run.jl all
# From a REPL (fastest):
#   include("dev/spherical1Delastic/run.jl")      # first time (runs the default: step3b)
#   spherical("step3a2", "--redo=eigen")          # afterwards
# Figures: dev/spherical1Delastic/outputs/<step>/
# --redo= stage names: cartesian, convergence, exactness, zones, eigen, systems, spectra,
#          gamma, shell, p6systems, p6spectra, rows, sources, coreMantle, p8spectra, all
# Code: lib/*.jl (helpers/sphericalTI_*.jl only forward there).

include(joinpath(@__DIR__, "..", "common", "devBootstrap.jl"))
SPH_STAGES = joinpath(@__DIR__, "stages")
spherical(args...) = devRun(joinpath(@__DIR__, "run.jl"), args...)

SPH_STEP_FILES = Dict(
    "step1"   => "step1_symbolic.jl",
    "step2"   => "step2_operators.jl",
    "step3a"  => "step3a_eigen.jl",
    "step3a1" => "step3a1_clipped.jl",
    "step3a2" => "step3a2_movingnu.jl",
    "step3b"  => "step3b_seismograms.jl",
    "p1"      => "p1_symbolic.jl",
    "p2"      => "p2_checks.jl",
    "p5"      => "p5_seismograms.jl",
    "p6"      => "p6_sources.jl",
    "p7"      => "p7_nearboundary.jl",
    "p8"      => "p8_convergence.jl",
    "p8t"     => "p8t_timing.jl",
)
DEV = parseDevArgs(@isdefined(DEV_ARGS) ? DEV_ARGS : ARGS;
    steps=["step1", "step2", "step3a", "step3a1", "step3a2", "step3b", "p1", "p2", "p5", "p6", "p7", "p8", "p8t"], defaultStep="p8")
DEV_OUTPUT_DIR[] = joinpath(@__DIR__, "outputs")
println("spherical1Delastic dev run: ", DEV)

include(joinpath(SPH_STAGES, "env.jl"))
include(joinpath(@__DIR__, "codeVersions.jl"))
include(joinpath(SPH_STAGES, "configs.jl"))
include(joinpath(SPH_STAGES, "plots.jl"))

for step in DEV.steps
    DEV_STEP[] = step
    println("\n================ ", step, " ================")
    @time step include(joinpath(SPH_STAGES, SPH_STEP_FILES[step]))
end
