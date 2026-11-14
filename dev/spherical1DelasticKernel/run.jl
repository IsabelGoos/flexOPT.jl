# dev/spherical1DelasticKernel/run.jl -- strain Green's tensor disks and Fréchet kernels
# (1-D spherical, frequency domain), built on dev/spherical1Delastic/lib.
#
# Logic of Zhao & Chevrot (2011, GJI 185, 922) and of the old DSM-Kernel
# SGTsh/SGTpsv + KernelMaker, re-implemented without hand-written azimuth
# tables: pole disks of the 20 (moment) / 10 (force) non-zero
# (source component, strain component) pairs, explicit 3×3 rotations for any
# source / receiver / scatterer position (see README.md).
#
# Steps:
#   k1   checks: strain formulas vs finite differences; reciprocity of the two
#        disk families for generic S, R; disk vs pointSourceSynthetics;
#        on-the-fly interpolation accuracy (PREM)                (stages/k1_checks.jl)
#   k2   PREM disks on the dev/spherical1Delastic grid (cached), truncation
#        diagnostics, storage, time snapshots (cf. ZC2011 Figs 9-10)  (stages/k2_disks.jl)
#   k3   Born kernels vs brute force (shell perturbations of a layered model:
#        δvs mantle, δρ mantle, δvp fluid core)                 (stages/k3_born.jl)
#   k4   first 3-D kernels: S traveltime (T component) on PREM, section + maps,
#        weights checked against brute-force cross-correlation   (stages/k4_traveltime.jl)
#   k5   time-domain waveform-kernel movie on the great-circle section, whole Earth
#        (PKP at 150°, T ≥ 10 s), per-frequency cache, run with -t auto  (stages/k5_movie.jl)
#
# From a terminal, in the repository root:
#   julia --project=. dev/spherical1DelasticKernel/run.jl k1
#   julia --project=. dev/spherical1DelasticKernel/run.jl k2 --redo=disks
# From a REPL:
#   include("dev/spherical1DelasticKernel/run.jl")    # default step k1
#   kernel1D("k3")
# --redo= stage names: systems, disks, exact, born, all

include(joinpath(@__DIR__, "..", "common", "devBootstrap.jl"))
KER_STAGES = joinpath(@__DIR__, "stages")
kernel1D(args...) = devRun(joinpath(@__DIR__, "run.jl"), args...)

KER_STEP_FILES = Dict("k1" => "k1_checks.jl", "k2" => "k2_disks.jl", "k3" => "k3_born.jl", "k4" => "k4_traveltime.jl", "k5" => "k5_movie.jl")
DEV = parseDevArgs(@isdefined(DEV_ARGS) ? DEV_ARGS : ARGS; steps=["k1", "k2", "k3", "k4", "k5"], defaultStep="k1")
DEV_OUTPUT_DIR[] = joinpath(@__DIR__, "outputs")
println("spherical1DelasticKernel dev run: ", DEV)

include(joinpath(KER_STAGES, "env.jl"))
include(joinpath(@__DIR__, "codeVersions.jl"))

for step in DEV.steps
    DEV_STEP[] = step
    println("\n================ ", step, " ================")
    @time step include(joinpath(KER_STAGES, KER_STEP_FILES[step]))
end
