# dev/raytheory/run.jl -- ray theory in flexOPT (new sub-project, Oct 2026)
#
# One Hamiltonian framework for P, S, qP, qS1, qS2, reflections, conversions:
#   layer 1  Hamiltonian propagation  (lib/rayHamiltonian.jl, lib/rayTracer.jl)
#   layer 2  dispersion relation       (isotropic; Christoffel: lib/rayChristoffel.jl)
#   layer 3  interface matching        (generalized Snell, interfaceSlowness)
#   phases   = event lists on interfaces (lib/rayPhases.jl), never new equations
#   fast 1-D path: τ-p integrals with the same event logic (lib/rayTauP.jl)
# Geometries: Spherical1D (ray parameter r sin i / v, s/rad), Cartesian2D and
# Cartesian3D (sin i / v, s/km, conserved only if laterally homogeneous).
# The core uses only the Julia standard library; src/ is not loaded.
#
# Steps:
#   check  every numerical check (analytic, tracer vs τ-p, TauP, anisotropy), no plots
#   s1     Cartesian isotropic: gradient, layer PmP/PmS, dipping reflector 2-D/3-D  (stages/s1_cartesian.jl)
#   s2     1-D spherical: tracer vs τ-p, PREM vs TauP, travel-time curves, ray paths (stages/s2_spherical.jl)
#          cases: PREM (default) ak135 iasp DWThot TAYAKr vprem
#   s3     anisotropy: Christoffel branches, VTI/TTI rays and wavefronts, anisotropic Snell (stages/s3_anisotropy.jl)
#   s4     acoustic rays in Marmousi (smoothed, C² B-spline) vs an eikonal solver (stages/s4_marmousi.jl)
#          cases: s100 (default, σ = 100 m) s50 s200
#
# From a terminal, in the repository root:
#   julia --project=. dev/raytheory/run.jl check
#   julia --project=. dev/raytheory/run.jl s2 PREM DWThot
#   julia --project=. dev/raytheory/run.jl all
# From a REPL:
#   include("dev/raytheory/run.jl")        # first time (runs the default: check)
#   raytheory("s2", "ak135")               # afterwards; lib edits: raytheory("s1", "--reload")
# Figures: dev/raytheory/outputs/<step>/
# Nothing is cached (everything runs in seconds to a minute).

include(joinpath(@__DIR__, "..", "common", "devBootstrap.jl"))
RT_STAGES = joinpath(@__DIR__, "stages")
raytheory(args...) = devRun(joinpath(@__DIR__, "run.jl"), args...)

RT_STEP_FILES = Dict(
    "check" => "check.jl",
    "s1" => "s1_cartesian.jl",
    "s2" => "s2_spherical.jl",
    "s3" => "s3_anisotropy.jl",
    "s4" => "s4_marmousi.jl",
)
DEV = parseDevArgs(@isdefined(DEV_ARGS) ? DEV_ARGS : ARGS;
    steps=["check", "s1", "s2", "s3", "s4"], defaultStep="check",
    cases=["PREM", "ak135", "iasp", "DWThot", "TAYAKr", "vprem", "s50", "s100", "s200"], defaultCases=String[])
DEV_OUTPUT_DIR[] = joinpath(@__DIR__, "outputs")
println("raytheory dev run: ", DEV)

include(joinpath(RT_STAGES, "env.jl"))
any(s -> s != "check", DEV.steps) && include(joinpath(RT_STAGES, "plots.jl"))

for step in DEV.steps
    DEV_STEP[] = step
    println("\n================ ", step, " ================")
    @time step include(joinpath(RT_STAGES, RT_STEP_FILES[step]))
end
