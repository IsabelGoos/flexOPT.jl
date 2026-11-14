# dev/spherical1Delastic/codeVersions.jl
#
# One version string per CACHED stage; each enters that stage's
# myProduceOrLoad parameters. A parameter change already gives a new cache
# file (new hash); a change of the CODE behind a stage (lib/*.jl,
# src/) is not seen by the hash -> bump the string here (or use --redo=<stage>).
# SH_EQUATION_REVISION (helpers/sphericalTI_SH_step2.jl) also enters every key.
CODE_VERSION = Dict(
    "step1_cartesian"   => "v1-2026-10-21",  # cartesianToroidalCheck (brute force, m ≠ 0)
    "step2_convergence" => "v1-2026-10-21",  # shConvergenceHomogeneous (uniform h and DSM re-grid)
    "step2_exactness"   => "v1-2026-10-21",  # shStencilExactness
    "step2_zones"       => "v1-2026-10-21",  # shZoneComparison (two-zone TI, PREM)
    "step3a_eigen"      => "v1-2026-10-21",  # shShellEigenBenchmark, Step 3a
    "step3a1_eigen"     => "v1-2026-10-21",  # Step 3a′ (clipped closures)
    "step3a2_eigen"     => "v1-2026-10-21",  # Step 3a″ (fallback vs moving ν)
    "step3b_systems"    => "v1-2026-10-21",  # shGlobalSystem, PREM
    "step3b_spectra"    => "v1-2026-10-21",  # shSynthetics
    "p1_cartesian"      => "v1-2026-10-26",  # cartesianSpheroidalCheck
    "p2_zones"          => "v1-2026-10-26",  # flexOPTZoneBasis vs dsmZoneBasis (lib/elastic1D.jl)
    "p2_coreMantle"     => "v1-2026-10-26",  # coreMantleEigenBenchmark (lib/elastic1D_checks.jl)
    "p5_systems"        => "v1-2026-10-26",  # elasticGlobalSystem (P-SV + SH)
    "p5_spectra"        => "v1-2026-10-26",  # elasticSynthetics
    "p6_shell"          => "v1-2026-10-27",  # shellSourceConvergence (lib/sources1D.jl)
    "p6_systems"        => "v1-2026-10-27",  # elasticGlobalSystem + anelasticSystems (P-SV + SH)
    "p6_spectra"        => "v1-2026-10-27",  # pointSourceSynthetics (moment tensor, φ = 90°)
    "p7_rows"           => "v2-2026-10-28",  # shellRowResiduals (lib/nearBoundary1D.jl)
    "p7_eigen"          => "v2-2026-10-28",  # shellEigenNearRows
    "p7_sources"        => "v2-2026-10-28",  # shellSourceConvergence per nearRows variant
    "p7_coreMantle"     => "v2-2026-10-28",  # coreMantleEigenBenchmark per nearRows variant
    "p8_spectra"        => "v1-2026-10-28",  # pointSourceSynthetics (split, threaded), one file per (scheme, h)
)
# myProduceOrLoad keeps these under DrWatson's datadir() = <repo>/data (git-ignored)
SPH_CACHE_DIR = joinpath("dev", "spherical1Delastic")
