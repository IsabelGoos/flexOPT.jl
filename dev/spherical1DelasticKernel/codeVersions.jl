# dev/spherical1DelasticKernel/codeVersions.jl -- one version string per cached stage
# (bump when the code behind it changes: lib/*.jl here or dev/spherical1Delastic/lib).
CODE_VERSION_K = Dict(
    "k_systems" => "v1-2026-10-28",   # elasticGlobalSystem PREM / layered (dev/spherical1Delastic/lib)
    "k1_disks"  => "v2-2026-10-28",   # v2: :force disk with Kᵀ (adjoint);   # computeSGTDisk for the checks
    "k2_disks"  => "v2-2026-10-28",   # v2: :force disk with Kᵀ (adjoint);   # computeSGTDisk PREM, grid radii
    "k4_disks"  => "v1-2026-10-28",   # k4: PREM disks for the traveltime kernels
    "k3_born"   => "v2-2026-10-28",   # v2: :force disk with Kᵀ (adjoint);   # bornBruteForceCheck
    "k5_freq"   => "v1-2026-11-13",   # k5: per-frequency δln α / δln β kernels on the section
)
KER_CACHE_DIR = joinpath("dev", "spherical1DelasticKernel")
