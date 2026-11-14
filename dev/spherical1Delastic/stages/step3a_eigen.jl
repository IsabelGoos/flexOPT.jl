# dev/spherical1Delastic/stages/step3a_eigen.jl
# Ported from notebooks/SphericalTI_SH_oct2026.ipynb, Step 3a (cells sh-step3-01 … sh-step3-03).
# Homogeneous shell, free surface + CMB, toroidal eigenfrequencies: cached.

eigL10Params = Dict{String,Any}(
    "stage" => "step3a_eigen", "codeVersion" => CODE_VERSION["step3a_eigen"],
    "equationRevision" => SH_EQUATION_REVISION, "configs" => shConfigs,
    "l" => 10, "nIntervals" => [50, 100, 200, 400],
    "closures" => Any[:dsmHalf, (:traction, 2), (:traction, 3), (:traction, 4)],
)
devRedo!("eigen", eigL10Params, SPH_CACHE_DIR, "step3a_eigL10")
eigL10 = myProduceOrLoad(eigL10Params, SPH_CACHE_DIR, "step3a_eigL10") do cfg
    Dict{String,Any}("rows" => shShellEigenBenchmark(cfg["configs"]; l=cfg["l"],
        nIntervals=cfg["nIntervals"], closures=cfg["closures"]))
end["rows"]

# lower angular order (stronger curvature terms relative to the stiffness)
eigL2Params = Dict{String,Any}(
    "stage" => "step3a_eigen", "codeVersion" => CODE_VERSION["step3a_eigen"],
    "equationRevision" => SH_EQUATION_REVISION, "configs" => shConfigs,
    "l" => 2, "nIntervals" => [50, 100, 200, 400],
    "closures" => Any[:dsmHalf, (:traction, 3), (:traction, 4)],
)
devRedo!("eigen", eigL2Params, SPH_CACHE_DIR, "step3a_eigL2")
eigL2 = myProduceOrLoad(eigL2Params, SPH_CACHE_DIR, "step3a_eigL2") do cfg
    Dict{String,Any}("rows" => shShellEigenBenchmark(cfg["configs"]; l=cfg["l"],
        nIntervals=cfg["nIntervals"], closures=cfg["closures"]))
end["rows"]

plotShellEigen(eigL10; title="(l=10)")
plotShellEigen(eigL2; title="(l=2)")
