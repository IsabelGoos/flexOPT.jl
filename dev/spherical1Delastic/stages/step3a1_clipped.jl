# dev/spherical1Delastic/stages/step3a1_clipped.jl
# Ported from notebooks/SphericalTI_SH_oct2026.ipynb, Step 3a′ (cells sh-step3p-01 … sh-step3p-03).
# Row-level check (cheap, rerun every time), then eigen benchmarks (cached).

# 1) row-level check: clipped boundary rows vs the conventional DSM half rows
let model = layeredSHModel(; R=6371.0, radii=[3480.0, 6371.0], ρ=[3.3], Vsv=[4.5], Vsh=[4.5]),
    seg = only(shRadialGrid(model; mode=:uniform, h=(6371.0 - 3480.0) / 100).segments),
    opt3 = only(filter(c -> c.name == "OPT3", shConfigs))
    for p in (2, 3, 5)
        try
            rec = shRecipe(shClippedBoundaryConfig(opt3; p=p))
            for r in compareClippedWithDSMHalf(rec, model, seg, 10, 100.0)
                @printf("clip p=%d  %-6s  s/h=% .6f  ‖row - s·DSMhalf‖/‖s·DSMhalf‖ = %.3e  support=%s\n",
                        p, r.side, r.scale_over_h, r.relDiff, string(r.support[[1, end]]))
            end
        catch err
            println("clip p=$p FAILED: ", first(sprint(showerror, err), 400))
        end
    end
end

eigL10pParams = Dict{String,Any}(
    "stage" => "step3a1_eigen", "codeVersion" => CODE_VERSION["step3a1_eigen"],
    "equationRevision" => SH_EQUATION_REVISION, "configs" => shConfigs,
    "l" => 10, "nIntervals" => [50, 100, 200, 400],
    "closures" => Any[:dsmHalf, (:traction, 4), (:traction, 5), (:traction, 6),
                      (:clipped, 2), (:clipped, 3), (:clipped, 5)],
)
devRedo!("eigen", eigL10pParams, SPH_CACHE_DIR, "step3a1_eigL10p")
eigL10p = myProduceOrLoad(eigL10pParams, SPH_CACHE_DIR, "step3a1_eigL10p") do cfg
    Dict{String,Any}("rows" => shShellEigenBenchmark(cfg["configs"]; l=cfg["l"],
        nIntervals=cfg["nIntervals"], closures=cfg["closures"]))
end["rows"]

eigL2pParams = Dict{String,Any}(
    "stage" => "step3a1_eigen", "codeVersion" => CODE_VERSION["step3a1_eigen"],
    "equationRevision" => SH_EQUATION_REVISION, "configs" => shConfigs,
    "l" => 2, "nIntervals" => [50, 100, 200, 400],
    "closures" => Any[:dsmHalf, (:traction, 4), (:traction, 6), (:clipped, 2), (:clipped, 3), (:clipped, 5)],
)
devRedo!("eigen", eigL2pParams, SPH_CACHE_DIR, "step3a1_eigL2p")
eigL2p = myProduceOrLoad(eigL2pParams, SPH_CACHE_DIR, "step3a1_eigL2p") do cfg
    Dict{String,Any}("rows" => shShellEigenBenchmark(cfg["configs"]; l=cfg["l"],
        nIntervals=cfg["nIntervals"], closures=cfg["closures"]))
end["rows"]

plotShellEigen(eigL10p; title="(l=10, 3a′)")
plotShellEigen(eigL2p; title="(l=2, 3a′)")
