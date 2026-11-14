# dev/spherical1Delastic/stages/step3a2_movingnu.jl
# Ported from notebooks/SphericalTI_SH_oct2026.ipynb, Step 3a″ (cells sh-step3pp-01 … sh-step3pp-03).

# row-level check again, both reference-node choices (cheap, rerun every time)
let model = layeredSHModel(; R=6371.0, radii=[3480.0, 6371.0], ρ=[3.3], Vsv=[4.5], Vsh=[4.5]),
    seg = only(shRadialGrid(model; mode=:uniform, h=(6371.0 - 3480.0) / 100).segments),
    opt3 = only(filter(c -> c.name == "OPT3", shConfigs))
    for rp in (:extended, :as2D), p in (2, 3, 5)
        try
            rec = shRecipe(shClippedBoundaryConfig(opt3; p=p, refPoints=rp))
            for r in compareClippedWithDSMHalf(rec, model, seg, 10, 100.0)
                @printf("clip p=%d %-9s %-6s  s/h=% .6f  ‖row - s·DSMhalf‖/‖s·DSMhalf‖ = %.3e  support=%s\n",
                        p, string(rp), r.side, r.scale_over_h, r.relDiff, string(r.support[[1, end]]))
            end
        catch err
            println("clip p=$p $rp FAILED: ", first(sprint(showerror, err), 400))
        end
    end
end

# l = 10: OPT3 fallback (as before) vs moving ν, with the best closures: cached
optOnlyConfigs = filter(c -> c.name in ("OPT3", "OPT4", "OPT5"), shConfigs)
eigL10ppResults = Dict{String,Any}()
for nearRows in (:fallback, :movingNu)
    params = Dict{String,Any}(
        "stage" => "step3a2_eigen", "codeVersion" => CODE_VERSION["step3a2_eigen"],
        "equationRevision" => SH_EQUATION_REVISION, "configs" => optOnlyConfigs,
        "l" => 10, "nIntervals" => [50, 100, 200, 400], "nearRows" => nearRows,
        "closures" => Any[(:traction, 6), (:clipped, 2), (:clipped, 3), (:clipped, 5)],
    )
    prefix = "step3a2_eigL10pp_$(nearRows)"
    devRedo!("eigen", params, SPH_CACHE_DIR, prefix)
    eigL10ppResults[string(nearRows)] = myProduceOrLoad(params, SPH_CACHE_DIR, prefix) do cfg
        Dict{String,Any}("rows" => shShellEigenBenchmark(cfg["configs"]; l=cfg["l"],
            nIntervals=cfg["nIntervals"], nearRows=cfg["nearRows"], closures=cfg["closures"]))
    end["rows"]
end
eigL10pp_fb = eigL10ppResults["fallback"]
eigL10pp_mv = eigL10ppResults["movingNu"]
plotShellEigen(vcat(eigL10pp_fb, eigL10pp_mv); title="(l=10, 3a″)")
