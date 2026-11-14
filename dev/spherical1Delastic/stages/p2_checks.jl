# dev/spherical1Delastic/stages/p2_checks.jl
# P-SV Steps P2–P4 in one go (operators + fluid + coupling):
#   (a) PREM zone bases: flexOPT (OPT3, OPT5) vs conventional DSM, per zone
#       and per basis matrix (T, H0, Hk, HΛ): scale ≈ 1 and small distance;
#   (b) exact eigenfrequencies of a fluid core + solid mantle (5×5 determinant)
#       with DSMconv, OPT3, OPT5: P-SV recipes, fluid recipe, native clipped
#       half rows, G&O coupling -- all at once.
# Both cached.

psvRecipes = Dict{Any,Any}()   # in-memory (name, equation) => recipe (files cached in data/semiSymbolic)

# ---- (a) zone-by-zone operator distances on PREM
p2ZonesParams = Dict{String,Any}(
    "stage" => "p2_zones", "codeVersion" => CODE_VERSION["p2_zones"],
    "configs" => shConfigs, "model" => "PREM", "h" => 25.0, "schemes" => ["OPT3", "OPT5"],
    "equationRevisions" => EQUATION_REVISIONS,
)
devRedo!("zones", p2ZonesParams, SPH_CACHE_DIR, "p2_zones")
p2Zones = myProduceOrLoad(p2ZonesParams, SPH_CACHE_DIR, "p2_zones") do cfg
    model = planet1D.getSet1Dmodel!(Symbol(cfg["model"]))
    grid = elasticRadialGrid(model; h=cfg["h"], region=:all)
    rows = NamedTuple[]
    for seg in grid.segments
        kind = seg.fluid ? :fluid : :PSV
        Bd = dsmZoneBasis(kind, model, seg)
        for sc in cfg["schemes"]
            config = only(filter(c -> c.name == sc, cfg["configs"]))
            Bf = flexOPTZoneBasis(kind, config, psvRecipes, model, seg; configs=cfg["configs"])
            d = operatorDistance(Bf, Bd, length(seg.ξ), _nfields(kind))
            push!(rows, (zone=seg.zone, kind=kind, rbot=seg.r[1], rtop=seg.r[end],
                         n=length(seg.ξ) - 1, h_km=seg.h_km, scheme=sc, dist=d))
        end
    end
    Dict{String,Any}("rows" => rows)
end
for r in p2Zones["rows"]
    @printf("zone %2d %-5s [%7.1f,%7.1f] n=%4d h=%6.2f km %-5s |", r.zone, string(r.kind),
            r.rbot, r.rtop, r.n, r.h_km, r.scheme)
    for key in (:T, :H0, :Hk, :HΛ)
        d = r.dist[key]
        isnan(d.scale) ? @printf(" %s: — ", key) : @printf(" %s: s=%.4f Δ=%.1e", key, d.scale, d.relDiff)
    end
    println()
end

# ---- (b) fluid core + solid mantle: exact eigenfrequencies
p2EigenParams = Dict{String,Any}(
    "stage" => "p2_coreMantle", "codeVersion" => CODE_VERSION["p2_coreMantle"],
    "configs" => shConfigs, "schemes" => ["DSMconv", "OPT3", "OPT5"], "l" => [2, 10],
    "nIntervals" => [30, 60, 120], "equationRevisions" => EQUATION_REVISIONS,
)
devRedo!("eigen", p2EigenParams, SPH_CACHE_DIR, "p2_coreMantle")
p2Eigen = myProduceOrLoad(p2EigenParams, SPH_CACHE_DIR, "p2_coreMantle") do cfg
    Dict{String,Any}(string(l) => coreMantleEigenBenchmark(cfg["schemes"]; configs=cfg["configs"],
                     l=l, nIntervals=cfg["nIntervals"], recipes=psvRecipes) for l in cfg["l"])
end
for l in p2EigenParams["l"]
    println("\nl = $l")
    for r in p2Eigen[string(l)]
        @printf("%-8s n=%4d  max|Im/Re|=%.1e | %s\n", r.scheme, r.n, r.maxImRatio,
                join((@sprintf("%+.2e", e) for e in r.relErr), " "))
    end
end
