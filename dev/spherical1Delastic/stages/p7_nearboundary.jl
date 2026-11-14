# dev/spherical1Delastic/stages/p7_nearboundary.jl
# Step P7: near-boundary rows for OPT > 3.
# Every row whose stencil leaves the segment gets the scheme's own off-centred
# recipe (nearBoundaryConfig, nearRows = :clippedFamily): row j from the end has
# its test centre j nodes inside, reference nodes all inside the segment
# (one-sided), full test support except the end row (half element, natural
# boundary condition), Taylor exactness of the interior (degree P-1).
# Compared with the current default (:fallback: OPT3 rows + OPT3 clipped end
# row) and moving ν (:movingNu), on the homogeneous traction-free SH shell:
#   (a) recipes of the family: built once, listed;
#   (b) row-by-row residuals of an exact eigenfunction (where is the order lost?);
#   (c) eigenfrequencies of the first 8 toroidal modes;
#   (d) δ / δ' sources (the p6 check) with each variant;
#   (e) fluid core + solid mantle P-SV eigenfrequencies (coupling + fluid rows).
# Cached: (b)-(e).

p7Recipes = Dict{Any,Any}()
p7Schemes = ["DSMconv", "OPT3", "OPT4", "OPT5"]
p7Variants = [:fallback, :movingNu, :clippedFamily]

println("\n---- (a) near-boundary recipe family")
for c in filter(c -> c.name in ("OPT3", "OPT4", "OPT5"), shConfigs)
    base = getRecipe!(p7Recipes, c, SH_EQUATION)
    offs = recipeOffsets(base)
    J = max(maximum(offs), -minimum(offs), 1)
    for j in 0:J-1
        nb = nearBoundaryConfig(c; j=j)
        try
            r = getRecipe!(p7Recipes, nb, SH_EQUATION)
            @printf("%-6s j=%d  %-28s nuCentre=%s refs=%s bounds=%s degree=%d  stencil offsets %s\n",
                    c.name, j, nb.name, string(nb.extra["nuCentre"]), string(nb.extra["trialFunctionRefPoints"]),
                    string(nb.extra["testIntegrationBounds"]), nb.extra["exactTaylorTotalDegree"],
                    string(recipeOffsets(r)))
        catch err
            @warn "near-boundary recipe failed" c.name j exception = (err, catch_backtrace())
        end
    end
end

# flexOPT's own right-hand-side stencils (Γ) for δ and -∂δ
p7Gamma = Dict{String,Any}()
for c in filter(c -> c.name in ("OPT3", "OPT4", "OPT5"), shConfigs)
    try
        p7Gamma[c.name] = sourceGammaStencils(c; verbose=false)
    catch err
        @warn "Γ extraction failed for $(c.name)" exception = (err, catch_backtrace())
    end
end

p7Params(stage, extra::Pair...) = Dict{String,Any}("stage" => stage, "codeVersion" => CODE_VERSION[stage],
    "configs" => shConfigs, "equationRevisions" => EQUATION_REVISIONS,
    "schemes" => p7Schemes, "variants" => string.(p7Variants), extra...)

println("\n---- (b) row residuals of the exact eigenfunction (mode 4, l = 10), |row·W| / max|Ω²T W|")
p7RowsParams = p7Params("p7_rows", "nIntervals" => [50, 100, 200, 400])
devRedo!("rows", p7RowsParams, SPH_CACHE_DIR, "p7_rows")
p7Rows = myProduceOrLoad(p7RowsParams, SPH_CACHE_DIR, "p7_rows") do cfg
    shellRowResiduals(cfg["schemes"], Symbol.(cfg["variants"]); configs=cfg["configs"],
                      nIntervals=cfg["nIntervals"], recipes=p7Recipes)
end
for i in 1:3
    nearRowOrders(p7Rows, r -> r.top[i]; label="top row n-$(i-1)  ")
    nearRowOrders(p7Rows, r -> r.bottom[i]; label="bottom row $(i)     ")
end
nearRowOrders(p7Rows, r -> r.interior; label="interior max      ")

println("\n---- (c) shell eigenfrequencies (l = 10, 8 modes)")
p7EigenParams = p7Params("p7_eigen", "nIntervals" => [50, 100, 200, 400])
devRedo!("eigen", p7EigenParams, SPH_CACHE_DIR, "p7_eigen")
p7Eigen = myProduceOrLoad(p7EigenParams, SPH_CACHE_DIR, "p7_eigen") do cfg
    shellEigenNearRows(cfg["schemes"], Symbol.(cfg["variants"]); configs=cfg["configs"],
                       nIntervals=cfg["nIntervals"], recipes=p7Recipes)
end
nearRowOrders(p7Eigen, r -> maximum(abs.(r.relErr)); label="max |δΩ/Ω|")
let fig = Figure(size=(900, 650))
    ax = Axis(fig[1, 1]; xscale=log10, yscale=log10, xlabel="h (km)", ylabel="max |δΩ/Ω| (8 modes)",
              title="SH shell eigenfrequencies: near-boundary rows")
    cols = Dict("DSMconv" => :black, "OPT3" => :dodgerblue, "OPT4" => :green, "OPT5" => :crimson)
    sty = Dict(:fallback => :dash, :movingNu => :dot, :clippedFamily => :solid)
    for key in unique([(r.scheme, r.variant) for r in p7Eigen])
        rr = sort(filter(r -> (r.scheme, r.variant) == key, p7Eigen); by=r -> r.h_km)
        lines!(ax, [r.h_km for r in rr], [max(maximum(abs.(r.relErr)), 1e-16) for r in rr];
               color=cols[key[1]], linestyle=sty[key[2]], label="$(key[1]) / $(key[2])")
    end
    axislegend(ax; position=:rb)
    devSave(fig, "shell_eigen_nearrows")
end

println("\n---- (d) δ / δ' sources (p6 check) per near-boundary variant")
p7SrcParams = p7Params("p7_sources", "nIntervals" => [50, 100, 200, 400, 800],
                       "gammaSchemes" => sort(collect(keys(p7Gamma))))
devRedo!("sources", p7SrcParams, SPH_CACHE_DIR, "p7_sources")
p7Src = myProduceOrLoad(p7SrcParams, SPH_CACHE_DIR, "p7_sources") do cfg
    out = NamedTuple[]
    for v in Symbol.(cfg["variants"])
        rr = shellSourceConvergence(cfg["schemes"]; configs=cfg["configs"],
                                    methods=[("lagrange own", :lagrange, nothing), ("Γ own", :gammaOwn, nothing)],
                                    nIntervals=cfg["nIntervals"], recipes=p7Recipes, nearRows=v,
                                    gammas=p7Gamma)
        append!(out, [merge(r, (method="$(r.method) / $(v)",)) for r in rr])
    end
    out
end
sourceConvergenceOrders(p7Src)

println("\n---- (e) fluid core + solid mantle, P-SV eigenfrequencies (l = 2, 3 modes)")
p7CMParams = p7Params("p7_coreMantle", "nIntervals" => [30, 60, 120, 240])
devRedo!("coreMantle", p7CMParams, SPH_CACHE_DIR, "p7_coreMantle")
p7CM = myProduceOrLoad(p7CMParams, SPH_CACHE_DIR, "p7_coreMantle") do cfg
    out = Dict{String,Any}()
    for v in Symbol.(cfg["variants"])
        println("nearRows = ", v)
        out[string(v)] = coreMantleEigenBenchmark(cfg["schemes"]; configs=cfg["configs"],
                                                  nIntervals=cfg["nIntervals"], recipes=p7Recipes,
                                                  nearRows=v)
    end
    out
end
for (v, rows) in sort(collect(filter(kv -> kv[2] isa AbstractVector, p7CM)); by=first)   # skip cache metadata
    for sc in p7Schemes
        rr = sort(filter(r -> r.scheme == sc, rows); by=r -> r.n)
        length(rr) < 2 && continue
        e = [maximum(abs.(r.relErr)) for r in rr]
        @printf("%-14s %-8s max|δΩ/Ω| %s  orders %s\n", v, sc, join((@sprintf("%.2e", x) for x in e), " "),
                join((@sprintf("%5.2f", log2(e[i] / e[i+1])) for i in 1:length(e)-1), " "))
    end
end
