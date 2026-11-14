# dev/spherical1Delastic/stages/p5_seismograms.jl
# P-SV Step P5: 3-component synthetics (Z = u_r, R = u_θ, T = u_φ) on PREM
# (inner core from r = 0, fluid outer core in Q = P/ω with the G&O coupling,
# mantle + crust), point forces at the north pole, 590 km depth, receivers on
# φ = 0 every 10°: f_z = 1 (P-SV only, m = 0) and f_y = 1 (P-SV m = ±1 + SH),
# flat in ω.  Global systems and spectra cached; traces/plots every time.

prem = planet1D.getSet1Dmodel!(:PREM)
Rprem = prem.averagedPlanetRadiusInKilometer
p5SrcDepth = 590.0
p5Forces = Dict("fz" => (1.0, 0.0, 0.0), "fy" => (0.0, 0.0, 1.0))   # (f_r, f_θ, f_φ) at the pole, φs = 0
p5Recipes = Dict{Any,Any}()

# ---- global systems (P-SV whole planet + SH outer solid), per scheme: cached
p5SystemsParams = Dict{String,Any}(
    "stage" => "p5_systems", "codeVersion" => CODE_VERSION["p5_systems"],
    "configs" => shConfigs, "model" => "PREM", "equationRevisions" => EQUATION_REVISIONS,
    # key => (scheme, uniform h in km)
    "systems" => Dict("DSMconv" => ("DSMconv", 10.0), "DSMconv_ref" => ("DSMconv", 5.0),
                      "OPT3" => ("OPT3", 10.0), "OPT5" => ("OPT5", 10.0)),
)
devRedo!("systems", p5SystemsParams, SPH_CACHE_DIR, "p5_systems")
p5Systems = myProduceOrLoad(p5SystemsParams, SPH_CACHE_DIR, "p5_systems") do cfg
    model = planet1D.getSet1Dmodel!(Symbol(cfg["model"]))
    out = Dict{String,Any}()
    for (key, (scheme, h)) in cfg["systems"]
        gP = elasticRadialGrid(model; h=h, region=:all)
        gS = elasticRadialGrid(model; h=h, region=:outerSolid)
        @time psv = elasticGlobalSystem(:PSV, model, gP; scheme=scheme, configs=cfg["configs"],
                                        recipes=p5Recipes, label=key)
        @time sh = elasticGlobalSystem(:SH, model, gS; scheme=scheme, configs=cfg["configs"],
                                       recipes=p5Recipes, label=key)
        out[key] = (psv=psv, sh=sh)
    end
    out
end
for k in ("DSMconv_ref", "DSMconv", "OPT3", "OPT5")
    s = p5Systems[k]
    @printf("%-12s P-SV: %6d dofs (nnz %7d), Cpl nnz %d | SH: %5d dofs\n", k, size(s.psv.T, 1),
            nnz(s.psv.H0), nnz(s.psv.Cpl), size(s.sh.T, 1))
end

# ---- spectra: cached
p5SpectraParams = Dict{String,Any}(
    "stage" => "p5_spectra", "codeVersion" => CODE_VERSION["p5_spectra"],
    "systemsCache" => basename(devCacheFile(p5SystemsParams, SPH_CACHE_DIR, "p5_systems")),
    "R" => Rprem, "rs_km" => Rprem - p5SrcDepth, "forces" => p5Forces,
    "receiversθdeg" => 10:10:180, "fmax" => 0.05,
    "keys" => ["DSMconv_ref", "DSMconv", "OPT3", "OPT5"],
)
devRedo!("spectra", p5SpectraParams, SPH_CACHE_DIR, "p5_spectra")
p5Spectra = myProduceOrLoad(p5SpectraParams, SPH_CACHE_DIR, "p5_spectra") do cfg
    out = Dict{String,Any}()
    for key in cfg["keys"]
        s = p5Systems[key]
        @time out[key] = elasticSynthetics(s.psv, s.sh, cfg["R"]; rs_km=cfg["rs_km"],
            forces=cfg["forces"], receiversθdeg=cfg["receiversθdeg"], fmax=cfg["fmax"])
    end
    out
end
# sanity: f_z at the pole → no transverse motion; f_y on φ = 0 → u_r, u_θ from P-SV, u_φ from both
for k in p5SpectraParams["keys"]
    sp = p5Spectra[k]["fz"]
    @printf("%-12s f_z: max|u_φ|/max|u_r| = %.2e\n", k, maximum(abs.(sp.uφ)) / maximum(abs.(sp.ur)))
end

# ---- time series, record sections (Z, R, T), misfit to the finer DSM
p5Traces = Dict((k, f) => elasticTimeSeries(p5Spectra[k][f]; np=4096, taper=0.2)
                for k in p5SpectraParams["keys"] for f in keys(p5Forces))
for f in sort(collect(keys(p5Forces)))
    plotRecordSection3C(Dict(k => p5Traces[(k, f)] for k in p5SpectraParams["keys"]),
                        collect(p5SpectraParams["receiversθdeg"]); force=f,
                        normalize=:trace, gain=2.0)   # each station: peak = 2 × 10°
    let ref = p5Traces[("DSMconv_ref", f)]
        for comp in (:Z, :R, :T)
            r = getfield(ref, comp); sel = r.t .<= 2000.0
            norm(r.traces[sel, :]) == 0 && continue
            for k in ("DSMconv", "OPT3", "OPT5")
                tr = getfield(p5Traces[(k, f)], comp)
                e = norm(tr.traces[sel, :] .- r.traces[sel, :]) / norm(r.traces[sel, :])
                @printf("%s %s: %-8s (h = 10 km) vs DSMconv (h = 5 km): relative L2 misfit = %.3e\n",
                        f, comp, k, e)
            end
        end
    end
end
