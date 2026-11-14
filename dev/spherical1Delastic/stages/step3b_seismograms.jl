# dev/spherical1Delastic/stages/step3b_seismograms.jl
# Ported from notebooks/SphericalTI_SH_oct2026.ipynb, Step 3b (cells sh-step3b-01 … sh-step3b-05).
# PREM SH synthetics, f_y point force at 590 km: global systems and spectra
# cached; time series, record section and misfits redone every time.

prem = planet1D.getSet1Dmodel!(:PREM)
Rprem = prem.averagedPlanetRadiusInKilometer
gridSeis = shRadialGrid(prem; mode=:uniform, h=10.0)
gridRef  = shRadialGrid(prem; mode=:uniform, h=2.5)
for s in gridSeis.segments
    @printf("zone %2d: %7.1f – %7.1f km, %4d intervals, h = %6.3f km\n", s.zone, s.r[1], s.r[end], length(s.r) - 1, s.h_km)
end
srcDepth = 590.0
rsSeis = Rprem - srcDepth

# ---- global radial systems (one per scheme): cached
seisSystemsParams = Dict{String,Any}(
    "stage" => "step3b_systems", "codeVersion" => CODE_VERSION["step3b_systems"],
    "equationRevision" => SH_EQUATION_REVISION, "configs" => shConfigs, "model" => "PREM",
    # key => (scheme, uniform h in km)
    "systems" => Dict("DSMconv" => ("DSMconv", 10.0), "DSMconv_ref" => ("DSMconv", 2.5),
                      "OPT3" => ("OPT3", 10.0), "OPT5" => ("OPT5", 10.0)),
)
devRedo!("systems", seisSystemsParams, SPH_CACHE_DIR, "step3b_systems")
seisSystems = myProduceOrLoad(seisSystemsParams, SPH_CACHE_DIR, "step3b_systems") do cfg
    model = planet1D.getSet1Dmodel!(Symbol(cfg["model"]))
    out = Dict{String,Any}()
    for (key, (scheme, h)) in cfg["systems"]
        grid = shRadialGrid(model; mode=:uniform, h=h)
        @time out[key] = scheme == "DSMconv" ? shGlobalSystem(model, grid; scheme=scheme) :
            shGlobalSystem(model, grid; scheme=scheme, configs=cfg["configs"])
    end
    out
end
for k in ("DSMconv", "DSMconv_ref", "OPT3", "OPT5")
    s = seisSystems[k]
    @printf("%-12s n = %5d nodes, nnz(T) = %d\n", k, length(s.ξ), nnz(s.T))
end

# ---- spectra: cached (depends on the systems cache via its file name)
seisSpectraParams = Dict{String,Any}(
    "stage" => "step3b_spectra", "codeVersion" => CODE_VERSION["step3b_spectra"],
    "systemsCache" => basename(devCacheFile(seisSystemsParams, SPH_CACHE_DIR, "step3b_systems")),
    "R" => Rprem, "rs_km" => rsSeis, "force" => (0.0, 0.0, 1.0),
    "receiversθdeg" => 10:10:180, "fmax" => 0.05,
    "keys" => ["DSMconv_ref", "DSMconv", "OPT3", "OPT5"],
)
devRedo!("spectra", seisSpectraParams, SPH_CACHE_DIR, "step3b_spectra")
seisSpectra = myProduceOrLoad(seisSpectraParams, SPH_CACHE_DIR, "step3b_spectra") do cfg
    out = Dict{String,Any}()
    for key in cfg["keys"]
        @time out[key] = shSynthetics(seisSystems[key], cfg["R"]; rs_km=cfg["rs_km"],
            force=cfg["force"], receiversθdeg=cfg["receiversθdeg"], fmax=cfg["fmax"])
    end
    out
end
# sanity: on φ = 0 an f_y source gives a purely transverse SH field
for k in seisSpectraParams["keys"]
    sp = seisSpectra[k]
    @printf("%-12s max|uθ|/max|uφ| = %.2e\n", k, maximum(abs.(sp.uθ)) / maximum(abs.(sp.uφ)))
end

# ---- time series, record section, misfit to the fine DSM reference
seisTraces = Dict(k => shTimeSeries(seisSpectra[k]; np=4096, taper=0.2) for k in seisSpectraParams["keys"])
plotRecordSection(seisTraces, seisSpectra["DSMconv_ref"].θdeg)
let ref = seisTraces["DSMconv_ref"], sel = ref.t .<= 2000.0
    for k in ("DSMconv", "OPT3", "OPT5")
        tr = seisTraces[k]
        e = norm(tr.traces[sel, :] .- ref.traces[sel, :]) / norm(ref.traces[sel, :])
        @printf("%-8s (h = 10 km) vs DSMconv (h = 2.5 km): relative L2 misfit = %.3e\n", k, e)
    end
end
