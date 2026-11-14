# dev/spherical1Delastic/stages/p6_sources.jl
# Step P6: moment tensors (double couples and beyond) and the classic DSM benchmark.
#   (a) Γ: flexOPT's own right-hand-side stencils for δ and -∂δ (convFD3, OPT3,
#       OPT4, OPT5) next to the Fornberg ones (p = 2, 4);
#   (b) exact check in a homogeneous traction-free SH shell: surface response to
#       δ (force) and δ' (moment) for every scheme × radial source treatment
#       (linear-spline "hat", Lagrange/Fornberg p, each scheme's own Γ);
#   (c) PREM, a general moment tensor at 590 km, receivers 10°-170° on φ = 90°
#       (stations due east of the source, as in the classic runs):
#       DSMconv / OPT3 / OPT5 vs classic DSM (tish + tipsv), elastic and
#       anelastic (complex moduli, Q(r) as in the model);
#   (d) wrap-around: tipsv with ω_I·tlen = ln 10², ln 10⁴, ln 10⁶ against runs
#       with period 2·tlen at the same ω_I (pure wrap-around error).
# Classic DSM spectra: <repo>/tmp/classicDSM_out/{elastic,anelastic,wrap/*}
# (computed outside the repository; tish 1.9.0 needs two fixes, see
# tmp/classicDSM_out/README.md).  (b), the systems and the spectra are cached.

p6Recipes = Dict{Any,Any}()
CLASSIC_DIR = joinpath(flexopt_root, "tmp", "classicDSM_out")
p6Depth = 590.0
p6θ = 10:10:170
# (Mrr, Mrθ, Mrφ, Mθθ, Mθφ, Mφφ) in 1e25 dyn·cm, every |m| ≤ 2 excited
p6Sources = Dict("MT" => pointSource(M=(1.0, 0.6, -0.8, -0.3, 0.4, -0.7)))
p6OmegaiTlen = 1e-4          # ω_I·tlen = ln 10⁴ (the classic runs use the same)

# ---------------------------------------------------------------------------
println("\n---- (a) right-hand-side stencils Γ (node units, source on the node at offset 0)")
p6Gamma = Dict{String,Any}()
for c in filter(c -> c.name in ("convFD3", "OPT3", "OPT4", "OPT5"), shConfigs)
    try
        p6Gamma[c.name] = sourceGammaStencils(c)
    catch err
        @warn "Γ extraction failed for $(c.name)" exception = (err, catch_backtrace())
    end
end
for p in (2, 4)
    st = fornbergSourceStencil(p)
    @printf("%-14s offsets %s\n   g0 = %s\n   g1 = %s\n", st.name, string(st.offsets),
            join((@sprintf("%+.6f", x) for x in st.g0), " "), join((@sprintf("%+.6f", x) for x in st.g1), " "))
end
if haskey(p6Gamma, "convFD3") && haskey(p6Gamma, "OPT3")
    a, b = p6Gamma["convFD3"], p6Gamma["OPT3"]
    if a.offsets == b.offsets
        @printf("FD3 vs OPT3: max|Δg0| = %.2e, max|Δg1| = %.2e\n", maximum(abs.(a.g0 .- b.g0)), maximum(abs.(a.g1 .- b.g1)))
    else
        println("FD3 vs OPT3: different supports ", a.offsets, " / ", b.offsets)
    end
end

# ---------------------------------------------------------------------------
println("\n---- (b) exact δ / δ' check, homogeneous traction-free SH shell (l = 10)")
p6ShellParams = Dict{String,Any}(
    "stage" => "p6_shell", "codeVersion" => CODE_VERSION["p6_shell"],
    "configs" => shConfigs, "equationRevisions" => EQUATION_REVISIONS,
    "schemes" => ["DSMconv", "convFD3", "OPT3", "OPT4", "OPT5"],
    "nIntervals" => [50, 100, 200, 400, 800], "gammaSchemes" => sort(collect(keys(p6Gamma))),
)
devRedo!("shell", p6ShellParams, SPH_CACHE_DIR, "p6_shell")
p6Shell = myProduceOrLoad(p6ShellParams, SPH_CACHE_DIR, "p6_shell") do cfg
    methods = [("hat", :hat, 1), ("lagrange p=2", :lagrange, 2), ("lagrange p=4", :lagrange, 4),
               ("Γ own", :gammaOwn, nothing)]
    shellSourceConvergence(cfg["schemes"]; configs=cfg["configs"], methods=methods,
                           nIntervals=cfg["nIntervals"], recipes=p6Recipes, gammas=p6Gamma)
end
sourceConvergenceOrders(p6Shell)
let fig = Figure(size=(1500, 650))
    for (ic, (fld, ttl)) in enumerate(((:errForce, "δ (point force)"), (:errDprime, "δ' (moment tensor)")))
        ax = Axis(fig[1, ic]; xscale=log10, yscale=log10, xlabel="h (km)", ylabel="relative error",
                  title="SH shell, surface response to $(ttl)")
        styles = Dict("hat" => :dot, "lagrange p=2" => :dash, "lagrange p=4" => :solid, "Γ own" => :dashdot)
        cols = Dict("DSMconv" => :black, "convFD3" => :orange, "OPT3" => :dodgerblue, "OPT4" => :green, "OPT5" => :crimson)
        for key in unique([(r.scheme, r.method) for r in p6Shell])
            rr = sort(filter(r -> (r.scheme, r.method) == key, p6Shell); by=r -> r.h_km)
            lines!(ax, [r.h_km for r in rr], [max(getfield(r, fld), 1e-16) for r in rr];
                   color=cols[key[1]], linestyle=styles[key[2]], label="$(key[1]) / $(key[2])")
        end
        ic == 2 && Legend(fig[1, 3], ax; nbanks=1, framevisible=false)
    end
    devSave(fig, "shell_source_convergence")
end

# ---------------------------------------------------------------------------
println("\n---- (c) PREM, moment tensor at $(p6Depth) km, receivers on φ = 90°, vs classic DSM")
prem = planet1D.getSet1Dmodel!(:PREM)
Rprem = prem.averagedPlanetRadiusInKilometer
p6SystemsParams = Dict{String,Any}(
    "stage" => "p6_systems", "codeVersion" => CODE_VERSION["p6_systems"],
    "configs" => shConfigs, "model" => "PREM", "equationRevisions" => EQUATION_REVISIONS,
    # key => (scheme, uniform h in km)
    "systems" => Dict("DSMconv" => ("DSMconv", 10.0), "DSMconv_ref" => ("DSMconv", 5.0),
                      "OPT3" => ("OPT3", 10.0), "OPT5" => ("OPT5", 10.0)),
    "anelastic" => true,
)
devRedo!("p6systems", p6SystemsParams, SPH_CACHE_DIR, "p6_systems")
p6Systems = myProduceOrLoad(p6SystemsParams, SPH_CACHE_DIR, "p6_systems") do cfg
    model = planet1D.getSet1Dmodel!(Symbol(cfg["model"]))
    out = Dict{String,Any}()
    for (key, (scheme, h)) in cfg["systems"]
        gP = elasticRadialGrid(model; h=h, region=:all)
        gS = elasticRadialGrid(model; h=h, region=:outerSolid)
        kw = (scheme=scheme, configs=cfg["configs"], recipes=p6Recipes, label=key)
        @time psv = elasticGlobalSystem(:PSV, model, gP; kw...)
        @time sh = elasticGlobalSystem(:SH, model, gS; kw...)
        q = nothing
        if cfg["anelastic"]
            @time q = (psv=anelasticSystems(:PSV, model, gP; kw...), sh=anelasticSystems(:SH, model, gS; kw...))
        end
        out[key] = (scheme=scheme, psv=psv, sh=sh, q=q)
    end
    out
end

p6SpectraParams = Dict{String,Any}(
    "stage" => "p6_spectra", "codeVersion" => CODE_VERSION["p6_spectra"],
    "systemsCache" => basename(devCacheFile(p6SystemsParams, SPH_CACHE_DIR, "p6_systems")),
    "R" => Rprem, "rs_km" => Rprem - p6Depth, "sources" => p6Sources,
    "receiversθdeg" => p6θ, "φr" => π / 2, "tlen" => 3276.8, "fmax" => 0.05, "omegaiTlen" => p6OmegaiTlen,
    # (system key, anelastic?)
    "runs" => [("DSMconv_ref", false), ("DSMconv", false), ("OPT3", false), ("OPT5", false),
               ("DSMconv_ref", true), ("DSMconv", true), ("OPT3", true), ("OPT5", true)],
)
devRedo!("p6spectra", p6SpectraParams, SPH_CACHE_DIR, "p6_spectra")
p6Spectra = myProduceOrLoad(p6SpectraParams, SPH_CACHE_DIR, "p6_spectra") do cfg
    out = Dict{String,Any}()
    for (key, an) in cfg["runs"]
        s = p6Systems[key]
        tag = an ? key * "+Q" : key
        println("spectra: ", tag, "  (radial source stencil p = ", sourceOrderOf(s.scheme), ")")
        @time out[tag] = pointSourceSynthetics(s.psv, s.sh, cfg["R"]; rs_km=cfg["rs_km"],
            sources=cfg["sources"], radial=(method=:lagrange, p=sourceOrderOf(s.scheme)),
            anelastic=an ? s.q : nothing, receiversθdeg=cfg["receiversθdeg"], φr=cfg["φr"],
            tlen=cfg["tlen"], fmax=cfg["fmax"], omegaiTlen=cfg["omegaiTlen"])["MT"]
    end
    out
end

p6Classic = Dict{String,Any}()
# "classic DSM": tish/tipsv with the zero frequency computed (dc_*, the reference);
# "classic DSM no ω=0": the codes as distributed (i = 0 skipped → late trend near the source)
for (case, tag) in (("dc_elastic", "classic DSM"), ("dc_anelastic", "classic DSM+Q"),
                    ("elastic", "classic DSM no ω=0"), ("anelastic", "classic DSM no ω=0+Q"))
    dir = joinpath(CLASSIC_DIR, case)
    if isdir(dir)
        p6Classic[tag] = classicSpectra(dir, p6θ)
    else
        println("no classic DSM spectra in ", dir, " -- comparisons with classic DSM skipped")
    end
end

# frequency-domain misfits (all frequencies 0-0.05 Hz, all stations)
for (suffix, ctag) in (("", "classic DSM"), ("+Q", "classic DSM+Q"))
    println("\nspectral misfit ‖A - B‖/‖B‖ (Z, R, T)", suffix == "" ? "  [elastic]" : "  [anelastic]")
    refs = [(ctag, get(p6Classic, ctag, nothing)), ("DSMconv_ref" * suffix, get(p6Spectra, "DSMconv_ref" * suffix, nothing))]
    for (rname, ref) in refs
        isnothing(ref) && continue
        for k in ("DSMconv", "OPT3", "OPT5", "DSMconv_ref", "classic DSM", "classic DSM no ω=0")
            key = k == "classic DSM" ? ctag : k * suffix
            (key == rname || !(haskey(p6Spectra, key) || haskey(p6Classic, key))) && continue
            A = haskey(p6Spectra, key) ? p6Spectra[key] : p6Classic[key]
            m = spectralMisfit(A, ref)
            @printf("  %-16s vs %-16s  Z %.3e  R %.3e  T %.3e\n", key, rname, m.Z, m.R, m.T)
        end
    end
end
if haskey(p6Spectra, "OPT5") && haskey(p6Spectra, "OPT5+Q")
    m = spectralMisfit(p6Spectra["OPT5+Q"], p6Spectra["OPT5"])
    @printf("\nsize of the attenuation effect: ‖OPT5+Q - OPT5‖/‖OPT5‖  Z %.3f  R %.3f  T %.3f\n", m.Z, m.R, m.T)
end

# time domain: record sections and misfits over 0-2000 s
# the cached Dict may also carry non-spectrum entries (e.g. metadata strings): keep spectra only
filter!(kv -> kv[2] isa NamedTuple && hasproperty(kv[2], :ur), p6Spectra)
p6Traces = Dict(k => elasticTimeSeries(v; np=4096, taper=0.2) for (k, v) in merge(p6Spectra, p6Classic))
for (suffix, ctag, ftag) in (("", "classic DSM", "MT_elastic"), ("+Q", "classic DSM+Q", "MT_anelastic"))
    keysPlot = Tuple(filter(k -> haskey(p6Traces, k), [ctag, "classic DSM no ω=0" * suffix,
                                                        "DSMconv" * suffix, "OPT3" * suffix, "OPT5" * suffix]))
    length(keysPlot) < 2 && continue
    for norm_ in (:trace, :global)
        plotRecordSection3C(Dict(k => p6Traces[k] for k in keysPlot), collect(p6θ); force="MT",
                            keys=keysPlot, normalize=norm_, gain=2.0,
                            colors=(:black, :gray60, :orange, :dodgerblue, :crimson),
                            sourceLabel="moment tensor at $(Int(p6Depth)) km" * (suffix == "" ? "" : ", Q"),
                            fileTag=ftag)
    end
    ref = p6Traces[first(keysPlot)]
    for comp in (:Z, :R, :T)
        r = getfield(ref, comp); sel = r.t .<= 2000.0
        nr = norm(r.traces[sel, :]); nr == 0 && continue
        for k in keysPlot[2:end]
            tr = getfield(p6Traces[k], comp)
            @printf("%s %s: %-12s vs %-14s: relative L2 misfit (0-2000 s) = %.3e\n", ftag, comp, k,
                    first(keysPlot), norm(tr.traces[sel, :] .- r.traces[sel, :]) / nr)
        end
    end
end
if haskey(p6Traces, "OPT5") && haskey(p6Traces, "OPT5+Q")
    plotRecordSection3C(Dict(k => p6Traces[k] for k in ("OPT5", "OPT5+Q")), collect(p6θ); force="MT",
                        keys=("OPT5", "OPT5+Q"), normalize=:trace, gain=2.0,
                        sourceLabel="moment tensor, elastic vs Q (OPT5)", fileTag="MT_Q_effect")
end

# ---------------------------------------------------------------------------
"Traces on t = 0:dt:tmax with a cosine taper fixed in frequency (f1 → f2), as shTimeSeries otherwise."
function p6TracesFreqTaper(spec, comp; dt=0.8, tmax=2000.0, f1=0.04, f2=0.05)
    U = getfield(spec, comp); f = spec.ω ./ 2π
    w = [x <= f1 ? 1.0 : x >= f2 ? 0.0 : 0.5 * (1 + cos(π * (x - f1) / (f2 - f1))) for x in f]
    t = collect(0:dt:tmax); Δω = 2π / spec.tlen
    out = zeros(length(t), size(U, 2))
    for ir in 1:size(U, 2), (it, tt) in enumerate(t)
        acc = 0.5 * real(w[1] * U[1, ir])
        for k in 2:length(f)
            acc += real(w[k] * U[k, ir] * cis(spec.ω[k] * tt))
        end
        out[it, ir] = acc * Δω / π * exp(spec.ωI * tt)
    end
    return out
end

println("\n---- (d) wrap-around (tipsv, Z): same ω_I, reference with period 2·tlen")
let wrapDir = joinpath(CLASSIC_DIR, "wrap"), θw = [40, 80, 120, 150, 170]
    if isdir(wrapDir)
        @printf("%-12s %s\n", "ω_I·tlen", join((@sprintf("%8d°", θ) for θ in θw), " "))
        for (short, long, lab) in (("w2", "La", "ln 10^2"), ("w4", "Lb", "ln 10^4"), ("w6", "Lc", "ln 10^6"))
            A = classicSpectra(joinpath(wrapDir, short), θw; parts=("psv",))
            B = classicSpectra(joinpath(wrapDir, long), θw; parts=("psv",))
            # the taper must be the SAME filter in frequency for both periods
            # (a taper defined as a fraction of the samples differs by ~1 %)
            ta = p6TracesFreqTaper(A, :ur); tb = p6TracesFreqTaper(B, :ur)
            e = [maximum(abs.(ta[:, i] .- tb[:, i])) / maximum(abs.(tb[:, i])) for i in eachindex(θw)]
            @printf("%-12s %s\n", lab, join((@sprintf("%9.1e", x) for x in e), " "))
        end
        println("(max |difference| / max |trace| over 0-2000 s; ln 10^2 is the classic DSM default)")
    else
        println("no wrap-around spectra in ", wrapDir)
    end
end
