# dev/spherical1Delastic/stages/p8_convergence.jl
# Step P8: convergence of PREM synthetics, P-SV and SH separately.
#   source: moment tensor at 500 km (middle of the 400-600 km zone: OPT5's
#           5-point source stencil Γ must not cross a discontinuity -- a
#           one-sided Γ at interfaces is still to be done), receivers 10°-170°
#           on φ = 90°, elastic, tlen = 3276.8 s, f ≤ 0.05 Hz, ω_I·tlen = ln 10⁴;
#   ours:   DSMconv (Lagrange p = 2 source), OPT3 and OPT5 with the near-boundary
#           family (nearRows = :clippedFamily) and flexOPT's own Γ for the source,
#           h = 20, 10, 5 km (+ OPT5 at 2.5 km), lmax = ωR/2.5 + 250;
#   reference: classic DSM with the zero frequency (tish fixed, tipsv), ratl = 1e-12,
#           re = 1e-5; classic re = 1e-2 … 1e-4 give its own convergence curve.
# Error = ‖A - B‖/‖B‖ over all frequencies, stations and components (Z, R, T).
# One cache file per (scheme, h): partial runs are kept.
# Run with threads:  julia --project=. -t auto dev/spherical1Delastic/run.jl p8

p8Recipes = Dict{Any,Any}()
p8Depth = 500.0
p8θ = 10:10:170
p8Sources = Dict("MT" => pointSource(M=(1.0, 0.6, -0.8, -0.3, 0.4, -0.7)))
p8ClassicDir = joinpath(flexopt_root, "tmp", "classicDSM_out", "conv500")
p8Runs = [("DSMconv", 20.0), ("DSMconv", 10.0), ("DSMconv", 5.0),
          ("OPT3", 20.0), ("OPT3", 10.0), ("OPT3", 5.0),
          ("OPT5", 20.0), ("OPT5", 10.0), ("OPT5", 5.0), ("OPT5", 2.5)]
p8Lmax(R) = ω -> ceil(Int, ω * R / 2.5) + 250
println("Julia threads: ", Threads.nthreads(), Threads.nthreads() == 1 ? "  (start with -t auto to use all cores)" : "")

println("\n---- Γ stencils (source term of the weak OPT schemes)")
p8Gamma = Dict{String,Any}()
for c in filter(c -> c.name in ("OPT3", "OPT5"), shConfigs)
    p8Gamma[c.name] = sourceGammaStencils(c)
end

prem = planet1D.getSet1Dmodel!(:PREM)
Rprem = prem.averagedPlanetRadiusInKilometer

"Spectra (split P-SV / SH) of one (scheme, h), cached."
function p8Spectra(scheme, h)
    params = Dict{String,Any}("stage" => "p8_spectra", "codeVersion" => CODE_VERSION["p8_spectra"],
        "configs" => shConfigs, "equationRevisions" => EQUATION_REVISIONS, "model" => "PREM",
        "scheme" => scheme, "h" => h, "nearRows" => scheme == "DSMconv" ? "fallback" : "clippedFamily",
        "source" => scheme == "DSMconv" ? "lagrange2" : "gammaOwn", "depth" => p8Depth,
        "sources" => p8Sources, "receiversθdeg" => p8θ, "φr" => π / 2, "tlen" => 3276.8,
        "fmax" => 0.05, "omegaiTlen" => 1e-4, "lmax" => "ωR/2.5+250")
    devRedo!("p8spectra", params, SPH_CACHE_DIR, "p8_spectra")
    res = myProduceOrLoad(params, SPH_CACHE_DIR, "p8_spectra") do cfg
        model = planet1D.getSet1Dmodel!(Symbol(cfg["model"]))
        nr = Symbol(cfg["nearRows"])
        kw = (scheme=scheme, configs=cfg["configs"], recipes=p8Recipes, label="$(scheme) h=$(h)", nearRows=nr)
        @time psv = elasticGlobalSystem(:PSV, model, elasticRadialGrid(model; h=h, region=:all); kw...)
        @time sh = elasticGlobalSystem(:SH, model, elasticRadialGrid(model; h=h, region=:outerSolid); kw...)
        radial = scheme == "DSMconv" ? (method=:lagrange, p=2) :
                 (method=:gamma, p=sourceOrderOf(scheme), gamma=p8Gamma[scheme])
        @time sp = pointSourceSynthetics(psv, sh, Rprem; rs_km=Rprem - cfg["depth"], sources=cfg["sources"],
            radial=radial, receiversθdeg=cfg["receiversθdeg"], φr=cfg["φr"], tlen=cfg["tlen"],
            fmax=cfg["fmax"], omegaiTlen=cfg["omegaiTlen"], lmaxOf=p8Lmax(Rprem), split=true)
        Dict{String,Any}("psv" => sp["MT|psv"], "sh" => sp["MT|sh"],
                         "nd_psv" => size(psv.T, 1), "nd_sh" => size(sh.T, 1))
    end
    return res
end

println("\n---- spectra (cached per scheme and h)")
p8Ours = Dict{Tuple{String,Float64},Any}()
for (sc, h) in p8Runs
    println(sc, "  h = ", h, " km")
    p8Ours[(sc, h)] = p8Spectra(sc, h)
end

"Relative misfit over all components, frequencies and stations."
function p8Misfit(A, B)
    num = 0.0; den = 0.0
    for c in (:ur, :uθ, :uφ)
        a = getfield(A, c); b = getfield(B, c); n = min(size(a, 1), size(b, 1))
        num += sum(abs2, a[1:n, :] .- b[1:n, :]); den += sum(abs2, b[1:n, :])
    end
    return sqrt(num / den)
end

# classic DSM: re → (spectra, mean radial spacing)
p8Classic = Dict{Tuple{String,String},Any}()
p8ClassicH = Dict{Tuple{String,String},Float64}()
for re in ("1e-2", "1e-3", "1e-4", "1e-5"), part in ("psv", "sh")
    dir = joinpath(p8ClassicDir, "re_$(re)_$(part)")
    isdir(dir) || continue
    p8Classic[(re, part)] = classicSpectra(dir, p8θ; parts=(part,))
    g = joinpath(dir, "grid.txt")
    isfile(g) && (p8ClassicH[(re, part)] = parse(Float64, strip(read(g, String))))
end
haskey(p8Classic, ("1e-5", "psv")) || error("no classic reference in $(p8ClassicDir)")

# dsmti-3.1.8 (Takeuchi), zero frequency added: 3 components together (P-SV + SH);
# file per station: header, then "i Re/Im(Z) Re/Im(R) Re/Im(T)".  ngrid_r = total
# radial points (coarsened over 4 frequency bands) -> mean spacing R/ngrid_r.
p8DsmtiDir = joinpath(flexopt_root, "tmp", "classicDSM_out", "dsmti500")
function p8ReadDsmti(dir, θdeg)
    rows = [readlines(joinpath(dir, @sprintf("S%03d.spc", round(Int, θ)))) for θ in θdeg]
    hdr = split(rows[1][1]); tlen = parse(Float64, hdr[1]); ωI = parse(Float64, hdr[3])
    nω = length(rows[1]) - 1
    U = [zeros(ComplexF64, nω, length(θdeg)) for _ in 1:3]
    for (ir, L) in enumerate(rows), k in 1:nω
        v = parse.(Float64, split(L[k+1])[2:7])
        for c in 1:3
            U[c][k, ir] = complex(v[2c-1], v[2c])
        end
    end
    ωs = 2π .* (0:nω-1) ./ tlen
    return (ω=collect(ωs), ωI=ωI, tlen=tlen, θdeg=collect(θdeg), ur=U[1], uθ=U[2], uφ=U[3])
end
p8Dsmti = Dict{Int,Any}()
for ng in (2000, 4000, 8000, 16000)
    dir = joinpath(p8DsmtiDir, "ng$(ng)")
    isdir(dir) && (p8Dsmti[ng] = p8ReadDsmti(dir, p8θ))
end
p8Sub(A, B) = (ω=A.ω, ωI=A.ωI, tlen=A.tlen, θdeg=A.θdeg, ur=A.ur .- B.ur, uθ=A.uθ .- B.uθ, uφ=A.uφ .- B.uφ)

# reference: P-SV -> tipsv re = 1e-5 (converges cleanly, ~h², good to ~1e-5);
#            SH   -> our OPT5 at h = 2.5 km: tish is not self-consistent below
#            re = 1e-3 (successive refinements differ more and more; a third tish
#            problem, not traced yet), so its points are plotted against OPT5.
p8RefName = Dict("psv" => "classic DSM (tipsv re=1e-5)", "sh" => "OPT5 h=2.5 km")
println("\n---- relative misfit (reference: P-SV ", p8RefName["psv"], ", SH ", p8RefName["sh"], ")")
p8Rows = NamedTuple[]
for part in ("psv", "sh")
    own = p8Ours[("OPT5", 2.5)][part]
    ref = part == "psv" ? p8Classic[("1e-5", part)] : own
    println(uppercase(part), ":")
    for (sc, h) in p8Runs
        A = p8Ours[(sc, h)][part]
        e = p8Misfit(A, ref); e2 = p8Misfit(A, own)
        push!(p8Rows, (part=part, scheme=sc, h=h, err=e, errOwn=e2, setting=string(h)))
        @printf("  %-8s h = %5.1f km   vs reference %.3e   vs OPT5(2.5 km) %.3e\n", sc, h, e, e2)
    end
    for re in (part == "psv" ? ("1e-2", "1e-3", "1e-4") : ("1e-2", "1e-3", "1e-4", "1e-5"))
        haskey(p8Classic, (re, part)) || continue
        e = p8Misfit(p8Classic[(re, part)], ref); e2 = p8Misfit(p8Classic[(re, part)], own)
        hh = get(p8ClassicH, (re, part), NaN)
        push!(p8Rows, (part=part, scheme="classic DSM", h=hh, err=e, errOwn=e2, setting="re=$(re)"))
        @printf("  classic re=%s (mean h %.2f km)   vs reference %.3e   vs OPT5(2.5 km) %.3e\n", re, hh, e, e2)
    end
    @printf("  classic(re=1e-5) vs OPT5(2.5 km): %.3e\n", p8Misfit(p8Classic[("1e-5", part)], own))
    # dsmti: its P-SV part = total - (OPT5 2.5 km SH); its SH part = total - (tipsv 1e-5 P-SV)
    other = part == "psv" ? p8Ours[("OPT5", 2.5)]["sh"] : p8Classic[("1e-5", "psv")]
    for ng in sort(collect(keys(p8Dsmti)))
        A = p8Sub(p8Dsmti[ng], other)
        e = p8Misfit(A, ref); hh = Rprem / ng
        push!(p8Rows, (part=part, scheme="dsmti", h=hh, err=e, errOwn=p8Misfit(A, own), setting="ngrid_r=$(ng)"))
        @printf("  dsmti ngrid_r=%5d (mean h %.2f km)   vs reference %.3e\n", ng, hh, e)
    end
end
for part in ("psv", "sh"), sc in ("DSMconv", "OPT3", "OPT5")
    rr = sort(filter(r -> r.part == part && r.scheme == sc && r.err > 0, p8Rows); by=r -> -r.h)
    length(rr) < 2 && continue
    @printf("%-4s %-8s orders vs classic: %s\n", part, sc,
            join((@sprintf("%5.2f", log2(rr[i].err / rr[i+1].err) / log2(rr[i].h / rr[i+1].h)) for i in 1:length(rr)-1), " "))
end

let fig = Figure(size=(1400, 600))
    cols = Dict("DSMconv" => :orange, "OPT3" => :dodgerblue, "OPT5" => :crimson, "classic DSM" => :black,
                "dsmti" => :gray40)
    for (ip, part) in enumerate(("psv", "sh"))
        ax = Axis(fig[1, ip]; xscale=log10, yscale=log10, xlabel="radial spacing h (km)",
                  ylabel="relative misfit (Z, R, T, all f and stations)",
                  title=(part == "psv" ? "P-SV" : "SH") * ": PREM, MT at 500 km, vs " * p8RefName[part])
        for sc in ("DSMconv", "OPT3", "OPT5", "classic DSM", "dsmti")
            rr = sort(filter(r -> r.part == part && r.scheme == sc && isfinite(r.h) && r.err > 0, p8Rows); by=r -> r.h)
            isempty(rr) && continue
            scatterlines!(ax, [r.h for r in rr], [r.err for r in rr]; color=cols[sc],
                          marker=sc == "classic DSM" ? :rect : sc == "dsmti" ? :diamond : :circle,
                          label=sc == "classic DSM" ? (part == "psv" ? "tipsv" : "tish") : sc)
        end
        hs = [2.5, 20.0]
        for (p, ls) in ((2, :dash), (4, :dot))
            y0 = minimum(r.err for r in p8Rows if r.part == part && r.err > 0)
            lines!(ax, hs, y0 .* (hs ./ hs[1]) .^ p; color=:gray, linestyle=ls, label="h^$(p)")
        end
        axislegend(ax; position=:rb)
    end
    devSave(fig, "convergence_psv_sh")
end

# record section at h = 10 km: classic (fine) vs DSMconv / OPT3 / OPT5
p8Total(d) = (ω=d["psv"].ω, ωI=d["psv"].ωI, tlen=d["psv"].tlen, θdeg=d["psv"].θdeg,
              ur=d["psv"].ur .+ d["sh"].ur, uθ=d["psv"].uθ .+ d["sh"].uθ, uφ=d["psv"].uφ .+ d["sh"].uφ)
let tot = p8Total
    c = Dict("psv" => p8Classic[("1e-5", "psv")], "sh" => p8Classic[("1e-5", "sh")])
    tr = Dict("classic DSM (re 1e-5)" => elasticTimeSeries(tot(c); np=4096, taper=0.2))
    for sc in ("DSMconv", "OPT3", "OPT5")
        tr["$(sc) h=10"] = elasticTimeSeries(tot(p8Ours[(sc, 10.0)]); np=4096, taper=0.2)
    end
    keysPlot = ("classic DSM (re 1e-5)", "DSMconv h=10", "OPT3 h=10", "OPT5 h=10")
    plotRecordSection3C(tr, collect(p8θ); force="MT", keys=keysPlot, normalize=:trace, gain=2.0,
                        sourceLabel="moment tensor at 500 km", fileTag="MT500_h10")
end

# ---------------------------------------------------------------------------
# accuracy vs CPU time (one core): outputs/p8/timing_ours.csv (stage p8t, -t 1) and
# tmp/classic_timing/timing.csv (bash run_timing.sh, same Mac).  Our time = systems +
# spectra of the 3 components (P-SV and SH together); tipsv = P-SV only, tish = SH
# only, dsmti = both.
let fo = joinpath(DEV_OUTPUT_DIR[], "p8", "timing_ours.csv"),
    fc = joinpath(flexopt_root, "tmp", "classic_timing", "timing.csv")
    if isfile(fo) && isfile(fc)
        tOurs = Dict{Tuple{String,String},Float64}()
        for l in readlines(fo)[2:end]
            sc, h, th, ts, tsp = split(l, ",")
            th == "1" && (tOurs[(sc, string(parse(Float64, h)))] = parse(Float64, ts) + parse(Float64, tsp))
        end
        tCl = Dict{Tuple{String,String},Float64}()
        for l in readlines(fc)[2:end]
            code, setting, sec = split(l, ",")
            tCl[(code, setting)] = parse(Float64, sec)
        end
        timeOf(r) = r.scheme in ("DSMconv", "OPT3", "OPT5") ? get(tOurs, (r.scheme, r.setting), NaN) :
                    r.scheme == "dsmti" ? get(tCl, ("dsmti", r.setting), NaN) :
                    get(tCl, (r.part == "psv" ? "tipsv" : "tish", r.setting), NaN)
        println("\n---- accuracy vs CPU time (one core)")
        fig = Figure(size=(1400, 600))
        cols = Dict("DSMconv" => :orange, "OPT3" => :dodgerblue, "OPT5" => :crimson, "classic DSM" => :black,
                    "dsmti" => :gray40)
        for (ip, part) in enumerate(("psv", "sh"))
            ax = Axis(fig[1, ip]; xscale=log10, yscale=log10, xlabel="CPU time, one core (s)",
                      ylabel="relative misfit", title=(part == "psv" ? "P-SV" : "SH") * ": accuracy vs cost (vs " * p8RefName[part] * ")")
            for sc in ("DSMconv", "OPT3", "OPT5", "classic DSM", "dsmti")
                rr = filter(r -> r.part == part && r.scheme == sc && r.err > 0 && isfinite(timeOf(r)), p8Rows)
                isempty(rr) && continue
                rr = sort(rr; by=timeOf)
                for r in rr
                    @printf("  %-4s %-12s %-14s  %8.1f s   misfit %.3e\n", part, sc, r.setting, timeOf(r), r.err)
                end
                scatterlines!(ax, timeOf.(rr), [r.err for r in rr]; color=cols[sc],
                              marker=sc == "classic DSM" ? :rect : sc == "dsmti" ? :diamond : :circle,
                              label=sc == "classic DSM" ? (part == "psv" ? "tipsv" : "tish") : sc)
            end
            axislegend(ax; position=:lb)
        end
        devSave(fig, "accuracy_vs_cputime")
    else
        println("\n(accuracy vs CPU time: run stage p8t with -t 1 and tmp/classic_timing/run_timing.sh first)")
    end
end
