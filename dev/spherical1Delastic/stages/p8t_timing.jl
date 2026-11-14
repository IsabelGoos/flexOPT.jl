# dev/spherical1Delastic/stages/p8t_timing.jl
# CPU time of our schemes for the p8 benchmark, NOT cached: global systems
# (P-SV + SH) and the 3-component spectra (P-SV and SH together), exactly as in p8.
# For the fair comparison with the classic codes run it on ONE core:
#     julia --project=. -t 1 dev/spherical1Delastic/run.jl p8t
# (and, if you like, again with -t auto: both are kept, keyed by the thread count).
# Writes outputs/p8/timing_ours.csv (scheme, h_km, threads, s_systems, s_spectra).

p8tRecipes = Dict{Any,Any}()
p8tRuns = [("DSMconv", 20.0), ("DSMconv", 10.0), ("DSMconv", 5.0),
           ("OPT3", 20.0), ("OPT3", 10.0), ("OPT3", 5.0),
           ("OPT5", 20.0), ("OPT5", 10.0), ("OPT5", 5.0), ("OPT5", 2.5)]
p8tSources = Dict("MT" => pointSource(M=(1.0, 0.6, -0.8, -0.3, 0.4, -0.7)))
p8tGamma = Dict(c.name => sourceGammaStencils(c; verbose=false)
                for c in filter(c -> c.name in ("OPT3", "OPT5"), shConfigs))
prem = planet1D.getSet1Dmodel!(:PREM)
Rprem = prem.averagedPlanetRadiusInKilometer
nthr = Threads.nthreads()
println("p8t: ", nthr, " thread(s)")

function p8tOne(scheme, h)
    nr = scheme == "DSMconv" ? :fallback : :clippedFamily
    kw = (scheme=scheme, configs=shConfigs, recipes=p8tRecipes, label="$(scheme) h=$(h)", nearRows=nr)
    t0 = time()
    psv = elasticGlobalSystem(:PSV, prem, elasticRadialGrid(prem; h=h, region=:all); kw...)
    sh = elasticGlobalSystem(:SH, prem, elasticRadialGrid(prem; h=h, region=:outerSolid); kw...)
    t1 = time()
    radial = scheme == "DSMconv" ? (method=:lagrange, p=2) :
             (method=:gamma, p=sourceOrderOf(scheme), gamma=p8tGamma[scheme])
    pointSourceSynthetics(psv, sh, Rprem; rs_km=Rprem - 500.0, sources=p8tSources, radial=radial,
                          receiversθdeg=10:10:170, φr=π / 2, tlen=3276.8, fmax=0.05, omegaiTlen=1e-4,
                          lmaxOf=ω -> ceil(Int, ω * Rprem / 2.5) + 250, verbose=false)
    t2 = time()
    return (t1 - t0, t2 - t1)
end

p8tOne("DSMconv", 50.0)            # warm-up: compilation is not part of the timing
csv = joinpath(DEV_OUTPUT_DIR[], "p8", "timing_ours.csv")
mkpath(dirname(csv))
rows = String[]
for (sc, h) in p8tRuns
    ts, tsp = p8tOne(sc, h)
    @printf("%-8s h = %5.1f km  threads %d   systems %7.2f s   spectra %8.2f s\n", sc, h, nthr, ts, tsp)
    push!(rows, @sprintf("%s,%.1f,%d,%.3f,%.3f", sc, h, nthr, ts, tsp))
end
keep = isfile(csv) ? [l for l in readlines(csv)[2:end] if split(l, ",")[3] != string(nthr)] : String[]
open(csv, "w") do io
    println(io, "scheme,h_km,threads,s_systems,s_spectra")
    foreach(l -> println(io, l), vcat(keep, rows))
end
println("-> ", csv)
