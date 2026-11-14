# application/Sonification/run.jl -- hear a signal: the output "sink" that turns any
# flexOPT trace (Vector + sampling interval) into a WAV file. Standalone demo, no recipe.
#
#   julia --project=. application/Sonification/run.jl                    # synthetic Earth hum (normal modes)
#   julia --project=. application/Sonification/run.jl days=2 acceleration=50000
#   julia --project=. application/Sonification/run.jl file=path/to/trace.csv dt=0.01 acceleration=100
#     (a CSV/text file with one sample per line)
#
# Output: application/Sonification/outputs/*.wav and a waveform PNG.
# Requires WAV.jl in the environment (Pkg.add("WAV") once if it is missing).

FLEXOPT_MODULES = [:sonification]
include(joinpath(@__DIR__, "..", "common", "appBootstrap.jl"))
using CairoMakie, DelimitedFiles
CairoMakie.activate!()

let A = appArgs(@isdefined(APP_ARGS) ? APP_ARGS : ARGS;
                days=4.0, dt=10.0, acceleration=1.0e5, file="", octave=0.0)
    dir = appOutputDir(@__DIR__)
    if isempty(A.file)
        # A few fundamental spheroidal modes of the Earth (frequencies in mHz, quality factors Q):
        # an illustrative "hum" after a large earthquake, not a computed synthetic.
        modes = [(name="0S2", f=0.3093, Q=510), (name="0S3", f=0.4686, Q=420),
                 (name="0S4", f=0.6470, Q=370), (name="0S0", f=0.8143, Q=5300),
                 (name="0S5", f=0.8400, Q=350), (name="1S0", f=1.6313, Q=1900)]
        dt = Float64(A.dt)
        t = collect(0:dt:Float64(A.days) * 86400)
        signal = zeros(length(t))
        for (k, m) in enumerate(modes)
            ω = 2π * m.f * 1e-3
            signal .+= (1 / k) .* cos.(ω .* t .+ k) .* exp.(-ω .* t ./ (2 * m.Q))
        end
        name = "earth_hum_$(A.days)days"
        println("synthetic normal-mode signal: ", join([m.name for m in modes], ", "))
    else
        signal = vec(Float64.(readdlm(A.file)))
        dt = Float64(A.dt)
        t = collect(0:length(signal)-1) .* dt
        name = splitext(basename(A.file))[1]
    end
    rendered = sonify_signal(signal, dt; acceleration=Float64(A.acceleration), octave=Float64(A.octave))
    wav = write_sonification(joinpath(dir, name * ".wav"), rendered)
    println("audio -> ", wav, "   (", round(rendered.audio_duration, digits=2), " s, playback ×",
            round(rendered.playback_factor, sigdigits=3), ")")
    fig = Figure(size=(800, 300))
    ax = Axis(fig[1, 1], xlabel="time (h)", ylabel="amplitude", title=name)
    lines!(ax, t ./ 3600, signal)
    appSave(fig, dir, name)
end
