module sonification

# Sonification of seismograms held as plain arrays (signal + sampling interval).
#
# Ported from SonifSismo.jl (examples/sonify_one_trace.jl), with the Seis.jl
# trace type removed: Seis pins DSP to 0.6/0.7, which conflicts with flexOPT
# (DSP 0.8). Reading miniSEED/JMA archives stays in SonifSismo; everything
# here works on flexOPT synthetics or any `Vector`/`Matrix` of samples.
#
# Dependencies: Statistics, Base64 (stdlib), DSP, Makie (already in flexOPT),
# WAV (add once with `Pkg.add("WAV")`).

using Statistics
using Base64
using DSP: conv, resample
using Makie
using WAV

export resample_linear, normalize_audio, fade_edges!, add_harmonics, set_volume
export suggest_acceleration, sonify_signal, sonify_stereo
export concatenate_audio, write_sonification, read_audio_mono
export convolve_with_audio, audio_html, waveform_video

const AUDIO_RATE = 44_100

# ---------------------------------------------------------------------------
# Elementary signal operations (columns = channels; vectors = mono)
# ---------------------------------------------------------------------------

"""
    resample_linear(signal, input_rate, output_rate)

Dependency-light linear interpolation between sampling rates. Adequate for
exploratory sonification; `resampler=:dsp` in `sonify_signal` uses a
band-limited polyphase resampler instead.
"""
function resample_linear(signal::AbstractVector, input_rate::Real, output_rate::Real)
    input_rate > 0 || throw(ArgumentError("input_rate must be positive"))
    output_rate > 0 || throw(ArgumentError("output_rate must be positive"))
    isempty(signal) && return Float64[]
    length(signal) == 1 && return [Float64(first(signal))]

    output_length = max(1, round(Int, length(signal) * output_rate / input_rate))
    locations = range(1.0, stop=Float64(length(signal)), length=output_length)
    output = Vector{Float64}(undef, output_length)
    for (i, location) in pairs(locations)
        left = floor(Int, location)
        right = min(left + 1, length(signal))
        weight = location - left
        output[i] = (1 - weight) * signal[left] + weight * signal[right]
    end
    return output
end

function _resample(signal::AbstractVector, input_rate::Real, output_rate::Real,
    resampler::Symbol)
    if resampler === :linear
        return resample_linear(signal, input_rate, output_rate)
    elseif resampler === :dsp
        return Float64.(resample(Float64.(signal), output_rate / input_rate))
    end
    throw(ArgumentError("resampler must be :linear or :dsp"))
end

"""
    normalize_audio(signal; peak=0.92)

Remove the mean of each channel and scale so that the largest absolute sample
over *all* channels equals `peak`. A common scale keeps relative amplitudes
between receivers.
"""
function normalize_audio(signal::AbstractVecOrMat; peak::Real=0.92)
    0 < peak <= 1 || throw(ArgumentError("peak must be between 0 and 1"))
    centered = Float64.(signal) .- mean(signal; dims=1)
    amplitude = maximum(abs, centered; init=0.0)
    amplitude == 0 && return centered
    return centered .* (peak / amplitude)
end

function fade_edges!(signal::AbstractVecOrMat, rate::Real; seconds::Real=0.03)
    seconds <= 0 && return signal
    nsamples = size(signal, 1)
    nfade = min(round(Int, seconds * rate), nsamples ÷ 2)
    nfade == 0 && return signal
    ramp = range(0.0, stop=1.0, length=nfade)
    for j in axes(signal, 2)
        view(signal, 1:nfade, j) .*= ramp
        view(signal, nsamples-nfade+1:nsamples, j) .*= reverse(ramp)
    end
    return signal
end

"""
    add_harmonics(signal; even=0, odd=0)

Add nonlinear timbral colour: `even` weights x² (second harmonic), `odd`
weights x³ (third harmonic). Duration is unchanged.
"""
function add_harmonics(signal::AbstractVecOrMat; even::Real=0.0, odd::Real=0.0)
    x = normalize_audio(signal; peak=1.0)
    second_harmonic = x .^ 2 .- mean(x .^ 2; dims=1)
    third_harmonic = x .^ 3
    return normalize_audio(x .+ even .* second_harmonic .+ odd .* third_harmonic)
end

"""
    set_volume(signal; level=0.92, gain=1, saturation=:none)

Set the exported peak `level`. With `saturation=:soft` (tanh) or `:hard`
(clip), `gain > 1` drives the normalized signal into saturation so that quiet
parts become louder relative to the peaks.
"""
function set_volume(
    signal::AbstractVecOrMat;
    level::Real=0.92,
    gain::Real=1.0,
    saturation::Symbol=:none,
)
    0 < level <= 1 || throw(ArgumentError("level must be between 0 and 1"))
    gain > 0 || throw(ArgumentError("gain must be positive"))
    x = normalize_audio(signal; peak=1.0)
    if saturation === :none
        gain == 1 || throw(ArgumentError(
            "use saturation=:soft or :hard when gain differs from 1",
        ))
        return level .* x
    elseif saturation === :soft
        return level .* tanh.(gain .* x)
    elseif saturation === :hard
        return level .* clamp.(gain .* x, -1.0, 1.0)
    end
    throw(ArgumentError("saturation must be :none, :soft, or :hard"))
end

# ---------------------------------------------------------------------------
# Sampling helpers
# ---------------------------------------------------------------------------

"""
    sampling_interval(times; rtol=1e-6)

Return the constant sampling interval of `times`, or throw if the samples are
not uniformly spaced.
"""
function sampling_interval(times::AbstractVector; rtol::Real=1e-6)
    length(times) >= 2 || throw(ArgumentError("need at least two time samples"))
    steps = diff(Float64.(times))
    dt = mean(steps)
    dt > 0 || throw(ArgumentError("times must be increasing"))
    maximum(abs.(steps .- dt)) <= rtol * abs(dt) * max(1, length(times)) ||
        throw(ArgumentError("times are not uniformly sampled; resample first"))
    return dt
end

"""
    suggest_acceleration(dominant_frequency; target_frequency=220)

Playback factor that maps a seismic `dominant_frequency` (Hz) onto
`target_frequency` (Hz, default A3). Pitch and duration are coupled: the audio
lasts `seismic_duration / acceleration`.
"""
suggest_acceleration(dominant_frequency::Real; target_frequency::Real=220.0) =
    target_frequency / dominant_frequency

# ---------------------------------------------------------------------------
# Rendering
# ---------------------------------------------------------------------------

function _render(
    channels::AbstractMatrix,
    dt::Real;
    acceleration::Real,
    octave::Real,
    harmonics::Tuple{<:Real,<:Real},
    audio_rate::Integer,
    resampler::Symbol,
    fade_seconds::Real,
)
    dt > 0 || throw(ArgumentError("dt must be positive"))
    acceleration > 0 || throw(ArgumentError("acceleration must be positive"))
    playback_factor = acceleration * 2.0^octave
    seismic_rate = inv(dt)
    virtual_rate = seismic_rate * playback_factor
    columns = [_resample(view(channels, :, j), virtual_rate, audio_rate, resampler)
               for j in axes(channels, 2)]
    nsamples = minimum(length, columns)
    audio = Matrix{Float64}(undef, nsamples, length(columns))
    for (j, c) in enumerate(columns)
        audio[:, j] .= view(c, 1:nsamples)
    end
    fade_edges!(audio, audio_rate; seconds=fade_seconds)
    audio = add_harmonics(audio; even=harmonics[1], odd=harmonics[2])
    return audio, playback_factor, size(channels, 1) / seismic_rate
end

function _rendered(audio, audio_rate, playback_factor, seismic_duration)
    return (
        signal=audio,
        sample_rate=audio_rate,
        playback_factor=playback_factor,
        seismic_duration=seismic_duration,
        audio_duration=size(audio, 1) / audio_rate,
    )
end

"""
    sonify_signal(signal, dt; acceleration=100, octave=0, harmonics=(0.0, 0.0),
                  level=0.92, gain=1, saturation=:none, audio_rate=44100,
                  resampler=:linear, fade_seconds=0.03)
    sonify_signal(times, signal; kwargs...)

Convert one seismogram (`signal` sampled every `dt` seconds) into mono audio.
`acceleration` compresses time and raises frequency; `octave=1` doubles that
playback factor (pitch and duration are not independent). `harmonics=(even,
odd)` adds timbral colour. `level` sets the output peak; `gain > 1` with
`saturation=:soft`/`:hard` makes quiet sections louder.

Returns a named tuple `(signal, sample_rate, playback_factor,
seismic_duration, audio_duration)`.
"""
function sonify_signal(
    signal::AbstractVector,
    dt::Real;
    acceleration::Real=100.0,
    octave::Real=0.0,
    harmonics::Tuple{<:Real,<:Real}=(0.0, 0.0),
    level::Real=0.92,
    gain::Real=1.0,
    saturation::Symbol=:none,
    audio_rate::Integer=AUDIO_RATE,
    resampler::Symbol=:linear,
    fade_seconds::Real=0.03,
)
    audio, playback_factor, seismic_duration = _render(
        reshape(Float64.(signal), :, 1), dt;
        acceleration, octave, harmonics, audio_rate, resampler, fade_seconds,
    )
    audio = vec(set_volume(audio; level, gain, saturation))
    return _rendered(audio, audio_rate, playback_factor, seismic_duration)
end

sonify_signal(times::AbstractVector, signal::AbstractVector; kwargs...) =
    sonify_signal(signal, sampling_interval(times); kwargs...)

"""
    sonify_stereo(signals, dt; positions=nothing, width=0.8, kwargs...)
    sonify_stereo(times, signals; kwargs...)

Render several receivers (`signals` is `nt × nreceivers`) into one stereo
mix, panning each receiver by its position. `positions` (e.g. receiver x in
metres) is mapped linearly onto `[-width, width]`; without it the receivers
are spread evenly. Relative amplitudes between receivers are preserved
(common normalization). Keywords are those of `sonify_signal`.
"""
function sonify_stereo(
    signals::AbstractMatrix,
    dt::Real;
    positions::Union{Nothing,AbstractVector}=nothing,
    width::Real=0.8,
    acceleration::Real=100.0,
    octave::Real=0.0,
    harmonics::Tuple{<:Real,<:Real}=(0.0, 0.0),
    level::Real=0.92,
    gain::Real=1.0,
    saturation::Symbol=:none,
    audio_rate::Integer=AUDIO_RATE,
    resampler::Symbol=:linear,
    fade_seconds::Real=0.03,
)
    0 <= width <= 1 || throw(ArgumentError("width must be between 0 and 1"))
    nreceivers = size(signals, 2)
    pan = if nreceivers == 1
        [0.0]
    elseif positions === nothing
        collect(range(-width, width; length=nreceivers))
    else
        length(positions) == nreceivers ||
            throw(ArgumentError("need one position per receiver (column)"))
        p = Float64.(positions)
        lo, hi = extrema(p)
        hi == lo ? zeros(nreceivers) : width .* (2 .* (p .- lo) ./ (hi - lo) .- 1)
    end

    audio, playback_factor, seismic_duration = _render(
        Float64.(signals), dt;
        acceleration, octave, harmonics, audio_rate, resampler, fade_seconds,
    )
    # Equal-power panning: θ ∈ [0, π/2], L = cos θ, R = sin θ.
    θ = (pan .+ 1) .* (π / 4)
    stereo = hcat(audio * cos.(θ), audio * sin.(θ))
    stereo = set_volume(stereo; level, gain, saturation)
    return merge(
        _rendered(stereo, audio_rate, playback_factor, seismic_duration),
        (pan=pan,),
    )
end

sonify_stereo(times::AbstractVector, signals::AbstractMatrix; kwargs...) =
    sonify_stereo(signals, sampling_interval(times); kwargs...)

"""
    concatenate_audio(rendered...; gap_seconds=0.25)

Play several renderings one after the other (same sample rate and channel
count), separated by `gap_seconds` of silence.
"""
function concatenate_audio(rendered...; gap_seconds::Real=0.25)
    isempty(rendered) && throw(ArgumentError("nothing to concatenate"))
    rate = first(rendered).sample_rate
    all(r -> r.sample_rate == rate, rendered) ||
        throw(ArgumentError("all renderings must share the same sample rate"))
    nchannels = size(first(rendered).signal, 2)
    all(r -> size(r.signal, 2) == nchannels, rendered) ||
        throw(ArgumentError("all renderings must have the same number of channels"))
    gap = zeros(round(Int, gap_seconds * rate), nchannels)
    pieces = Any[]
    for (i, r) in enumerate(rendered)
        i > 1 && push!(pieces, gap)
        push!(pieces, reshape(r.signal, :, nchannels))
    end
    signal = reduce(vcat, pieces)
    nchannels == 1 && (signal = vec(signal))
    return (
        signal=signal,
        sample_rate=rate,
        playback_factor=first(rendered).playback_factor,
        seismic_duration=sum(r.seismic_duration for r in rendered),
        audio_duration=size(signal, 1) / rate,
    )
end

# ---------------------------------------------------------------------------
# I/O
# ---------------------------------------------------------------------------

function write_sonification(path::AbstractString, rendered)
    mkpath(dirname(abspath(path)))
    wavwrite(
        rendered.signal,
        path;
        Fs=rendered.sample_rate,
        nbits=16,
        compression=WAVE_FORMAT_PCM,
    )
    return path
end

"""
    read_audio_mono(path) -> (signal, sample_rate)

Read a WAV file and average its channels.
"""
function read_audio_mono(path::AbstractString)
    samples, rate, _, _ = wavread(path)
    signal = ndims(samples) == 1 ? vec(Float64.(samples)) :
             vec(mean(Float64.(samples); dims=2))
    return signal, Float64(rate)
end

"""
    audio_html(rendered)

Inline audio player for Jupyter/Pluto (16-bit PCM WAV embedded as base64).
"""
function audio_html(rendered)
    io = IOBuffer()
    wavwrite(rendered.signal, io; Fs=rendered.sample_rate, nbits=16,
        compression=WAVE_FORMAT_PCM)
    payload = base64encode(take!(io))
    return HTML("<audio controls src=\"data:audio/wav;base64,$payload\"></audio>")
end

"""
    convolve_with_audio(rendered, music_path; kernel_seconds=2, mix=0.35)

Use the opening part of a WAV file as a convolution kernel. A measured room
impulse response gives conventional reverberation; a music excerpt gives a
smeared hybrid texture. Works for mono and stereo renderings.
"""
function convolve_with_audio(
    rendered,
    music_path::AbstractString;
    kernel_seconds::Real=2.0,
    mix::Real=0.35,
)
    0 <= mix <= 1 || throw(ArgumentError("mix must be between 0 and 1"))
    kernel_seconds > 0 || throw(ArgumentError("kernel_seconds must be positive"))
    source, source_rate = read_audio_mono(music_path)
    kernel = resample_linear(source, source_rate, rendered.sample_rate)
    kernel = kernel[1:min(length(kernel), round(Int, kernel_seconds * rendered.sample_rate))]
    fade_edges!(kernel, rendered.sample_rate; seconds=min(0.02, kernel_seconds / 4))
    kernel = vec(normalize_audio(kernel; peak=1.0))

    dry_in = reshape(Float64.(rendered.signal), size(rendered.signal, 1), :)
    nwet = size(dry_in, 1) + length(kernel) - 1
    wet = Matrix{Float64}(undef, nwet, size(dry_in, 2))
    for j in axes(dry_in, 2)
        wet[:, j] .= conv(dry_in[:, j], kernel)
    end
    wet = normalize_audio(wet)
    dry = vcat(dry_in, zeros(size(wet, 1) - size(dry_in, 1), size(dry_in, 2)))
    mixed = normalize_audio((1 - mix) .* dry .+ mix .* wet)
    ndims(rendered.signal) == 1 && (mixed = vec(mixed))
    return merge(rendered, (signal=mixed, audio_duration=size(mixed, 1) / rendered.sample_rate))
end

# ---------------------------------------------------------------------------
# Video
# ---------------------------------------------------------------------------

"""
    waveform_video(times, signals, rendered, path; labels=nothing, framerate=30,
                   window_seconds=nothing, title="Sonified waveform",
                   with_audio=true)

Animate the seismograms (`signals`: vector, or `nt × nreceivers` matrix, one
panel per column) with a playhead synchronized to `rendered`, which should be
a sonification of the same time window. The MP4 gets the rendered audio
attached with Makie's bundled ffmpeg. `window_seconds=nothing` shows the whole
record. Works with any Makie backend (e.g. CairoMakie).
"""
function waveform_video(
    times::AbstractVector,
    signals::AbstractVecOrMat,
    rendered,
    path::AbstractString;
    labels=nothing,
    framerate::Integer=30,
    window_seconds::Union{Nothing,Real}=nothing,
    title::AbstractString="Sonified waveform",
    with_audio::Bool=true,
)
    traces = reshape(Float64.(signals), size(signals, 1), :)
    size(traces, 1) == length(times) ||
        throw(ArgumentError("times and signals must have the same number of samples"))
    framerate > 0 || throw(ArgumentError("framerate must be positive"))
    ntraces = size(traces, 2)
    labels = labels === nothing ? ["trace $j" for j in 1:ntraces] : string.(labels)

    seconds = Float64.(times) .- first(times)
    seismic_seconds = last(seconds)
    window = window_seconds === nothing ? seismic_seconds : Float64(window_seconds)
    window > 0 || throw(ArgumentError("window_seconds must be positive"))
    video_seconds = rendered.audio_duration
    # A short held-final-frame tail prevents H.264 muxing from ending before
    # the WAV soundtrack on very short clips.
    nframes = max(2, ceil(Int, video_seconds * framerate) + ceil(Int, 0.25framerate))
    cursor = Observable([0.0])
    scale = maximum(abs, traces .- mean(traces; dims=1); init=0.0)
    scale = scale == 0 ? 1.0 : scale

    figure = Figure(size=(1100, max(300, 180ntraces)))
    axes_list = Axis[]
    for j in 1:ntraces
        axis = Axis(
            figure[j, 1];
            ylabel=labels[j],
            xlabel=j == ntraces ? "Seconds in seismic window" : "",
            title=j == 1 ? title : "",
        )
        lines!(axis, seconds, (traces[:, j] .- mean(traces[:, j])) ./ scale;
            color=:black, linewidth=1)
        vlines!(axis, cursor; color=:darkorange, linewidth=2)
        ylims!(axis, -1.05, 1.05)
        j < ntraces && hidexdecorations!(axis; grid=false)
        push!(axes_list, axis)
    end
    length(axes_list) > 1 && linkxaxes!(axes_list...)
    xlims!(first(axes_list), 0, max(window, eps()))

    mkpath(dirname(abspath(path)))
    silent_path = with_audio ? tempname() * ".mp4" : path
    record(figure, silent_path, 0:(nframes-1); framerate=framerate) do frame
        elapsed_audio = min(frame / framerate, video_seconds)
        elapsed_seismic = min(elapsed_audio * rendered.playback_factor, seismic_seconds)
        left = max(0.0, elapsed_seismic - window / 2)
        right = min(seismic_seconds, left + window)
        left = max(0.0, right - window)
        cursor[] = [elapsed_seismic]
        xlims!(first(axes_list), left, max(right, left + eps()))
    end

    if with_audio
        audio_path = tempname() * ".wav"
        try
            write_sonification(audio_path, rendered)
            ffmpeg = Makie.FFMPEG_jll.ffmpeg()
            run(`$ffmpeg -y -loglevel error -i $silent_path -i $audio_path -c:v libx264 -pix_fmt yuv420p -c:a aac -shortest $path`)
        finally
            rm(silent_path; force=true)
            rm(audio_path; force=true)
        end
    end
    return path
end

waveform_video(dt::Real, signals::AbstractVecOrMat, rendered, path::AbstractString; kwargs...) =
    waveform_video((0:size(signals, 1)-1) .* dt, signals, rendered, path; kwargs...)

end # module sonification
