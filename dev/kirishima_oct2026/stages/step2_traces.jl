# dev/kirishima_oct2026/stages/step2_traces.jl
# Ported from notebooks/Elastic2DtopoBenchmark_oct2026.ipynb, cell step2-*-traces.
waveformQuantity = @isdefined(waveformQuantity) ? waveformQuantity : :velocity
function displacement_traces(history, times, xaxis, zindex, receivers)
    traces = Matrix{Float64}(undef, length(times), length(receivers))
    for (j, receiver) in enumerate(receivers)
        ix = argmin(abs.(xaxis .- receiver))
        traces[:, j] .= Float64.(history[ix, zindex, :])
    end
    return (time=Float64.(times), values=traces)
end

function velocity_traces(history, times, xaxis, zindex, receivers)
    traces = Matrix{Float64}(undef, length(times) - 1, length(receivers))
    for (j, receiver) in enumerate(receivers)
        ix = argmin(abs.(xaxis .- receiver))
        traces[:, j] .= diff(Float64.(history[ix, zindex, :])) ./ diff(times)
    end
    return (time=(times[1:end-1] .+ times[2:end]) ./ 2, values=traces)
end
function velocity_traces_at_points(history, times, xaxis, zaxis, receivers;
    sampling=:nearest)
    sampling === :nearest || error("only :nearest is enabled in this cheap benchmark path")
    traces = Matrix{Float64}(undef, length(times) - 1, length(receivers))
    for (j, receiver) in enumerate(receivers)
        ix = argmin(abs.(xaxis .- receiver.x))
        iz = argmin(abs.(zaxis .- receiver.z))
        traces[:, j] .= diff(Float64.(history[ix, iz, :])) ./ diff(times)
    end
    (time=(times[1:end-1] .+ times[2:end]) ./ 2, values=traces,
        receivers=collect(receivers), sampling)
end
fdReceiverZIndex = argmin(abs.(fdCoordinates.z .- receiverZ))
optReceiverZIndex = argmin(abs.(zOPT .- receiverZ))
fdDisplacementX = displacement_traces(uxFD, fdTimes, fdCoordinates.x,
    fdReceiverZIndex, receiverX)
fdDisplacementZ = displacement_traces(uzFD, fdTimes, fdCoordinates.x,
    fdReceiverZIndex, receiverX)
optDisplacementX = displacement_traces(uxOPT, optTimes, xOPT,
    optReceiverZIndex, receiverX)
optDisplacementZ = displacement_traces(uzOPT, optTimes, xOPT,
    optReceiverZIndex, receiverX)
convFDDisplacementX = displacement_traces(uxConvFD, convFDTimes, xOPT,
    optReceiverZIndex, receiverX)
convFDDisplacementZ = displacement_traces(uzConvFD, convFDTimes, xOPT,
    optReceiverZIndex, receiverX)
optNoGammaDisplacementX = runOPTWithoutGamma ? displacement_traces(
    uxOPTWithoutGamma, optTimes, xOPT, optReceiverZIndex, receiverX) : nothing
optNoGammaDisplacementZ = runOPTWithoutGamma ? displacement_traces(
    uzOPTWithoutGamma, optTimes, xOPT, optReceiverZIndex, receiverX) : nothing
opt4DisplacementX = runOPT4 ? displacement_traces(
    uxOPT4, opt4Times, xOPT, optReceiverZIndex, receiverX) : nothing
opt4DisplacementZ = runOPT4 ? displacement_traces(
    uzOPT4, opt4Times, xOPT, optReceiverZIndex, receiverX) : nothing
opt5DisplacementX = runOPT5 ? displacement_traces(
    uxOPT5, opt5Times, xOPT, optReceiverZIndex, receiverX) : nothing
opt5DisplacementZ = runOPT5 ? displacement_traces(
    uzOPT5, opt5Times, xOPT, optReceiverZIndex, receiverX) : nothing

fdVelocityX = velocity_traces(uxFD, fdTimes, fdCoordinates.x,
    fdReceiverZIndex, receiverX)
fdVelocityZ = velocity_traces(uzFD, fdTimes, fdCoordinates.x,
    fdReceiverZIndex, receiverX)
optVelocityX = velocity_traces(uxOPT, optTimes, xOPT,
    optReceiverZIndex, receiverX)
optVelocityZ = velocity_traces(uzOPT, optTimes, xOPT,
    optReceiverZIndex, receiverX)
convFDVelocityX = velocity_traces(uxConvFD, convFDTimes, xOPT,
    optReceiverZIndex, receiverX)
convFDVelocityZ = velocity_traces(uzConvFD, convFDTimes, xOPT,
    optReceiverZIndex, receiverX)
optNoGammaVelocityX = runOPTWithoutGamma ? velocity_traces(
    uxOPTWithoutGamma, optTimes, xOPT, optReceiverZIndex, receiverX) : nothing
optNoGammaVelocityZ = runOPTWithoutGamma ? velocity_traces(
    uzOPTWithoutGamma, optTimes, xOPT, optReceiverZIndex, receiverX) : nothing
opt4VelocityX = runOPT4 ? velocity_traces(uxOPT4, opt4Times, xOPT,
    optReceiverZIndex, receiverX) : nothing
opt4VelocityZ = runOPT4 ? velocity_traces(uzOPT4, opt4Times, xOPT,
    optReceiverZIndex, receiverX) : nothing
opt5VelocityX = runOPT5 ? velocity_traces(uxOPT5, opt5Times, xOPT,
    optReceiverZIndex, receiverX) : nothing
opt5VelocityZ = runOPT5 ? velocity_traces(uzOPT5, opt5Times, xOPT,
    optReceiverZIndex, receiverX) : nothing
if waveformQuantity === :velocity
    fdWaveformX, fdWaveformZ = fdVelocityX, fdVelocityZ
    optWaveformX, optWaveformZ = optVelocityX, optVelocityZ
    convFDWaveformX, convFDWaveformZ = convFDVelocityX, convFDVelocityZ
    optNoGammaWaveformX, optNoGammaWaveformZ =
        optNoGammaVelocityX, optNoGammaVelocityZ
    opt4WaveformX, opt4WaveformZ = opt4VelocityX, opt4VelocityZ
    opt5WaveformX, opt5WaveformZ = opt5VelocityX, opt5VelocityZ
    waveformUnitLabel, waveformComponentPrefix = "m/s", "v"
else
    fdWaveformX, fdWaveformZ = fdDisplacementX, fdDisplacementZ
    optWaveformX, optWaveformZ = optDisplacementX, optDisplacementZ
    convFDWaveformX, convFDWaveformZ = convFDDisplacementX, convFDDisplacementZ
    optNoGammaWaveformX, optNoGammaWaveformZ =
        optNoGammaDisplacementX, optNoGammaDisplacementZ
    opt4WaveformX, opt4WaveformZ = opt4DisplacementX, opt4DisplacementZ
    opt5WaveformX, opt5WaveformZ = opt5DisplacementX, opt5DisplacementZ
    waveformUnitLabel, waveformComponentPrefix = "m", "u"
end
fdVelocity, optVelocity = fdVelocityZ, optVelocityZ # physical velocity aliases
receiver = cld(length(receiverX), 2)
fdTrace = (time=fdWaveformZ.time, values=fdWaveformZ.values[:, receiver])
optTrace = (time=optWaveformZ.time, values=optWaveformZ.values[:, receiver])
convFDTrace = (time=convFDWaveformZ.time,
    values=convFDWaveformZ.values[:, receiver])
optNoGammaTrace = runOPTWithoutGamma ? (time=optNoGammaWaveformZ.time,
    values=optNoGammaWaveformZ.values[:, receiver]) : nothing
opt4Trace = runOPT4 ? (time=opt4WaveformZ.time,
    values=opt4WaveformZ.values[:, receiver]) : nothing
opt5Trace = runOPT5 ? (time=opt5WaveformZ.time,
    values=opt5WaveformZ.values[:, receiver]) : nothing
preSPECFEMTraces = (explicitFD3=fdTrace,
    convFD3=convFDTrace, OPT3=optTrace)
runOPTWithoutGamma && (preSPECFEMTraces = merge(preSPECFEMTraces,
    (OPT3_no_Gamma=optNoGammaTrace,)))
runOPT4 && (preSPECFEMTraces = merge(preSPECFEMTraces, (OPT4=opt4Trace,)))
runOPT5 && (preSPECFEMTraces = merge(preSPECFEMTraces, (OPT5=opt5Trace,)))
preSPECFEM = plot_solver_benchmark(preSPECFEMTraces;
    title="$(caseTitle) interior model — normalized $(waveformComponentPrefix)_z")
devSave(preSPECFEM.figure, "traces_before_SPECFEM_" * string(waveformComponentPrefix) * "z")
nothing

