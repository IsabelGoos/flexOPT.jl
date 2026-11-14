# dev/kirishima_oct2026/stages/step2_compare.jl
# Ported from notebooks/Elastic2DtopoBenchmark_oct2026.ipynb, cell step2-*-comparison.
function compact_waveform_metrics(reference, candidate)
    full = waveform_metrics(reference, candidate; samples=1001)
    (correlation=full.correlation, relative_error=full.relative_error,
     relative_rmse=full.relative_rmse,
     optimal_amplitude=full.optimal_amplitude,
     bias=full.bias, error_variance=full.error_variance)
end

function solver_traces_for_station(station)
    traces = (
        FD3=(time=fdWaveformZ.time, values=fdWaveformZ.values[:, station]),
        convFD3=(time=convFDWaveformZ.time,
            values=convFDWaveformZ.values[:, station]),
        OPT3=(time=optWaveformZ.time, values=optWaveformZ.values[:, station]),
    )
    runOPT4 && (traces = merge(traces, (OPT4=(time=opt4WaveformZ.time,
        values=opt4WaveformZ.values[:, station]),)))
    runOPT5 && (traces = merge(traces, (OPT5=(time=opt5WaveformZ.time,
        values=opt5WaveformZ.values[:, station]),)))
    traces
end

stationErrorMetrics = map(eachindex(receiverX)) do station
    reference = specfemWaveformsZ[station]
    candidates = solver_traces_for_station(station)
    metrics = map(trace -> compact_waveform_metrics(reference, trace), candidates)
    (station=station, x_km=receiverX[station]/1e3, metrics=metrics)
end
foreach(display, stationErrorMetrics)

methodNames = keys(first(stationErrorMetrics).metrics)
globalErrorMetrics = NamedTuple{methodNames}(map(methodNames) do method
    rows = [getproperty(item.metrics, method) for item in stationErrorMetrics]
    (mean_correlation=mean(getproperty.(rows, :correlation)),
     rms_relative_error=sqrt(mean(abs2, getproperty.(rows, :relative_error))),
     maximum_relative_error=maximum(getproperty.(rows, :relative_error)),
     mean_raw_relative_rmse=mean(getproperty.(rows, :relative_rmse)))
end)
display(globalErrorMetrics)

errorFigure = Figure(size=(1050, 430))
correlationAxis = Axis(errorFigure[1, 1]; xlabel="receiver x (km)",
    ylabel="correlation with SPECFEM2D")
relativeAxis = Axis(errorFigure[1, 2]; xlabel="receiver x (km)",
    ylabel="shape-relative error")
for method in methodNames
    rows = [getproperty(item.metrics, method) for item in stationErrorMetrics]
    xs = getproperty.(stationErrorMetrics, :x_km)
    scatterlines!(correlationAxis, xs, getproperty.(rows, :correlation);
        marker=:circle, label=string(method))
    scatterlines!(relativeAxis, xs, getproperty.(rows, :relative_error);
        marker=:circle, label=string(method))
end
axislegend(correlationAxis); axislegend(relativeAxis)
devSave(errorFigure, "errors_vs_SPECFEM2D")

# keep the numbers next to the figures
metricsPath = joinpath(DEV_OUTPUT_DIR[], DEV_STEP[], "errors_vs_SPECFEM2D.jld2")
jldsave(metricsPath; stationErrorMetrics, globalErrorMetrics,
    receiverX, optCacheFile, fdCacheFile, specfemCacheFile)
println("metrics -> ", metricsPath)

