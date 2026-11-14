#!/usr/bin/env julia

using CairoMakie
using JLD2
using Statistics
using LinearAlgebra

include(joinpath(@__DIR__, "..", "src", "elasticGreens2D.jl"))
using .elasticGreens2D

const ROOT = normpath(joinpath(@__DIR__, ".."))
const DATA_DIRECTORY = joinpath(ROOT, "data", "elastic2d_convergence",
    "homogeneous_flat_free_surface")

CairoMakie.activate!(type="png")

function linear_sample(trace, times)
    x, y = Float64.(trace.time), Float64.(trace.values)
    result = zeros(length(times))
    for (i, time) in pairs(times)
        (time < first(x) || time > last(x)) && continue
        upper = searchsortedfirst(x, time)
        if upper <= 1
            result[i] = y[1]
        elseif upper > length(x)
            result[i] = y[end]
        else
            fraction = (time - x[upper - 1]) / (x[upper] - x[upper - 1])
            result[i] = muladd(fraction, y[upper] - y[upper - 1], y[upper - 1])
        end
    end
    result
end

function trace(result, method, component, station)
    block = getproperty(result, method)
    values = getproperty(block, component)
    method === :SPECFEM2D ? values[station] :
        (time=values.time, values=values.values[:, station])
end

function normalized(values)
    scale = maximum(abs, values)
    scale > eps(Float64) ? values ./ scale : zero(values)
end

function relative_rms(reference, candidate)
    denominator = norm(reference)
    denominator > eps(Float64) || return NaN
    norm(candidate - reference) / denominator
end

"""Stationary travel time of an incident `v1` and reflected `v2` ray."""
function reflected_arrival(source, receiver, v1, v2; iterations=100)
    # Reflection point q=(x,0). The convex travel-time function is minimized
    # with golden-section search; v1 != v2 naturally enforces Snell's law.
    span = abs(receiver.x - source.x) + abs(source.z) + abs(receiver.z)
    left = min(source.x, receiver.x) - span
    right = max(source.x, receiver.x) + span
    phi = (sqrt(5) - 1) / 2
    travel(q) = hypot(q - source.x, source.z) / v1 +
        hypot(receiver.x - q, receiver.z) / v2
    c, d = right - phi * (right - left), left + phi * (right - left)
    fc, fd = travel(c), travel(d)
    for _ in 1:iterations
        if fc < fd
            right, d, fd = d, c, fc
            c = right - phi * (right - left); fc = travel(c)
        else
            left, c, fc = c, d, fd
            d = left + phi * (right - left); fd = travel(d)
        end
    end
    q = (left + right) / 2
    (time=travel(q), reflection_x=q)
end

function arrivals(source, receiver, vp, vs, source_delay)
    distance = hypot(receiver.x - source.x, receiver.z - source.z)
    (
        P=source_delay + distance / vp,
        S=source_delay + distance / vs,
        PP=source_delay + reflected_arrival(source, receiver, vp, vp).time,
        PS=source_delay + reflected_arrival(source, receiver, vp, vs).time,
        SP=source_delay + reflected_arrival(source, receiver, vs, vp).time,
        SS=source_delay + reflected_arrival(source, receiver, vs, vs).time,
    )
end

function load_results()
    files = isdir(DATA_DIRECTORY) ? sort(filter(
        file -> endswith(file, ".jld2"), readdir(DATA_DIRECTORY; join=true))) : String[]
    isempty(files) && error("no convergence data in $DATA_DIRECTORY")
    sort([load(file, "result") for file in files]; by=r -> r.spacing_m, rev=true)
end

function station_label(receiver)
    "x=$(receiver.x / 1e3) km, z=$(receiver.z / 1e3) km"
end

function plot_specfem_waveforms(results; component=:z, normalized_amplitude=false)
    finest = last(results)
    receivers = finest.receivers
    figure = Figure(size=(1200, 230 * length(receivers)))
    colors = Makie.wong_colors()
    for station in eachindex(receivers)
        receiver = receivers[station]
        axis = Axis(figure[station, 1]; xlabel="time (s)",
            ylabel=normalized_amplitude ? "normalized" : "displacement (m)",
            title="SPECFEM2D $(component), station $station: $(station_label(receiver))")
        for (resolution, result) in enumerate(results)
            candidate = trace(result, :SPECFEM2D, component, station)
            values = normalized_amplitude ? normalized(candidate.values) : candidate.values
            lines!(axis, candidate.time, values; color=colors[resolution],
                label="Δ≈$(round(Int, result.spacing_m)) m")
        end
        travel = arrivals(finest.source.position, receiver, 6000.0, 3464.0,
            finest.source.delay_s)
        for (name, time) in pairs(travel)
            vlines!(axis, [time]; color=name in (:P, :S) ? :black : :gray50,
                linestyle=name in (:P, :S) ? :solid : :dash, linewidth=0.8)
        end
        station == 1 && axislegend(axis; position=:rt, nbanks=3)
    end
    figure
end

function plot_specfem_station_convergence(results)
    finest = last(results)
    spacings = getproperty.(results[1:end-1], :spacing_m)
    figure = Figure(size=(1100, 950))
    for (column, component) in enumerate((:x, :z))
        axis = Axis(figure[1, column]; xscale=log10, yscale=log10,
            xlabel="mean GLL spacing (m)", ylabel="relative RMS to finest",
            title="SPECFEM2D u$(component), nonzero physical traces")
        null_axis = Axis(figure[2, column]; xscale=log10, yscale=log10,
            xlabel="mean GLL spacing (m)", ylabel="absolute RMS displacement (m)",
            title="SPECFEM2D u$(component), symmetry/null traces")
        reference_norms = [norm(trace(finest, :SPECFEM2D, component,
            station).values) for station in eachindex(finest.receivers)]
        null_threshold = maximum(reference_norms) * 1e-6
        for station in eachindex(finest.receivers)
            reference = trace(finest, :SPECFEM2D, component, station)
            common_times = collect(range(maximum(first(trace(result,
                :SPECFEM2D, component, station).time) for result in results),
                minimum(last(trace(result, :SPECFEM2D, component, station).time)
                    for result in results); length=1501))
            ref_values = linear_sample(reference, common_times)
            if reference_norms[station] > null_threshold
                errors = [relative_rms(ref_values, linear_sample(trace(result,
                    :SPECFEM2D, component, station), common_times))
                    for result in results[1:end-1]]
                valid = isfinite.(errors) .& (errors .> 0)
                any(valid) && scatterlines!(axis, spacings[valid], errors[valid];
                    label="station $station", markersize=7)
            else
                errors = [sqrt(mean(abs2, linear_sample(trace(result,
                    :SPECFEM2D, component, station), common_times)))
                    for result in results[1:end-1]]
                valid = isfinite.(errors) .& (errors .> 0)
                any(valid) && scatterlines!(null_axis, spacings[valid], errors[valid];
                    label="station $station", markersize=7)
            end
        end
        axislegend(axis; position=:rb, nbanks=2, labelsize=9)
        !isempty(null_axis.scene.plots) && axislegend(null_axis;
            position=:rb, nbanks=2, labelsize=9)
    end
    figure
end

function plot_cross_method_finest(results)
    result = last(results)
    figure = Figure(size=(1200, 230 * length(result.receivers)))
    colors = Dict(:FD3 => :dodgerblue, :OPT3 => :darkorange,
        :SPECFEM2D => :black)
    for (station, receiver) in pairs(result.receivers)
        axis = Axis(figure[station, 1]; xlabel="time (s)",
            ylabel="normalized displacement",
            title="finest grid, u_z, station $station: $(station_label(receiver))")
        for method in (:FD3, :OPT3, :SPECFEM2D)
            candidate = trace(result, method, :z, station)
            lines!(axis, candidate.time, normalized(candidate.values);
                color=colors[method], label=String(method))
        end
        travel = arrivals(result.source.position, receiver, 6000.0, 3464.0,
            result.source.delay_s)
        for (name, time) in pairs(travel)
            vlines!(axis, [time]; color=name in (:P, :S) ? :red : :gray55,
                linestyle=name in (:P, :S) ? :solid : :dash, linewidth=0.8)
        end
        station == 1 && axislegend(axis; position=:rt)
    end
    figure
end

function plot_absolute_all_methods(results; component=:z)
    receivers = last(results).receivers
    methods = (:FD3, :OPT3, :SPECFEM2D)
    colors = Dict(:FD3 => :dodgerblue, :OPT3 => :darkorange,
        :SPECFEM2D => :black)
    figure = Figure(size=(1450, 225 * length(receivers)))
    for (row, receiver) in pairs(receivers), (column, result) in pairs(results)
        axis = Axis(figure[row, column]; xlabel="time (s)",
            ylabel="u$(component) (m)",
            title="station $row, $(round(Int, result.spacing_m)) m: " *
                station_label(receiver))
        for method in methods
            candidate = trace(result, method, component, row)
            lines!(axis, candidate.time, candidate.values;
                color=colors[method], label=String(method))
        end
        travel = arrivals(result.source.position, receiver, 6000.0, 3464.0,
            result.source.delay_s)
        for (name, time) in pairs(travel)
            vlines!(axis, [time]; color=name in (:P, :S) ? :red : :gray55,
                linestyle=name in (:P, :S) ? :solid : :dash, linewidth=0.7)
        end
        row == 1 && column == 1 && axislegend(axis; position=:rt)
    end
    figure
end

function plot_late_absolute_finest(results; component=:z)
    result = last(results)
    methods = (:FD3, :OPT3, :SPECFEM2D)
    colors = Dict(:FD3 => :dodgerblue, :OPT3 => :darkorange,
        :SPECFEM2D => :black)
    figure = Figure(size=(1200, 225 * length(result.receivers)))
    for (station, receiver) in pairs(result.receivers)
        travel = arrivals(result.source.position, receiver, 6000.0, 3464.0,
            result.source.delay_s)
        late_start = max(0.0, min(travel.PP, travel.PS, travel.SP, travel.SS) -
            0.5 / result.source.frequency_hz)
        axis = Axis(figure[station, 1]; xlabel="time (s)",
            ylabel="u$(component) (m)",
            title="late phases, station $station: $(station_label(receiver))",
            limits=(late_start, result.duration_s, nothing, nothing))
        for method in methods
            candidate = trace(result, method, component, station)
            lines!(axis, candidate.time, candidate.values;
                color=colors[method], label=String(method))
        end
        for (name, time) in pairs(travel)
            name in (:PP, :PS, :SP, :SS) || continue
            vlines!(axis, [time]; color=:gray45, linestyle=:dash, linewidth=0.8)
        end
        station == 1 && axislegend(axis; position=:rt)
    end
    figure
end

function global_self_error(reference_result, candidate_result, method, component)
    reference_values = Float64[]
    candidate_values = Float64[]
    for station in eachindex(reference_result.receivers)
        reference = trace(reference_result, method, component, station)
        candidate = trace(candidate_result, method, component, station)
        common_times = collect(range(max(first(reference.time), first(candidate.time)),
            min(last(reference.time), last(candidate.time)); length=1501))
        append!(reference_values, linear_sample(reference, common_times))
        append!(candidate_values, linear_sample(candidate, common_times))
    end
    relative_rms(reference_values, candidate_values)
end

function plot_coupled_space_time_convergence(results)
    finest = last(results)
    figure = Figure(size=(1200, 720))
    for (column, component) in enumerate((:x, :z))
        axis = Axis(figure[1, column]; xscale=log10, yscale=log10,
            xlabel="spatial spacing Δx or mean GLL interval (m)",
            ylabel="global relative RMS to finest",
            title="coupled Δx–Δt convergence, u$(component)")
        for (method, color, marker) in ((:FD3, :dodgerblue, :circle),
                (:OPT3, :darkorange, :rect), (:SPECFEM2D, :black, :diamond))
            candidates = results[1:end-1]
            spacings = [method === :SPECFEM2D ?
                result.SPECFEM2D.mean_gll_interval_m.x : result.spacing_m
                for result in candidates]
            errors = [global_self_error(finest, result, method, component)
                for result in candidates]
            scatterlines!(axis, spacings, errors; color, marker,
                label=String(method))
            for (spacing, error, result) in zip(spacings, errors, candidates)
                dt = getproperty(result, method).dt
                text!(axis, spacing, error;
                    text="  Δt=$(round(dt; sigdigits=3)) s",
                    color, fontsize=10, align=(:left, :bottom))
            end
        end
        axislegend(axis; position=:rb)
    end
    Label(figure[0, :],
        "Δx and Δt are refined together in the saved runs (not an isolated temporal test)",
        fontsize=16)
    figure
end

function plot_fullspace_green_window(results)
    result = last(results)
    selected = findall(receiver -> receiver.z < 0 &&
        hypot(receiver.x-result.source.position.x,
            receiver.z-result.source.position.z) >= 5e3, result.receivers)
    receivers = result.receivers[selected]
    times = collect(0.0:0.02:result.duration_s)
    frequency = result.source.frequency_hz
    delay = result.source.delay_s
    ricker(time) = begin
        a = π * frequency * (time - delay)
        (1 - 2a^2) * exp(-a^2)
    end
    green = aki_richards_line_force_2d(times, receivers;
        source=result.source.position, vp=6000.0, vs=3464.0, rho=2700.0,
        force_time_function=ricker, force_amplitude=result.source.force,
        source_frequency=frequency, out_of_plane_step=250.0)
    figure = Figure(size=(1200, 230 * length(selected)))
    for (row, station) in pairs(selected)
        receiver = result.receivers[station]
        axis = Axis(figure[row, 1]; xlabel="time (s)",
            ylabel="normalized displacement",
            title="pre-surface Green2D gate, station $station: $(station_label(receiver))")
        reflected = arrivals(result.source.position, receiver, 6000.0, 3464.0,
            delay).PP
        for (component_index, component, color) in ((1, :x, :dodgerblue),
                (2, :z, :darkorange))
            numerical = trace(result, :SPECFEM2D, component, station)
            mask = numerical.time .<= reflected
            lines!(axis, numerical.time[mask], normalized(numerical.values[mask]);
                color, label="SPECFEM u$(component)")
            analytic_values = green.displacement[:, row, component_index]
            analytic_mask = green.time .<= reflected
            lines!(axis, green.time[analytic_mask], normalized(
                analytic_values[analytic_mask]); color, linestyle=:dash,
                label="Green2D u$(component)")
        end
        vlines!(axis, [reflected]; color=:black, linestyle=:dot,
            label="earliest PP")
        row == 1 && axislegend(axis; position=:rt, nbanks=3)
    end
    figure
end

function main()
    mkpath(DATA_DIRECTORY)
    results = load_results()
    products = (
        specfem_ux_absolute=plot_specfem_waveforms(results;
            component=:x, normalized_amplitude=false),
        specfem_uz_absolute=plot_specfem_waveforms(results;
            component=:z, normalized_amplitude=false),
        specfem_ux_normalized=plot_specfem_waveforms(results;
            component=:x, normalized_amplitude=true),
        specfem_uz_normalized=plot_specfem_waveforms(results;
            component=:z, normalized_amplitude=true),
        specfem_station_convergence=plot_specfem_station_convergence(results),
        finest_cross_method=plot_cross_method_finest(results),
        fullspace_green_gate=plot_fullspace_green_window(results),
        all_methods_ux_absolute=plot_absolute_all_methods(results; component=:x),
        all_methods_uz_absolute=plot_absolute_all_methods(results; component=:z),
        late_ux_absolute=plot_late_absolute_finest(results; component=:x),
        late_uz_absolute=plot_late_absolute_finest(results; component=:z),
        coupled_space_time_convergence=plot_coupled_space_time_convergence(results),
    )
    paths = String[]
    for (name, figure) in pairs(products)
        path = joinpath(DATA_DIRECTORY, "$(name).png")
        save(path, figure)
        push!(paths, path)
    end
    println.(paths)
    products
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
