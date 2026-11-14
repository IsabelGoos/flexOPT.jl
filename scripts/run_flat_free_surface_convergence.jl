#!/usr/bin/env julia

using JSON
using JLD2
using Dates
using Statistics

const ROOT = normpath(joinpath(@__DIR__, ".."))
const NOTEBOOK = joinpath(ROOT, "notebooks",
    "HomogeneousElastic2DBenchmark_FreeSurface.ipynb")
const OUTPUT_DIRECTORY = joinpath(ROOT, "data",
    "elastic2d_convergence", "homogeneous_flat_free_surface")
const LONG_OPT5_OUTPUT_DIRECTORY = joinpath(ROOT, "data",
    "elastic2d_convergence", "flat_free_surface_long_opt5")

function parse_spacings(arguments)
    # 1000 m violates the notebook's minimum S-wavelength sampling for the
    # fixed 0.6495 Hz source. Keep the source unchanged and use only
    # physically admissible resolutions spanning a factor of two.
    isempty(arguments) && return [500.0, 375.0, 250.0]
    values = parse.(Float64, arguments)
    all(>(0), values) || error("all spacings must be positive")
    sort!(unique(values); rev=true)
end

function notebook_cell(notebook, id)
    index = findfirst(cell -> get(cell, "id", "") == id, notebook["cells"])
    isnothing(index) && error("notebook cell '$id' was not found")
    join(notebook["cells"][index]["source"])
end

function compact_trace(trace)
    (time=Float64.(trace.time), values=Float64.(trace.values))
end

function run_worker(spacing; long_opt5=false)
    ENV["FLEXOPT_BENCHMARK_DX"] = string(spacing)
    if long_opt5
        ENV["FLEXOPT_BENCHMARK_QUICK"] = "false"
        ENV["FLEXOPT_BENCHMARK_DURATION"] = "14.0"
        ENV["FLEXOPT_DOMAIN_HALF_WIDTH"] = "90000.0"
        ENV["FLEXOPT_BUILD_HIGHER_ORDER"] = "true"
        ENV["FLEXOPT_RUN_OPT5"] = "true"
        ENV["FLEXOPT_CASE_SUFFIX"] = "_long14s_opt5"
    end
    ENV["FLEXOPT_RUN_TAKEUCHI_A90"] = "false"
    ENV["FLEXOPT_LOAD_SPECFEM_WAVEFIELD"] = "false"
    case_prefix = "homogeneous_flat_free_surface"
    case_spacing = long_opt5 ? spacing / 2 : spacing
    case_suffix = long_opt5 ? "_long14s_opt5" : ""
    specfem_case = joinpath(ROOT, "data", "specfem2d_benchmarks",
        "$(case_prefix)_dx$(round(Int, case_spacing))m$(case_suffix)")
    # A completed case is deterministic for this fixed benchmark and can be
    # reused after an extraction/viewer failure.
    ENV["FLEXOPT_RUN_SPECFEM2D"] =
        isfile(joinpath(specfem_case, "solver.log")) ? "false" : "true"

    notebook = JSON.parsefile(NOTEBOOK)
    # These are the computational cells only. Plotting/error/video cells are
    # deliberately excluded; the saved product is a small viewer data set.
    cell_ids = (
        "bootstrap", "configuration", "fd-propagation", "opt-recipes",
        "opt-operator", "opt-propagation", "source-plot", "traces",
    )
    for id in cell_ids
        @info "executing convergence cell" spacing id
        source = notebook_cell(notebook, id)
        # The worker already starts with --project=ROOT. Re-activating from
        # the notebook writes Pkg's global usage log and is both redundant
        # and unfriendly to batch/sandbox execution.
        id == "bootstrap" && (source = replace(source,
            "Pkg.activate(flexopt_root)" => "nothing # project set at startup"))
        Base.include_string(Main, source,
            "$(basename(NOTEBOOK)):$id")
    end

    @info "executing SPECFEM2D reference" spacing
    specfem_wall_time = @elapsed Base.include_string(Main,
        notebook_cell(notebook, "specfem"), "$(basename(NOTEBOOK)):specfem")

    # Julia 1.12 gives definitions introduced by include_string a newer world
    # age than this driver function. Build the collector in that new world and
    # enter it once through invokelatest.
    collector = Core.eval(Main, quote
        function _collect_flat_surface_convergence(specfem_wall_time)
            sampler = waveformQuantity === :velocity ?
                velocity_traces_at_points : displacement_traces_at_points
            fd_x = sampler(uxFD, fdTimes,
                fdCoordinates.x, fdCoordinates.z, receiverGrid)
            fd_z = sampler(uzFD, fdTimes,
                fdCoordinates.x, fdCoordinates.z, receiverGrid)
            opt_x = sampler(uxOPT, optTimes, xOPT, zOPT, receiverGrid)
            opt_z = sampler(uzOPT, optTimes, xOPT, zOPT, receiverGrid)
            shifted_x = sampler(uxOPTShifted, optShiftedTimes,
                xOPT, zOPT, receiverGrid)
            shifted_z = sampler(uzOPTShifted, optShiftedTimes,
                xOPT, zOPT, receiverGrid)
            result = (
                schema_version=1,
                created_at=string(now()),
                spacing_m=Float64(dx),
                source=(position=sourcePosition, frequency_hz=sourceFrequency,
                    delay_s=sourceDelay, force=sourceForce,
                    spatial_sigma_m=sourceSpatialSigma),
                receivers=collect(receiverGrid),
                duration_s=duration,
                waveform_quantity=waveformQuantity,
                FD3=(x=compact_trace(fd_x), z=compact_trace(fd_z),
                    dt=Float64(fd.dt)),
                OPT3=(x=compact_trace(opt_x), z=compact_trace(opt_z),
                    dt=Float64(dtOPT), timing=propagationOPT.timing),
                OPT3_shifted_full=(x=compact_trace(shifted_x),
                    z=compact_trace(shifted_z), dt=Float64(dtOPT),
                    timing=propagationOPTShifted.timing),
                SPECFEM2D=(
                    x=[compact_trace(trace) for trace in specfemGridWaveformsX],
                    z=[compact_trace(trace) for trace in specfemGridWaveformsZ],
                    dt=Float64(dtSPECFEM),
                    mean_gll_interval_m=
                        discretizationSummary.SPECFEM2D.mean_GLL_interval,
                    reported_cfl=specfemCFL,
                    wall_time_s=specfem_wall_time,
                    case_directory=specfemCase.case_directory,
                ),
            )
            if runOPT5
                opt5_x = sampler(uxOPT5, opt5Times, xOPT, zOPT, receiverGrid)
                opt5_z = sampler(uzOPT5, opt5Times, xOPT, zOPT, receiverGrid)
                result = merge(result, (OPT5=(x=compact_trace(opt5_x),
                    z=compact_trace(opt5_z), dt=Float64(dtOPT),
                    timing=propagationOPT5.timing),))
            end
            return result
        end
    end)
    result = Base.invokelatest(collector, specfem_wall_time)

    output_directory = long_opt5 ? LONG_OPT5_OUTPUT_DIRECTORY : OUTPUT_DIRECTORY
    mkpath(output_directory)
    output = joinpath(output_directory,
        "convergence_dx$(round(Int, result.spacing_m))m.jld2")
    jldsave(output; result)
    @info "saved convergence result" output
    return output
end

function velocity_to_displacement(trace)
    times = Float64.(trace.time)
    values = Float64.(trace.values)
    length(times) >= 2 || error("velocity trace needs at least two samples")
    dt = median(diff(times))
    edge_times = [max(0.0, first(times) - dt / 2); times .+ dt / 2]
    displacement = vcat(zeros(1, size(values, 2)),
        cumsum(values .* dt; dims=1))
    return (time=edge_times, values=displacement)
end

function repair_existing_results()
    files = isdir(OUTPUT_DIRECTORY) ? filter(
        file -> endswith(file, ".jld2"),
        readdir(OUTPUT_DIRECTORY; join=true)) : String[]
    isempty(files) && error("no convergence files found in $OUTPUT_DIRECTORY")
    for file in files
        result = load(file, "result")
        result.waveform_quantity === :velocity || continue
        repaired = merge(result, (
            waveform_quantity=:displacement,
            FD3=merge(result.FD3, (
                x=velocity_to_displacement(result.FD3.x),
                z=velocity_to_displacement(result.FD3.z),
            )),
            OPT3=merge(result.OPT3, (
                x=velocity_to_displacement(result.OPT3.x),
                z=velocity_to_displacement(result.OPT3.z),
            )),
        ))
        jldsave(file; result=repaired)
        @info "repaired velocity/displacement mismatch" file
    end
end

function launch(spacings)
    mkpath(OUTPUT_DIRECTORY)
    for spacing in spacings
        command = `$(Base.julia_cmd()) --project=$(ROOT) --startup-file=no --threads=8 $(@__FILE__) --worker $(spacing)`
        @info "launching isolated convergence worker" spacing command
        run(command)
    end
    println("Convergence data are ready in $OUTPUT_DIRECTORY")
end

if !isempty(ARGS) && first(ARGS) == "--worker"
    length(ARGS) == 2 || error("usage: --worker spacing_m")
    run_worker(parse(Float64, ARGS[2]))
elseif !isempty(ARGS) && first(ARGS) == "--long-opt5-worker"
    length(ARGS) == 2 || error("usage: --long-opt5-worker nominal_spacing_m")
    run_worker(parse(Float64, ARGS[2]); long_opt5=true)
elseif !isempty(ARGS) && first(ARGS) == "--long-opt5"
    spacings = isempty(ARGS[2:end]) ? [1000.0, 750.0, 500.0] :
        parse.(Float64, ARGS[2:end])
    mkpath(LONG_OPT5_OUTPUT_DIRECTORY)
    for spacing in sort(unique(spacings); rev=true)
        command = `$(Base.julia_cmd()) --project=$(ROOT) --startup-file=no --threads=8 $(@__FILE__) --long-opt5-worker $(spacing)`
        @info "launching long OPT3/OPT5 worker" spacing command
        run(command)
    end
elseif ARGS == ["--repair-existing"]
    repair_existing_results()
else
    launch(parse_spacings(ARGS))
end
