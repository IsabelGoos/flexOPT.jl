module specfemBenchmark

using CairoMakie

export specfem2d_root, specfem2d_status
export write_specfem2d_tomography, write_specfem2d_interfaces
export prepare_specfem2d_case
export run_specfem2d_case, read_specfem2d_trace
export find_specfem2d_traces
export waveform_metrics, plot_solver_benchmark

const DEFAULT_SPECFEM2D_ROOT =
    normpath(joinpath(@__DIR__, "..", "..", "othersPackages", "specfem2d"))

"""
    specfem2d_root(; root=nothing)

Resolve SPECFEM2D without making it a flexOPT dependency. An explicit `root`
wins, followed by `ENV["SPECFEM2D_ROOT"]`, then the conventional sibling
checkout `Github/othersPackages/specfem2d`.
"""
function specfem2d_root(; root=nothing)
    candidate = if !isnothing(root)
        String(root)
    elseif haskey(ENV, "SPECFEM2D_ROOT")
        ENV["SPECFEM2D_ROOT"]
    else
        DEFAULT_SPECFEM2D_ROOT
    end
    path = abspath(expanduser(candidate))
    isdir(path) || throw(ArgumentError(
        "SPECFEM2D was not found at $path. Set ENV[\"SPECFEM2D_ROOT\"].",
    ))
    path
end

function specfem2d_status(; root=nothing)
    path = specfem2d_root(; root)
    mesher = joinpath(path, "bin", "xmeshfem2D")
    solver = joinpath(path, "bin", "xspecfem2D")
    (
        root=path,
        mesher=mesher,
        solver=solver,
        mesher_ready=isfile(mesher) && isexecutable(mesher),
        solver_ready=isfile(solver) && isexecutable(solver),
    )
end

function _regular_axis(values, name)
    axis = Float64.(collect(values))
    length(axis) >= 2 || throw(ArgumentError("$name needs at least two points"))
    all(diff(axis) .> 0) || throw(ArgumentError("$name must increase"))
    spacing = diff(axis)
    isapprox(extrema(spacing)...; rtol=1e-8, atol=eps(maximum(abs, axis))) ||
        throw(ArgumentError("$name must be regularly spaced for SPECFEM tomography"))
    axis
end

"""
    write_specfem2d_tomography(path, x, z, vp, vs, rho)

Write the documented SPECFEM2D ASCII tomography format. Arrays use flexOPT's
`(x,z)` order and SI units: coordinates in metres, velocities in m/s and
density in kg/m³.
"""
function write_specfem2d_tomography(path, x, z, vp, vs, rho)
    xaxis = _regular_axis(x, "x")
    zaxis = _regular_axis(z, "z")
    expected = (length(xaxis), length(zaxis))
    for (name, field) in ((:vp, vp), (:vs, vs), (:rho, rho))
        size(field) == expected ||
            throw(DimensionMismatch("$name has size $(size(field)); expected $expected"))
        all(isfinite, field) || throw(ArgumentError("$name contains non-finite values"))
    end
    minimum(vp) > 0 || throw(ArgumentError("vp must be positive"))
    minimum(vs) > 0 || throw(ArgumentError("vs must be positive in an elastic case"))
    minimum(rho) > 0 || throw(ArgumentError("rho must be positive"))

    target = abspath(path)
    mkpath(dirname(target))
    open(target, "w") do io
        println(io, "# flexOPT → SPECFEM2D tomography")
        println(io, first(xaxis), " ", first(zaxis), " ",
                last(xaxis), " ", last(zaxis))
        println(io, xaxis[2] - xaxis[1], " ", zaxis[2] - zaxis[1])
        println(io, length(xaxis), " ", length(zaxis))
        println(io, minimum(vp), " ", maximum(vp), " ",
                minimum(vs), " ", maximum(vs), " ",
                minimum(rho), " ", maximum(rho))
        # SPECFEM2D expects x to vary fastest.
        for iz in eachindex(zaxis), ix in eachindex(xaxis)
            println(io, xaxis[ix], " ", zaxis[iz], " ",
                    vp[ix, iz], " ", vs[ix, iz], " ", rho[ix, iz])
        end
    end
    target
end

"""
    write_specfem2d_interfaces(path, x, surface_z; bottom_z, nz_elements)

Write a two-interface mesh (bottom and free surface). `surface_z` may be a
scalar for the flat benchmark or one value per x point for Kirishima.
"""
function write_specfem2d_interfaces(
    path,
    x,
    surface_z;
    bottom_z,
    nz_elements::Integer,
)
    xaxis = Float64.(collect(x))
    surface = surface_z isa Real ?
        fill(Float64(surface_z), length(xaxis)) :
        Float64.(collect(surface_z))
    length(surface) == length(xaxis) ||
        throw(DimensionMismatch("surface_z must have one value per x point"))
    all(surface .> bottom_z) ||
        throw(ArgumentError("the free surface must lie above bottom_z"))
    nz_elements > 0 || throw(ArgumentError("nz_elements must be positive"))

    target = abspath(path)
    mkpath(dirname(target))
    open(target, "w") do io
        println(io, "2")
        println(io, "2")
        println(io, first(xaxis), " ", bottom_z)
        println(io, last(xaxis), " ", bottom_z)
        println(io, length(xaxis))
        for i in eachindex(xaxis)
            println(io, xaxis[i], " ", surface[i])
        end
        println(io, nz_elements)
    end
    target
end

function _set_parameter(text, name, value)
    pattern = Regex("(?m)^\\s*" * name * "\\s*=.*\$")
    occursin(pattern, text) ||
        throw(ArgumentError("parameter $name was not found in SPECFEM Par_file"))
    replace(text, pattern => "$name = $value")
end

"""
    prepare_specfem2d_case(case_directory, x, z, vp, vs, rho, surface_z; ...)

Create an isolated serial SPECFEM2D elastic case from flexOPT arrays. The
official checkout is used only as a source of executables and a validated
`Par_file` template; it is never modified.
"""
function prepare_specfem2d_case(
    case_directory,
    x,
    z,
    vp,
    vs,
    rho,
    surface_z;
    root=nothing,
    source=(x=0.0, z=-2_000.0),
    receivers=range(first(x), last(x); length=5),
    duration=30.0,
    dt=0.001,
    f0=1.0,
    source_factor=1.0e10,
    source_angle=0.0,
    nx_elements=max(4, cld(length(x) - 1, 4)),
    nz_elements=max(4, cld(length(z) - 1, 4)),
)
    status = specfem2d_status(; root)
    template = joinpath(
        status.root, "EXAMPLES", "benchmarks",
        "semi_infinite_homogeneous", "DATA",
    )
    case_path = abspath(case_directory)
    data_path = joinpath(case_path, "DATA")
    output_path = joinpath(case_path, "OUTPUT_FILES")
    mkpath(data_path)
    mkpath(output_path)

    tomography_file = write_specfem2d_tomography(
        joinpath(data_path, "kirishima_tomography.xyz"),
        x, z, vp, vs, rho,
    )
    interfaces_file = write_specfem2d_interfaces(
        joinpath(data_path, "interfaces.dat"),
        x, surface_z; bottom_z=first(z), nz_elements,
    )

    par = read(joinpath(template, "Par_file"), String)
    settings = (
        ("title", "flexOPT FD3 OPT3 SPECFEM2D benchmark"),
        ("NPROC", "1"),
        ("NSTEP", string(ceil(Int, duration / dt))),
        ("DT", string(dt)),
        ("MODEL", "default"),
        ("TOMOGRAPHY_FILE", "./DATA/$(basename(tomography_file))"),
        ("seismotype", "2"),
        ("NTSTEP_BETWEEN_OUTPUT_SAMPLE", "1"),
        ("USER_T0", "0.d0"),
        ("nreceiversets", "1"),
        ("nrec", string(length(receivers))),
        ("xdeb", string(first(receivers))),
        ("zdeb", string(maximum(surface_z isa Real ? [surface_z] : surface_z))),
        ("xfin", string(last(receivers))),
        ("zfin", string(maximum(surface_z isa Real ? [surface_z] : surface_z))),
        ("record_at_surface_same_vertical", ".true."),
        ("interfacesfile", basename(interfaces_file)),
        ("xmin", string(first(x))),
        ("xmax", string(last(x))),
        ("nx", string(nx_elements)),
        ("absorbbottom", ".true."),
        ("absorbright", ".true."),
        ("absorbtop", ".false."),
        ("absorbleft", ".true."),
        ("nbmodels", "1"),
        ("nbregions", "1"),
    )
    for (name, value) in settings
        par = _set_parameter(par, name, value)
    end
    par = replace(
        par,
        r"(?m)^1 1 2700\.d0 3000\.d0 1732\.05d0.*$" =>
            # The fifth value is a positive Vs marker: it tells the mesher
            # that the tomography region is elastic rather than acoustic.
            "1 -1 0 0 1 0 0 0 0 0 0 0 0 0 0",
        r"(?m)^1 50 1\s+50 1\s*$" =>
            "1 $nx_elements 1 $nz_elements 1",
    )
    write(joinpath(data_path, "Par_file"), par)

    source_text = read(joinpath(template, "SOURCE"), String)
    source_text = _set_parameter(source_text, "xs", source.x)
    source_text = _set_parameter(source_text, "zs", source.z)
    source_text = _set_parameter(source_text, "f0", f0)
    source_text = _set_parameter(source_text, "source_type", 1)
    source_text = _set_parameter(source_text, "anglesource", source_angle)
    source_text = _set_parameter(source_text, "factor", source_factor)
    write(joinpath(data_path, "SOURCE"), source_text)

    (
        case_directory=case_path,
        data=data_path,
        output=output_path,
        tomography=tomography_file,
        interfaces=interfaces_file,
        par_file=joinpath(data_path, "Par_file"),
        source_file=joinpath(data_path, "SOURCE"),
        source_factor=Float64(source_factor),
        source_angle=Float64(source_angle),
        surface=surface_z,
        nx_elements,
        nz_elements,
    )
end

"""
    run_specfem2d_case(case_directory; root=nothing)

Run the mesher and solver inside an already prepared SPECFEM2D case directory.
The external installation is never modified.
"""
function run_specfem2d_case(case_directory; root=nothing)
    status = specfem2d_status(; root)
    status.mesher_ready && status.solver_ready ||
        error("SPECFEM2D is present but not compiled: $(status.root)")
    case_path = abspath(case_directory)
    isfile(joinpath(case_path, "DATA", "Par_file")) ||
        throw(ArgumentError("$case_path does not contain DATA/Par_file"))
    mesher_log = joinpath(case_path, "mesh.log")
    solver_log = joinpath(case_path, "solver.log")
    open(mesher_log, "w") do io
        command = Cmd(Cmd([status.mesher]); dir=case_path)
        run(pipeline(command, stdout=io, stderr=io))
    end
    open(solver_log, "w") do io
        command = Cmd(Cmd([status.solver]); dir=case_path)
        run(pipeline(command, stdout=io, stderr=io))
    end
    (; case_directory=case_path, output=joinpath(case_path, "OUTPUT_FILES"),
       mesher_log, solver_log)
end

function read_specfem2d_trace(path)
    rows = readlines(path)
    values = Tuple{Float64,Float64}[]
    for row in rows
        stripped = strip(row)
        (isempty(stripped) || startswith(stripped, "#")) && continue
        columns = split(stripped)
        length(columns) >= 2 || continue
        push!(values, (parse(Float64, columns[1]), parse(Float64, columns[2])))
    end
    isempty(values) && throw(ArgumentError("no trace samples found in $path"))
    (; time=first.(values), values=last.(values), path=abspath(path))
end

"""
    find_specfem2d_traces(output_directory; component=:z)

Return one SPECFEM2D ASCII velocity trace per station. SPECFEM uses different
channel prefixes depending on the output convention (`CXZ`, `FXZ`, ...), so
selection is based on the physical final component rather than a fixed prefix.
"""
function find_specfem2d_traces(output_directory; component::Symbol=:z)
    component in (:x, :z) ||
        throw(ArgumentError("SPECFEM2D component must be :x or :z"))
    suffix = component === :x ? "XX.semv" : "XZ.semv"
    paths = filter(
        path -> endswith(basename(path), suffix),
        readdir(output_directory; join=true),
    )
    sort!(paths; by=path -> begin
        fields = split(basename(path), '.')
        length(fields) >= 4 ? fields[2] : basename(path)
    end)
    paths
end

function _linear_sample(time, values, query)
    query <= first(time) && return first(values)
    query >= last(time) && return last(values)
    i = searchsortedlast(time, query)
    α = (query - time[i]) / (time[i + 1] - time[i])
    (1 - α) * values[i] + α * values[i + 1]
end

function waveform_metrics(reference, candidate; samples=2001)
    start_time = max(first(reference.time), first(candidate.time))
    end_time = min(last(reference.time), last(candidate.time))
    start_time < end_time || throw(ArgumentError("traces do not overlap in time"))
    time = collect(range(start_time, end_time; length=samples))
    a = [_linear_sample(reference.time, reference.values, t) for t in time]
    b = [_linear_sample(candidate.time, candidate.values, t) for t in time]
    a .-= sum(a) / length(a)
    b .-= sum(b) / length(b)
    denom = sqrt(sum(abs2, a) * sum(abs2, b))
    correlation = denom == 0 ? NaN : sum(a .* b) / denom
    amplitude = sum(abs2, b) == 0 ? NaN : sum(a .* b) / sum(abs2, b)
    relative_error = sum(abs2, a) == 0 ? NaN :
        sqrt(sum(abs2, a .- amplitude .* b) / sum(abs2, a))
    (; correlation, relative_error, optimal_amplitude=amplitude, time)
end

function plot_solver_benchmark(traces; normalize=true, title="Elastic 2D benchmark")
    figure = Figure(size=(1100, 650))
    axis = Axis(figure[1, 1]; xlabel="time (s)",
                ylabel=normalize ? "normalized amplitude" : "amplitude",
                title)
    colors = Makie.wong_colors()
    for (i, (name, trace)) in enumerate(pairs(traces))
        values = Float64.(trace.values)
        scale = normalize ? max(maximum(abs, values), eps(Float64)) : 1.0
        lines!(axis, trace.time, values ./ scale;
               label=String(name), color=colors[mod1(i, length(colors))])
    end
    axislegend(axis)
    (; figure, axis)
end

end
