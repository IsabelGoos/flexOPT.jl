module specfemBenchmark

using CairoMakie

export specfem2d_root, specfem2d_status
export write_specfem2d_tomography, write_specfem2d_interfaces
export prepare_specfem2d_case
export run_specfem2d_case, read_specfem2d_trace
export find_specfem2d_traces
export read_specfem2d_wavefield_dumps
export make_specfem2d_snapshot_video
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
    sources=nothing,
    receivers=range(first(x), last(x); length=5),
    receiver_points=nothing,
    duration=30.0,
    dt=0.001,
    f0=1.0,
    source_factor=1.0e10,
    source_angle=0.0,
    source_type::Integer=1,
    moment_tensor=(Mxx=1.0, Mzz=-1.0, Mxz=0.0),
    source_time_function=nothing,
    receiver_z=maximum(surface_z isa Real ? [surface_z] : surface_z),
    free_surface=true,
    record_at_surface_same_vertical=free_surface,
    nx_elements=max(4, cld(length(x) - 1, 4)),
    nz_elements=max(4, cld(length(z) - 1, 4)),
    snapshot_interval_steps=100,
    snapshot_image_type=5,
    output_wavefield_dumps=false,
    wavefield_dump_type=1,
    binary_wavefield_dumps=true,
)
    source_type in (1, 2) ||
        throw(ArgumentError("source_type must be 1 (force) or 2 (moment tensor)"))
    all(name -> hasproperty(moment_tensor, name), (:Mxx, :Mzz, :Mxz)) ||
        throw(ArgumentError("moment_tensor needs Mxx, Mzz and Mxz fields"))
    all(isfinite, Float64.((moment_tensor.Mxx, moment_tensor.Mzz,
        moment_tensor.Mxz))) ||
        throw(ArgumentError("moment-tensor components must be finite"))
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
    receiver_specs = isnothing(receiver_points) ?
        [(x=Float64(xr), z=Float64(receiver_z)) for xr in receivers] :
        [(x=Float64(point.x), z=Float64(point.z)) for point in receiver_points]
    isempty(receiver_specs) && throw(ArgumentError("receiver list cannot be empty"))
    all(point -> isfinite(point.x) && isfinite(point.z), receiver_specs) ||
        throw(ArgumentError("receiver coordinates must be finite"))
    stations_file = if isnothing(receiver_points)
        nothing
    else
        path = joinpath(data_path, "STATIONS")
        open(path, "w") do io
            for (index, point) in enumerate(receiver_specs)
                println(io, "S", lpad(index, 4, '0'), " FX ",
                    point.x, " ", point.z, " 0.0 0.0")
            end
        end
        path
    end
    settings = (
        ("title", "flexOPT FD3 OPT3 SPECFEM2D benchmark"),
        ("NPROC", "1"),
        ("NSTEP", string(ceil(Int, duration / dt))),
        ("NSOURCES", string(isnothing(sources) ? 1 : length(sources))),
        ("DT", string(dt)),
        ("MODEL", "default"),
        ("TOMOGRAPHY_FILE", "./DATA/$(basename(tomography_file))"),
        ("seismotype", "2"),
        ("NTSTEP_BETWEEN_OUTPUT_SAMPLE", "1"),
        ("NTSTEP_BETWEEN_OUTPUT_IMAGES", string(snapshot_interval_steps)),
        ("output_color_image", ".true."),
        ("imagetype_JPEG", string(snapshot_image_type)),
        ("output_wavefield_dumps", output_wavefield_dumps ? ".true." : ".false."),
        ("imagetype_wavefield_dumps", string(wavefield_dump_type)),
        ("use_binary_for_wavefield_dumps", binary_wavefield_dumps ? ".true." : ".false."),
        ("USE_SNAPSHOT_NUMBER_IN_FILENAME", ".true."),
        ("output_postscript_snapshot", ".false."),
        ("USER_T0", "0.d0"),
        ("nreceiversets", "1"),
        ("nrec", string(length(receiver_specs))),
        ("xdeb", string(first(receiver_specs).x)),
        ("zdeb", string(receiver_z)),
        ("xfin", string(last(receiver_specs).x)),
        ("zfin", string(receiver_z)),
        ("use_existing_STATIONS", isnothing(stations_file) ? ".false." : ".true."),
        ("record_at_surface_same_vertical",
            record_at_surface_same_vertical ? ".true." : ".false."),
        ("interfacesfile", basename(interfaces_file)),
        ("xmin", string(first(x))),
        ("xmax", string(last(x))),
        ("nx", string(nx_elements)),
        ("absorbbottom", ".true."),
        ("absorbright", ".true."),
        ("absorbtop", free_surface ? ".false." : ".true."),
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

    source_template = read(joinpath(template, "SOURCE"), String)
    source_specs = isnothing(sources) ?
        [(x=source.x, z=source.z, weight=1.0)] : collect(sources)
    isempty(source_specs) && throw(ArgumentError("sources cannot be empty"))
    all(s -> hasproperty(s, :x) && hasproperty(s, :z) &&
        hasproperty(s, :weight), source_specs) || throw(ArgumentError(
        "each distributed source needs x, z and weight fields",
    ))
    all(s -> isfinite(s.x) && isfinite(s.z) && isfinite(s.weight),
        source_specs) || throw(ArgumentError(
        "distributed source coordinates and weights must be finite",
    ))
    isapprox(sum(s.weight for s in source_specs), 1.0; atol=1e-12) ||
        throw(ArgumentError("distributed source weights must sum to one"))
    source_time_function_file = nothing
    # SPECFEM labels seismograms from -1.2/f0 for an external STF as well.
    # The external array itself starts at physical time zero, so consumers
    # must add this offset to compare it with FD/OPT time axes.
    time_axis_shift = isnothing(source_time_function) ? 0.0 : 1.2 / f0
    if !isnothing(source_time_function)
        nstep = ceil(Int, duration / dt)
        times = (0:nstep-1) .* dt
        values = source_time_function isa Function ?
            Float64.(source_time_function.(times)) :
            Float64.(collect(source_time_function))
        length(values) == nstep || throw(DimensionMismatch(
            "external SPECFEM source needs $nstep samples, got $(length(values))",
        ))
        all(isfinite, values) || throw(ArgumentError(
            "external SPECFEM source contains non-finite values",
        ))
        source_time_function_file = joinpath(data_path, "source_time_function.txt")
        open(source_time_function_file, "w") do io
            for (time, value) in zip(times, values)
                println(io, time, " ", value)
            end
        end
    end
    source_blocks = map(source_specs) do spec
        block = _set_parameter(source_template, "xs", spec.x)
        block = _set_parameter(block, "zs", spec.z)
        block = _set_parameter(block, "f0", f0)
        block = _set_parameter(block, "source_type", source_type)
        block = _set_parameter(block, "anglesource", source_angle)
        block = _set_parameter(block, "Mxx", moment_tensor.Mxx)
        block = _set_parameter(block, "Mzz", moment_tensor.Mzz)
        block = _set_parameter(block, "Mxz", moment_tensor.Mxz)
        block = _set_parameter(block, "factor", source_factor * spec.weight)
        if !isnothing(source_time_function)
            block = _set_parameter(block, "time_function_type", 8)
            block = _set_parameter(
                block, "name_of_source_file",
                "./DATA/$(basename(source_time_function_file))",
            )
        end
        block
    end
    source_text = join(source_blocks, "\n")
    write(joinpath(data_path, "SOURCE"), source_text)

    (
        case_directory=case_path,
        data=data_path,
        output=output_path,
        tomography=tomography_file,
        interfaces=interfaces_file,
        par_file=joinpath(data_path, "Par_file"),
        source_file=joinpath(data_path, "SOURCE"),
        source_time_function_file,
        time_axis_shift=Float64(time_axis_shift),
        source_factor=Float64(source_factor),
        source_specs,
        source_angle=Float64(source_angle),
        receiver_points=receiver_specs,
        stations_file,
        receiver_z=Float64(receiver_z),
        free_surface=Bool(free_surface),
        surface=surface_z,
        nx_elements,
        nz_elements,
        snapshot_interval_steps=Int(snapshot_interval_steps),
        snapshot_image_type=Int(snapshot_image_type),
        output_wavefield_dumps=Bool(output_wavefield_dumps),
        wavefield_dump_type=Int(wavefield_dump_type),
        binary_wavefield_dumps=Bool(binary_wavefield_dumps),
    )
end

"""
    read_specfem2d_wavefield_dumps(output; dt, time_axis_shift=0)

Read binary P-SV vector dumps produced with `imagetype_wavefield_dumps=1`.
The result has `ux` and `uz` arrays ordered `(x, z, time)`. This helper is
restricted to the structured rectangular meshes generated by
`prepare_specfem2d_case`.
"""
function read_specfem2d_wavefield_dumps(
    output;
    dt::Real,
    time_axis_shift::Real=0.0,
)
    output_path = abspath(output)
    grid_file = joinpath(output_path, "wavefield_grid_for_dumps.bin")
    isfile(grid_file) || throw(ArgumentError(
        "SPECFEM wavefield grid is missing: $grid_file. " *
        "Rerun with output_wavefield_dumps=true.",
    ))
    grid_values = reinterpret(Float32, read(grid_file))
    iseven(length(grid_values)) || error("invalid SPECFEM wavefield grid")
    grid = reshape(grid_values, 2, :)
    # Neighboring elements can write the same conceptual GLL coordinate with
    # slightly different Float32 roundoff. Cluster those values before testing
    # the tensor-product structure.
    coordinate_tolerance = max(
        Float64(maximum(grid[1, :]) - minimum(grid[1, :])),
        Float64(maximum(grid[2, :]) - minimum(grid[2, :])),
    ) * 1e-6
    function clustered_axis(values)
        ordered = sort!(Float64.(collect(values)))
        clusters = Vector{Vector{Float64}}()
        for value in ordered
            if isempty(clusters) || value - last(last(clusters)) > coordinate_tolerance
                push!(clusters, [value])
            else
                push!(last(clusters), value)
            end
        end
        [sum(cluster) / length(cluster) for cluster in clusters]
    end
    x = clustered_axis(grid[1, :])
    z = clustered_axis(grid[2, :])
    length(x) * length(z) == size(grid, 2) || throw(ArgumentError(
        "wavefield dump is not a complete tensor grid",
    ))
    nearest_axis_index(axis, value) = begin
        upper = searchsortedfirst(axis, Float64(value))
        upper <= 1 && return 1
        upper > length(axis) && return length(axis)
        abs(axis[upper] - value) < abs(axis[upper - 1] - value) ? upper : upper - 1
    end
    grid_indices = [(nearest_axis_index(x, grid[1, point]),
        nearest_axis_index(z, grid[2, point])) for point in axes(grid, 2)]
    length(unique(grid_indices)) == size(grid, 2) || throw(ArgumentError(
        "clustered SPECFEM GLL coordinates still contain duplicate points",
    ))
    files = sort(filter(path -> occursin(r"wavefield\d+_\d+\.bin$", basename(path)),
        readdir(output_path; join=true)))
    isempty(files) && throw(ArgumentError(
        "no binary SPECFEM wavefield dumps were found in $output_path",
    ))
    ux = Array{Float32}(undef, length(x), length(z), length(files))
    uz = similar(ux)
    steps = Vector{Int}(undef, length(files))
    for (frame, file) in enumerate(files)
        values = reinterpret(Float32, read(file))
        length(values) == 2size(grid, 2) || error(
            "unexpected vector count in $file",
        )
        vector_field = reshape(values, 2, :)
        for point in axes(vector_field, 2)
            ix, iz = grid_indices[point]
            ux[ix, iz, frame] = vector_field[1, point]
            uz[ix, iz, frame] = vector_field[2, point]
        end
        matched = match(r"wavefield(\d+)_\d+\.bin$", basename(file))
        isnothing(matched) && error("cannot parse SPECFEM dump step from $file")
        steps[frame] = parse(Int, matched.captures[1])
    end
    times = (steps .- 1) .* Float64(dt) .+ Float64(time_axis_shift)
    return (; x, z, time=times, ux, uz, files)
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
    # SPECFEM does not remove seismograms from a preceding run.  Keeping those
    # files makes a subsequent receiver configuration look as if it had extra
    # stations, and stale wavefield dumps may belong to another mesh. Remove
    # only generated time-series/snapshot products; mesh databases are left
    # untouched.
    output_directory = joinpath(case_path, "OUTPUT_FILES")
    if isdir(output_directory)
        for path in readdir(output_directory; join=true)
            name = basename(path)
            generated = occursin(
                r"^[^.]+\.[^.]+\.[^.]+\.sem[avd]$", name,
            ) || occursin(r"^wavefield[0-9]+_[0-9]+\.bin$", name) ||
                occursin(r"^forward_(?:image|img)[0-9]+\.jpg$", name) ||
                occursin(r"^SPECFEM2D_.*\.(?:mp4|gif)$", name)
            generated && rm(path)
        end
    end
    mesher_wall_time_s = @elapsed open(mesher_log, "w") do io
        command = Cmd(Cmd([status.mesher]); dir=case_path)
        run(pipeline(command, stdout=io, stderr=io))
    end
    solver_wall_time_s = @elapsed open(solver_log, "w") do io
        command = Cmd(Cmd([status.solver]); dir=case_path)
        run(pipeline(command, stdout=io, stderr=io))
    end
    (; case_directory=case_path, output=joinpath(case_path, "OUTPUT_FILES"),
       mesher_log, solver_log, mesher_wall_time_s, solver_wall_time_s,
       total_wall_time_s=mesher_wall_time_s + solver_wall_time_s)
end

function read_specfem2d_trace(path; time_shift=0.0)
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
    (; time=first.(values) .+ Float64(time_shift),
       values=last.(values), path=abspath(path),
       time_shift=Float64(time_shift))
end

"""
    make_specfem2d_snapshot_video(output_directory; output_path, framerate=20)

Assemble SPECFEM2D's genuine `forward_image*.jpg` or `forward_img*.jpg`
wavefield snapshots into an
MP4. The JPEGs are generated by SPECFEM itself; this function only encodes
them and therefore never reconstructs a field from receiver seismograms.
"""
function make_specfem2d_snapshot_video(
    output_directory;
    output_path=joinpath(output_directory, "SPECFEM2D_wavefield.mp4"),
    framerate::Integer=20,
    overwrite::Bool=true,
)
    directory = abspath(output_directory)
    all_images = readdir(directory; join=true)
    # Some old cases contain both naming conventions from different runs.
    # Never concatenate the two sequences: that creates a non-chronological
    # video. Prefer SPECFEM's current `forward_image` sequence and fall back
    # to the legacy `forward_img` sequence only when necessary.
    current_images = sort(filter(path -> occursin(
        r"^forward_image[0-9]+\.jpg$", basename(path)), all_images))
    legacy_images = sort(filter(path -> occursin(
        r"^forward_img[0-9]+\.jpg$", basename(path)), all_images))
    images = isempty(current_images) ? legacy_images : current_images
    isempty(images) && throw(ArgumentError(
        "no SPECFEM2D forward_image*.jpg or forward_img*.jpg snapshots found in $directory",
    ))
    ffmpeg = Sys.which("ffmpeg")
    isnothing(ffmpeg) && throw(ArgumentError(
        "ffmpeg is required to encode SPECFEM2D snapshots",
    ))
    destination = abspath(output_path)
    mkpath(dirname(destination))
    overwrite_flag = overwrite ? "-y" : "-n"
    image_prefix = startswith(basename(first(images)), "forward_img") ?
        "forward_img" : "forward_image"
    pattern = joinpath(directory, image_prefix * "*.jpg")
    video_filter = "pad=ceil(iw/2)*2:ceil(ih/2)*2"
    ffmpeg_command = `$ffmpeg $overwrite_flag -loglevel error \
        -framerate $framerate -pattern_type glob -i $pattern \
        -vf $video_filter -c:v libx264 -pix_fmt yuv420p $destination`
    try
        run(pipeline(ffmpeg_command; stdout=devnull, stderr=devnull))
        return (; path=destination, frames=length(images), images,
            encoder=:ffmpeg)
    catch ffmpeg_error
        # Homebrew ffmpeg can temporarily be unusable after a codec-library
        # upgrade. An animated GIF still exposes the genuine SPECFEM snapshots
        # inline and requires no reconstruction of the wavefield.
        magick = Sys.which("magick")
        isnothing(magick) && rethrow(ffmpeg_error)
        gif_path = replace(destination, r"\.[^.]+$" => ".gif")
        delay_centiseconds = max(1, round(Int, 100 / framerate))
        arguments = String[magick, "-delay", string(delay_centiseconds),
            "-loop", "0"]
        append!(arguments, images)
        push!(arguments, gif_path)
        run(Cmd(arguments))
        return (; path=gif_path, frames=length(images), images,
            encoder=:imagemagick_gif, ffmpeg_error)
    end
end

"""
    find_specfem2d_traces(output_directory; component=:z)

Return one SPECFEM2D ASCII velocity trace per station. SPECFEM uses different
channel prefixes depending on the output convention (`CXZ`, `FXZ`, ...), so
selection is based on the physical final component rather than a fixed prefix.
"""
function find_specfem2d_traces(
    output_directory;
    component::Symbol=:z,
    network::Union{Nothing,AbstractString}=nothing,
    stations::Union{Nothing,AbstractVector{<:AbstractString}}=nothing,
)
    component in (:x, :z) ||
        throw(ArgumentError("SPECFEM2D component must be :x or :z"))
    suffix = component === :x ? "XX.semv" : "XZ.semv"
    paths = filter(
        path -> begin
            fields = split(basename(path), '.')
            endswith(basename(path), suffix) &&
                (isnothing(network) || (!isempty(fields) && fields[1] == network)) &&
                (isnothing(stations) || (length(fields) >= 2 && fields[2] in stations))
        end,
        readdir(output_directory; join=true),
    )
    # SPECFEM changes the channel prefix when seismotype/output conventions
    # change (for example HXZ -> CXZ). Old files are not always removed by a
    # new run. Keep only the channel family written most recently.
    families = Dict{String,Vector{String}}()
    for path in paths
        fields = split(basename(path), '.')
        channel = length(fields) >= 4 ? fields[3] : basename(path)
        push!(get!(families, channel, String[]), path)
    end
    if length(families) > 1
        channels = collect(keys(families))
        newest_channel = argmax(
            channel -> maximum(mtime, families[channel]), channels,
        )
        paths = families[newest_channel]
    end
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
    raw_residual = b .- a
    bias = sum(raw_residual) / length(raw_residual)
    mse = sum(abs2, raw_residual) / length(raw_residual)
    error_variance = sum(abs2, raw_residual .- bias) / length(raw_residual)
    reference_mean = sum(a) / length(a)
    reference_variance = sum(abs2, a .- reference_mean) / length(a)
    relative_rmse = reference_variance == 0 ? NaN : sqrt(mse / reference_variance)
    a .-= sum(a) / length(a)
    b .-= sum(b) / length(b)
    denom = sqrt(sum(abs2, a) * sum(abs2, b))
    correlation = denom == 0 ? NaN : sum(a .* b) / denom
    amplitude = sum(abs2, b) == 0 ? NaN : sum(a .* b) / sum(abs2, b)
    relative_error = sum(abs2, a) == 0 ? NaN :
        sqrt(sum(abs2, a .- amplitude .* b) / sum(abs2, a))
    (; correlation, relative_error, optimal_amplitude=amplitude,
       bias, mse, error_variance, reference_variance, relative_rmse, time)
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
