const DEFAULT_NIED_VELOCITY_PAGE = Ref(
    "https://www.hinet.bosai.go.jp/topics/sokudo_kozo/alljpn_download.php?LANG=ja",
)

const DEFAULT_NIED_VELOCITY_SOURCE = Ref(
    "https://www.hinet.bosai.go.jp/topics/sokudo_kozo/dlDialogue.php?ID=ALJ",
)

const DEFAULT_NIED_VELOCITY_CACHE = Ref(
    joinpath(first(DEPOT_PATH), "flexOPT", "niedVelocityModels"),
)

"""
    set_default_nied_velocity_source!(source; cache_dir=nothing)

Change the default NIED velocity archive URL or local ZIP path. The source is
used only when `velocity_model=:NIED2023`. A per-call `nied_source` keyword can
override it without changing global state.
"""
function set_default_nied_velocity_source!(
    source::AbstractString;
    cache_dir::Union{Nothing,AbstractString}=nothing,
)
    isempty(strip(source)) &&
        throw(ArgumentError("NIED velocity source cannot be empty"))
    DEFAULT_NIED_VELOCITY_SOURCE[] = String(source)
    if !isnothing(cache_dir)
        DEFAULT_NIED_VELOCITY_CACHE[] = abspath(expanduser(cache_dir))
    end
    return (
        source=DEFAULT_NIED_VELOCITY_SOURCE[],
        cache_dir=DEFAULT_NIED_VELOCITY_CACHE[],
    )
end

struct NIEDVelocityModel
    longitude::Vector{Float64}
    latitude::Vector{Float64}
    depth_km::Vector{Float64}
    Vp::Array{Float32,3}
    Vs::Array{Float32,3}
    P_confidence::Array{Float32,3}
    S_confidence::Array{Float32,3}
    source::String
end

_is_url(source::AbstractString) =
    startswith(lowercase(source), "https://") ||
    startswith(lowercase(source), "http://")

function _validate_zip(path::AbstractString)
    magic = open(path, "r") do input
        read(input, 4)
    end
    length(magic) == 4 && magic[1:2] == UInt8[0x50, 0x4b] ||
        throw(ArgumentError(
            "NIED source is not a ZIP archive: $path. Use the archive " *
            "download endpoint (for example dlDialogue.php?ID=ALJ), not " *
            "the HTML information page.",
        ))
    return path
end

function _source_key(source::AbstractString)
    if _is_url(source)
        return bytes2hex(SHA.sha1(source))[1:12]
    end
    path = abspath(expanduser(source))
    isfile(path) || throw(ArgumentError("NIED archive does not exist: $path"))
    return bytes2hex(SHA.sha1(read(path)))[1:12]
end

nied_velocity_source_key(
    source::AbstractString=DEFAULT_NIED_VELOCITY_SOURCE[],
) = _source_key(source)

function _nied_archive_path(
    source::AbstractString,
    cache_dir::AbstractString;
    force_download::Bool=false,
)
    mkpath(cache_dir)
    if !_is_url(source)
        path = abspath(expanduser(source))
        isfile(path) || throw(ArgumentError("NIED archive does not exist: $path"))
        return _validate_zip(path)
    end

    archive = joinpath(cache_dir, "NIED_velocity_$(_source_key(source)).zip")
    if force_download || !isfile(archive)
        temporary = archive * ".download"
        last_error = nothing
        for attempt in 1:3
            try
                Downloads.download(source, temporary)
                filesize(temporary) > 0 ||
                    error("downloaded NIED archive is empty")
                _validate_zip(temporary)
                mv(temporary, archive; force=true)
                last_error = nothing
                break
            catch error
                last_error = error
                isfile(temporary) && rm(temporary; force=true)
                attempt < 3 && sleep(2.0^attempt)
            end
        end
        if !isnothing(last_error)
            isfile(temporary) && rm(temporary; force=true)
            throw(last_error)
        end
    end
    return _validate_zip(archive)
end

function _zip_entry(reader, suffix::AbstractString)
    index = findfirst(file -> endswith(file.name, suffix), reader.files)
    isnothing(index) &&
        throw(ArgumentError("NIED archive does not contain $suffix"))
    return reader.files[index]
end

function _read_zip_vector(reader, suffix)
    entry = _zip_entry(reader, suffix)
    return parse.(Float64, split(String(read(entry))))
end

function _extract_velocity_table(archive, cache_dir, source_key)
    table_path = joinpath(cache_dir, "VEL_$source_key.dat")
    isfile(table_path) && return table_path
    temporary = table_path * ".extracting"
    reader = ZipFile.Reader(archive)
    try
        entry = _zip_entry(reader, "/VEL.dat")
        open(temporary, "w") do output
            while !eof(entry)
                write(output, read(entry, 1 << 20))
            end
        end
        mv(temporary, table_path; force=true)
    finally
        close(reader)
        isfile(temporary) && rm(temporary; force=true)
    end
    return table_path
end

function _parse_nied_archive(
    archive::AbstractString,
    source::AbstractString,
    cache_dir::AbstractString,
)
    reader = ZipFile.Reader(archive)
    longitude = latitude = depth_km = nothing
    try
        longitude = _read_zip_vector(reader, "/LON.lst")
        latitude = _read_zip_vector(reader, "/LAT.lst")
        depth_km = _read_zip_vector(reader, "/DEP.lst")
    finally
        close(reader)
    end

    nlon, nlat, ndepth =
        length(longitude), length(latitude), length(depth_km)
    expected_rows = nlon * nlat * ndepth
    Vp = Array{Float32}(undef, nlon, nlat, ndepth)
    Vs = similar(Vp)
    P_confidence = similar(Vp)
    S_confidence = similar(Vp)

    source_key = _source_key(source)
    table_path = _extract_velocity_table(archive, cache_dir, source_key)
    row_count = 0
    for row in CSV.File(
        table_path;
        header=false,
        delim=' ',
        ignorerepeated=true,
        types=Float32,
        silencewarnings=true,
    )
        row_count += 1
        zero_based = row_count - 1
        j = mod(zero_based, nlat) + 1
        i = mod(div(zero_based, nlat), nlon) + 1
        k = div(zero_based, nlat * nlon) + 1
        k <= ndepth ||
            throw(ArgumentError("NIED VEL.dat has more than $expected_rows rows"))

        if row_count == 1 || row_count == expected_rows
            isapprox(Float64(row.Column1), longitude[i]; atol=5e-4) &&
            isapprox(Float64(row.Column2), latitude[j]; atol=5e-4) &&
            isapprox(Float64(row.Column3), depth_km[k]; atol=5e-3) ||
                throw(ArgumentError("NIED coordinate lists do not match VEL.dat"))
        end
        Vp[i, j, k] = row.Column4
        Vs[i, j, k] = row.Column5
        P_confidence[i, j, k] = row.Column6
        S_confidence[i, j, k] = row.Column7
    end
    row_count == expected_rows ||
        throw(ArgumentError(
            "NIED VEL.dat has $row_count rows; expected $expected_rows",
        ))
    rm(table_path; force=true)

    return NIEDVelocityModel(
        longitude,
        latitude,
        depth_km,
        Vp,
        Vs,
        P_confidence,
        S_confidence,
        String(source),
    )
end

"""
    load_nied_velocity_model(; source=DEFAULT_NIED_VELOCITY_SOURCE[],
        cache_dir=DEFAULT_NIED_VELOCITY_CACHE[], force_download=false,
        force_parse=false)

Download or open an ALJ-format NIED archive, then load a compact parsed JLD2
cache. URL and local ZIP sources are both supported.
"""
function load_nied_velocity_model(;
    source::AbstractString=DEFAULT_NIED_VELOCITY_SOURCE[],
    cache_dir::AbstractString=DEFAULT_NIED_VELOCITY_CACHE[],
    force_download::Bool=false,
    force_parse::Bool=false,
)
    cache_path = abspath(expanduser(cache_dir))
    mkpath(cache_path)
    source_key = _source_key(source)
    parsed_path = joinpath(cache_path, "NIED_velocity_$source_key.jld2")
    if isfile(parsed_path) && !force_parse && !force_download
        return JLD2.load(parsed_path, "model")
    end

    archive = _nied_archive_path(
        source,
        cache_path;
        force_download,
    )
    model = _parse_nied_archive(archive, source, cache_path)
    JLD2.jldsave(parsed_path; model)
    return model
end

function _bracketing_index(coordinates, value)
    coordinates[1] <= value <= coordinates[end] || return nothing
    index = searchsortedlast(coordinates, value)
    index == length(coordinates) && (index -= 1)
    return index
end

function _trilinear(array, i, j, k, tx, ty, tz)
    c00 = muladd(tx, array[i + 1, j, k] - array[i, j, k], array[i, j, k])
    c10 = muladd(
        tx,
        array[i + 1, j + 1, k] - array[i, j + 1, k],
        array[i, j + 1, k],
    )
    c01 = muladd(
        tx,
        array[i + 1, j, k + 1] - array[i, j, k + 1],
        array[i, j, k + 1],
    )
    c11 = muladd(
        tx,
        array[i + 1, j + 1, k + 1] - array[i, j + 1, k + 1],
        array[i, j + 1, k + 1],
    )
    c0 = muladd(ty, c10 - c00, c00)
    c1 = muladd(ty, c11 - c01, c01)
    return Float64(muladd(tz, c1 - c0, c0))
end

function _trilinear_slowness(velocity, i, j, k, tx, ty, tz)
    c000 = inv(Float64(velocity[i, j, k]))
    c100 = inv(Float64(velocity[i + 1, j, k]))
    c010 = inv(Float64(velocity[i, j + 1, k]))
    c110 = inv(Float64(velocity[i + 1, j + 1, k]))
    c001 = inv(Float64(velocity[i, j, k + 1]))
    c101 = inv(Float64(velocity[i + 1, j, k + 1]))
    c011 = inv(Float64(velocity[i, j + 1, k + 1]))
    c111 = inv(Float64(velocity[i + 1, j + 1, k + 1]))
    c00 = muladd(tx, c100 - c000, c000)
    c10 = muladd(tx, c110 - c010, c010)
    c01 = muladd(tx, c101 - c001, c001)
    c11 = muladd(tx, c111 - c011, c011)
    c0 = muladd(ty, c10 - c00, c00)
    c1 = muladd(ty, c11 - c01, c01)
    return muladd(tz, c1 - c0, c0)
end

function _sample_nied(model::NIEDVelocityModel, longitude, latitude, depth_km)
    i = _bracketing_index(model.longitude, longitude)
    j = _bracketing_index(model.latitude, latitude)
    k = _bracketing_index(model.depth_km, depth_km)
    any(isnothing, (i, j, k)) && return nothing

    tx = (longitude - model.longitude[i]) /
         (model.longitude[i + 1] - model.longitude[i])
    ty = (latitude - model.latitude[j]) /
         (model.latitude[j + 1] - model.latitude[j])
    tz = (depth_km - model.depth_km[k]) /
         (model.depth_km[k + 1] - model.depth_km[k])

    # Reproduce NIED's reference C/Fortran implementation: interpolate
    # slowness and invert, rather than interpolating velocity directly.
    inv_Vp = _trilinear_slowness(model.Vp, i, j, k, tx, ty, tz)
    inv_Vs = _trilinear_slowness(model.Vs, i, j, k, tx, ty, tz)
    return (
        Vp=inv(inv_Vp),
        Vs=inv(inv_Vs),
        P_confidence=_trilinear(
            model.P_confidence,
            i,
            j,
            k,
            tx,
            ty,
            tz,
        ),
        S_confidence=_trilinear(
            model.S_confidence,
            i,
            j,
            k,
            tx,
            ty,
            tz,
        ),
    )
end

function apply_nied_velocity_model!(
    seismic_model,
    all_grids_in_geopoints;
    source::AbstractString=DEFAULT_NIED_VELOCITY_SOURCE[],
    cache_dir::AbstractString=DEFAULT_NIED_VELOCITY_CACHE[],
    confidence_max::Union{Nothing,Real}=0.8,
    outside::Symbol=:planet1D,
    low_confidence::Symbol=:planet1D,
    force_download::Bool=false,
    force_parse::Bool=false,
)
    outside in (:planet1D, :error) ||
        throw(ArgumentError("outside must be :planet1D or :error"))
    low_confidence in (:planet1D, :accept, :error) ||
        throw(ArgumentError(
            "low_confidence must be :planet1D, :accept, or :error",
        ))
    isnothing(confidence_max) || confidence_max >= 0 ||
        throw(ArgumentError("confidence_max must be nonnegative or nothing"))

    model = load_nied_velocity_model(;
        source,
        cache_dir,
        force_download,
        force_parse,
    )
    P_confidence = fill(NaN, size(all_grids_in_geopoints))
    S_confidence = fill(NaN, size(all_grids_in_geopoints))
    nied_P_mask = falses(size(all_grids_in_geopoints))
    nied_S_mask = falses(size(all_grids_in_geopoints))

    Threads.@threads for index in eachindex(all_grids_in_geopoints)
        point = all_grids_in_geopoints[index]
        sample = _sample_nied(
            model,
            point.lon,
            point.lat,
            -point.alt * 1e-3,
        )
        if isnothing(sample)
            outside === :error &&
                error("point $(point.lon), $(point.lat), $(point.alt) is outside NIED")
            continue
        end
        P_confidence[index] = sample.P_confidence
        S_confidence[index] = sample.S_confidence
        P_reliable = isnothing(confidence_max) ||
                     sample.P_confidence <= confidence_max
        S_reliable = isnothing(confidence_max) ||
                     sample.S_confidence <= confidence_max
        if !(P_reliable && S_reliable) && low_confidence === :error
            error("NIED confidence exceeds threshold at index $index")
        end
        if low_confidence === :accept
            P_reliable = true
            S_reliable = true
        end
        if !P_reliable && !S_reliable
            continue
        end

        if P_reliable
            seismic_model.Vpv[index] = sample.Vp
            seismic_model.Vph[index] = sample.Vp
            nied_P_mask[index] = true
        end
        if S_reliable
            seismic_model.Vsv[index] = sample.Vs
            seismic_model.Vsh[index] = sample.Vs
            nied_S_mask[index] = true
        end
    end

    return (
        nied_P_confidence=P_confidence,
        nied_S_confidence=S_confidence,
        nied_P_mask=nied_P_mask,
        nied_S_mask=nied_S_mask,
        nied_mask=nied_P_mask .& nied_S_mask,
        nied_source=String(source),
    )
end
