module seismicStations

using CSV
using Downloads
using Makie
using StaticArrays
using ..GeoPoints

export SeismicStation, StationBounds
export fetch_nied_stations, fetch_fdsn_stations
export station_local_coordinates, stations_in_box, plot_stations!

"""
    SeismicStation

Provider-neutral seismic station/channel metadata. Elevations and depths are in
metres. `ground_elevation_m` is the ground or station elevation,
`sensor_depth_m` is positive downward from that level, and
`sensor_elevation_m` is the sensor elevation above mean sea level.
"""
Base.@kwdef struct SeismicStation
    provider::Symbol
    network::String
    station::String
    location::String = ""
    channel::String = ""
    name::String = ""
    latitude::Float64
    longitude::Float64
    ground_elevation_m::Float64
    sensor_depth_m::Float64 = 0.0
    sensor_elevation_m::Float64 =
        ground_elevation_m - sensor_depth_m
    active::Bool = true
end

"""
    StationBounds(minlatitude, maxlatitude, minlongitude, maxlongitude)

Geographic rectangle in decimal degrees.
"""
struct StationBounds
    minlatitude::Float64
    maxlatitude::Float64
    minlongitude::Float64
    maxlongitude::Float64
    function StationBounds(minlatitude, maxlatitude, minlongitude, maxlongitude)
        minlatitude <= maxlatitude ||
            throw(ArgumentError("latitude bounds must be increasing"))
        minlongitude <= maxlongitude ||
            throw(ArgumentError("longitude bounds must be increasing"))
        new(minlatitude, maxlatitude, minlongitude, maxlongitude)
    end
end

function StationBounds(points::AbstractArray{<:GeoPoint})
    latitudes = (point.lat for point in points)
    longitudes = (point.lon for point in points)
    minlatitude, maxlatitude = extrema(latitudes)
    minlongitude, maxlongitude = extrema(longitudes)
    return StationBounds(
        minlatitude,
        maxlatitude,
        minlongitude,
        maxlongitude,
    )
end

_inside(station::SeismicStation, bounds::StationBounds) =
    bounds.minlatitude <= station.latitude <= bounds.maxlatitude &&
    bounds.minlongitude <= station.longitude <= bounds.maxlongitude

function _float(row, name)
    value = row[Symbol(name)]
    value === missing &&
        throw(ArgumentError("missing $name in station metadata"))
    return value isa Real ? Float64(value) : parse(Float64, strip(string(value), ['\'', '"']))
end

function _string(row, name; default="")
    value = row[Symbol(name)]
    value === missing && return default
    return strip(string(value), ['\'', '"', ' '])
end

"""
    fetch_nied_stations(; bounds=nothing, network=:HiNet, active_only=true,
        url=NIED_STATION_CSV_URL)

Download the current official NIED Hi-net/F-net station CSV and return normalized
station metadata. This metadata endpoint is public; downloading Hi-net waveform
data is a separate authenticated operation.
"""
const NIED_STATION_CSV_URL =
    "https://www.hinet.bosai.go.jp/st_info/detail/dlDialogue.php?f=CSV&LANG=en"

function fetch_nied_stations(;
    bounds::Union{Nothing,StationBounds}=nothing,
    network::Symbol=:HiNet,
    active_only::Bool=true,
    url::AbstractString=NIED_STATION_CSV_URL,
)
    network in (:HiNet, :FNet, :all) ||
        throw(ArgumentError("network must be :HiNet, :FNet, or :all"))
    buffer = IOBuffer()
    Downloads.request(url; output=buffer)
    seekstart(buffer)

    stations = SeismicStation[]
    for row in CSV.File(buffer; normalizenames=false, silencewarnings=true)
        station_code = _string(row, "station_cd")
        is_hinet = endswith(station_code, "H")
        is_fnet = endswith(station_code, "F")
        network === :HiNet && !is_hinet && continue
        network === :FNet && !is_fnet && continue
        active = _string(row, "repeal_station(1 = repeal)") != "1"
        active_only && !active && continue

        required = (
            Symbol("latitude"),
            Symbol("longitude"),
            Symbol("height(m)"),
            Symbol("sensor_height(m)"),
            Symbol("hole_depth(m)"),
        )
        any(name -> row[name] === missing, required) && continue

        ground = _float(row, "height(m)")
        sensor_elevation = _float(row, "sensor_height(m)")
        depth = _float(row, "hole_depth(m)")
        station = SeismicStation(
            provider=:NIED,
            network=is_hinet ? "Hi-net" : is_fnet ? "F-net" : _string(row, "network_id"),
            station=station_code,
            name=_string(row, "station_name"),
            latitude=_float(row, "latitude"),
            longitude=_float(row, "longitude"),
            ground_elevation_m=ground,
            sensor_depth_m=depth,
            sensor_elevation_m=sensor_elevation,
            active=active,
        )
        isnothing(bounds) || _inside(station, bounds) || continue
        push!(stations, station)
    end
    return stations
end

function _fdsn_url(
    service::AbstractString,
    bounds::StationBounds;
    network::AbstractString="*",
    station::AbstractString="*",
    location::AbstractString="*",
    channel::AbstractString="*Z",
    starttime=nothing,
    endtime=nothing,
)
    parameters = [
        "format=text",
        "level=channel",
        "nodata=404",
        "network=$(network)",
        "station=$(station)",
        "location=$(location)",
        "channel=$(channel)",
        "minlatitude=$(bounds.minlatitude)",
        "maxlatitude=$(bounds.maxlatitude)",
        "minlongitude=$(bounds.minlongitude)",
        "maxlongitude=$(bounds.maxlongitude)",
    ]
    isnothing(starttime) || push!(parameters, "starttime=$(starttime)")
    isnothing(endtime) || push!(parameters, "endtime=$(endtime)")
    return rstrip(service, '/') * "/query?" * join(parameters, '&')
end

"""
    fetch_fdsn_stations(bounds; service=EARTHSCOPE_STATION_SERVICE, ...)

Query an FDSN station service at channel level so that sensor depth is retained.
The default is the current EarthScope (formerly IRIS) endpoint. One entry is
returned per channel epoch.
"""
const EARTHSCOPE_STATION_SERVICE =
    "https://service.earthscope.org/fdsnws/station/1"

function fetch_fdsn_stations(
    bounds::StationBounds;
    service::AbstractString=EARTHSCOPE_STATION_SERVICE,
    network::AbstractString="*",
    station::AbstractString="*",
    location::AbstractString="*",
    channel::AbstractString="*Z",
    starttime=nothing,
    endtime=nothing,
)
    url = _fdsn_url(
        service,
        bounds;
        network,
        station,
        location,
        channel,
        starttime,
        endtime,
    )
    buffer = IOBuffer()
    Downloads.request(url; output=buffer)
    seekstart(buffer)

    stations = SeismicStation[]
    for line in eachline(buffer)
        startswith(line, '#') && continue
        isempty(strip(line)) && continue
        fields = split(line, '|'; keepempty=true)
        length(fields) >= 17 ||
            throw(ArgumentError("unexpected FDSN channel row: $line"))
        ground = parse(Float64, fields[7])
        depth = parse(Float64, fields[8])
        push!(
            stations,
            SeismicStation(
                provider=:FDSN,
                network=fields[1],
                station=fields[2],
                location=fields[3],
                channel=fields[4],
                latitude=parse(Float64, fields[5]),
                longitude=parse(Float64, fields[6]),
                ground_elevation_m=ground,
                sensor_depth_m=depth,
                sensor_elevation_m=ground - depth,
                name=fields[11],
                active=isempty(fields[17]),
            ),
        )
    end
    return stations
end

function _local_point(station::SeismicStation, box; sensor::Bool)
    altitude = sensor ?
               station.sensor_elevation_m :
               station.ground_elevation_m
    point = GeoPoint(station.latitude, station.longitude; alt=altitude)
    return p_ECEF_to_local(
        point.ecef,
        box.pOriginECEF,
        box.rotationMatrix,
    )
end

"""
    station_local_coordinates(stations, box; position=:ground)

Convert station metadata to the Cartesian frame of a `constructLocalBox`
result. Coordinates are returned in metres.
"""
function station_local_coordinates(
    stations::AbstractVector{SeismicStation},
    box;
    position::Symbol=:ground,
)
    position in (:ground, :sensor) ||
        throw(ArgumentError("position must be :ground or :sensor"))
    sensor = position === :sensor
    points = [_local_point(station, box; sensor) for station in stations]
    return (
        x=[point[1] for point in points],
        y=[point[2] for point in points],
        z=[point[3] for point in points],
    )
end

function stations_in_box(
    stations::AbstractVector{SeismicStation},
    box,
)
    xlimits = extrema(point.xyz[1] for point in box.allGridsInCartesian)
    ylimits = extrema(point.xyz[2] for point in box.allGridsInCartesian)
    ground = station_local_coordinates(stations, box; position=:ground)
    indices = [
        index for index in eachindex(stations) if
        xlimits[1] <= ground.x[index] <= xlimits[2] &&
        ylimits[1] <= ground.y[index] <= ylimits[2]
    ]
    return stations[indices]
end

"""
    plot_stations!(axis, stations, box; units=:km, show_boreholes=true,
        labels=true)

Overlay station ground locations and, optionally, buried sensors on a Makie
`Axis3`. Surface triangles and sensor circles share the local box coordinates.
"""
function plot_stations!(
    axis,
    stations::AbstractVector{SeismicStation},
    box;
    units::Symbol=:km,
    show_boreholes::Bool=true,
    labels::Bool=true,
    ground_color=:red,
    sensor_color=:white,
)
    units in (:m, :km) || throw(ArgumentError("units must be :m or :km"))
    scale = units === :km ? 1e-3 : 1.0
    selected = stations_in_box(stations, box)
    ground = station_local_coordinates(selected, box; position=:ground)
    sensor = station_local_coordinates(selected, box; position=:sensor)

    ground_points = [
        Point3f(
            ground.x[index] * scale,
            ground.y[index] * scale,
            ground.z[index] * scale,
        ) for index in eachindex(selected)
    ]
    sensor_points = [
        Point3f(
            sensor.x[index] * scale,
            sensor.y[index] * scale,
            sensor.z[index] * scale,
        ) for index in eachindex(selected)
    ]

    ground_plot = scatter!(
        axis,
        ground_points;
        marker=:utriangle,
        color=ground_color,
        markersize=18,
    )
    sensor_plot = nothing
    if show_boreholes
        for index in eachindex(selected)
            lines!(
                axis,
                [ground_points[index], sensor_points[index]];
                color=(:black, 0.7),
                linewidth=2,
            )
        end
        sensor_plot = scatter!(
            axis,
            sensor_points;
            color=sensor_color,
            strokecolor=:black,
            strokewidth=1,
            markersize=10,
        )
    end
    if labels && !isempty(selected)
        text!(
            axis,
            [station.station for station in selected];
            position=ground_points,
            space=:data,
            fontsize=12,
            offset=(6, 6),
        )
    end
    return (; stations=selected, ground, sensor, ground_plot, sensor_plot)
end

end
