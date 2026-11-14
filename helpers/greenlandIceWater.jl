# helpers/greenlandIceWater.jl
#
# Claude-authored helpers for dev/greenland_nov2026 (3-D elastic waves through
# Greenland, ice + ocean + PREM, single force near Nuuk). Nothing in src/ is
# changed; these functions only combine src/ results:
#
#   * ETOPO 2022 ice-surface and bedrock grids (NOAA NCEI) -> ice thickness and
#     water depth (GMT's @earth_relief has the ice SURFACE only, no bed);
#   * sub-cell mixing of air / ice / water / rock in every grid cell
#     (a 3 km ice sheet is far thinner than a 20-50 km cell);
#   * a 3-D free-surface closure on the prepared OPT linear system: air rows
#     pinned (u = 0) and the Dietrich mirror moves on the solid rows beside air,
#     the direct 3-D generalisation of src/numSolvers/linearSolution.jl's
#     `_pinned_void_rows_2d` + `_apply_dietrich_surface_2d` (prepareLinearSystem
#     only applies those in 2-D; in 3-D the mask is ignored);
#   * the landslide force history and local ZNE bases at stations.

module GreenlandIceWater

using Downloads, Interpolations, LinearAlgebra, SparseArrays, StaticArrays
using Printf
import GMT

export etopo_region, fetch_etopo, etopo_interpolators
export Phase, mix_cell, greenland_material, DepthProfile, depth_profile, average_rock, rock_top
export whole_material_mask, apply_free_surface_3d
export landslide_force_shape, zne_basis_local, local_direction
export estimate_operator_memory_gb

# ---------------------------------------------------------------------------
# ETOPO 2022 (60 arc-second), ice surface and bedrock, read on the user's machine
# ---------------------------------------------------------------------------

const ETOPO_BASE = "https://www.ngdc.noaa.gov/thredds"
etopo_path(kind) = "global/ETOPO2022/60s/60s_$(kind)_elev_netcdf/ETOPO_2022_v1_60s_N90W180_$(kind).nc"

"""
    etopo_region(points; margin=1.0)

[west, east, south, north] in degrees covering an array of `GeoPoint`s.
"""
function etopo_region(points; margin=1.0)
    lats = [p.lat for p in points]
    lons = [p.lon for p in points]
    return [max(-180.0, floor(minimum(lons) - margin)), min(180.0, ceil(maximum(lons) + margin)),
            max(-90.0, floor(minimum(lats) - margin)), min(90.0, ceil(maximum(lats) + margin))]
end

# GMTgrid -> (lon, lat, Z[lon, lat]) with ascending axes (same handling of the
# pixel-registration edge vectors as getParamsAndTopo in src/GeoPoints).
function _grid_to_arrays(G)
    z = Float64.(G.z)
    x = collect(Float64, G.x)
    y = collect(Float64, G.y)
    length(x) == size(z, 2) + 1 && (x = (x[1:end-1] .+ x[2:end]) ./ 2)
    length(y) == size(z, 1) + 1 && (y = (y[1:end-1] .+ y[2:end]) ./ 2)
    Z = permutedims(z)                       # (lon, lat)
    if !issorted(x)
        p = sortperm(x); x = x[p]; Z = Z[p, :]
    end
    if !issorted(y)
        p = sortperm(y); y = y[p]; Z = Z[:, p]
    end
    return x, y, Z
end

function _read_ncss(kind, region, stride, workdir)
    w, e, s, n = region
    url = "$(ETOPO_BASE)/ncss/grid/$(etopo_path(kind))?var=z&north=$(n)&west=$(w)&east=$(e)&south=$(s)" *
          "&horizStride=$(stride)&accept=netcdf3"
    file = joinpath(workdir, "etopo2022_60s_$(kind)_ncss.nc")
    println("ETOPO $(kind): NCSS ", url)
    Downloads.download(url, file)
    return _grid_to_arrays(GMT.gmtread(file; grid=true))..., url
end

# OPeNDAP ASCII fallback: no netCDF library needed. ETOPO 2022 60 s is
# cell-centred: lat_i = -90 + (i + 1/2)/60, lon_j = -180 + (j + 1/2)/60.
function _read_opendap_ascii(kind, region, stride, workdir)
    w, e, s, n = region
    i0 = clamp(floor(Int, (s + 90) * 60), 0, 10799)
    i1 = clamp(floor(Int, (n + 90) * 60), 0, 10799)
    j0 = clamp(floor(Int, (w + 180) * 60), 0, 21599)
    j1 = clamp(floor(Int, (e + 180) * 60), 0, 21599)
    url = "$(ETOPO_BASE)/dodsC/$(etopo_path(kind)).ascii?z[$(i0):$(stride):$(i1)][$(j0):$(stride):$(j1)]"
    file = joinpath(workdir, "etopo2022_60s_$(kind)_opendap.txt")
    println("ETOPO $(kind): OPeNDAP ", url)
    Downloads.download(url, file)
    rows = Vector{Vector{Float64}}()
    inData = false
    for line in eachline(file)
        if !inData
            # "z[nlat][nlon]" (plain array) or "z.z[nlat][nlon]" (grid)
            occursin(r"^z(\.z)?\[\d+\]\[\d+\]", strip(line)) && (inData = true)
            continue
        end
        m = match(r"^\[(\d+)\],(.*)$", strip(line))
        isnothing(m) && break
        push!(rows, parse.(Float64, split(m.captures[2], ',')))
    end
    isempty(rows) && error("no data rows found in $file")
    Zlatlon = reduce(vcat, permutedims.(rows))          # (nlat, nlon)
    lat = [-90 + (i + 0.5) / 60 for i in i0:stride:i1]
    lon = [-180 + (j + 0.5) / 60 for j in j0:stride:j1]
    size(Zlatlon) == (length(lat), length(lon)) ||
        error("OPeNDAP grid size $(size(Zlatlon)) differs from the request ($(length(lat)), $(length(lon)))")
    return lon, lat, permutedims(Zlatlon), url
end

function _read_local(kind, region, localDir)
    isdir(localDir) || error("no local directory $localDir")
    files = filter(f -> occursin(kind, lowercase(f)) &&
                        any(endswith(lowercase(f), ext) for ext in (".nc", ".grd", ".tif", ".tiff")),
                   readdir(localDir))
    isempty(files) && error("no *$(kind)*.nc/.grd/.tif in $localDir")
    file = joinpath(localDir, first(files))
    w, e, s, n = region
    println("ETOPO $(kind): local file ", file)
    G = GMT.grdcut(file; region=(w, e, s, n))
    return _grid_to_arrays(G)..., file
end

"""
    fetch_etopo(params) -> Dict

For `myProduceOrLoad`. `params` holds "region" ([w,e,s,n] in degrees),
"strideArcmin" (1 = native 60 s grid) and "localDir". Each of the two grids
(surface, bed) is read from the first source that works:

1. NOAA THREDDS NetCDF Subset Service (read with GMT);
2. NOAA THREDDS OPeNDAP, ASCII (no netCDF reader needed);
3. files you put yourself in `localDir`, any name containing "surface" / "bed"
   (.nc, .grd or GeoTIFF), e.g. downloaded from the NCEI ETOPO 2022 Grid Extract
   page (choose "Ice surface" and "Bedrock").

The bed grid is interpolated on the surface grid if their nodes differ.
"""
function fetch_etopo(params)
    region = Float64.(params["region"])
    stride = Int(params["strideArcmin"])
    localDir = String(params["localDir"])
    workdir = mktempdir()
    out = Dict{String,Any}()
    for kind in ("surface", "bed")
        result = nothing
        errors = String[]
        for reader in (() -> _read_ncss(kind, region, stride, workdir),
                       () -> _read_opendap_ascii(kind, region, stride, workdir),
                       () -> _read_local(kind, region, localDir))
            try
                result = reader()
                break
            catch err
                push!(errors, sprint(showerror, err)[1:min(end, 300)])
            end
        end
        isnothing(result) && error("""
            Could not read the ETOPO 2022 $(kind) grid. Tried:
              $(join(errors, "\n  "))
            Download it by hand from NCEI (ETOPO 2022 'Grid Extract', $(kind == "bed" ? "Bedrock" : "Ice surface"),
            region $(region)) and put the file in $(localDir) with '$(kind)' in its name.""")
        lon, lat, Z, source = result
        out[kind * "_lon"], out[kind * "_lat"], out[kind], out[kind * "_source"] = lon, lat, Z, source
    end
    if out["surface_lon"] != out["bed_lon"] || out["surface_lat"] != out["bed_lat"]
        itp = linear_interpolation((out["bed_lon"], out["bed_lat"]), out["bed"]; extrapolation_bc=Flat())
        out["bed"] = [itp(x, y) for x in out["surface_lon"], y in out["surface_lat"]]
    end
    out["lon"], out["lat"] = out["surface_lon"], out["surface_lat"]
    thickness = max.(out["surface"] .- out["bed"], 0.0)
    @printf("ETOPO 2022: %d x %d nodes, max ice thickness %.0f m, min bed %.0f m\n",
        size(out["surface"])..., maximum(thickness), minimum(out["bed"]))
    return out
end

"Bilinear (lon, lat) -> metres interpolators of the ETOPO surface and bed."
function etopo_interpolators(etopo)
    s = linear_interpolation((etopo["lon"], etopo["lat"]), etopo["surface"]; extrapolation_bc=Flat())
    b = linear_interpolation((etopo["lon"], etopo["lat"]), etopo["bed"]; extrapolation_bc=Flat())
    return s, b
end

# ---------------------------------------------------------------------------
# Sub-cell mixing of air / ice / water / rock
# ---------------------------------------------------------------------------

"Isotropic phase: density (g/cm³), Vp, Vs (km/s), the units of getParamsAndTopo."
struct Phase
    ρ::Float64
    vp::Float64
    vs::Float64
end

_overlap(a, b, zbot, ztop) = max(0.0, min(b, ztop) - max(a, zbot)) / (ztop - zbot)

"""
    DepthProfile / depth_profile(depth_m, ρ, vp, vs)

1-D rock profile (e.g. planet1D PREM) as a function of depth below the local
rock surface, stored as cumulative integrals of ρ, ρVp², ρVs² so that
`average_rock(profile, d1, d2)` returns the exact Voigt average over any
depth interval (discontinuities such as the Moho included).
"""
struct DepthProfile
    d::Vector{Float64}
    Cρ::Vector{Float64}
    CM::Vector{Float64}
    Cμ::Vector{Float64}
end
function depth_profile(depth, ρ, vp, vs)
    p = sortperm(depth)   # stable: duplicated depths (discontinuities) keep below/above order
    d = Float64.(depth[p]); r = Float64.(ρ[p]); M = r .* Float64.(vp[p]) .^ 2; μ = r .* Float64.(vs[p]) .^ 2
    cum(f) = vcat(0.0, cumsum([(d[i+1] - d[i]) * (f[i] + f[i+1]) / 2 for i in 1:length(d)-1]))
    return DepthProfile(d, cum(r), cum(M), cum(μ))
end
function _cumulative(prof::DepthProfile, C, x)
    x <= prof.d[1] && return C[1] + (x - prof.d[1]) * (C[2] - C[1]) / max(prof.d[2] - prof.d[1], eps())
    i = clamp(searchsortedlast(prof.d, x), 1, length(prof.d) - 1)
    while i < length(prof.d) - 1 && prof.d[i+1] == prof.d[i]; i += 1; end
    w = prof.d[i+1] - prof.d[i]
    return w == 0 ? C[i] : C[i] + (x - prof.d[i]) / w * (C[i+1] - C[i])
end
function average_rock(prof::DepthProfile, d1, d2)
    d1, d2 = max(min(d1, d2), 0.0), max(d1, d2, 0.0)
    if d2 - d1 < 1.0
        d2 = d1 + 1.0
    end
    ρ = (_cumulative(prof, prof.Cρ, d2) - _cumulative(prof, prof.Cρ, d1)) / (d2 - d1)
    M = (_cumulative(prof, prof.CM, d2) - _cumulative(prof, prof.CM, d1)) / (d2 - d1)
    μ = (_cumulative(prof, prof.Cμ, d2) - _cumulative(prof, prof.Cμ, d1)) / (d2 - d1)
    return Phase(ρ, sqrt(M / ρ), sqrt(max(μ, 0.0) / ρ))
end

"Top of the rock (m above sea level): bed under ice or ocean, surface on bare land."
function rock_top(surface, bed; iceThreshold=1.0)
    isIce = (surface - bed > iceThreshold) && surface > 0
    return (isIce || surface <= 0) ? bed : surface
end

"""
    mix_cell(zbot, ztop, surface, bed, rock, ice, water; iceThreshold=1.0)

Volume fractions of air, ice, water and rock in the vertical interval
[zbot, ztop] (m above sea level) of a column whose ETOPO ice/land surface and
bed are `surface`, `bed` (m), and the Voigt average of the non-air part
(ρ, M = ρVp², μ = ρVs² averaged by volume).

* ice where surface - bed > iceThreshold and surface > 0 (assumed grounded:
  ice from bed to surface, rock below bed);
* ocean where surface <= 0: water from bed to sea level, air above 0;
* bare land otherwise: rock up to the surface.
"""
function mix_cell(zbot, ztop, surface, bed, rock::Phase, ice::Phase, water::Phase;
                  iceThreshold=1.0)
    isIce = (surface - bed > iceThreshold) && surface > 0
    isOcean = !isIce && surface <= 0
    fAir = isOcean ? _overlap(0.0, Inf, zbot, ztop) : _overlap(surface, Inf, zbot, ztop)
    fIce = isIce ? _overlap(bed, surface, zbot, ztop) : 0.0
    fWater = isOcean ? _overlap(bed, 0.0, zbot, ztop) : 0.0
    rockTop = (isIce || isOcean) ? bed : surface
    fRock = _overlap(-Inf, rockTop, zbot, ztop)
    fSolidish = fIce + fWater + fRock
    if fSolidish <= 1e-12
        return (ρ=0.0, vp=0.0, vs=0.0, fAir=1.0, fIce=0.0, fWater=0.0, fRock=0.0)
    end
    ρ = (fIce * ice.ρ + fWater * water.ρ + fRock * rock.ρ) / fSolidish
    M = (fIce * ice.ρ * ice.vp^2 + fWater * water.ρ * water.vp^2 + fRock * rock.ρ * rock.vp^2) / fSolidish
    μ = (fIce * ice.ρ * ice.vs^2 + fWater * water.ρ * water.vs^2 + fRock * rock.ρ * rock.vs^2) / fSolidish
    return (ρ=ρ, vp=sqrt(M / ρ), vs=sqrt(μ / ρ), fAir=fAir, fIce=fIce, fWater=fWater, fRock=fRock)
end

"""
    greenland_material(rockProfile, points, dz, surf, bed; ice, water, air, airThreshold=0.5)

Cell-by-cell model on the GeoPoints box. The rock of a cell is the Voigt
average of `rockProfile` (planet1D as a function of depth BELOW THE LOCAL
ROCK SURFACE, i.e. the bed under ice and ocean) over the rock part of that
cell, so the crust follows the curved Earth surface and a crust thinner than
one cell is averaged rather than hit-or-missed. `points` are the box
GeoPoints, `dz` the vertical spacing (m), `surf`/`bed` the ETOPO
interpolators. A cell is material when its air fraction is < `airThreshold`;
the others get the `air` phase (pinned later by the free-surface closure).
Returns the fields (g/cm³, km/s), fractions and the material mask.
"""
function greenland_material(rockProfile::DepthProfile, points, dz, surf, bed;
                            ice::Phase, water::Phase, air::Phase,
                            airThreshold=0.5, iceThreshold=1.0)
    sz = size(points)
    ρ = zeros(sz); vp = zeros(sz); vs = zeros(sz)
    fAir = zeros(sz); fIce = zeros(sz); fWater = zeros(sz); fRock = zeros(sz)
    material = falses(sz)
    Threads.@threads for k in 1:sz[3]
        for j in 1:sz[2], i in 1:sz[1]
            p = points[i, j, k]
            s = surf(p.lon, p.lat); b = bed(p.lon, p.lat)
            zbot, ztop = p.alt - dz / 2, p.alt + dz / 2
            rt = rock_top(s, b; iceThreshold)
            r = average_rock(rockProfile, rt - min(ztop, rt), rt - zbot)
            c = mix_cell(zbot, ztop, s, b, r, ice, water; iceThreshold)
            fAir[i, j, k], fIce[i, j, k], fWater[i, j, k], fRock[i, j, k] = c.fAir, c.fIce, c.fWater, c.fRock
            if c.fAir < airThreshold
                material[i, j, k] = true
                ρ[i, j, k], vp[i, j, k], vs[i, j, k] = c.ρ, c.vp, c.vs
            else
                ρ[i, j, k], vp[i, j, k], vs[i, j, k] = air.ρ, air.vp, air.vs
            end
        end
    end
    # column maps (top layer lat/lon) for figures
    iceThickness = [max(surf(points[i, j, end].lon, points[i, j, end].lat) -
                        bed(points[i, j, end].lon, points[i, j, end].lat), 0.0) *
                    (surf(points[i, j, end].lon, points[i, j, end].lat) > 0)
                    for i in 1:sz[1], j in 1:sz[2]]
    waterDepth = [let s = surf(points[i, j, end].lon, points[i, j, end].lat)
                      s <= 0 ? -bed(points[i, j, end].lon, points[i, j, end].lat) : 0.0
                  end for i in 1:sz[1], j in 1:sz[2]]
    surfaceElevation = [surf(points[i, j, end].lon, points[i, j, end].lat) for i in 1:sz[1], j in 1:sz[2]]
    return (ρ=ρ, Vp=vp, Vs=vs, material=material, fAir=fAir, fIce=fIce, fWater=fWater, fRock=fRock,
            iceThickness=iceThickness, waterDepth=waterDepth, surfaceElevation=surfaceElevation)
end

# ---------------------------------------------------------------------------
# 3-D free surface on the prepared OPT system (pinned void + Dietrich mirror)
# ---------------------------------------------------------------------------

"Material mask on the whole (Cerjan-padded) grid, padding = nearest model point."
function whole_material_mask(prepared, materialMask::AbstractArray{Bool,3})
    conv = prepared.left_operator.geometry.conv
    shape = prepared.spaceShape
    whole = falses(shape)
    msz = size(materialMask)
    for wp in CartesianIndices(shape)
        raw = conv.whole2model(wp)
        mp = CartesianIndex(ntuple(d -> clamp(raw[d], 1, msz[d]), 3))
        whole[wp] = materialMask[mp]
    end
    return whole
end

@inline function _move!(M, row, oldc, newc, sign)
    v = M[row, oldc]
    iszero(v) && return
    M[row, oldc] = 0.0
    M[row, newc] += sign * v
    return
end

"""
    apply_free_surface_3d(prepared, materialMask; prepared_with_matrices)

Same closure as prepareLinearSystem applies in 2-D with
`free_surface_mode=:zero_traction_flux` and a material mask, generalised to 3-D:

* for every material point and every axis direction d (±x, ±y, ±z) whose
  neighbour is air and whose mirror point (point - d) is material, the
  coefficients multiplying the air neighbour are moved onto the mirror point,
  with sign +1 for the field normal to that face and -1 for the two
  tangential fields (exactly the sign rule of `_apply_dietrich_surface_2d`),
  in the future matrix A and in every known-time block of L;
* every air point gets the identity row (u = 0) and empty L/R rows
  (`_pinned_void_rows_2d`).

Row index = field + NField*(point-1); column = point + Npoints*(field-1)
(+ time-slot offsets in L), the layout of src/numSolvers/linearSolution.jl.
`prepared_with_matrices` is `flexOPT._prepared_with_matrices`.
"""
function apply_free_surface_3d(prepared, materialMask::AbstractArray{Bool,3}; prepared_with_matrices)
    NField = prepared.NField
    NField == 3 || throw(DimensionMismatch("3-D closure needs three displacement fields"))
    whole = whole_material_mask(prepared, materialMask)
    Np = length(whole)
    Np == prepared.NpointsSpace || throw(DimensionMismatch("mask and operator sizes differ"))
    LI = LinearIndices(whole)
    A = copy(prepared.A_unknown)
    L = copy(prepared.L_known)
    R = copy(prepared.R_force)
    nKnownTime = div(size(L, 2), Np * NField)
    directions = (CartesianIndex(-1, 0, 0), CartesianIndex(1, 0, 0),
                  CartesianIndex(0, -1, 0), CartesianIndex(0, 1, 0),
                  CartesianIndex(0, 0, -1), CartesianIndex(0, 0, 1))
    surfaceRows = Int[]
    nMoves = 0
    for point in CartesianIndices(whole)
        whole[point] || continue
        p = LI[point]
        touched = false
        for d in directions
            airPt = point + d
            checkbounds(Bool, whole, airPt) || continue
            whole[airPt] && continue
            mirror = point - d
            (checkbounds(Bool, whole, mirror) && whole[mirror]) || continue
            touched = true
            ap, mp = LI[airPt], LI[mirror]
            normalField = d[1] != 0 ? 1 : (d[2] != 0 ? 2 : 3)
            for eq in 1:NField, field in 1:NField
                row = eq + NField * (p - 1)
                sign = field == normalField ? 1.0 : -1.0
                oldc = ap + (field - 1) * Np
                newc = mp + (field - 1) * Np
                _move!(A, row, oldc, newc, sign)
                for slot in 1:nKnownTime
                    off = (slot - 1) * Np * NField
                    _move!(L, row, oldc + off, newc + off, sign)
                end
                nMoves += 1
            end
        end
        touched && append!(surfaceRows, (NField * (p - 1) + 1):(NField * p))
    end
    # pin the air (void) points
    keep = ones(size(A, 1))
    voidRows = Int[]; voidCols = Int[]
    for point in CartesianIndices(whole)
        whole[point] && continue
        p = LI[point]
        for field in 1:NField
            row = field + NField * (p - 1)
            keep[row] = 0.0
            push!(voidRows, row); push!(voidCols, p + (field - 1) * Np)
        end
    end
    D = Diagonal(keep)
    A = D * A + sparse(voidRows, voidCols, 1.0, size(A)...)
    L = D * L
    R = D * R
    dropzeros!(A); dropzeros!(L); dropzeros!(R)
    @printf("3-D free surface: %d material points touch air (%d coefficient groups moved), %d air points pinned\n",
        div(length(surfaceRows), NField), nMoves, div(length(voidRows), NField))
    return prepared_with_matrices(prepared, A, L, R;
        free_surface_rows=sort!(unique(surfaceRows)), void_rows=voidRows,
        whole_material=whole)
end

# ---------------------------------------------------------------------------
# Source and stations
# ---------------------------------------------------------------------------

"""
    landslide_force_shape(t; onset, duration)

Dimensionless force history (max |g| = 1) of a mass that accelerates then
stops over `duration` s, starting at `onset`:
  a(τ) ∝ sin(2πτ) sin²(πτ),  τ = (t - onset)/duration ∈ [0, 1].
Smooth at both ends, zero net impulse (the mass ends at rest) and positive
run-out. The force ON THE EARTH is -m a(t): opposite to the slide direction
while the mass accelerates, along it while it decelerates.
"""
function landslide_force_shape(t; onset, duration)
    τ = (t - onset) / duration
    (0 <= τ <= 1) || return 0.0
    return sin(2π * τ) * sin(π * τ)^2 / 0.6495190528383291   # max of sin(2πτ)sin²(πτ)
end

"""
    zne_basis_local(g, R) -> (z, n, e)

Unit vertical (ellipsoid normal), north and east vectors at GeoPoint `g`,
expressed in the local box frame of rotation matrix `R` (local = R' * ecef,
as GeoPoints.p_ECEF_to_local).
"""
function zne_basis_local(g, R)
    φ, λ = deg2rad(g.lat), deg2rad(g.lon)
    up = SVector(cos(φ) * cos(λ), cos(φ) * sin(λ), sin(φ))
    east = SVector(-sin(λ), cos(λ), 0.0)
    north = cross(up, east)
    return R' * up, R' * north, R' * east
end

"""
    local_direction(g, R; azimuth, plunge=0)

Unit vector, in the box frame, pointing to `azimuth` (degrees clockwise from
north) and `plunge` degrees below the horizontal at GeoPoint `g`.
"""
function local_direction(g, R; azimuth, plunge=0.0)
    z, n, e = zne_basis_local(g, R)
    a, p = deg2rad(azimuth), deg2rad(plunge)
    v = cos(p) * (cos(a) * n + sin(a) * e) - sin(p) * z
    return v / norm(v)
end

"""
    estimate_operator_memory_gb(nWholePoints; nField=3, stencil=27, timeLevels=3)

Rough size of A (future), L (2 past levels) and R (3 force levels) as sparse
Float64 matrices, plus a copy during the free-surface closure.
"""
function estimate_operator_memory_gb(nWholePoints; nField=3, stencil=27, timeLevels=3)
    rows = nWholePoints * nField
    nnzRow = stencil * nField
    blocks = 1 + (timeLevels - 1) + timeLevels
    return 2 * rows * nnzRow * blocks * 12 / 1e9
end

end # module
