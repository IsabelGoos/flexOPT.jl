# earth2D_StagYY_neutrino_seismo.jl
#
# 2-D Earth disk built from a StagYY mantle-convection snapshot:
#   * spherically averaged (1-D) and heterogeneous (2-D) models,
#   * long-period P-SV waveforms (FD3: 3 points in space and time, traction-free
#     surface) for a vertical point force at 600 km depth, receivers every 5°,
#   * 3-flavour ν and ν̄ oscillations through the same two disks, folded with the
#     Daemonflux atmospheric flux and CC cross sections (perfect detector),
#   * one 2×3 summary figure (δT, δC, δwater, wavefield, δN_νe events, δu_r).
#
# Derived from notebooks/neutrinoPropagation.ipynb (StagYY reading, Z/A, nₑ,
# creationPaths, Pνν) and notebooks/HomogeneousElastic2DBenchmark_FreeSurface.ipynb
# (elasticWave2D FD3 operators, force source, receiver sampling).
#
# Run from the repository root:
#     julia --project=. --threads=auto scripts/earth2D_StagYY_neutrino_seismo.jl
# Optional environment overrides:
#     FLEXOPT_STAGYY_DIR, FLEXOPT_STAGYY_SNAPSHOT, FLEXOPT_EARTH2D_DX (m),
#     FLEXOPT_EARTH2D_DURATION (s), FLEXOPT_EARTH2D_PERIOD (s)
#
# NB on velocities: the StagYY directory has no seismic velocities (test_vp* is the
# velocity/pressure vector field, test_cs* is surface topography). δlnVs is
# therefore obtained from δT, δbasalt and δwater with simple, editable
# sensitivities (see `seismicScaling` below); δlnρ comes directly from test_rho*.

import Pkg
function find_flexopt_root(start=@__DIR__)
    candidates = haskey(ENV, "FLEXOPT_ROOT") ? [ENV["FLEXOPT_ROOT"]] : String[]
    directory = abspath(start)
    while true
        push!(candidates, directory)
        parent = dirname(directory)
        parent == directory && break
        directory = parent
    end
    for candidate in unique(abspath.(expanduser.(candidates)))
        isfile(joinpath(candidate, "src", "commonBatchs.jl")) && return candidate
    end
    error("Cannot locate flexOPT; set ENV[\"FLEXOPT_ROOT\"]")
end
flexopt_root = find_flexopt_root()
Pkg.activate(flexopt_root)

include(joinpath(flexopt_root, "src", "commonBatchs.jl"))
include(joinpath(flexopt_root, "src", "planet1D.jl"))
planet1D.configure_input!()
include(joinpath(flexopt_root, "src", "neutrinoOscillation.jl"))
include(joinpath(flexopt_root, "src", "GeoPoints.jl"))           # creationPaths
include(joinpath(flexopt_root, "src", "elasticWave2D.jl"))       # FD3 state + source
include(joinpath(flexopt_root, "src", "neutrinoFlux.jl"))
include(joinpath(flexopt_root, "src", "neutrinoCrossSections.jl"))
include(joinpath(flexopt_root, "src", "neutrinoDetector.jl"))
using .commonBatchs, .planet1D, .neutrinoOscillation, .GeoPoints, .elasticWave2D
using .neutrinoFlux, .Neutrino_Cross_Sections, .neutrinoDetector
using Interpolations, Statistics, LinearAlgebra, JLD2, CairoMakie, Printf
using Base.Threads
CairoMakie.activate!(type="png")
@show VERSION Threads.nthreads()

# ─────────────────────────────────────────────────────────────────────────────
# 0. User parameters
# ─────────────────────────────────────────────────────────────────────────────
stagyyDirectory = get(ENV, "FLEXOPT_STAGYY_DIR",
    "/Users/nobuaki/Documents/MantleConvectionTakashi/op_old_full_mars_2025")
stagyySnapshot = parse(Int, get(ENV, "FLEXOPT_STAGYY_SNAPSHOT", "199"))
stagyyFile(field) = joinpath(stagyyDirectory, @sprintf("test_%s%05d", field, stagyySnapshot))
stagyyFiles = (temperature=stagyyFile("t"), basalt=stagyyFile("bs"),
    water=stagyyFile("wtr"), density=stagyyFile("rho"))
all(isfile, stagyyFiles) || error("missing StagYY file(s): $(filter(!isfile, collect(stagyyFiles)))")
stagyyRotationDeg = 0.0        # rigid in-plane rotation of the convection pattern
seamBlendAngle = deg2rad(2.0)  # removes the StagYY θ=0 seam (as in the neutrino notebook)

dx = parse(Float64, get(ENV, "FLEXOPT_EARTH2D_DX", "10e3"))       # m, square grid
duration = parse(Float64, get(ENV, "FLEXOPT_EARTH2D_DURATION", "2400.0"))  # s
dominantPeriod = parse(Float64, get(ENV, "FLEXOPT_EARTH2D_PERIOD", "30.0")) # s
sourceDepth = 600e3                       # m
sourceForce = 1.0e15                      # N/m (line force, unit-thickness slice)
receiverSpacingDeg = 5.0
receiverDepth = 2dx                       # avoid the staircase surface row itself
outputSampling = 1.0                      # s between stored receiver samples
cfl = 0.38                                # same as ElasticThreePointConfig2D
surfaceDamping = 2e-3                     # per-step damping on the outermost material cells only

# LID values replace crust + ocean above this radius (long-period cheat: no fluid
# layer at the free surface, no thin low-velocity crust to resolve at 10 km).
mohoRadius = 6346.6e3

# δlnVs = aT(depth)·δT + aC·δbasalt + aW·δwater ; δlnVp = δlnVs / RSP
seismicScaling = (aT_top=-1.2e-4, aT_cmb=-0.5e-4,   # 1/K, linear in depth
    aC=0.015, aW=0.0, RSP=2.0, clip=0.10)

# Neutrinos: same binning and flux as neutrinoPropagation.ipynb, so the cached
# Daemonflux table in data/neutrinoFlux is reused (no Python call needed).
energies = logrange(1.0, 100.0, 200)      # GeV
cosθgrid = range(-1.0, 0.0, 400)
binning = (energies=energies, cosθgrid=cosθgrid)
Δbaseline = 50e3                          # m
detectorDepth = 2.5e3                     # m below the surface, above the source

outputDirectory = joinpath(flexopt_root, "data", "earth2D_stagyy",
    @sprintf("snap%05d_dx%dkm_T%ds", stagyySnapshot, round(Int, dx / 1e3), round(Int, dominantPeriod)))
mkpath(outputDirectory)

# ─────────────────────────────────────────────────────────────────────────────
# 1. Regular Cartesian disk centred on the planet (spherical Earth)
# ─────────────────────────────────────────────────────────────────────────────
planetRadius = planet1D.my1DDSMmodel.averagedPlanetRadiusInKilometer * 1e3
cmbRadius = planet1D.my1DDSMmodel.averagedPlanetCMBInKilometer * 1e3
icbRadius = 1_221_500.0
halfCells = ceil(Int, planetRadius / dx) + 3
xs = collect((-halfCells:halfCells) .* dx)
zs = copy(xs)
X = [x for x in xs, z in zs]
Z = [z for x in xs, z in zs]
radius = hypot.(X, Z)
material = radius .<= planetRadius
mantle = material .& (radius .>= cmbRadius)   # StagYY anomalies live only here
@show size(X) count(material) count(mantle)

# ─────────────────────────────────────────────────────────────────────────────
# 2. 1-D reference (planet1D, PREM by default) on the disk
# ─────────────────────────────────────────────────────────────────────────────
radialTableKm = collect(0.0:1.0:planetRadius / 1e3)
tableRadii, tableParams = planet1D.compute1DseismicParamtersFromPolynomialCoefficientsWithGivenRadiiArray(
    planet1D.my1DDSMmodel, radialTableKm, "below")
radialItp(values) = LinearInterpolation(tableRadii .* 1e3, values; extrapolation_bc=Flat())
ρItp, vpItp, vsItp = radialItp(tableParams.ρ), radialItp(tableParams.Vpv), radialItp(tableParams.Vsv)
lookupRadius = min.(radius, mohoRadius - 1.0)
ρ1D = ifelse.(material, ρItp.(lookupRadius), 0.0)    # g cm⁻³
vp1D = ifelse.(material, vpItp.(lookupRadius), 0.0)  # km s⁻¹
vs1D = ifelse.(material, vsItp.(lookupRadius), 0.0)  # km s⁻¹ (0 in the outer core)

# ─────────────────────────────────────────────────────────────────────────────
# 3. StagYY fields: periodic polar sampling, CMB→surface radial mapping
# ─────────────────────────────────────────────────────────────────────────────
rotationAngles = (θshift=deg2rad(stagyyRotationDeg), ϕshift=0.0)
sampleStagYY(file) = getCartesianField(file, X, Z;
    rotationAngles=rotationAngles, interpolation_method=:polar,
    target_cmb_radius=cmbRadius, target_surface_radius=planetRadius,
    clamp_to_surface=true, seam_blend_angle=seamBlendAngle)
TStag = sampleStagYY(stagyyFiles.temperature)
CStag = sampleStagYY(stagyyFiles.basalt)
WStag = sampleStagYY(stagyyFiles.water)
ρStag = sampleStagYY(stagyyFiles.density)
mantleOnly(field) = ifelse.(mantle, field, 0.0)
δT = mantleOnly(TStag.diffField)          # K, wrt the angular (spherical) mean
δC = mantleOnly(CStag.diffField)          # basalt fraction
δW = mantleOnly(WStag.diffField)          # water fraction
waterFraction = mantleOnly(clamp.(WStag.field, 0.0, 1.0))
δlnρ = mantleOnly(ρStag.diffField ./ max.(ρStag.avField, 1.0))

# ─────────────────────────────────────────────────────────────────────────────
# 4. Heterogeneous seismic model. Inner/outer core (ρ, Vp, Vs) stay at the 1-D
#    reference values: every perturbation below is multiplied by the mantle mask.
# ─────────────────────────────────────────────────────────────────────────────
depthFraction = clamp.((planetRadius .- radius) ./ (planetRadius - cmbRadius), 0.0, 1.0)
aT = seismicScaling.aT_top .+ (seismicScaling.aT_cmb - seismicScaling.aT_top) .* depthFraction
δlnVs = mantleOnly(clamp.(aT .* δT .+ seismicScaling.aC .* δC .+ seismicScaling.aW .* δW,
    -seismicScaling.clip, seismicScaling.clip))
δlnVp = δlnVs ./ seismicScaling.RSP
ρ2D = ρ1D .* (1 .+ δlnρ)
vp2D = vp1D .* (1 .+ δlnVp)
vs2D = vs1D .* (1 .+ δlnVs)
core = material .& .!mantle
@assert ρ2D[core] == ρ1D[core] && vp2D[core] == vp1D[core] && vs2D[core] == vs1D[core]
@show extrema(δT[mantle]) extrema(δC[mantle]) extrema(δW[mantle])
@show extrema(δlnVs[mantle]) extrema(δlnρ[mantle])

# ─────────────────────────────────────────────────────────────────────────────
# 5. FD3 (elasticWave2D operators) on a closed curved free surface.
#    Same stress / zero-stress-flux scheme as elasticWave2D.step_elastic_wave_2d!,
#    except that strain is taken one-sided next to a void node. The centred
#    version (fine on the flat benchmark) grows without bound on the circular
#    staircase within a few hundred seconds. A tiny per-step damping on the
#    outermost material cells (`surfaceDamping`) removes the remaining slow
#    growth at a few staircase corners; body waves are unchanged by it, the
#    late surface-wave coda is attenuated.
# ─────────────────────────────────────────────────────────────────────────────
@inline function dpx(u, m, i, k, h)
    p, q = m[i+1, k], m[i-1, k]
    p && q && return (u[i+1, k] - u[i-1, k]) / (2h)
    p && return (u[i+1, k] - u[i, k]) / h
    q && return (u[i, k] - u[i-1, k]) / h
    return zero(eltype(u))
end
@inline function dpz(u, m, i, k, h)
    p, q = m[i, k+1], m[i, k-1]
    p && q && return (u[i, k+1] - u[i, k-1]) / (2h)
    p && return (u[i, k+1] - u[i, k]) / h
    q && return (u[i, k] - u[i, k-1]) / h
    return zero(eltype(u))
end

function step_closed_surface!(state::ElasticWaveState2D{T}, edgeFactor) where T
    (; ux, uz, ux_previous, uz_previous, ux_next, uz_next,
       ρ, λ, μ, σxx, σzz, σxz, material) = state
    dx, dz = state.spacing
    nx, nz = size(ux)
    @threads for k in 2:nz-1
        for i in 2:nx-1
            if !material[i, k]
                σxx[i, k] = σzz[i, k] = σxz[i, k] = zero(T)
                continue
            end
            ux_x = dpx(ux, material, i, k, dx); uz_x = dpx(uz, material, i, k, dx)
            ux_z = dpz(ux, material, i, k, dz); uz_z = dpz(uz, material, i, k, dz)
            σxx[i, k] = (λ[i, k] + 2μ[i, k]) * ux_x + λ[i, k] * uz_z
            σzz[i, k] = λ[i, k] * ux_x + (λ[i, k] + 2μ[i, k]) * uz_z
            σxz[i, k] = μ[i, k] * (ux_z + uz_x)
        end
    end
    dt2 = state.dt^2
    @threads for k in 2:nz-1
        for i in 2:nx-1
            if !material[i, k]
                ux_next[i, k] = uz_next[i, k] = zero(T)
                continue
            end
            σxx_r = material[i+1, k] ? (σxx[i, k] + σxx[i+1, k]) / 2 : zero(T)
            σxx_l = material[i-1, k] ? (σxx[i, k] + σxx[i-1, k]) / 2 : zero(T)
            σxz_t = material[i, k+1] ? (σxz[i, k] + σxz[i, k+1]) / 2 : zero(T)
            σxz_b = material[i, k-1] ? (σxz[i, k] + σxz[i, k-1]) / 2 : zero(T)
            σxz_r = material[i+1, k] ? (σxz[i, k] + σxz[i+1, k]) / 2 : zero(T)
            σxz_l = material[i-1, k] ? (σxz[i, k] + σxz[i-1, k]) / 2 : zero(T)
            σzz_t = material[i, k+1] ? (σzz[i, k] + σzz[i, k+1]) / 2 : zero(T)
            σzz_b = material[i, k-1] ? (σzz[i, k] + σzz[i, k-1]) / 2 : zero(T)
            ax = ((σxx_r - σxx_l) / dx + (σxz_t - σxz_b) / dz) / ρ[i, k]
            az = ((σxz_r - σxz_l) / dx + (σzz_t - σzz_b) / dz) / ρ[i, k]
            f = edgeFactor[i, k]
            ux_next[i, k] = f * (2ux[i, k] - ux_previous[i, k] + dt2 * ax)
            uz_next[i, k] = f * (2uz[i, k] - uz_previous[i, k] + dt2 * az)
        end
    end
    state.ux_previous, state.ux, state.ux_next = state.ux, state.ux_next, state.ux_previous
    state.uz_previous, state.uz, state.uz_next = state.uz, state.uz_next, state.uz_previous
    state.step += 1
    state.time += state.dt
    state
end

# Source: vertical point force under the "north pole".
sourceIndex = CartesianIndex(argmin(abs.(xs)), argmin(abs.(zs .- (planetRadius - sourceDepth))))
@assert mantle[sourceIndex]
sourceFrequency = 1 / dominantPeriod
sourceDelay = 1.5 / sourceFrequency

# Receivers every 5° of epicentral distance, on the −x half (same side as the
# counter-clockwise neutrino rays of creationPaths).
epicentralDistances = collect(0.0:receiverSpacingDeg:180.0)
receiverIndices = map(epicentralDistances) do Δ
    r = planetRadius - receiverDepth
    while true
        idx = CartesianIndex(argmin(abs.(xs .+ r * sind(Δ))), argmin(abs.(zs .- r * cosd(Δ))))
        material[idx] && return idx
        r -= dx / 2
    end
end
receiverRadialUnit = [(xs[I[1]], zs[I[2]]) ./ hypot(xs[I[1]], zs[I[2]]) for I in receiverIndices]

# One common Δt for both models (identical time axes for the δ seismograms).
vpMaximum = 1e3 * max(maximum(vp1D[material]), maximum(vp2D[material]))
dtCommon = cfl / (vpMaximum * sqrt(2 / dx^2))

function run_fd3(ρ, vp, vs; label)
    model = (ρ=ρ, Vpv=vp, Vsv=vs)   # g cm⁻³, km s⁻¹ (prepare_elastic_wave_2d converts)
    bc = (cerjan=nothing,)   # no absorbing layer: the whole planet is modelled
    state = prepare_elastic_wave_2d(model, (dx, dx); material_mask=material,
        boundary_conditions=bc, dt=dtCommon)
    m = state.material
    edgeFactor = ones(Float32, size(m))
    for k in 2:size(m, 2)-1, i in 2:size(m, 1)-1
        m[i, k] && !(m[i+1, k] && m[i-1, k] && m[i, k+1] && m[i, k-1]) &&
            (edgeFactor[i, k] = 1 - Float32(surfaceDamping))
    end
    steps = ceil(Int, duration / state.dt)
    stride = max(1, round(Int, outputSampling / state.dt))
    nstore = steps ÷ stride + 1
    ux = zeros(nstore, length(receiverIndices)); uz = similar(ux)
    times = zeros(nstore)
    snapshotTimes = (300.0, 600.0, 900.0)
    snapshots = Dict{Float64,Matrix{Float32}}()
    stored = 1
    elapsed = @elapsed for step in 1:steps
        step_closed_surface!(state, edgeFactor)
        add_ricker_source!(state, sourceIndex; f0=sourceFrequency, t0=sourceDelay,
            amplitude=sourceForce, component=:z, source_kind=:force)
        if step % stride == 0 && stored < nstore
            stored += 1
            times[stored] = state.time
            for (j, I) in enumerate(receiverIndices)
                ux[stored, j] = state.ux[I]; uz[stored, j] = state.uz[I]
            end
        end
        for ts in snapshotTimes
            abs(state.time - ts) < state.dt / 2 && (snapshots[ts] = Float32.(hypot.(state.ux, state.uz)))
        end
        if step % 2000 == 0
            amplitude = maximum(abs, state.uz)
            isfinite(amplitude) || error("$label FD3 diverged at t=$(state.time) s")
            @info "$label" t = round(state.time) max_uz = amplitude
        end
    end
    ur = similar(ux)
    for j in eachindex(receiverIndices)
        ex, ez = receiverRadialUnit[j]
        ur[:, j] .= ux[:, j] .* ex .+ uz[:, j] .* ez
    end
    @info "$label done" dt = state.dt steps elapsed
    return (time=times[1:stored], ur=ur[1:stored, :], ux=ux[1:stored, :],
        uz=uz[1:stored, :], snapshots, dt=state.dt, elapsed)
end

seismo1D = run_fd3(ρ1D, vp1D, vs1D; label="1D")
seismo2D = run_fd3(ρ2D, vp2D, vs2D; label="2D")

# ─────────────────────────────────────────────────────────────────────────────
# 6. Neutrinos: Z/A with mantle water, nₑ, radial (1-D) average, P(νe→νe)
# ─────────────────────────────────────────────────────────────────────────────
zOverALayers = makeZOverALayers(icb_radius=icbRadius, cmb_radius=cmbRadius,
    surface_radius=planetRadius + dx, inner_core=0.466, outer_core=0.466, mantle=0.496)
zOverA = layeredZOverA(radius; layers=zOverALayers, water_fraction=waterFraction,
    water_z_over_a=5 / 9, material_mask=material, outside=0.0).mixed
nₑ2D = electronDensity(ρ2D, zOverA)
nₑStats = radialAverageAnomaly(nₑ2D, radius; bin_width=dx, mask=material, outside=0.0)
nₑ1D = nₑStats.radial_mean

detectorLocal = SVector(0.0, planetRadius - detectorDepth)
osc = set_oscillation_parameters()
function oscillation_probabilities(nₑ)
    sampling = creationPaths(nₑ, xs, zs, detectorLocal, binning;
        Δbaseline=Δbaseline, material_mask=material, modifiedLongitude2Disk=0.0,
        center=(0.0, 0.0), boundary_search_step=dx / 2, boundary_tolerance=10.0)
    conversion = pathsFromElectronDensity(sampling; reference_z_over_a=0.5)
    ν = Pνν(osc, energies, conversion.paths;
        zoa=conversion.reference_z_over_a, roundU_and_H=false)
    antiν = Pνν(osc, energies, conversion.paths;
        zoa=conversion.reference_z_over_a, anti=true, roundU_and_H=false)
    return ν, antiν, sampling                      # (energy, cosθ, initial, final)
end
P1D, antiP1D, sampling1D = oscillation_probabilities(nₑ1D)
P2D, antiP2D, _ = oscillation_probabilities(nₑ2D)
Pee1D = Array(parent(P1D))[:, :, 1, 1]
Pee2D = Array(parent(P2D))[:, :, 1, 1]

# Perfect detector (20 Mt × 1 yr) with the atmospheric flux: νe + ν̄e CC events,
# including the νμ → νe appearance channel, as in plot_interacting_events.
set_default_neutrino_flux!(:Daemonflux; flux_mode=:total, return_uncertainties=false)
neutrinoFluxTable = produce_or_load_neutrino_flux(binning)
neutrinoCrossSectionTable = evaluate_neutrino_cross_sections(neutrinoFluxTable.energies)
perfectDetector = PerfectDetector(mass_kg=20.0e9, livetime_s=365.25 * 24.0 * 3600.0)
events1D = interacting_events(neutrinoFluxTable, neutrinoCrossSectionTable, P1D, antiP1D;
    detector=perfectDetector)
events2D = interacting_events(neutrinoFluxTable, neutrinoCrossSectionTable, P2D, antiP2D;
    detector=perfectDetector)
Nνe1D = events1D.combined.νe
Nνe2D = events2D.combined.νe
δNνe = Nνe2D .- Nνe1D
@show sum(Nνe1D) sum(δNνe) extrema(δNνe)

# ─────────────────────────────────────────────────────────────────────────────
# 7. Save products
# ─────────────────────────────────────────────────────────────────────────────
jldsave(joinpath(outputDirectory, "earth2D_products.jld2");
    xs, zs, material, mantle, ρ1D, vp1D, vs1D, ρ2D, vp2D, vs2D,
    δT, δC, δW, δlnVs, δlnρ, nₑ1D, nₑ2D,
    epicentralDistances, sourceIndex, receiverIndices,
    time=seismo1D.time, ur1D=seismo1D.ur, ur2D=seismo2D.ur,
    energies=collect(energies), cosθgrid=collect(cosθgrid), Pee1D, Pee2D, Nνe1D, Nνe2D,
    seismicScaling, stagyyFiles, stagyyRotationDeg, dx, dominantPeriod, sourceDepth)

# ─────────────────────────────────────────────────────────────────────────────
# 8. ANR summary figure
# ─────────────────────────────────────────────────────────────────────────────
km = xs ./ 1e3
masked(field) = ifelse.(mantle, field, NaN)
symmetricLimit(field; q=0.99) = quantile(abs.(filter(isfinite, vec(field))), q)
rayAngles = (-1.0, -0.8, -0.5, -0.2)   # a few neutrino trajectories for the map

function disk_panel!(position, field, title, label, colormap; scale=1.0)
    layout = GridLayout(position)
    ax = Axis(layout[1, 1]; aspect=DataAspect(), title=title,
        xticksvisible=false, yticksvisible=false,
        xticklabelsvisible=false, yticklabelsvisible=false,
        leftspinevisible=false, rightspinevisible=false,
        topspinevisible=false, bottomspinevisible=false)
    lim = symmetricLimit(masked(field)) * scale
    hm = heatmap!(ax, km, km, masked(field) .* scale; colormap=colormap,
        colorrange=(-lim, lim), nan_color=:transparent)
    θ = range(0, 2π, 361)
    for r in (planetRadius, cmbRadius, icbRadius)
        lines!(ax, r / 1e3 .* cos.(θ), r / 1e3 .* sin.(θ); color=:gray35, linewidth=0.6)
    end
    Colorbar(layout[1, 2], hm; label=label, height=Relative(0.7))
    return ax
end

figure = Figure(size=(1500, 980), fontsize=13)
axT = disk_panel!(figure[1, 1], δT, "δT (wrt spherical mean)", "K", Reverse(:RdBu))
axC = disk_panel!(figure[1, 2], δC, "δ composition (basalt fraction)", "", :PuOr)
axW = disk_panel!(figure[1, 3], δW, "δ water", "wt %", :BrBG; scale=100)
axV = disk_panel!(figure[2, 1], δlnVs, "δlnVs, source ★, receivers ▼, ν rays", "%", :RdBu; scale=100)
scatter!(axV, [xs[sourceIndex[1]] / 1e3], [zs[sourceIndex[2]] / 1e3];
    marker=:star5, markersize=18, color=:gold, strokecolor=:black, strokewidth=1)
scatter!(axV, [xs[I[1]] / 1e3 for I in receiverIndices], [zs[I[2]] / 1e3 for I in receiverIndices];
    marker=:dtriangle, markersize=7, color=:black)
for cθ in rayAngles
    profile = sampling1D.profiles[argmin(abs.(sampling1D.cosθgrid .- cθ))]
    lines!(axV, [profile.source[1], profile.detector[1]] ./ 1e3,
        [profile.source[2], profile.detector[2]] ./ 1e3; color=(:magenta, 0.8), linewidth=1.2)
end

layoutP = GridLayout(figure[2, 2])
axP = Axis(layoutP[1, 1]; xscale=log10, xlabel="E (GeV)", ylabel="cos θ (zenith)",
    title="δN(νe+ν̄e CC) = N₂D − N₁D, perfect detector 20 Mt·yr")
limP = symmetricLimit(δNνe; q=0.995)
hmP = heatmap!(axP, collect(energies), collect(cosθgrid), δNνe;
    colormap=:balance, colorrange=(-limP, limP))
Colorbar(layoutP[1, 2], hmP; label="events / bin", height=Relative(0.7))

axS = Axis(figure[2, 3]; xlabel="time (s)", ylabel="epicentral distance (°)",
    title="δ seismogram: u_r 2-D (grey), δu_r = 2-D − 1-D (red)")
spacing = 0.48 * receiverSpacingDeg
bodyWaveWindow = 1500.0
n = min(size(seismo1D.ur, 1), size(seismo2D.ur, 1))
for j in eachindex(epicentralDistances)
    reference = maximum(abs, seismo2D.ur[seismo2D.time .< bodyWaveWindow, j])
    reference > 0 || continue
    lines!(axS, seismo2D.time[1:n], epicentralDistances[j] .+
        spacing .* clamp.(seismo2D.ur[1:n, j] ./ reference, -1.6, 1.6); color=:gray72, linewidth=0.8)
    lines!(axS, seismo2D.time[1:n], epicentralDistances[j] .+
        spacing .* clamp.((seismo2D.ur[1:n, j] .- seismo1D.ur[1:n, j]) ./ reference, -1.6, 1.6);
        color=:firebrick, linewidth=0.7)
end
xlims!(axS, 0, min(duration, bodyWaveWindow)); ylims!(axS, -4, 184)

Label(figure[0, :], @sprintf("StagYY snapshot %05d → 2-D Earth: geodynamics · neutrinos · seismic waves (T = %d s, source %d km)",
    stagyySnapshot, round(Int, dominantPeriod), round(Int, sourceDepth / 1e3)); fontsize=17)
save(joinpath(outputDirectory, "earth2D_panel.png"), figure; px_per_unit=2)
save(joinpath(outputDirectory, "earth2D_panel.pdf"), figure)
@info "figure written" outputDirectory
figure
