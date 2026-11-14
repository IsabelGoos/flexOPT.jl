# dev/greenland_nov2026/stages/config.jl
#
# Every physical / numerical choice of the Greenland project, in one place.
# `greenlandCase` (quick | standard | full) only changes the grid spacing and
# the Cerjan layers; everything else is shared.

# ---------------------------------------------------------------- stations
# Toyokuni et al. (2018, PEPI 277, 81-98), Fig. 1 and Table 1: the 16 GLISN
# stations (lat, lon in degrees, elevation in m, sensor depth in m).
GLISN_STATIONS = [
    (code="DY2G",  lat=66.480, lon=-46.309, elevation=2132.0, depth=0.0,  onIce=true),
    (code="ICESG", lat=69.092, lon=-39.647, elevation=2932.0, depth=0.0,  onIce=true),
    (code="NEEM",  lat=77.445, lon=-51.074, elevation=2513.0, depth=0.0,  onIce=true),
    (code="SUMG",  lat=72.574, lon=-38.462, elevation=3200.0, depth=2.0,  onIce=true),
    (code="ANGG",  lat=65.616, lon=-37.637, elevation=9.0,    depth=0.0,  onIce=false),
    (code="BORG",  lat=64.747, lon=-21.327, elevation=110.0,  depth=95.0, onIce=false),
    (code="DAG",   lat=76.771, lon=-18.655, elevation=23.0,   depth=0.0,  onIce=false),
    (code="DBG",   lat=74.307, lon=-20.219, elevation=1.0,    depth=0.0,  onIce=false),
    (code="IVI",   lat=61.206, lon=-48.171, elevation=15.0,   depth=0.0,  onIce=false),
    (code="KBS",   lat=78.915, lon=11.939,  elevation=87.0,   depth=3.0,  onIce=false),
    (code="NOR",   lat=81.605, lon=-16.661, elevation=36.0,   depth=0.0,  onIce=false),
    (code="NRS",   lat=61.160, lon=-45.419, elevation=80.0,   depth=0.0,  onIce=false),
    (code="NUUK",  lat=64.184, lon=-51.668, elevation=110.0,  depth=0.0,  onIce=false),
    (code="SFJD",  lat=66.996, lon=-50.621, elevation=330.0,  depth=0.0,  onIce=false),
    (code="SOEG",  lat=68.204, lon=-31.376, elevation=179.0,  depth=0.0,  onIce=false),
    (code="TULEG", lat=76.537, lon=-68.824, elevation=38.0,   depth=0.0,  onIce=false),
]

# ---------------------------------------------------------------- source
# Ice slide on the coast near Nuuk, modelled as a single force (Kawakatsu 1989;
# Ekström & Stark 2013). The force on the Earth is -m a(t): first opposite to
# the slide direction, then along it (helpers: landslide_force_shape).
SOURCE = (
    lat=64.18, lon=-51.70,      # Nuuk coast; the force sits on the top material node of that column
    slideAzimuth=245.0,         # direction the mass moves (deg from N): seaward, out of the fjord
    slidePlunge=0.0,            # 0 = horizontal slide; > 0 adds a downward component
    peakForce=1.0e11,           # N, illustrative (large landslide); everything is linear in it
    onset=20.0,                 # s
    duration=100.0,             # s, acceleration + deceleration
)

# ---------------------------------------------------------------- materials
# Ice from Toyokuni et al. (2018), section 3.5: (0.920 g/cm³, 3.81, 1.94 km/s).
ICE   = Phase(0.920, 3.81, 1.94)
WATER = Phase(1.03, 1.45, 0.0)          # sea water (the paper's meltwater: 1.00, 1.44)
AIR   = Phase(0.001, 0.314, 0.0)        # getParamsAndTopo defaults; air rows are pinned anyway
AIR_THRESHOLD = 0.5                     # cell = material if its air fraction < 0.5
ICE_THRESHOLD = 1.0                     # m of surface - bed to call a column ice-covered
ETOPO_LOCAL_DIR = joinpath(flexopt_root, "dataInput", "greenland")   # manual fallback for the grids

# ---------------------------------------------------------------- box
# Local Cartesian box of GeoPoints.constructLocalBox. Centre and x-axis azimuth
# minimise the box that contains the 16 stations (x points ~south, y ~east):
# the half-width is set in step1_box.jl from the station positions
# (farthest station + MIN_STATION_EDGE_DISTANCE), at least BOX_HALF_WIDTH_MIN.
BOX_CENTRE = (lat=72.230, lon=-32.421)
BOX_AXIS_AZIMUTH = 170.0
BOX_HALF_WIDTH_MIN = 1200e3        # m, physical part (Cerjan layers are added outside)
DEPTH_BELOW_SURFACE = 400e3        # m of Earth kept below the surface at the box corners
MIN_STATION_EDGE_DISTANCE = 150e3  # m: stations this far from the absorbing layers

# ---------------------------------------------------------------- grid presets
GRID_PRESETS = Dict(
    #             dx (m)   Cerjan cells   shortest period reasonably resolved (OPT3, ~4 pts/λ_S in crust)
    "quick"    => (dx=50e3, cerjan=8,  minPeriod=60.0),
    "standard" => (dx=30e3, cerjan=12, minPeriod=40.0),
    "full"     => (dx=20e3, cerjan=16, minPeriod=25.0),
)
preset = GRID_PRESETS[greenlandCase]
dx = preset.dx                      # isotropic grid, Δx = Δy = Δz (recipe in scaled coordinates)
cerjanCells = preset.cerjan
cerjanDamping = 3.0 / cerjanCells^2 # edge factor exp(-3), the Kirishima value (0.0053 for 24 cells)
minResolvedPeriod = preset.minPeriod

# ---------------------------------------------------------------- time marching
duration = 1200.0          # s: KBS is ~2550 km from Nuuk (Rayleigh ~730 s) + source
courant = 0.25             # dt = courant * dx / (sqrt(3) * Vp_max)
optTimeMarchingScheme = :predictor_corrector
optPredictorMode = :lumped_future
optCorrectorIterations = 8   # 3 left accepted residuals up to 4e-2 (median 3e-4) on the quick grid
optCorrectorTolerance = 1e-6
optCorrectorRelaxation = 1.0
blowupLimit = 1e3          # m: stop if the displacement gets absurd
snapshotEvery = 5.0        # s, surface maps of the wavefield (movie frames; 25 s was too coarse)
traceSampling = 1.0        # s, output traces (the run itself keeps every dt)

# ---------------------------------------------------------------- OPT recipe
optSupplementaryOrder = 2  # OPT3 as in the Kirishima benchmarks; 0 = cheaper recipe
recipeOnGPU = false        # CPU Float64 recipes (index note: GPU only if 4-D recipes get too slow)
freeSurfaceClosure = :dietrich3d   # :dietrich3d (pinned air + mirror moves) | :none (air kept as material)

# ---------------------------------------------------------------- output band
filterBand = (20.0, 1500.0)   # s, the requested band; see README for what the box can carry

GREENLAND_CACHE_DIR = joinpath("dev", "greenland_nov2026")
println("Greenland case = ", greenlandCase, ": dx = ", dx / 1e3, " km, Cerjan ", cerjanCells,
    " cells, periods ≳ ", minResolvedPeriod, " s are resolved")
