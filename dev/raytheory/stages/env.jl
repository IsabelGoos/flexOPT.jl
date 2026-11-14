# dev/raytheory/stages/env.jl
# The ray-theory core (lib/) uses only the Julia standard library: src/ is not
# loaded at all. CairoMakie is loaded only by the plotting steps (plots.jl).
# lib/*.jl and checks.jl are re-read with --reload.
rtLoad = !@isdefined(RAYTHEORY_ENV) || DEV.reload
@show VERSION Threads.nthreads() Base.active_project()
using LinearAlgebra, Printf
RT_DIR = joinpath(flexopt_root, "dev", "raytheory")
rtLoad && include(joinpath(RT_DIR, "lib", "rayTheory.jl"))
rtLoad && include(joinpath(RT_DIR, "stages", "checks.jl"))

# 1-D models from dataInput/model (DSM .psv.model format), by case name
RT_MODELS = Dict(
    "PREM"   => joinpath("earth", "PREM.psv.model"),
    "ak135"  => joinpath("earth", "ak135.psv.model"),
    "iasp"   => joinpath("earth", "iasp.psv.model"),
    "DWThot" => joinpath("mars", "DWThot.psv.model"),
    "TAYAKr" => joinpath("mars", "TAYAKr.psv.model"),
    "vprem"  => joinpath("moon", "vprem.psv.model"),
)
rtModel(name; isotropic=:voigt) = readDSMModel(joinpath(flexopt_root, "dataInput", "model", RT_MODELS[name]); isotropic=isotropic, name=name)
RT_TAUP_CSV = joinpath(RT_DIR, "reference", "taup_premDSM.csv")
# Marmousi (step s4, and `check` when the file exists): the model the Marmousi
# notebooks build with src/imageReader; 802 × 676 samples over 9.2 km × 3.0 km,
# surface in the LAST column. JLD2 comes with devBootstrap.
RT_MARMOUSI_FILE = joinpath(flexopt_root, "notebooks", "tmp", "seismicModelMarmousi.jld2")
"Marmousi vp (km/s) as V[ix, iz], iz = 1 at the surface, and dx, dz (km)"
function loadMarmousiVp(file=RT_MARMOUSI_FILE; width=9.2, depth=3.0)
    m = load(file, "output")
    V = Matrix{Float64}(getfield(m, :Vpv))[:, end:-1:1]
    nx, nz = size(V)
    return V, width / (nx - 1), depth / (nz - 1)
end
RAYTHEORY_ENV = true
