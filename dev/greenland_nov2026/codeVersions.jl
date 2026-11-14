# dev/greenland_nov2026/codeVersions.jl
# One version string per cached stage (enters its myProduceOrLoad key).
# Changing a parameter already gives a new cache file; changing the CODE of a
# stage (or helpers/greenlandIceWater.jl, or src/ behind it) is not seen ->
# bump the string here, or run once with --redo=<stage>.
CODE_VERSION = Dict(
    "etopo" => "v1-2026-11-13",   # ETOPO 2022 surface + bed download/reading
    "model" => "v2-2026-11-13",   # PREM by depth below bed, cell-averaged + ice/water/air mixing
    "opt"   => "v1-2026-11-13",   # OPT3 3-D recipe, operators, free surface, P-C marching
)
