# dev/kirishima_oct2026/codeVersions.jl
#
# One version string per CACHED stage. Each enters that stage's
# myProduceOrLoad parameters, so:
#   * changing a parameter (dx, duration, source, model, flags...) already
#     gives a new cache file automatically (new hash);
#   * changing the CODE of a stage (stages/*.jl, or src/ behind it) is not
#     seen by the hash -> bump the string here, and only that stage reruns.
#   * one-off recomputation without editing anything:  --redo=opt,fd,specfem
CODE_VERSION = Dict(
    "step2_opt"     => "v1-2026-10-21",  # recipes + operators + convFD3/OPT3/OPT3-noΓ/OPT4/OPT5 propagation
    "step2_fd"      => "v1-2026-10-21",  # explicit elasticWave2D FD3
    "step2_specfem" => "v1-2026-10-21",  # SPECFEM2D reference run + traces
)
# where myProduceOrLoad keeps the Step-2 results (under DrWatson's datadir(),
# i.e. <repo>/data/, which is git-ignored)
KIRISHIMA_CACHE_DIR = joinpath("dev", "kirishima_oct2026")
# Step-1/2 velocity models: lazyProduceOrLoad by name, in the folder the
# notebook already uses (notebooks/tmp), so the existing files are reused.
KIRISHIMA_MODEL_CACHE = joinpath(flexopt_root, "notebooks", "tmp")
