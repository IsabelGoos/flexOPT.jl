# application/Poisson_1D/schemes.jl -- the numerical schemes this application offers.
#
# One NamedTuple per scheme, exactly the parameters makeOPTsemiSymbolic needs.
# `status` says how far it is benchmarked for THIS equation:
#   :validated     reproduces the analytical solution at the expected order
#   :experimental  under benchmark; may change, used only when asked by name
# `from` points to where the configuration was taken / benchmarked.
#
# When an operator is updated during a benchmark, edit it here and bump
# `version`: the version is part of the recipe-cache key (see app.jl), so
# the old recipe is not reused by mistake.

POISSON1D_SCHEMES = Dict{String,NamedTuple}(
    "convFD3" => (
        status = :validated, version = "v1-2026-10-29",
        from = "notebooks/poisson1Dbenchmark_small.ipynb",
        orderBtime = 1, orderBspace = -1, pointsInSpace = 3, pointsInTime = 1, supplementaryOrder = 0,
        fieldItpl = (ptsSpace = 1, ptsTime = 1, offsetSpace = 1.0, offsetTime = 1, YorderBspace = -1, YorderBtime = 1),
        materItpl = (ptsSpace = 1, ptsTime = 1, offsetSpace = 1.0, offsetTime = 1, YorderBspace = -1, YorderBtime = 1),
    ),
    "convFD5" => (
        status = :validated, version = "v1-2026-10-29",
        from = "notebooks/poisson1Dbenchmark_small.ipynb",
        orderBtime = 1, orderBspace = -1, pointsInSpace = 5, pointsInTime = 1, supplementaryOrder = 0,
        fieldItpl = (ptsSpace = 1, ptsTime = 1, offsetSpace = 2.0, offsetTime = 1, YorderBspace = -1, YorderBtime = 1),
        materItpl = (ptsSpace = 1, ptsTime = 1, offsetSpace = 2.0, offsetTime = 1, YorderBspace = -1, YorderBtime = 1),
    ),
    "OPT3" => (
        status = :validated, version = "v1-2026-10-29",
        from = "notebooks/poisson1Dbenchmark_small.ipynb",
        orderBtime = 1, orderBspace = 1, pointsInSpace = 3, pointsInTime = 1, supplementaryOrder = 2,
        fieldItpl = (ptsSpace = 1, ptsTime = 1, offsetSpace = 1.0, offsetTime = 1, YorderBspace = -1, YorderBtime = 1),
        materItpl = (ptsSpace = 1, ptsTime = 1, offsetSpace = 1.0, offsetTime = 1, YorderBspace = -1, YorderBtime = 1),
    ),
    "OPT3_stag" => (
        status = :experimental, version = "v1-2026-10-29",
        from = "scripts/runPoisson1DBenchmark_mini.jl",
        orderBtime = 1, orderBspace = 1, pointsInSpace = 3, pointsInTime = 1, supplementaryOrder = 2,
        fieldItpl = (ptsSpace = 5, ptsTime = 1, offsetSpace = 0, offsetTime = 1, YorderBspace = 1, YorderBtime = 1),
        materItpl = (ptsSpace = 3, ptsTime = 1, offsetSpace = 1/2, offsetTime = 1, YorderBspace = 1, YorderBtime = 1),
    ),
)

"Scheme parameters by name; refuses experimental ones unless `experimental=true`."
function poisson1DScheme(name::AbstractString; experimental::Bool=false)
    haskey(POISSON1D_SCHEMES, name) ||
        error("unknown scheme '$name'; available: $(join(sort(collect(keys(POISSON1D_SCHEMES))), ", "))")
    s = POISSON1D_SCHEMES[name]
    if s.status != :validated
        experimental || error("scheme '$name' is $(s.status) for Poisson_1D; pass experimental=true to use it")
        @warn "scheme '$name' is $(s.status) for Poisson_1D (from $(s.from))"
    end
    return s
end
