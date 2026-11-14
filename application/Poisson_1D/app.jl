# application/Poisson_1D/app.jl -- 1-D heterogeneous Poisson equation
#
#     d/dx( κ(x) dT/dx ) = f(x)        on [0, L],  T given at both ends
#
# Manufactured solution T(x) = cos(x); f is computed symbolically from κ, so
# every run has an exact answer to compare with. Equation in src/:
# famousEquation "1DpoissonHetero" (src/motorsOPT/famousEquations.jl).
#
# The whole flexOPT pipeline in one function (poisson1DForward):
#   1. model on the grid          -> modelFam  (models, modelName, Δ ...)
#   2. OPT recipe (cached)        -> makeOPTsemiSymbolic via myProduceOrLoad
#   3. numerical operators        -> numericalOperatorConstruction
#   4. linear system + solve      -> prepareLinearSystem, timeMarchingSchemeLinear
#   5. "observe"                  -> field on the grid + error vs exact
#
# Requires: FLEXOPT_MODULES = [:flexOPT], and Symbolics in Main.

using Symbolics, LinearAlgebra

include(joinpath(@__DIR__, "schemes.jl"))

# conductivity models κ(ξ) (ξ is a Symbolics variable); same cases as the notebook
POISSON1D_CASES = Dict{String,Function}(
    "homogeneous"       => ξ -> 1.0 + 0 * ξ,
    "sameλ"             => ξ -> sin(ξ) + 2,
    "twiceλ"            => ξ -> sin(ξ / 2) + 2,
    "sameλ_shifted_π_3" => ξ -> sin(ξ + π / 3) + 2,
    "λ_2"               => ξ -> cos(ξ)^2 + 1,
    "quadratic"         => ξ -> ξ^2 + 1,
)

"""
    poisson1DModel(case; h, L=10π) -> (; X, κ, f, exact, modelFam)

Grid, conductivity, source term and exact solution for one case.
`h` is the requested spacing; the actual Δx divides L exactly.
"""
function poisson1DModel(case::AbstractString; h::Real, L::Real=10π)
    haskey(POISSON1D_CASES, case) ||
        error("unknown case '$case'; available: $(join(sort(collect(keys(POISSON1D_CASES))), ", "))")
    @variables ξ
    D = Differential(ξ)
    κsym = POISSON1D_CASES[case](ξ)
    Tsym = cos(ξ)
    fsym = expand_derivatives(D(κsym * D(Tsym)))

    Nx = Int(L ÷ h) + 1
    Δx = L / (Nx - 1)
    X = [Δx * (i - 1) for i in 1:Nx]
    at(expr, xi) = Float64(Symbolics.value(substitute(expr, Dict(ξ => xi))))
    κ = [at(κsym, xi) for xi in X]
    f = [at(fsym, xi) for xi in X]
    exact = cos.(X)

    modelFam = (models = [κ], modelName = "poisson1D_$(case)_$(Nx)", modelPoints = (Nx),
                Δ = (Δx), symbols = (Nx = Nx, Δx = Δx, X = X))
    return (; X, κ, f, exact, modelFam, Nx, Δx)
end

"""
    poisson1DRecipe(scheme, Δx; experimental=false)

The OPT recipe for this scheme and spacing, built once and then loaded from
data/semiSymbolic/ (DrWatson cache). The scheme `version` is in the key.
"""
function poisson1DRecipe(scheme::AbstractString, Δx::Real; experimental::Bool=false)
    s = poisson1DScheme(scheme; experimental)
    famousEquationType = "1DpoissonHetero"
    @unpack orderBtime, orderBspace, pointsInSpace, pointsInTime, supplementaryOrder, fieldItpl, materItpl = s
    params = @strdict famousEquationType Δ=Δx orderBtime orderBspace pointsInSpace pointsInTime supplementaryOrder fieldItpl materItpl
    params["applicationSchemeVersion"] = s.version
    # A box-car test function (orderBspace = -1, convFD*) is only allowed in
    # the strong form; OPT schemes keep the package default (DEFAULT_FORM).
    params["variationalForm"] = orderBspace == -1 ? :strong : DEFAULT_FORM[]
    return myProduceOrLoad(makeOPTsemiSymbolic, params, "semiSymbolic")
end

"""
    poisson1DForward(; scheme="OPT3", case="sameλ", h=exp(-2.0), L=10π,
                       precision=Float64, experimental=false) -> NamedTuple

Solve once and return
`(X, T, exact, error, κ, f, scheme, case, Nx, Δx, seconds)` where `error` is
‖T − exact‖₂ / Nx (the metric of notebooks/poisson1Dbenchmark_small.ipynb).
"""
function poisson1DForward(; scheme::AbstractString="OPT3", case::AbstractString="sameλ",
                          h::Real=exp(-2.0), L::Real=10π, precision=Float64, experimental::Bool=false)
    t0 = time()
    m = poisson1DModel(case; h, L)
    optRec = poisson1DRecipe(scheme, m.Δx; experimental)
    numOpt = numericalOperatorConstruction(
        @strdict optRec=optRec modelFam=m.modelFam absorbingBoundaries=nothing maskedRegionInSpace=nothing)
    preparedLin = prepareLinearSystem(numOpt["numOperators"]; T=precision)
    out = timeMarchingSchemeLinear(preparedLin, 1, m.modelFam.Δ, m.modelFam.modelName;
        videoMode=false, sourceType="Explicit", sourceFull=reshape(m.f, m.Nx, 1, 1),
        iExperiment=nothing, boundaryConditionForced=true)
    T = Float64.(vec(reduce(vcat, out)))
    err = norm(T - m.exact) / m.Nx
    return (; X=m.X, T, exact=m.exact, error=err, κ=m.κ, f=m.f, scheme, case,
              Nx=m.Nx, Δx=m.Δx, seconds=time() - t0)
end

"""
    poisson1DConvergence(; schemes=["convFD3","OPT3"], case="sameλ",
                           logInvH=0:0.5:3, kwargs...) -> Vector of NamedTuples

Error versus spacing h = exp(-logInvH) for several schemes, with the
observed order between successive spacings.
"""
function poisson1DConvergence(; schemes=["convFD3", "OPT3"], case::AbstractString="sameλ",
                              logInvH=0:0.5:3, kwargs...)
    rows = NamedTuple[]
    for s in schemes, ℓ in logInvH
        r = poisson1DForward(; scheme=s, case, h=exp(-ℓ), kwargs...)
        push!(rows, (scheme=s, case=case, Δx=r.Δx, Nx=r.Nx, error=r.error, seconds=r.seconds))
        println(rpad(s, 12), " Nx=", lpad(r.Nx, 6), "  Δx=", round(r.Δx, sigdigits=4),
                "  error=", round(r.error, sigdigits=4), "  (", round(r.seconds, digits=2), " s)")
    end
    return rows
end

"Observed orders log(e₁/e₂)/log(h₁/h₂) for one scheme's rows (sorted by Δx)."
function poisson1DOrders(rows, scheme)
    r = sort(filter(x -> x.scheme == scheme, rows); by=x -> -x.Δx)
    return [log(r[i].error / r[i+1].error) / log(r[i].Δx / r[i+1].Δx) for i in 1:length(r)-1]
end
