using KernelAbstractions

opt_integral_order = :ln_lc


function _opt_paramget(params::Dict, name::Symbol, default)
    return get(params, String(name), get(params, name, default))
end

function _normalise_recipe_backend(backend_spec)
    if backend_spec === nothing || backend_spec === :auto || backend_spec == "auto"
        return backend
    elseif backend_spec === :cpu || backend_spec == "cpu"
        return CPU()
    elseif backend_spec === :gpu || backend_spec == "gpu"
        return backend
    elseif backend_spec === :metal || backend_spec == "metal"
        if isdefined(Main, :Metal)
            return getproperty(Main.Metal, :MetalBackend)()
        end
        @warn "recipe_backend=:metal requested, but Metal is not loaded in Main; using current flexOPT.backend" current_backend=typeof(backend)
        return backend
    elseif backend_spec === :cuda || backend_spec == "cuda"
        if isdefined(Main, :CUDA) && isdefined(Main, :CUDABackend)
            return Main.CUDABackend()
        end
        @warn "recipe_backend=:cuda requested, but CUDA/CUDABackend is not loaded in Main; using current flexOPT.backend" current_backend=typeof(backend)
        return backend
    else
        return backend_spec
    end
end

"""
Resolve `recipe_backend = :auto` once the equation is known.

Recipe coefficients are computed in Float32 on a GPU backend (see
`coefficientFloat` below), which leaves a ~1e-7 relative inexactness in every
stencil; it shows up as an accuracy floor ~1e-7/(kh)^2 on fine grids.  So
`:auto` means the CPU (Float64) for equations with up to 3 coordinates, where a
recipe takes seconds, and keeps the current (GPU) backend only for 4
coordinates (3D + t), where the recipe is expensive -- with a warning.  An
explicit backend (:cpu, :metal, :cuda, ...) is honoured unchanged.
"""
function _resolve_auto_recipe_backend(backend_spec, ncoordinates::Integer)
    is_auto = backend_spec === nothing || backend_spec === :auto || backend_spec == "auto"
    is_auto && ncoordinates <= 3 && return CPU()
    selected = _normalise_recipe_backend(backend_spec)
    if is_auto && !(selected isa KernelAbstractions.CPU)
        @warn "recipe_backend=:auto with $(ncoordinates) coordinates: recipe built on $(nameof(typeof(selected))) in Float32 (pass recipe_backend=:cpu for Float64 coefficients)" maxlog=1
    end
    return selected
end

function _recipe_backend_name(recipe_backend)
    recipe_backend isa KernelAbstractions.CPU && return :cpu
    return Symbol(lowercase(string(nameof(typeof(recipe_backend)))))
end

function TaylorOptions(itplParams,supplementaryOrder)
    options=(YorderBspace=itplParams.YorderBspace,YorderBtime=itplParams.YorderBtime,supplementaryOrder=supplementaryOrder,pointsμInSpace=itplParams.ptsSpace,pointsμInTime=itplParams.ptsTime,offsetμInΔyInSpace=itplParams.offsetSpace,offsetμInΔyInTime=itplParams.offsetTime)
    return options
end

function _hierarchical_test_orders(orderBspline)
    maximum_order = maximum(orderBspline)
    maximum_order <= 1 && return [maximum_order]
    return collect(1:maximum_order)
end

function _orders_at_test_level(ordersForSplines, level)
    orderBspline = map(ordersForSplines.orderBspline) do order
        order > 0 ? min(level, order) : order
    end
    return merge(ordersForSplines, (; orderBspline))
end

function _construct_test_blocks(
    equationCharacteristics, numbersOfTheSystem,
    ordersForSplinesμ, configsTaylorμ,
    ordersForSplinesμᶜ, configsTaylorμᶜ,
    Δnum, bigα, varM;
    hierarchical=false,
    even_order_half_shift_mode=:none,
    kwargs...,
)
    levels = hierarchical ?
        _hierarchical_test_orders(ordersForSplinesμ.orderBspline) :
        [maximum(ordersForSplinesμ.orderBspline)]
    results = map(levels) do level
        test_nu_offset =
            iseven(level) && even_order_half_shift_mode != :none ? 0.5 : 0.0
        interpolation_offset =
            iseven(level) && even_order_half_shift_mode == :all ? 0.5 : 0.0
        constructAmatrix(
            equationCharacteristics, numbersOfTheSystem,
            _orders_at_test_level(ordersForSplinesμ, level), configsTaylorμ,
            _orders_at_test_level(ordersForSplinesμᶜ, level), configsTaylorμᶜ,
            Δnum, bigα, varM;
            test_nu_offset=test_nu_offset,
            interpolation_offset=interpolation_offset,
            kwargs...,
        )
    end
    operators = length(results) == 1 ? results[1][1] :
        cat((result[1] for result in results)...; dims=5)
    offsets = [
        iseven(level) && even_order_half_shift_mode != :none ? 0.5 : 0.0
        for level in levels
    ]
    return operators, results[1][2], results[1][3], levels, offsets
end

function _normalise_variational_form(form)
    symbol = Symbol(lowercase(String(form)))
    symbol === :weak && return :natural_weak
    symbol === :strong && return :weighted_residual
    symbol in (:natural_weak, :weighted_residual) && return symbol
    throw(ArgumentError(
        "variationalForm must be :weak, :strong, :natural_weak, or " *
        ":weighted_residual",
    ))
end

_public_variational_form(form::Symbol) =
    form === :natural_weak ? :weak : :strong

function makeOPTsemiSymbolic(params::Dict)
    @unpack famousEquationType, Δ, orderBtime, orderBspace, pointsInSpace, pointsInTime, supplementaryOrder, fieldItpl, materItpl = params
    # resolved below, once the number of coordinates of the equation is known
    recipe_backend_spec = _opt_paramget(params, :recipe_backend, _opt_paramget(params, :coefficient_backend, :auto))
    nuGeometryMode = Symbol(_opt_paramget(params, :nuGeometryMode, :middle))
    nuCentre = _opt_paramget(params, :nuCentre, nothing)
    exactTaylorTotalDegree = _opt_paramget(
        params, :exactTaylorTotalDegree, nothing,
    )
    taylorInverseMode = Symbol(_opt_paramget(
        params, :taylorInverseMode,
        _opt_paramget(params, :taylor_inverse_mode, :hierarchical_constrained),
    ))
    trialFunctionRefPoints = _opt_paramget(params, :trialFunctionRefPoints, nothing)
    testIntegrationBounds = _opt_paramget(params, :testIntegrationBounds, nothing)
    hierarchicalTestFunctions = Bool(_opt_paramget(params, :hierarchicalTestFunctions, false))
    testDerivativeOrders = _opt_paramget(params, :testDerivativeOrders, nothing)
    variationalForm = _normalise_variational_form(
        _opt_paramget(params, :variationalForm, DEFAULT_FORM[]),
    )
    form = _public_variational_form(variationalForm)
    if variationalForm === :natural_weak &&
       (orderBspace == -1 || orderBtime == -1)
        throw(ArgumentError(
            "orderBspace/orderBtime = -1 is available only with " *
            "form=:strong (set_default_form!(:strong) or " *
            "variationalForm=:strong). YorderBspace/YorderBtime = -1 " *
            "remains valid because Y is not the differentiated test function.",
        ))
    end
    variationalForm === :natural_weak && !isnothing(testDerivativeOrders) &&
        throw(ArgumentError(
            "testDerivativeOrders is inferred term by term in :natural_weak mode",
        ))
    evenOrderHalfShiftMode =
        Symbol(_opt_paramget(params, :evenOrderHalfShiftMode, :none))
    evenOrderHalfShiftMode in (:none, :w_only, :all) ||
        throw(ArgumentError("evenOrderHalfShiftMode must be :none, :w_only, or :all"))
    nuGeometryMode in (:middle, :all) ||
        throw(ArgumentError("nuGeometryMode must be :middle or :all"))
    Δnum = SVector(Δ)
    # construction of NamedTuples
    trialFunctionsCharacteristics=(orderBtime=orderBtime,orderBspace=orderBspace,
        pointsInSpace=pointsInSpace,pointsInTime=pointsInTime,
        nuGeometryMode=nuGeometryMode, nuCentre=nuCentre,
        exactTaylorTotalDegree=exactTaylorTotalDegree)

    # here we can compute the different interpolated Taylor expansion options
    TaylorOptionsμ=TaylorOptions(fieldItpl,supplementaryOrder)
    TaylorOptionsμᶜ=TaylorOptions(materItpl,supplementaryOrder)

    equationCharacteristics,equationCharacteristicsForce=famousEquations(famousEquationType)
    recipe_backend = _resolve_auto_recipe_backend(
        recipe_backend_spec, length(equationCharacteristics.coordinates),
    )
    fields=equationCharacteristics.fields
    extfields=equationCharacteristicsForce.fields

    # Compact coefficients for the l.h.s.  In natural-weak mode, preserve the
    # public strong equation and split only its explicit outer spatial
    # divergences into groups carrying different derivatives of W.
    equationCharacteristicsStrong = equationCharacteristics
    weakForm = variationalForm === :natural_weak ?
        naturalWeakForm(equationCharacteristicsStrong) : nothing
    groupedExpressions = if variationalForm === :natural_weak
        weakTermGroups(weakForm)
    else
        orders = isnothing(testDerivativeOrders) ?
            ntuple(_ -> 0, length(equationCharacteristicsStrong.coordinates)) :
            Tuple(Int.(collect(testDerivativeOrders)))
        Dict(orders => collect(equationCharacteristicsStrong.exprs isa Tuple ?
            equationCharacteristicsStrong.exprs :
            (equationCharacteristicsStrong.exprs,)))
    end
    derivativeGroups = sort!(collect(keys(groupedExpressions)); by=identity)

    lhsResults = map(derivativeGroups) do derivativeOrders
        all(>=(0), derivativeOrders) ||
            throw(ArgumentError("test derivative orders must be non-negative"))
        groupedCharacteristics = merge(equationCharacteristicsStrong, (
            exprs=Tuple(groupedExpressions[derivativeOrders]),
        ))
        numbers = numbersOfTheExpression(
            groupedCharacteristics, trialFunctionsCharacteristics,
            TaylorOptionsμ, TaylorOptionsμᶜ,
        )
        _, ordersμ, configsμ, ordersμc, configsμc = investigateDependencies(
            groupedCharacteristics, numbers, trialFunctionsCharacteristics,
            TaylorOptionsμ, TaylorOptionsμᶜ,
        )
        α, materialVariables, dependencies = bigαFinder(
            groupedCharacteristics, numbers, ordersμ,
        )
        matrix, localFields, taylorCoefficients, testOrders, testOffsets =
            _construct_test_blocks(
                groupedCharacteristics, numbers,
                ordersμ, configsμ, ordersμc, configsμc,
                Δnum, α, materialVariables;
                hierarchical=hierarchicalTestFunctions,
                even_order_half_shift_mode=evenOrderHalfShiftMode,
                recipe_backend=recipe_backend,
                taylor_inverse_mode=taylorInverseMode,
                trial_function_ref_points=trialFunctionRefPoints,
                test_integration_bounds=testIntegrationBounds,
                test_derivative_orders=collect(derivativeOrders),
            )
        return (
            matrix=matrix, localFields=localFields,
            taylorCoefficients=taylorCoefficients,
            testOrders=testOrders, testOffsets=testOffsets,
            numbers=numbers, materialVariables=materialVariables,
            dependencies=dependencies, configs=configsμ,
        )
    end
    isempty(lhsResults) && error("the weak form produced no volume term")
    referenceLHS = first(lhsResults)
    all(result -> result.testOrders == referenceLHS.testOrders,
        lhsResults) || error("weak-term test orders differ")
    all(result -> result.testOffsets == referenceLHS.testOffsets,
        lhsResults) || error("weak-term test offsets differ")
    all(result -> size(result.matrix) == size(referenceLHS.matrix),
        lhsResults) || error("weak-term recipe matrices differ in size")

    Ajiννᶜ = copy(referenceLHS.matrix)
    for result in Iterators.drop(lhsResults, 1)
        Ajiννᶜ .+= result.matrix
    end
    Ulocal = referenceLHS.localFields
    Cˡη = referenceLHS.taylorCoefficients
    testFunctionOrders = referenceLHS.testOrders
    testFunctionOffsets = referenceLHS.testOffsets
    numbersOfTheSystemL = referenceLHS.numbers
    varM = referenceLHS.materialVariables
    CartesianDependencies = referenceLHS.dependencies
    configsTaylorμ = referenceLHS.configs
    lhs=(Ajiννᶜ=Ajiννᶜ,Ulocal=Ulocal,varM=varM,CartesianDependencies=CartesianDependencies)
    lhsTestDerivativeOrders = derivativeGroups

    # compact coefficients for r.h.s. of the equation
    equationCharacteristics=equationCharacteristicsForce
    numbersOfTheSystemR=numbersOfTheExpression(equationCharacteristicsForce,trialFunctionsCharacteristics,TaylorOptionsμ,TaylorOptionsμᶜ)
    numbersOfTheSystem = numbersOfTheSystemR
    _,ordersForSplinesμ,configsTaylorμ,ordersForSplinesμᶜ,configsTaylorμᶜ=investigateDependencies(equationCharacteristics,numbersOfTheSystem,trialFunctionsCharacteristics,TaylorOptionsμ,TaylorOptionsμᶜ)
    bigα,varM,CartesianDependencies=bigαFinder(equationCharacteristics,numbersOfTheSystem,ordersForSplinesμ)
    Γjiννᶜ,Flocal,CˡηForce,testFunctionOrdersForce,testFunctionOffsetsForce =_construct_test_blocks(equationCharacteristics,numbersOfTheSystem,ordersForSplinesμ,configsTaylorμ,ordersForSplinesμᶜ,configsTaylorμᶜ,Δnum,bigα,varM; hierarchical=hierarchicalTestFunctions, even_order_half_shift_mode=evenOrderHalfShiftMode, recipe_backend=recipe_backend, taylor_inverse_mode=taylorInverseMode, trial_function_ref_points=trialFunctionRefPoints, test_integration_bounds=testIntegrationBounds, test_derivative_orders=zeros(Int, length(equationCharacteristics.coordinates)))
    testFunctionOrders == testFunctionOrdersForce ||
        error("left and right hierarchical test orders differ")
    testFunctionOffsets == testFunctionOffsetsForce ||
        error("left and right hierarchical test offsets differ")
    rhs=(Γjiννᶜ=Γjiννᶜ,Flocal=Flocal,varF=varM,CartesianDependencies=CartesianDependencies)

    #
    nodes=configsTaylorμ.availablePointsConfigurations
    centresIndices=configsTaylorμ.centrePointConfigurations
    nConfigurations=size(nodes)[1]
    numbersOfTheSystem=(numbersOfTheSystemL=numbersOfTheSystemL,numbersOfTheSystemR=numbersOfTheSystemR,nConfigurations=nConfigurations)
    fieldNames=(fields=fields, extfields=extfields)
    recette=(lhs=lhs,rhs=rhs,nodes=nodes,centresIndices=centresIndices,numbersOfTheSystem=numbersOfTheSystem,fieldNames=fieldNames,Cˡη=Cˡη,CˡηForce=CˡηForce, testFunctionOrders=testFunctionOrders, testFunctionOffsets=testFunctionOffsets, testDerivativeOrders=lhsTestDerivativeOrders, form=form, variationalForm=variationalForm, weakForm=weakForm, hierarchicalTestFunctions=hierarchicalTestFunctions, evenOrderHalfShiftMode=evenOrderHalfShiftMode, taylorInverseMode=taylorInverseMode, recipe_backend=_recipe_backend_name(recipe_backend), recipe_backend_type=string(typeof(recipe_backend)))
    return @strdict(recette)

end

# Cache key of a recipe: record the backend actually used when it is the CPU
# (see _resolve_auto_recipe_backend), so that recipes built on a GPU in
# Float32 under the same parameters are never reloaded as if they were Float64.
# GPU-resolved recipes keep their historical key, so the existing (expensive)
# 4D caches remain valid.
function commonBatchs.cacheKeyExtras(::typeof(makeOPTsemiSymbolic), params)
    spec = _opt_paramget(params, :recipe_backend,
        _opt_paramget(params, :coefficient_backend, :auto))
    equation = _opt_paramget(params, :famousEquationType, nothing)
    isnothing(equation) && return nothing
    characteristics, _ = famousEquations(equation)
    selected = _resolve_auto_recipe_backend(spec, length(characteristics.coordinates))
    selected isa KernelAbstractions.CPU || return nothing
    return Dict{String,Any}("recipe_backend_resolved" => "cpu")
end


function constructAmatrix(equationCharacteristics,numbersOfTheSystem,ordersForSplinesμ,configsTaylorμ,Δnum,bigα,varM;ImakeReport=true, recipe_backend=backend, taylor_inverse_mode=:hierarchical_constrained)
    return constructAmatrix(equationCharacteristics,numbersOfTheSystem,ordersForSplinesμ,configsTaylorμ,ordersForSplinesμ,configsTaylorμ,Δnum,bigα,varM;ImakeReport=ImakeReport, recipe_backend=recipe_backend, taylor_inverse_mode=taylor_inverse_mode)
end

function constructAmatrix(equationCharacteristics,numbersOfTheSystem,ordersForSplinesμ,configsTaylorμ,ordersForSplinesμᶜ,configsTaylorμᶜ,Δnum,bigα,varM;ImakeReport=true, recipe_backend=backend, taylor_inverse_mode=:hierarchical_constrained, trial_function_ref_points=nothing, test_integration_bounds=nothing, test_nu_offset=0.0, interpolation_offset=0.0, test_derivative_orders=nothing)
    numberGeometries = configsTaylorμ.numberGeometries
    numberGeometries == configsTaylorμᶜ.numberGeometries ||
        throw(ArgumentError("field and material Taylor grids must expose the same ν geometries"))

    results = map(1:numberGeometries) do iGeometry
        configμ = merge(configsTaylorμ, (
            numberGeometries=1,
            availablePointsConfigurations=[configsTaylorμ.availablePointsConfigurations[iGeometry]],
            centrePointConfigurations=[configsTaylorμ.centrePointConfigurations[iGeometry]],
            availableμPoints=[configsTaylorμ.availableμPoints[iGeometry]],
            availableμaxes=[configsTaylorμ.availableμaxes[iGeometry]],
        ))
        configμᶜ = merge(configsTaylorμᶜ, (
            numberGeometries=1,
            availablePointsConfigurations=[configsTaylorμᶜ.availablePointsConfigurations[iGeometry]],
            centrePointConfigurations=[configsTaylorμᶜ.centrePointConfigurations[iGeometry]],
            availableμPoints=[configsTaylorμᶜ.availableμPoints[iGeometry]],
            availableμaxes=[configsTaylorμᶜ.availableμaxes[iGeometry]],
        ))
        _constructAmatrix_single(
            equationCharacteristics, numbersOfTheSystem,
            ordersForSplinesμ, configμ, ordersForSplinesμᶜ, configμᶜ,
            Δnum, bigα, varM;
            ImakeReport=ImakeReport, recipe_backend=recipe_backend,
            taylor_inverse_mode=taylor_inverse_mode,
            trial_function_ref_points=trial_function_ref_points,
            test_integration_bounds=test_integration_bounds,
            test_nu_offset=test_nu_offset,
            interpolation_offset=interpolation_offset,
            test_derivative_orders=test_derivative_orders,
        )
    end

    Ajiννᶜ = cat((result[1] for result in results)...; dims=4)
    return Ajiννᶜ, results[1][2], results[1][3]
end

function _constructAmatrix_single(equationCharacteristics,numbersOfTheSystem,ordersForSplinesμ,configsTaylorμ,ordersForSplinesμᶜ,configsTaylorμᶜ,Δnum,bigα,varM;ImakeReport=true, recipe_backend=backend, taylor_inverse_mode=:hierarchical_constrained, trial_function_ref_points=nothing, test_integration_bounds=nothing, test_nu_offset=0.0, interpolation_offset=0.0, test_derivative_orders=nothing)
    
    # for the future develpments: ν can move but it's already more or less coded! look at pointν and nGeometry

    # This function is for one iConfigGeometry

    #region preparation
    @unpack fields,vars = equationCharacteristics
    @unpack nCoordinates,NtypeofExpr,NtypeofExpr,NtypeofFields = numbersOfTheSystem
    @unpack multiOrdersIndices,availablePointsConfigurations,centrePointConfigurations,availableμPoints,availableμaxes, numberGeometries = configsTaylorμ
    exactTaylorTotalDegree = configsTaylorμ.exactTaylorTotalDegree
    availableμᶜPoints = configsTaylorμᶜ.availableμPoints
    availableμᶜaxes = configsTaylorμᶜ.availableμaxes
    
    #@show availableμᶜPoints, availableμPoints

    orderBspline = ordersForSplinesμ.orderBspline
    YorderBsplineμ = ordersForSplinesμ.YorderBspline
    YorderBsplineμᶜ= ordersForSplinesμᶜ.YorderBspline

    nGeometry = 1
    iConfigGeometry = 1
    testDerivativeOrders = isnothing(test_derivative_orders) ?
        zeros(Int, nCoordinates) : Int.(collect(test_derivative_orders))
    length(testDerivativeOrders) == nCoordinates ||
        throw(DimensionMismatch(
            "test_derivative_orders must contain one order per coordinate",
        ))

    @show pointsIndices=availablePointsConfigurations[iConfigGeometry] # CartesianIndices
    @show middleLinearν=centrePointConfigurations[iConfigGeometry] # scalar
    @show μPoints = availableμPoints[iConfigGeometry] # Array(SVector)
    @show μᶜPoints = availableμᶜPoints[iConfigGeometry]
    @show μaxes = availableμaxes[iConfigGeometry]
    @show μᶜaxes = availableμᶜaxes[iConfigGeometry]
    @show size(μPoints)
    @show pointν = pointsIndices[middleLinearν] # SVector

    # In :all mode shift both interpolation centres and the centres used to
    # acquire Cˡη.  This keeps Y and K mutually consistent.
    if !iszero(interpolation_offset)
        μPoints = map(p -> p .+ interpolation_offset, μPoints)
        μᶜPoints = map(p -> p .+ interpolation_offset, μᶜPoints)
        μaxes = map(axis -> axis .+ interpolation_offset, μaxes)
        μᶜaxes = map(axis -> axis .+ interpolation_offset, μᶜaxes)
    end
    
    # this is fully GPU optimised version of ASymbolic 
    
    nPoints = length(pointsIndices)
    nLs = length(multiOrdersIndices)

    # 

    L_MINUS_N = multiOrdersIndices
    L_MINUS_N = L_MINUS_N .-L_MINUS_N[1] 
    # here L_MINUS_N is truly \mathbf{l}-\mahtbf{n} ∈ \mathbb{Z}_{≥0}

    #endregion

    #region we compute the integral of WYYKK in 1D domain 
    
    # look at the debug1DKernelIntegral.ipynb    
    
    coefWYYKK = Array{Any,1}(undef,nCoordinates) # Array of (l_n_max+1,lᶜ_nᶜ_max+1,length(μs),length(μᶜs),length(ν) ) times nCoordinates
    axis_option(option, iCoord) = begin
        isnothing(option) && return nothing
        if option isa Tuple && length(option) == nCoordinates &&
           all(value -> value isa Tuple || value isa AbstractVector, option)
            return option[iCoord]
        end
        return option
    end
    for iCoord ∈ 1:nCoordinates
        maxNodes = pointsIndices[end][iCoord]
        nodesFromOne = [1,2,3] # ∈ Z like [1,2,3], an array of integers collect(1:1:Npoints) (nothing else!!)
        ν = pointν[iCoord]
        lᶜ_nᶜ_max = L_MINUS_N[end][iCoord] # variable
        l_n_max = L_MINUS_N[end][iCoord] # field
        referencePointsHere = axis_option(trial_function_ref_points, iCoord)
        νRefPoints = referencePointsHere === nothing ?
            collect(1:maxNodes) :
            collect(referencePointsHere)
        integrationBounds = axis_option(test_integration_bounds, iCoord)
        # WνOffset is deliberately applied inside WYYKK only to the test
        # family Wν.  The Yμ/Yμᶜ interpolation families and Kμ/Kμᶜ Taylor
        # kernels remain on their original field/material μ grids.
        params = (orderBspline1D=orderBspline[iCoord], YorderBspline1Dμᶜ=YorderBsplineμᶜ[iCoord], YorderBspline1Dμ=YorderBsplineμ[iCoord], μᶜs=μᶜaxes[iCoord], μs=μaxes[iCoord], maxNode = pointsIndices[end][iCoord], ν=ν, νRefPoints=νRefPoints, WνOffset=test_nu_offset, WDerivativeOrder=testDerivativeOrders[iCoord], integrationBounds=integrationBounds, lᶜ_nᶜ_max=lᶜ_nᶜ_max, l_n_max=l_n_max,  Δ=Δnum[iCoord],ImakeReport=ImakeReport)
        coefWYYKK[iCoord] = WYYKKIntegralNumerical(params) 
    end

    #endregion

    #region get CˡηGlobal (for ν)

    @show typeof(μPoints), μPoints[1], typeof(pointsIndices)

    coefInversionDict = Dict{String,Any}(@strdict multiOrdersIndices pointsIndices μpointsIndices=μPoints Δ=Δnum)
    coefInversionDict["pinv_version"] = "$(taylor_inverse_mode)_hierarchical_v1"
    coefInversionDict["taylor_inverse_mode"] = taylor_inverse_mode
    coefInversionDict["exact_taylor_total_degree"] = exactTaylorTotalDegree
    output=myProduceOrLoad(TaylorCoefInversion,coefInversionDict,"taylorCoefInv")
    Cˡη=output["CˡηGlobal"]

    coefInversionDict = Dict{String,Any}(@strdict multiOrdersIndices pointsIndices μpointsIndices=μᶜPoints Δ=Δnum)
    coefInversionDict["pinv_version"] = "$(taylor_inverse_mode)_hierarchical_v1"
    coefInversionDict["taylor_inverse_mode"] = taylor_inverse_mode
    coefInversionDict["exact_taylor_total_degree"] = exactTaylorTotalDegree
    output=myProduceOrLoad(TaylorCoefInversion,coefInversionDict,"taylorCoefInv")
    Cˡηᶜ=output["CˡηGlobal"]

    #endregion 

    #region 
    @show nTotalSmallα = sum(length(bigα[iExpr, iField]) for iField ∈ 1:NtypeofFields, iExpr ∈ 1:NtypeofExpr)

    # but this will already include bigα so the coefficients for each α_{nn'ji} should be given here
    #endregion


    #region useful LinearIndices conversion functions
    
    LI_points = LinearIndices(pointsIndices)
    LI_L_MINUS_N_plus_1 = LinearIndices(L_MINUS_N.+vec2car(ones(Int,nCoordinates)))
  
    #endregion


    #region make the table for each (x',x,n',n) (x=η+μ)
    
    tableForLoop = Array{Int32,3}(undef,2+nCoordinates*2,nLs*nLs,nTotalSmallα)

    fill!(tableForLoop, 0)
    indexLinearα = 1
    for iExpr ∈ 1:NtypeofExpr,iField ∈ 1:NtypeofFields
        α = bigα[iExpr,iField]
        for eachα ∈ α
            nᶜ = eachα.nᶜ - vec2car(ones(Int,nCoordinates))
            n = eachα.n - vec2car(ones(Int,nCoordinates))
            #nodeValue = eachα.node # not important at this point
            # Available indices
            Lᶜ_avail = (nᶜ .+ L_MINUS_N) ∩ L_MINUS_N
            L_avail = (n .+ L_MINUS_N) ∩ L_MINUS_N
            Lᶜ_Nᶜ_avail = Lᶜ_avail .- nᶜ
            L_N_avail = L_avail .- n
            
            iL=1
            
            for l ∈ L_avail, lᶜ ∈ Lᶜ_avail
                
                tableForLoop[1,iL,indexLinearα]= LI_L_MINUS_N_plus_1[lᶜ+vec2car(ones(Int,nCoordinates))]
                tableForLoop[2,iL,indexLinearα]= LI_L_MINUS_N_plus_1[l+vec2car(ones(Int,nCoordinates))]
                #tableForLoop[3,iL,indexLinearα]= LI_L_MINUS_N_plus_1[lᶜ-nᶜ+vec2car(ones(Int,nCoordinates))]
                #tableForLoop[4,iL,indexLinearα]= LI_L_MINUS_N_plus_1[l-n+vec2car(ones(Int,nCoordinates))]
                tmplᶜ_nᶜ = lᶜ-nᶜ+vec2car(ones(Int,nCoordinates))
                tmpl_n = l-n+vec2car(ones(Int,nCoordinates))
                for iCoord ∈ 1:nCoordinates
                    tableForLoop[2+iCoord,iL,indexLinearα] = tmplᶜ_nᶜ[iCoord]
                    tableForLoop[2+iCoord+nCoordinates,iL,indexLinearα] = tmpl_n[iCoord]
                end

                iL += 1
            end
       
            indexLinearα += 1
        end
    end

    #endregion


    #region make a dictionary for μ ∈ μPoints and its linearised version

    tableForμPoints = Array{Int32,2}(undef,nCoordinates,length(μPoints))
    linearμPoints = LinearIndices(μPoints)
    for iμ ∈ CartesianIndices(μPoints), iCoord ∈ 1:nCoordinates
        tableForμPoints[iCoord, linearμPoints[iμ]] = iμ[iCoord]
    end

    tableForμᶜPoints = Array{Int32,2}(undef,nCoordinates,length(μᶜPoints))
    linearμᶜPoints = LinearIndices(μᶜPoints)
    for iμᶜ ∈ CartesianIndices(μᶜPoints), iCoord ∈ 1:nCoordinates        
        tableForμᶜPoints[iCoord, linearμᶜPoints[iμᶜ]] = iμᶜ[iCoord]
    end
    


    #endregion


    
    selected_backend = _normalise_recipe_backend(recipe_backend)
    coefficientFloat = selected_backend isa KernelAbstractions.CPU ? Float64 : Float32

    #region adapt the arrays to the GPU backend
    tableForμPoints_gpu = Adapt.adapt(selected_backend,tableForμPoints)
    tableForμᶜPoints_gpu = Adapt.adapt(selected_backend,tableForμᶜPoints)
    tableForLoop_gpu = Adapt.adapt(selected_backend,tableForLoop)
    C_gpu = Adapt.adapt(selected_backend,coefficientFloat.(Cˡη))
    Cᶜ_gpu = Adapt.adapt(selected_backend,coefficientFloat.(Cˡηᶜ))

    #endregion

    #region preparation for GPU computation
    # Collect the size of each array
    @show all_sizes = collect.(size.(coefWYYKK))  # vector of tuples

    # Element-wise maximum across all dimensions
    max_size = map((xs...) -> maximum(xs), all_sizes...)

    int_total = Array{coefficientFloat,6}(undef, nCoordinates, max_size...)
    int_total .= zero(coefficientFloat)

    for iCoord ∈ 1:nCoordinates
        small_size = size(coefWYYKK[iCoord])
        tmpMatrix = coefficientFloat.(coefWYYKK[iCoord])
        tmpRange = CartesianIndices(tmpMatrix)
        int_total[iCoord,tmpRange] = tmpMatrix
    end
    int_gpu = Adapt.adapt(selected_backend,int_total)

    output_gpu = Adapt.adapt(selected_backend, zeros(coefficientFloat, nPoints, nPoints, nTotalSmallα, nGeometry))

    #endregion



    #region launch GPU computation

    # scalars in Int32
    P      = Int32(nPoints)
    Pμᶜ    = Int32(length(μᶜPoints))
    Pμ     = Int32(length(μPoints))
    L      = Int32(nLs)
    nDim   = Int32(nCoordinates)
    nα     = Int32(nTotalSmallα)
    nGeometry_gpu=Int32(nGeometry)
    int1max = Int32(max_size[1])
    int2max = Int32(max_size[2])

    @show typeof(output_gpu)         # must be MtlArray{Float32,3}
    #@show typeof(C_gpu)              # must be MtlArray{Float32,3}
    #@show typeof(int_gpu)            # must be MtlArray{Float32,5}
    #@show typeof(tableForLoop_gpu)   # must be MtlArray{Int32,3}
    #@show typeof(tableForPoints_gpu) # must be MtlMatrix{Int32}
    #@show typeof(P) typeof(L) typeof(nDim) typeof(nα) typeof(int1max) typeof(int2max)

    integralOrderFlag = (isdefined(@__MODULE__, :opt_integral_order) && opt_integral_order == :lc_ln) ? Int32(2) : Int32(1)
    kernel! = windowContraction!(selected_backend,(8,8,8))#,128,size(output_gpu))
    kernel!(output_gpu, C_gpu, Cᶜ_gpu, int_gpu, tableForLoop_gpu,tableForμPoints_gpu,tableForμᶜPoints_gpu,
       P,Pμᶜ,Pμ,L,nDim,nα,int1max,int2max,nGeometry_gpu,integralOrderFlag,ndrange=(P,P,nα*nGeometry))
    KernelAbstractions.synchronize(selected_backend)

    # Here output_gpu[x',x,eachα] ∑μ' ∑μ ∑l' ∑ l C[x',l',μ'] C[x,l,μ] ∏_iCoord K[iCoord][l-n,l'-n',μ,μ']
    # with (n',n) depends on eachα and x=η+μ0

    @show "GPU computation of Ajiννᶜ: done"
    newCoef = Adapt.adapt(CPU(), output_gpu)
    #endregion


    iConfigGeometry = 1

    #region contruct Ulocal

    # the order is: (νᶜ,) ν, i, j  here

    Ulocal = Array{Num,2}(undef,length(pointsIndices),NtypeofFields)
    for iField in eachindex(fields)
        newstring=split(string(fields[iField]),"(")[1]
        Ulocal[:,iField]=Symbolics.variables(Symbol(newstring),1:length(pointsIndices))
    end

    #endregion

    #region make Ajiννᶜ and AjiννᶜU symbolically (which we will soon remove!)
    Ajiννᶜ = Array{Num,4}(undef,length(pointsIndices),NtypeofFields,NtypeofExpr,nGeometry)
    #AjiννᶜU = Array{Num,2}(undef,NtypeofExpr,nGeometry)
    
    # this is the cost function for ν point so the number of elements is just the number of expressions (governing equations)
    Ajiννᶜ .= 0
    #AjiννᶜU .= 0
    indexLinearα = 1
   
    for iExpr ∈ 1:NtypeofExpr,iField ∈ 1:NtypeofFields
        α = bigα[iExpr,iField]
        for eachα ∈ α
            @show nodeValue=eachα.node
            for x ∈ pointsIndices
                xLinear = LI_points[svec2car(x)]
                
                localmapηᶜ=Dict()
                for iVar ∈ eachindex(vars)
                    localmapηᶜ[vars[iVar]]=varM[iVar,xLinear][]
                end
                for xᶜ ∈ pointsIndices
                    xᶜLinear = LI_points[svec2car(xᶜ)]
                    U_HERE = Ulocal[xᶜLinear,iField]                    
                    substitutedValue = substitute(nodeValue, localmapηᶜ)
                    # Equation (54): x is the material point (µ′+η′) and xᶜ is the final field column ν′.
                    # windowContraction! stores axes as (material point, field point, α, geometry).
                    Ajiννᶜ[xᶜLinear,iField,iExpr,iConfigGeometry] += newCoef[xLinear,xᶜLinear,indexLinearα,iConfigGeometry] * substitutedValue
                    #AjiννᶜU[iExpr,iConfigGeometry]+= Ajiννᶜ[xᶜLinear,iField,iExpr,iConfigGeometry] * U_HERE
                end

            end
            indexLinearα += 1
        end

    end

    #return coefWYYKK,Cˡη,tableForμPoints, newCoef
    return Ajiννᶜ,Ulocal,Cˡη
end






@kernel function windowContraction!(
    output::AbstractArray{T,4},            # (P, P, nα, nGeometry)
    C::AbstractArray{T,3},                 # (P, L, Pμ)
    Cc::AbstractArray{T,3},                # (P, L, Pμᶜ)
    int::AbstractArray{T,6},               # (nDim, int1max, int1max, int2max, int2max, nGeometry)
    table::AbstractArray{Int32,3},         # (2 + 2*nDim, nLoop, nα)
    tableμPoints::AbstractArray{Int32,2},  # (nDim, Pμ)
    tableμᶜPoints::AbstractArray{Int32,2}, # (nDim, Pμᶜ)
    P::Int32, Pμᶜ::Int32, Pμ::Int32, L::Int32,
    nDim::Int32, nα::Int32, int1max::Int32, int2max::Int32, nGeometry::Int32,
    integralOrderFlag::Int32
) where {T}
    (xᶜ, x, ag) = @index(Global, NTuple)

    if xᶜ <= P && x <= P && ag <= nα * nGeometry
        α = Int32(((ag - 1) % nα) + 1)
        iGeometry = Int32(((ag - 1) ÷ nα) + 1)

        acc = zero(T)

        @inbounds for idx in 1:size(table, 2)
            lᶜ = table[1, idx, α]
            l  = table[2, idx, α]

            if lᶜ > 0 && l > 0
                for μᶜ in 1:Pμᶜ
                    for μ in 1:Pμ
                        prod_int = one(T)

                        for iDim in 1:nDim
                            k = table[2 + iDim, idx, α]
                            lp = table[2 + nDim + iDim, idx, α]

                            if k > 0 && lp > 0
                                mᶜ = tableμᶜPoints[iDim, μᶜ]
                                m  = tableμPoints[iDim, μ]

                                # coefWYYKK is indexed as (l_n, lᶜ_nᶜ, μ, μᶜ, ν).
                                # For debugging the OPT convention, allow checking the transposed
                                # l/lᶜ read without changing the rest of the assembly.
                                if integralOrderFlag == Int32(2)
                                    prod_int *= int[iDim, k, lp, m, mᶜ, iGeometry]
                                else
                                    prod_int *= int[iDim, lp, k, m, mᶜ, iGeometry]
                                end
                            else
                                prod_int = zero(T)
                                break
                            end
                        end

                        acc += Cc[xᶜ, lᶜ, μᶜ] * C[x, l, μ] * prod_int
                    end
                end
            end
        end

        output[xᶜ, x, α, iGeometry] = acc
    end
end
