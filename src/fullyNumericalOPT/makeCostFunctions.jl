using SparseArrays

function CerjanBoundaryCondition(distance2; damping=0.0053)
    return exp(-damping * Float64(distance2))
end

struct CouplingTable{T}
    rows::Vector{Int32}
    cols::Vector{Int32}
    vals::Vector{T}
end

struct MatrixFreeNumericalOperator{T}
    table::CouplingTable{T}
    size::Tuple{Int,Int}
    side::String
    backend::Symbol
    representation::Symbol
    geometry
end

struct SparseNumericalOperator{T,S}
    matrix::S
    table::CouplingTable{T}
    side::String
    backend::Symbol
    representation::Symbol
    geometry
end

function applyNumericalOperator!(y, op::MatrixFreeNumericalOperator, x)
    fill!(y, zero(eltype(y)))
    rows = op.table.rows
    cols = op.table.cols
    vals = op.table.vals
    @inbounds for k in eachindex(vals)
        y[rows[k]] += vals[k] * x[cols[k]]
    end
    return y
end

function Base.:-(left::MatrixFreeNumericalOperator{TL}, right::MatrixFreeNumericalOperator{TR}) where {TL,TR}
    left.size == right.size || throw(DimensionMismatch("matrix-free operators must have the same size"))
    T = promote_type(TL, TR)
    rows = vcat(left.table.rows, right.table.rows)
    cols = vcat(left.table.cols, right.table.cols)
    vals = vcat(T.(left.table.vals), .-T.(right.table.vals))
    table = CouplingTable(rows, cols, vals)
    return MatrixFreeNumericalOperator(
        table,
        left.size,
        "$(left.side)-$(right.side)",
        left.backend,
        left.representation,
        left.geometry,
    )
end

function applyNumericalOperator!(y, op::SparseNumericalOperator, x)
    mul!(y, op.matrix, x)
    return y
end

function _residual_numerical_operator(left, right)
    left.size[1] == right.size[1] ||
        throw(DimensionMismatch("left and right numerical operators must have the same number of rows"))
    T = promote_type(eltype(left.table.vals), eltype(right.table.vals))
    rows = vcat(left.table.rows, right.table.rows)
    cols = vcat(left.table.cols, Int32.(right.table.cols .+ left.size[2]))
    vals = vcat(T.(left.table.vals), .-T.(right.table.vals))
    table = CouplingTable(rows, cols, vals)
    return MatrixFreeNumericalOperator(
        table,
        (left.size[1], left.size[2] + right.size[2]),
        "$(left.side)-$(right.side)",
        left.backend,
        :matrixfree,
        (left=left.geometry, right=right.geometry),
    )
end

function _normalise_residual_rows(operator::MatrixFreeNumericalOperator)
    row_norm2 = zeros(Float64, operator.size[1])
    for (row, value) in zip(operator.table.rows, operator.table.vals)
        row_norm2[row] += abs2(value)
    end
    row_scales = map(row_norm2) do norm2
        norm2 > eps(Float64) ? inv(sqrt(norm2)) : 1.0
    end
    vals = similar(operator.table.vals)
    for index in eachindex(vals)
        vals[index] = operator.table.vals[index] * row_scales[operator.table.rows[index]]
    end
    table = CouplingTable(copy(operator.table.rows), copy(operator.table.cols), vals)
    normalised = MatrixFreeNumericalOperator(
        table, operator.size, operator.side, operator.backend,
        operator.representation, operator.geometry,
    )
    return normalised, row_scales
end

_as_symbol(x::Symbol) = x
_as_symbol(x::AbstractString) = Symbol(x)

_field_symbol(fieldSyms::AbstractArray, i) = fieldSyms[i]
_field_symbol(fieldSyms, i) = fieldSyms

function _numerical_value(x)
    v = Symbolics.value(x)
    return v isa Number ? v : error("non-numerical coefficient remains after material substitution: $v")
end

function _coefficient_vector(coefficient_type)
    type_sym = _as_symbol(coefficient_type)
    if type_sym == :auto || type_sym == :float || type_sym == :float64
        return Float64[]
    elseif type_sym == :complex || type_sym == :complex64 || type_sym == :complexf64
        return ComplexF64[]
    else
        throw(ArgumentError("coefficient_type must be :auto, :float64, or :complexf64; got $coefficient_type"))
    end
end

function _push_coupling!(rows, cols, vals::Vector{Float64}, row, col, val)
    v = _numerical_value(val)
    iszero(v) && return vals

    if v isa Complex
        if iszero(imag(v))
            v = real(v)
        else
            vals_complex = ComplexF64.(vals)
            push!(rows, Int32(row))
            push!(cols, Int32(col))
            push!(vals_complex, ComplexF64(v))
            return vals_complex
        end
    end

    push!(rows, Int32(row))
    push!(cols, Int32(col))
    push!(vals, Float64(v))
    return vals
end

function _push_coupling!(rows, cols, vals::Vector{ComplexF64}, row, col, val)
    v = _numerical_value(val)
    iszero(v) && return vals
    push!(rows, Int32(row))
    push!(cols, Int32(col))
    push!(vals, ComplexF64(v))
    return vals
end

function _finalize_numerical_operator(table::CouplingTable{T}, nrow::Int, ncol::Int, side, backend, representation, geometry) where T
    backend_sym = _as_symbol(backend)
    repr_sym = representation === :auto ? :sparse : _as_symbol(representation)

    if backend_sym == :mpi
        repr_sym = repr_sym == :auto ? :matrixfree : repr_sym
        return MatrixFreeNumericalOperator(table, (nrow, ncol), side, backend_sym, repr_sym, geometry)
    elseif repr_sym == :matrixfree
        return MatrixFreeNumericalOperator(table, (nrow, ncol), side, backend_sym, repr_sym, geometry)
    elseif repr_sym == :sparse || repr_sym == :blocksparse
        A = sparse(Int.(table.rows), Int.(table.cols), table.vals, nrow, ncol)
        if backend_sym == :cuda
            if isdefined(Main, :CUDA)
                A = Main.CUDA.CUSPARSE.CuSparseMatrixCSR(A)
            else
                @warn "backend=:cuda requested, but CUDA is not loaded in Main; returning CPU sparse matrix"
                backend_sym = :cpu
            end
        elseif backend_sym != :cpu
            throw(ArgumentError("backend must be :cpu, :cuda, or :mpi; got $backend"))
        end
        return SparseNumericalOperator(A, table, side, backend_sym, repr_sym, geometry)
    else
        throw(ArgumentError("representation must be :auto, :matrixfree, :sparse, or :blocksparse; got $representation"))
    end
end

struct CoefficientRecipe{F,V}
    vars::V
    f::F
    constant::Any
end

struct ComplexCoefficientRecipe{R,I}
    real_recipe::R
    imag_recipe::I
end

function _compile_coefficient_recipe(expr::Complex)
    return ComplexCoefficientRecipe(
        _compile_coefficient_recipe(real(expr)),
        _compile_coefficient_recipe(imag(expr)),
    )
end

function _compile_coefficient_recipe(expr)
    vars = collect(Symbolics.get_variables(expr))
    if isempty(vars)
        return CoefficientRecipe(vars, nothing, _numerical_value(expr))
    end
    f = Symbolics.build_function(expr, vars; expression=Val(false))
    return CoefficientRecipe(vars, f, nothing)
end

function _compile_coefficient_recipes(symA)
    cache = IdDict{Any,Any}()
    recipes = Array{Any}(undef, size(symA))
    for idx in eachindex(symA)
        expr = symA[idx]
        recipes[idx] = get!(cache, expr) do
            _compile_coefficient_recipe(expr)
        end
    end
    return recipes
end

function _evaluate_recipe(recipe::CoefficientRecipe, materialMapping::Dict)
    if isempty(recipe.vars)
        return recipe.constant
    end
    return recipe.f([materialMapping[v] for v in recipe.vars])
end

function _evaluate_recipe(recipe::ComplexCoefficientRecipe, materialMapping::Dict)
    return complex(
        _evaluate_recipe(recipe.real_recipe, materialMapping),
        _evaluate_recipe(recipe.imag_recipe, materialMapping),
    )
end

function _make_field_arrays(fieldSyms, NtypeofFields, activeTimePoints, νWhole, wholeRegionPointsSpace)
    fieldArrays = Array{Any,2}(undef, NtypeofFields, activeTimePoints)
    for iT in 1:activeTimePoints
        for iField in 1:NtypeofFields
            fieldArrays[iField, iT] = reshape(
                [_field_symbol(fieldSyms, iField) for _ in 1:length(νWhole)],
                Tuple(wholeRegionPointsSpace),
            )
        end
    end
    return fieldArrays
end


function numericalOperatorConstruction(params::Dict)
    @unpack optRec,modelFam,absorbingBoundaries,maskedRegionInSpace=params
    boundaryConditions = _paramget(params, :boundaryConditions, nothing)
    backend = _paramget(params, :backend, :cpu)
    representation = _paramget(params, :representation, :auto)
    coefficient_type = _paramget(params, :coefficient_type, :auto)
    compatibility_outputs = _paramget(params, :compatibility_outputs, false)
    if !hasproperty(modelFam, :modelName) || modelFam.modelName === nothing
        modelFam = merge(modelFam, (; modelName = "model_" * Dates.format(now(), "yyyymmdd_HHMMSS")))
    end
    freeSurfaceBoundary = nothing
    cerjanDamping = 0.0053
    if !isnothing(boundaryConditions)
        freeSurfaceBoundary = boundaryConditions.free_surface
        if !isnothing(boundaryConditions.cerjan)
            isnothing(absorbingBoundaries) ||
                throw(ArgumentError(
                    "give Cerjan padding through boundaryConditions or " *
                    "absorbingBoundaries, not both",
                ))
            absorbingBoundaries = cerjan_padding(boundaryConditions.cerjan)
            cerjanDamping = boundaryConditions.cerjan.damping
        end
    end
    costLHS=numericalOperatorConstruction(optRec,modelFam,"left";absorbingBoundaries=absorbingBoundaries,maskedRegionInSpace=maskedRegionInSpace,freeSurfaceBoundary=freeSurfaceBoundary,cerjanDamping=cerjanDamping,backend=backend,representation=representation,coefficient_type=coefficient_type,compatibility_outputs=compatibility_outputs)
    costRHS=numericalOperatorConstruction(optRec,modelFam,"right";absorbingBoundaries=absorbingBoundaries,maskedRegionInSpace=maskedRegionInSpace,freeSurfaceBoundary=freeSurfaceBoundary,cerjanDamping=cerjanDamping,backend=backend,representation=representation,coefficient_type=coefficient_type,compatibility_outputs=compatibility_outputs)
    costfunctions=_residual_numerical_operator(costLHS.operator,costRHS.operator)
    hierarchical = hasproperty(optRec["recette"], :hierarchicalTestFunctions) &&
        optRec["recette"].hierarchicalTestFunctions
    rowNormalisation = ones(Float64, costfunctions.size[1])
    if hierarchical
        costfunctions, rowNormalisation = _normalise_residual_rows(costfunctions)
    end
    fieldLHS=costLHS.場
    fieldRHS=costRHS.場
    champsLimité=costRHS.champsLimité
    numericalOperators=(left=costLHS.operator,right=costRHS.operator,residual=costfunctions)
    numOperators=(costfunctions=costfunctions,fieldLHS=fieldLHS,fieldRHS=fieldRHS,champsLimité=champsLimité,numericalOperators=numericalOperators,rowNormalisation=rowNormalisation,boundaryConditions=boundaryConditions)
    return @strdict(numOperators)
end

function numericalOperatorConstruction(optRec,modelFam,side;absorbingBoundaries=nothing,maskedRegionInSpace=nothing,
    freeSurfaceBoundary=nothing,
    cerjanDamping=0.0053,
    backend=:cpu,        # :cpu, :cuda, :mpi
    representation=:auto, # :matrixfree, :sparse, :blocksparse
    coefficient_type=:auto, # :auto, :float64, :complexf64
    compatibility_outputs=false
)

    numericalGeometry = prepareNumericalOperatorGeometry(
        optRec, modelFam, side;
        absorbingBoundaries, maskedRegionInSpace, freeSurfaceBoundary,
    )

    @unpack νWhole, νGeometry, νRelative, localPointsIndices,
        ModelPoints, Models, maskingField, conv,
        timePointsUsedForOneStep, activeTimePoints, wholeMin, wholeMax,
        modelMin, modelDomainMax, wholeRegionPointsSpace = numericalGeometry

    @unpack lhs, rhs, numbersOfTheSystem, fieldNames = optRec["recette"]
    @unpack fields, extfields = fieldNames

    if side == "left"
        numbersSide = numbersOfTheSystem.numbersOfTheSystemL
        symA = lhs.Ajiννᶜ
        varM = lhs.varM
        fieldSyms = fields
    elseif side == "right"
        numbersSide = numbersOfTheSystem.numbersOfTheSystemR
        symA = rhs.Γjiννᶜ
        varM = rhs.varF
        fieldSyms = extfields
    else
        throw(ArgumentError("side must be \"left\" or \"right\", got $side"))
    end

    @unpack NtypeofExpr, NtypeofFields = numbersSide
    NtypeofCoefficientVariables = size(varM, 1)
    nTestOrderBlocks = ndims(symA) >= 5 ? size(symA, 5) : 1

    coefficientRecipes = _compile_coefficient_recipes(symA)
    fieldLinearIndices = LinearIndices(
        (NtypeofFields, activeTimePoints, Tuple(wholeRegionPointsSpace)...)
    )
    residualLinearIndices = LinearIndices((NtypeofExpr, length(νWhole), nTestOrderBlocks))

    rows = Int32[]
    cols = Int32[]
    vals = _coefficient_vector(coefficient_type)

    for iTestFunctions in eachindex(νWhole)
        νtmpWhole = νWhole[iTestFunctions]
        iGeometry = νGeometry[iTestFunctions]
        localPointsHere = localPointsIndices[iGeometry]
        middlepointHere = νRelative[iTestFunctions]
        localPointsSpaceHere = carDropDim.(localPointsHere)
        localPointsSpaceIndicesHere = CartesianIndices(Tuple(localPointsSpaceHere[end]))
        middlepointSpaceHere = svec2car(carDropDim(middlepointHere))
        iTimeMax = timePointsUsedForOneStep[iGeometry]

        νᶜtmpWhole = localPointsSpaceIndicesHere .+ (νtmpWhole .- middlepointSpaceHere)
        νᶜtmpModel = conv.whole2model.(νᶜtmpWhole)
        localLinearIndices = LinearIndices(Tuple(localPointsHere[end]))

        materialMapping = Dict{Any,Any}()
        for iT in 1:iTimeMax
            for iVar in 1:min(NtypeofCoefficientVariables, length(Models))
                spaceModelBouncedPoints = ModelPoints[1:end-1, iVar]
                iiT = ModelPoints[end, iVar] > 1 ? iT : 1
                νᶜtmpModelTruncated = BouncingCoordinates.(νᶜtmpModel, Ref(spaceModelBouncedPoints))

                for jPoint in νᶜtmpWhole
                    jPointLocal = jPoint - νtmpWhole + middlepointSpaceHere
                    jPointTLocal = carAddDim(jPointLocal, iT)
                    linearjPointTLocal = localLinearIndices[jPointTLocal]
                    materialMapping[varM[iVar, linearjPointTLocal]] =
                        Models[iVar][carAddDim(νᶜtmpModelTruncated[jPointLocal], iiT)]
                end
            end
        end

        for iTestOrderBlock in 1:nTestOrderBlocks, iExpr in 1:NtypeofExpr
            row = residualLinearIndices[iExpr, iTestFunctions, iTestOrderBlock]

            for iT in 1:iTimeMax
                for jPoint in νᶜtmpWhole
                    jPointLocal = jPoint - νtmpWhole + middlepointSpaceHere
                    jPointTLocal = carAddDim(jPointLocal, iT)
                    linearjPointTLocal = localLinearIndices[jPointTLocal]

                    jPointInWhole = is_all_less_than_or_equal(wholeMin, jPoint) &&
                                    is_all_less_than_or_equal(jPoint, wholeMax)
                    boundaryWeight = 0.0

                    if jPointInWhole
                        jPointModel = conv.whole2model(jPoint)
                        jPointInModel =
                            is_all_less_than_or_equal(modelMin, jPointModel) &&
                            is_all_less_than_or_equal(jPointModel, modelDomainMax)

                        if jPointInModel || iT == iTimeMax
                            boundaryWeight = 1.0
                        else
                            distance2 = distance2_point_to_box(
                                jPointModel,
                                modelMin,
                                modelDomainMax,
                            )
                            boundaryWeight = CerjanBoundaryCondition(
                                distance2;
                                damping=cerjanDamping,
                            )
                        end
                    end

                    weight = jPointInWhole ? maskingField[jPoint] * boundaryWeight : 0.0
                    iszero(weight) && continue

                    for iField in 1:NtypeofFields
                        col = fieldLinearIndices[iField, iT, Tuple(jPoint)...]
                        recipe = nTestOrderBlocks == 1 && ndims(symA) == 4 ?
                            coefficientRecipes[linearjPointTLocal, iField, iExpr, iGeometry] :
                            coefficientRecipes[linearjPointTLocal, iField, iExpr, iGeometry, iTestOrderBlock]
                        coef = weight * _evaluate_recipe(recipe, materialMapping)
                        vals = _push_coupling!(rows, cols, vals, row, col, coef)
                    end
                end
            end
        end
    end

    table = CouplingTable(rows, cols, vals)
    nrow = NtypeofExpr * length(νWhole) * nTestOrderBlocks
    ncol = NtypeofFields * activeTimePoints * length(νWhole)
    operator = _finalize_numerical_operator(table, nrow, ncol, side, backend, representation, numericalGeometry)
    costFunctions = operator isa SparseNumericalOperator ? operator.matrix : operator

    fieldArrays = compatibility_outputs ?
        _make_field_arrays(fieldSyms, NtypeofFields, activeTimePoints, νWhole, wholeRegionPointsSpace) :
        nothing

    return (
        costFunctions = costFunctions,
        operator = operator,
        couplingTable = table,
        fields = fieldArrays,
        場 = fieldArrays,
        champsLimité = nothing,
        geometry = numericalGeometry,
        coefficientRecipes = coefficientRecipes,
    )
end


function numericalOperaotrConstruction_slow_due_to_substitute_iteration(optRec,modelFam,side;absorbingBoundaries=nothing,maskedRegionInSpace=nothing,
    backend=:cpu,        # :cpu, :cuda, :mpi
    representation=:auto, # :matrixfree, :sparse, :blocksparse
    coefficient_type=:auto # :auto, :float64, :complexf64
)

    numericalGeometry = prepareNumericalOperatorGeometry(optRec,modelFam,side; absorbingBoundaries, maskedRegionInSpace)

    @unpack νWhole, νGeometry, νRelative, localPointsIndices,
        ModelPoints, Models, maskingField, conv,
        timePointsUsedForOneStep, activeTimePoints, wholeMin, wholeMax,
        modelMin, modelDomainMax, wholeRegionPointsSpace = numericalGeometry

    @unpack lhs, rhs, numbersOfTheSystem, fieldNames = optRec["recette"]
    @unpack fields, extfields = fieldNames

    if side == "left"
        numbersSide = numbersOfTheSystem.numbersOfTheSystemL
        symA = lhs.Ajiννᶜ
        varM = lhs.varM
        fieldSyms = fields
    elseif side == "right"
        numbersSide = numbersOfTheSystem.numbersOfTheSystemR
        symA = rhs.Γjiννᶜ
        varM = rhs.varF
        fieldSyms = extfields
    else
        throw(ArgumentError("side must be \"left\" or \"right\", got $side"))
    end

    @unpack NtypeofExpr, NtypeofFields = numbersSide
    NtypeofCoefficientVariables = size(varM, 1)

    fieldLinearIndices = LinearIndices(
        (NtypeofFields, activeTimePoints, Tuple(wholeRegionPointsSpace)...)
    )
    residualLinearIndices = LinearIndices((NtypeofExpr, length(νWhole)))

    rows = Int32[]
    cols = Int32[]
    vals = _coefficient_vector(coefficient_type)

    for iTestFunctions in eachindex(νWhole)
        νtmpWhole = νWhole[iTestFunctions]
        iGeometry = νGeometry[iTestFunctions]
        localPointsHere = localPointsIndices[iGeometry]
        middlepointHere = νRelative[iTestFunctions]
        localPointsSpaceHere = carDropDim.(localPointsHere)
        localPointsSpaceIndicesHere = CartesianIndices(Tuple(localPointsSpaceHere[end]))
        middlepointSpaceHere = svec2car(carDropDim(middlepointHere))
        iTimeMax = timePointsUsedForOneStep[iGeometry]

        νᶜtmpWhole = localPointsSpaceIndicesHere .+ (νtmpWhole .- middlepointSpaceHere)
        νᶜtmpModel = conv.whole2model.(νᶜtmpWhole)
        localLinearIndices = LinearIndices(Tuple(localPointsHere[end]))

        materialMapping = Dict()
        for iT in 1:iTimeMax
            for iVar in 1:NtypeofMaterialVariables
                spaceModelBouncedPoints = ModelPoints[1:end-1, iVar]
                iiT = ModelPoints[end, iVar] > 1 ? iT : 1
                νᶜtmpModelTruncated = BouncingCoordinates.(νᶜtmpModel, Ref(spaceModelBouncedPoints))

                for jPoint in νᶜtmpWhole
                    jPointLocal = jPoint - νtmpWhole + middlepointSpaceHere
                    jPointTLocal = carAddDim(jPointLocal, iT)
                    linearjPointTLocal = localLinearIndices[jPointTLocal]
                    materialMapping[varM[iVar, linearjPointTLocal]] =
                        Models[iVar][carAddDim(νᶜtmpModelTruncated[jPointLocal], iiT)]
                end
            end
        end

        for iExpr in 1:NtypeofExpr
            row = residualLinearIndices[iExpr, iTestFunctions]

            for iT in 1:iTimeMax
                for jPoint in νᶜtmpWhole
                    jPointLocal = jPoint - νtmpWhole + middlepointSpaceHere
                    jPointTLocal = carAddDim(jPointLocal, iT)
                    linearjPointTLocal = localLinearIndices[jPointTLocal]

                    jPointInWhole = is_all_less_than_or_equal(wholeMin, jPoint) &&
                                    is_all_less_than_or_equal(jPoint, wholeMax)
                    boundaryWeight = 0.0

                    if jPointInWhole
                        jPointModel = conv.whole2model(jPoint)
                        jPointInModel =
                            is_all_less_than_or_equal(modelMin, jPointModel) &&
                            is_all_less_than_or_equal(jPointModel, modelDomainMax)

                        if jPointInModel || iT == iTimeMax
                            boundaryWeight = 1.0
                        else
                            distance2 = distance2_point_to_box(
                                jPointModel,
                                modelMin,
                                modelDomainMax,
                            )
                            boundaryWeight = CerjanBoundaryCondition(distance2)
                        end
                    end

                    weight = jPointInWhole ? maskingField[jPoint] * boundaryWeight : 0.0
                    iszero(weight) && continue

                    for iField in 1:NtypeofFields
                        col = fieldLinearIndices[iField, iT, Tuple(jPoint)...]
                        coef = weight * substitute(symA[linearjPointTLocal, iField, iExpr, iGeometry], materialMapping)
                        vals = _push_coupling!(rows, cols, vals, row, col, coef)
                    end
                end
            end
        end
    end

    table = CouplingTable(rows, cols, vals)
    nrow = NtypeofExpr * length(νWhole)
    ncol = NtypeofFields * activeTimePoints * length(νWhole)
    operator = _finalize_numerical_operator(table, nrow, ncol, side, backend, representation, numericalGeometry)

    fieldArrays = Array{Any,2}(undef, NtypeofFields, activeTimePoints)
    for iT in 1:activeTimePoints
        for iField in 1:NtypeofFields
            fieldArrays[iField, iT] = reshape(
                [_field_symbol(fieldSyms, iField) for _ in 1:length(νWhole)],
                Tuple(wholeRegionPointsSpace),
            )
        end
    end

    costFunctions = operator isa SparseNumericalOperator ? operator.matrix : operator
    return (
        costFunctions = costFunctions,
        operator = operator,
        couplingTable = table,
        場 = fieldArrays,
        champsLimité = nothing,
        geometry = numericalGeometry,
    )
end


function numericalOperatorConstruction_too_heavy(optRec,modelFam,side;absorbingBoundaries=nothing,maskedRegionInSpace=nothing)

    #region general introduction
    #
    # 01/05/2026 Nobuaki Fuji
    #
    # after the construction of local (semi-)symbolic expressions with linearised operators#
    # here we will read the model parameters and construct the numerical operators
    #
    # Nobuaki Fuji @ IPGP/UPC/IUF since 2024
    #
    # 
    # encouraged by Thibault Duretz @ U. Frankfurt Goethe, Kurama Okubo @ NIED
    #
    # Julia hackathon October 2024, March 2025
    #
    #

    # coordinates: Model: the real model domain; Whole: computation domain with absorbing boundaries; 
    #              Empty: Whole + some more points to avoid missing reference to the field and material (they should be just zeros)
    #

    # intermediate presentations: IPGP-CIA workshop October 2024; IPGP-ERI workshop November 2024; lighttalk @ systemI December 2024
    #             EGU @ Vienna May 2025
    #    
    #     Fuji & Duretz in preparation
    #
    #
    #
    #
    #endregion

    #region

    numericalGeometry = prepareNumericalOperatorGeometry(optRec,modelFam,side; absorbingBoundaries, maskedRegionInSpace)

    @unpack νWhole, νGeometry, νRelative, localPointsIndices, middlepoints,
        ModelPoints, Models, maskingField, conv,
        timePointsUsedForOneStep, wholeMin, wholeMax, modelMin, modelDomainMax,
        wholeRegionPoints, wholeRegionPointsSpace = numericalGeometry

    @unpack models, modelName, modelPoints = modelFam

    @unpack lhs, rhs, nodes, centresIndices, numbersOfTheSystem,fieldNames=optRec["recette"]

    @unpack fields, extfields = fieldNames
    
    numbersL = numbersOfTheSystem.numbersOfTheSystemL
    numbersSide = side == "left" ? numbersOfTheSystem.numbersOfTheSystemL : numbersOfTheSystem.numbersOfTheSystemR
    @unpack timeMarching,nCoordinates,NtypeofExpr,NtypeofMaterialVariables = numbersL
    NtypeofFields = numbersSide.NtypeofFields
    
    nConfigurations=numbersOfTheSystem.nConfigurations

    if side == "left"
        symA = lhs.Ajiννᶜ
        varM = lhs.varM
        Ulocal = lhs.Ulocal
        fieldSyms = fields
    elseif side == "right"
        symA = rhs.Γjiννᶜ
        varM = rhs.varF
        Ulocal = rhs.Flocal
        fieldSyms = extfields
    else
        throw(ArgumentError("side must be \"left\" or \"right\", got $side"))
    end




    #region
    # one operator per test point in space, for now
    NtestfunctionsInSpace = length(νWhole)
    costFunctions = Array{Any,2}(undef, NtypeofExpr, NtestfunctionsInSpace)
    #costFunctionsCoefs = Array{Number,6}(undef, NtestfunctionsInSpace, NtypeofFields, NtypeofExpr, iTimeMax, NtestfunctionsInSpace)
    costFunctionsCoefs .= 0
    Threads.@threads for iTestFunctions in eachindex(νWhole)

        localCost = Array{Any,1}(undef,NtypeofExpr)
        localCost .= 0.0
        iPoint = iTestFunctions
        νtmpWhole = νWhole[iPoint]

        # preferred geometry for this test point
        iGeometry = νGeometry[iPoint]
        localPointsHere = localPointsIndices[iGeometry]
        middlepointHere = νRelative[iPoint]              # for now = middlepoints[iGeometry][1]
        localPointsSpaceHere = carDropDim.(localPointsHere)
        localPointsSpaceIndicesHere = CartesianIndices(Tuple(localPointsSpaceHere[end]))
        middlepointSpaceHere = svec2car(carDropDim(middlepointHere))
        
        # if time is appended, the active time depth depends on the geometry
        iTimeMax = timePointsUsedForOneStep[iGeometry]

        # shift local stencil to the current whole-space test point
     
        νᶜtmpWhole = localPointsSpaceIndicesHere .+ (νtmpWhole .- svec2car(carDropDim(middlepointHere)))
        νᶜtmpModel = conv.whole2model.(νᶜtmpWhole)

        # linear indexing on the chosen local stencil
        localLinearIndices = LinearIndices(Tuple(localPointsHere[end]))


        for iExpr in 1:NtypeofExpr
            tmpMapping = Dict()

            for iT in 1:iTimeMax
                # material variables
                for iVar in 1:NtypeofMaterialVariables
                    spaceModelBouncedPoints = ModelPoints[1:end-1, iVar]
                    iiT = ModelPoints[end, iVar] > 1 ? iT : 1

                    
                    νᶜtmpModelTruncated = BouncingCoordinates.(νᶜtmpModel, Ref(spaceModelBouncedPoints))

                    for jPoint in νᶜtmpWhole
                        jPointLocal = jPoint - νtmpWhole + svec2car(carDropDim(middlepointHere))
                        jPointTLocal = carAddDim(jPointLocal, iT)
                        linearjPointTLocal = localLinearIndices[jPointTLocal]

                        tmpMapping[varM[iVar, linearjPointTLocal]] =
                            Models[iVar][carAddDim(νᶜtmpModelTruncated[jPointLocal], iiT)]
                    end


                end

                # fields
                for jPoint in νᶜtmpWhole
                    jPointLocal = jPoint - νtmpWhole + middlepointSpaceHere
                    jPointTLocal = carAddDim(jPointLocal, iT)
                    linearjPointTLocal = localLinearIndices[jPointTLocal]

                    jPointInWhole = is_all_less_than_or_equal(wholeMin, jPoint) &&
                                    is_all_less_than_or_equal(jPoint, wholeMax)

                    if !jPointInWhole
                        for iField in 1:NtypeofFields
                            tmpMapping[Ulocal[linearjPointTLocal, iField]] = 0.0
                        end
                        continue
                    end

                    jPointModel = conv.whole2model(jPoint)

                    for iField in 1:NtypeofFields
                        if is_all_less_than_or_equal(modelMin, jPointModel) &&
                        is_all_less_than_or_equal(jPointModel, modelDomainMax)

                            #tmpMapping[Ulocal[linearjPointTLocal, iField]] =
                            #    場[iField, iT][jPoint] * maskingField[jPoint]
                            localCost[iExpr] += 
                                場[iField, iT][jPoint] * maskingField[jPoint] *
                                substitute(symA[linearjPointTLocal,iField,iExpr,iGeometry],tmpMapping)

                            #costFunctionsCoefs[jPoint] += NtypeofFields,  NtypeofExpr,iTimeMax, NtestfunctionsInSpace)
                        else
                            if iT == iTimeMax
                                #tmpMapping[Ulocal[linearjPointTLocal, iField]] =
                                #    場[iField, iT][jPoint] * maskingField[jPoint]
                                localCost[iExpr] += 
                                    場[iField, iT][jPoint] * maskingField[jPoint] *
                                    substitute(symA[linearjPointTLocal,iField,iExpr,iGeometry],tmpMapping)
                            else
                                distance2 = distance2_point_to_box(
                                    jPointModel,
                                    modelMin,
                                    modelDomainMax,
                                )
                                #tmpMapping[Ulocal[linearjPointTLocal, iField]] =
                                #    場[iField, iT][jPoint] *
                                #    maskingField[jPoint] *
                                #    CerjanBoundaryCondition(distance2)
                                localCost[iExpr] += 
                                    場[iField, iT][jPoint] * maskingField[jPoint] *
                                    CerjanBoundaryCondition(distance2) *
                                    substitute(symA[linearjPointTLocal,iField,iExpr,iGeometry],tmpMapping)
                            end
                        end
                    end
                end

            end

            # centre-only address: always use the geometry centre now
            #costFunctions[iExpr, iTestFunctions] +=
            #    substitute(Av[iExpr,iGeometry], tmpMapping)
        end
        costFunctions[:, iTestFunctions] .= localCost
    end



    #endregion

    return (costFunctions=costFunctions,場=場,champsLimité=champsLimité)

end

function prepareNumericalOperatorGeometry(
    optRec,
    modelFam,
    side;
    absorbingBoundaries,
    maskedRegionInSpace,
    freeSurfaceBoundary=nothing,
)


    @unpack models, modelName, modelPoints = modelFam


    @unpack lhs, rhs, nodes, centresIndices, numbersOfTheSystem,fieldNames=optRec["recette"]

    @unpack fields, extfields = fieldNames

    
    if side=="left"
        symA=lhs.Ajiννᶜ
        varM=lhs.varM
        Ulocal=lhs.Ulocal
        fields=fields
    elseif side=="right"
        fields = nothing
        symA=rhs.Γjiννᶜ
        varM=rhs.varF
        Ulocal=rhs.Flocal
        fields=extfields
    end

    

    numbersL = numbersOfTheSystem.numbersOfTheSystemL
    numbersSide = side == "left" ? numbersOfTheSystem.numbersOfTheSystemL : numbersOfTheSystem.numbersOfTheSystemR
    @unpack timeMarching,nCoordinates,NtypeofExpr,NtypeofMaterialVariables = numbersL
    NtypeofFields = numbersSide.NtypeofFields
    
    nConfigurations=numbersOfTheSystem.nConfigurations

    # normally the geometry configurations should be proposed in the preferred order
    nGeometry=nConfigurations
    Ndimension=nCoordinates

    

    # the last coordinate should be cosidered as time

    localPointsIndices = Vector{Any}(undef, nGeometry)
    middlepoints = Vector{Any}(undef, nGeometry)
    lhs_CartesianDependencies = lhs.CartesianDependencies
    rhs_CartesianDependencies = rhs.CartesianDependencies

    newD=Ndimension

    if !timeMarching
        newD = Ndimension + 1
        modelPoints = (modelPoints..., 1)

        lhs_CartesianDependencies = vcat(lhs.CartesianDependencies,zeros(Int, 1, NtypeofMaterialVariables),)
        rhs_CartesianDependencies = vcat(rhs.CartesianDependencies,zeros(Int, 1, NtypeofFields),)

        for iGeometry in 1:nGeometry
            localPointsIndices[iGeometry] = [SVector{newD,Int}(p..., 1) for p in nodes[iGeometry]]
            selected = nodes[iGeometry][centresIndices[iGeometry]]
            middlepoints[iGeometry] = SVector{newD,Int}(selected..., 1)
        end
    else
        
        for iGeometry in 1:nGeometry
            localPointsIndices[iGeometry] = nodes[iGeometry]
            selected = nodes[iGeometry][centresIndices[iGeometry]]
            middlepoints[iGeometry] = selected
        end
    end



    if length(models) !== NtypeofMaterialVariables 
        @error "Each material has to have its own model"
    end
    
    Models=Array{Any,1}(undef,NtypeofMaterialVariables)
    ModelPoints = Array{Int,2}(undef, newD, NtypeofMaterialVariables)
    
    for iVar ∈ 1:NtypeofMaterialVariables
        CartesianDependency=lhs_CartesianDependencies[:,iVar]
        if ndims(models[iVar]) !== sum(CartesianDependency)
            @error "model parameter dimension is not what you declared in the equation!"
        end
        if sum(CartesianDependency) == 0 # when it is a constant
            tmpModel = Array{Any,newD}(undef, (ones(Int, newD)...)...)
            ModelPoints[:, iVar] = ones(Int, newD)
            tmpModel[vec2car(ones(Int, newD))] = models[iVar]
            Models[iVar]=tmpModel
        else
            size(models[iVar])
            newCoords=expandVectors(size(models[iVar]),CartesianDependency)
            ModelPoints[:,iVar] = newCoords
            tmpModel=reshape(models[iVar],newCoords...)
            Models[iVar]=tmpModel

            for iCoord in eachindex(newCoords)
                if newCoords[iCoord]!== modelPoints[iCoord] && newCoords[iCoord] !== 1
                    @error "the model should have the same dimension! (or constant)"
                end
            end
        end
    end
    
    #endregion
    
    
    #region construction of the fields

    wholeRegionPoints = nothing

    if absorbingBoundaries === nothing
        wholeRegionPoints=modelPoints
        absorbingBoundaries = zeros(Int,2, newD)
    elseif absorbingBoundaries === "CerjanBoundary"
        throw(ArgumentError(
            "the string \"CerjanBoundary\" is ambiguous; use " *
            "CerjanBoundarySpec(lower, upper) to choose every side explicitly",
        ))
    else
        # Rows are lower/upper sides; columns are coordinates.
        if size(absorbingBoundaries)[1] !== 2
            throw(DimensionMismatch(
                "absorbingBoundaries must have two rows (lower and upper)",
            ))
        elseif size(absorbingBoundaries)[2] !== length(modelPoints) && !timeMarching
            throw(DimensionMismatch(
                "absorbingBoundaries needs one column per coordinate",
            ))
        elseif size(absorbingBoundaries)[2] === length(modelPoints)-1 && timeMarching
            absorbingBoundaries = hcat(
                absorbingBoundaries,
                zeros(Int, 2, 1),
            )
        elseif size(absorbingBoundaries)[2] !== length(modelPoints)
            throw(DimensionMismatch(
                "absorbingBoundaries needs one column per spatial coordinate",
            ))
        end
        wholeRegionPoints=modelPoints.+ sum(absorbingBoundaries, dims=1)[:]
    end
    wholeRegionPointsSpace=wholeRegionPoints[1:end-1]

    #endregion


    #region 

    spacePointsUsed=Vector{Any}(undef, nGeometry)
    timePointsUsedForOneStep=Vector{Any}(undef, nGeometry)
    for iGeometry ∈ 1:nGeometry
        localPointVecs = localPointsIndices[iGeometry]
        spacePointsUsed[iGeometry] = localPointVecs[end][1:end-1]
        timePointsUsedForOneStep[iGeometry] = localPointVecs[end][end]
    end

    # Every moving-ν geometry has the same time depth.
    activeTimePoints = timePointsUsedForOneStep[1]

    # Prefer the historical middle ν in the interior, then move ν just far
    # enough inward for the complete local stencil to remain in the domain.
    # The geometry index is a single extra lookup dimension, so this remains
    # compatible with the existing CPU/GPU contraction and with 2D--4D grids.
    geometryPreference = fill(1, Tuple(wholeRegionPointsSpace))
    if nGeometry > 1
        geometryByCentre = Dict(
            Tuple(middlepoints[iGeometry][1:end-1]) => iGeometry
            for iGeometry in 1:nGeometry
        )
        stencilSize = collect(localPointsIndices[1][end][1:end-1])
        preferredCentre = ((stencilSize .- 1) .÷ 2) .+ 1
        if timeMarching && localPointsIndices[1][end][end] > 1
            # Only spatial ν moves; candidate generation fixes the time level.
            preferredCentre = collect(preferredCentre)
        end

        for point in CartesianIndices(Tuple(wholeRegionPointsSpace))
            globalPoint = collect(Tuple(point))
            lower = max.(1, globalPoint .+ stencilSize .- collect(wholeRegionPointsSpace))
            upper = min.(stencilSize, globalPoint)
            selectedCentre = clamp.(preferredCentre, lower, upper)
            geometryPreference[point] = geometryByCentre[Tuple(selectedCentre)]
        end
    end

   


    #since everything is super clumsy, here we make several useful functions to change one coordinate to another
    
    conv=spaceCoordinatesConversionfunctions(absorbingBoundaries[:,1:end-1], newD-1)
    freeSurfaceBoundaryWhole = if isnothing(freeSurfaceBoundary)
        nothing
    else
        expectedDimension = newD - 1
        all(length(point) == expectedDimension for point in freeSurfaceBoundary.points) ||
            throw(DimensionMismatch(
                "free-surface point dimension differs from the PDE",
            ))
        FreeSurfacePointSet(
            conv.model2whole.(freeSurfaceBoundary.points),
            copy(freeSurfaceBoundary.normals),
        )
    end
    #endregion 

    #region

    # Useful point lists
    PointsSpace = CartesianIndices(Tuple(wholeRegionPointsSpace))
    NpointsSpace = length(PointsSpace)

    # For now, test functions are still identified with spatial points
    νWhole = collect(PointsSpace)
    wholeMin = νWhole[1]
    wholeMax = νWhole[end]
    modelMin = CartesianIndex(ones(Int, newD-1)...)
    modelDomainMax = vec2car(collect(modelPoints[1:end-1]))



    # Pointwise preferred geometry and corresponding relative centre
    νGeometry = Vector{Int}(undef, NpointsSpace)
    νRelative = Vector{Any}(undef, NpointsSpace)

    for iPoint in eachindex(νWhole)
        νtmpWhole = νWhole[iPoint]
        iGeometry = geometryPreference[νtmpWhole]
        νGeometry[iPoint] = iGeometry

        νRelative[iPoint] = middlepoints[iGeometry]
    end

 

    #endregion


    #region making a maskingField (for limited source areas, boundary conditions, etc.)

    maskingField=Array{Any,newD-1}(undef,Tuple(wholeRegionPointsSpace)) # maskingField is defined only for whole domain
    champsLimité = nothing
    if maskedRegionInSpace === nothing
        maskingField .= 1.0
    elseif typeof(maskedRegionInSpace) === Array{CartesianIndex,1}
        champsLimité = Array{Any,2}(undef,NtypeofFields,timePointsUsedForOneStep)
        for it in 1:timePointsUsedForOneStep
            for iField in eachindex(fields)
                newstring=split(string(fields[iField]),"(")[1]*"_mod_limited"*"_t="*string(it)
                champsLimité[iField,it] = Array{Any,1}(undef,length(maskedRegionInSpace))
            end
        end
        maskingField .= 0.0
        tmpIndex=1
        for iSpace in maskedRegionInSpace
            jSpace = conv.model2whole(iSpace)
            maskingField[jSpace] =1.0
            for it in 1:timePointsUsedForOneStep
                for iField in eachindex(fields)
                    
                    #tmpChampsLimitéContents= (jSpace,場[iField,it][jSpace])
                    champsLimité[iField,it][tmpIndex]=場[iField,it][jSpace]
                end
            end
            tmpIndex += 1
        end
    else
        @error "maskedRegionInSpace should be a 1D array of CartesianIndex (if it is CartesianIndices, you need to collect(Tuple()))"
    end

    #endregion

    fieldLinearIndices = LinearIndices(
        (NtypeofFields, activeTimePoints, Tuple(wholeRegionPointsSpace)...)
    )

    residualLinearIndices = LinearIndices(
        (NtypeofExpr, length(νWhole))
    )



    return (modelPoints=modelPoints,
    wholeRegionPoints=wholeRegionPoints,
    absorbingBoundaries=absorbingBoundaries,
    νWhole=νWhole,
    νGeometry=νGeometry,
    νRelative=νRelative,
    localPointsIndices=localPointsIndices,
    middlepoints=middlepoints,
    ModelPoints=ModelPoints,
    geometryPreference=geometryPreference,
    maskingField=maskingField,
    conv=conv,
    Models = Models,
    timePointsUsedForOneStep = timePointsUsedForOneStep,
    activeTimePoints = activeTimePoints,
    wholeRegionPointsSpace = wholeRegionPointsSpace,
    wholeMin = wholeMin,
    wholeMax = wholeMax,
    modelMin = modelMin,
    modelDomainMax = modelDomainMax,
    freeSurfaceBoundary = freeSurfaceBoundaryWhole,
    )
end





carDropDim(v::SVector{N,T}) where {N,T} = SVector{N-1,T}(ntuple(i -> v[i], N-1))

function spaceCoordinatesConversionfunctions(absorbingBoundaries, NdimensionMinusTime)
    offset_model = vec2car(absorbingBoundaries[1, 1:NdimensionMinusTime])
    #offset_empty = vec2car(spacePointsUsed)

    model2whole(a::CartesianIndex) = a + offset_model
    whole2model(a::CartesianIndex) = a - offset_model
    #whole2empty(a::CartesianIndex) = a + offset_empty
    #empty2whole(a::CartesianIndex) = a - offset_empty
    #model2empty(a::CartesianIndex) = whole2empty(model2whole(a))
    #empty2model(a::CartesianIndex) = whole2model(empty2whole(a))
    return(; model2whole, whole2model)
    #return (; model2whole, whole2model, whole2empty, empty2whole, model2empty, empty2model)
end


function BouncingCoordinates(a::CartesianIndex,PointsUsed)
    #
    # this will bounce the boundary inside the PointsUsed vector
    #
    # i.e. get the nearby coordinates inside the domain to fake the continuity

    if length(a) !== length(PointsUsed)
        @error "cannot bound this CartesianIndex due to the dimension mismatch"
    end
    avector=car2vec(a)
    for iCoord in eachindex(avector)
        if avector[iCoord] < 1
            avector[iCoord] = 1
        elseif avector[iCoord] > PointsUsed[iCoord]
            avector[iCoord] = PointsUsed[iCoord]
        end
    end
    a=vec2car(avector)
    return a
end
