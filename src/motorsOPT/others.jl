
timeDimensionString="t" 

_as_variable_tuple(vars) = vars === nothing ? () : vars

"A volume contribution paired with a derivative of the OPT test function."
struct WeakTerm{E,N}
    expression::E
    test_derivative_orders::NTuple{N,Int}
end

"A physical flux on the right-hand side of the variational boundary integral."
struct BoundaryFlux{E,N}
    expression::E
    normal_coordinate::Int
    test_derivative_orders::NTuple{N,Int}
end

_weak_operation(expression) = try
    Symbolics.operation(Symbolics.value(expression))
catch
    nothing
end

_weak_arguments(expression) = try
    Symbolics.arguments(Symbolics.value(expression))
catch
    ()
end

function _weak_additive_terms(expression)
    _weak_operation(expression) === (+) || return Any[expression]
    return reduce(vcat, (_weak_additive_terms(argument) for argument in
                         _weak_arguments(expression)); init=Any[])
end

function _numeric_factor_and_outer_derivative(term)
    operation = _weak_operation(term)
    if operation isa Differential
        arguments = _weak_arguments(term)
        length(arguments) == 1 || return nothing
        return (coefficient=1, differential=operation, flux=arguments[1])
    elseif operation === (*)
        arguments = collect(_weak_arguments(term))
        derivative_positions = findall(argument ->
            _weak_operation(argument) isa Differential, arguments)
        length(derivative_positions) == 1 || return nothing
        position = only(derivative_positions)
        coefficient_factors = deleteat!(copy(arguments), position)
        all(factor -> Symbolics.value(factor) isa Number,
            coefficient_factors) || return nothing
        derivative_term = arguments[position]
        derivative_arguments = _weak_arguments(derivative_term)
        length(derivative_arguments) == 1 || return nothing
        coefficient = prod(Symbolics.value(factor) for factor in
                           coefficient_factors; init=1)
        return (
            coefficient=coefficient,
            differential=_weak_operation(derivative_term),
            flux=derivative_arguments[1],
        )
    end
    return nothing
end

function _coordinate_index(variable, coordinates)
    target = Symbolics.value(variable)
    return findfirst(coordinate ->
        isequal(Symbolics.value(coordinate), target), coordinates)
end

"""
    naturalWeakForm(equationCharacteristics)

Convert only explicit *outer spatial divergences* in the symbolic strong
residual into volume terms containing derivatives of the test function.  Time
derivatives remain on the unknown field.  Boundary fluxes are reported with
the sign they have on the right-hand side of the variational equation.

The transformation is deliberately conservative: an expression such as
`a * Dx(q)` is transformed only when `a` is a numerical constant.  A
spatially varying factor outside `Dx` is not silently interpreted as part of
the flux.
"""
function naturalWeakForm(equationCharacteristics)
    (; exprs, coordinates) = equationCharacteristics
    expressions = exprs isa Tuple ? exprs : (exprs,)
    coordinate_tuple = Tuple(coordinates)
    ncoordinates = length(coordinate_tuple)
    time_index = findfirst(coordinate ->
        string(coordinate) == timeDimensionString, coordinate_tuple)
    zero_orders = ntuple(_ -> 0, ncoordinates)

    volume_terms = map(expressions) do expression
        terms = WeakTerm[]
        for additive_term in _weak_additive_terms(expression)
            outer = _numeric_factor_and_outer_derivative(additive_term)
            coordinate_index = isnothing(outer) ? nothing :
                _coordinate_index(outer.differential.x, coordinate_tuple)
            if isnothing(outer) || isnothing(coordinate_index) ||
               coordinate_index == time_index
                push!(terms, WeakTerm(additive_term, zero_orders))
            else
                derivative_orders = ntuple(d -> d == coordinate_index ? 1 : 0,
                                            ncoordinates)
                # ∫ W c ∂j(q) = -∫ (∂jW) c q + ∮ W c q n_j.
                push!(terms, WeakTerm(
                    -outer.coefficient * outer.flux,
                    derivative_orders,
                ))
            end
        end
        terms
    end

    boundary_fluxes = map(expressions) do expression
        fluxes = BoundaryFlux[]
        for additive_term in _weak_additive_terms(expression)
            outer = _numeric_factor_and_outer_derivative(additive_term)
            isnothing(outer) && continue
            coordinate_index = _coordinate_index(
                outer.differential.x, coordinate_tuple)
            (isnothing(coordinate_index) || coordinate_index == time_index) &&
                continue
            # Move the integration-by-parts boundary contribution to the RHS.
            push!(fluxes, BoundaryFlux(
                -outer.coefficient * outer.flux,
                coordinate_index,
                zero_orders,
            ))
        end
        fluxes
    end

    return (
        volume_terms=volume_terms,
        boundary_fluxes=boundary_fluxes,
        coordinates=coordinate_tuple,
        time_coordinate=time_index,
    )
end

function weakTermGroups(weak_form)
    n_equations = length(weak_form.volume_terms)
    groups = Dict{Any,Vector{Any}}()
    for (equation, terms) in enumerate(weak_form.volume_terms), term in terms
        expressions = get!(groups, term.test_derivative_orders) do
            Any[0 for _ in 1:n_equations]
        end
        expressions[equation] += term.expression
    end
    return groups
end
# if the user does not want to use "t" for the time marching scheme, it should be changed
# and this "t" should be the last element in coordinates

function findCartesianDependency(expression,Ndimension,∂)
    expressionDependency=ones(Int,Ndimension)
    for iDimension in 1:Ndimension
        if typeof(expand_derivatives(∂[iDimension](expression))==0) <: Bool 
            if expand_derivatives(∂[iDimension](expression))==0
                expressionDependency[iDimension] = 0
            end
        end
    end
    return expressionDependency
end

function makeMixPartials(orders,coordinates;field=identity)
    # this function will give a matrix of mixed partial derivative operators 
    # coordinates should be an array of Symbolics variables 
    # orders is a matrix 


    Ndimension=length(coordinates)
    if length(orders)!== Ndimension
        @error "the highest orders array has not the same dimension as that of the coordinates"
    end

    ∂ = []
    for iDim in 1:Ndimension
        ∂ = push!(∂,Differential(coordinates[iDim]))
    end

    ∇ = Array{Any,Ndimension}(undef, Tuple(orders))
    R=CartesianIndices(∇)
        
    ∇ .= field
    for I in R
        for iDim in 1:Ndimension
            ∇[I] = (∂[iDim]^(I[iDim]-1))(∇[I])
        end
    end
   
    return ∇
end

function _unique_derivative_orders(indices::Vector{Int})
    isempty(indices) && return [Int[]]
    result = Vector{Vector{Int}}()
    for coordinate in unique(indices)
        position = findfirst(==(coordinate), indices)
        remaining = copy(indices)
        deleteat!(remaining, position)
        for suffix in _unique_derivative_orders(remaining)
            push!(result, [coordinate; suffix])
        end
    end
    return result
end

const _commuting_terms_cache = Dict{Any,Any}()
const _commuting_coefficient_cache = Dict{Any,Any}()
const _commuting_cache_lock = ReentrantLock()

_cartesian_orders(index) = Tuple(index[d] - 1 for d in 1:length(index))

function _symbolic_operation(expression)
    value = Symbolics.value(expression)
    return try
        Symbolics.operation(value)
    catch
        nothing
    end
end

function _symbolic_arguments(expression)
    value = Symbolics.value(expression)
    return try
        Symbolics.arguments(value)
    catch
        ()
    end
end

"Maximum consecutive symbolic Differential depth applied to `field`."
function _maximum_derivative_depth(expression, field, depth=0)
    operation = _symbolic_operation(expression)
    field_operation = _symbolic_operation(field)
    operation == field_operation && return depth
    arguments = _symbolic_arguments(expression)
    isempty(arguments) && return -1
    next_depth = operation isa Differential ? depth + 1 : 0
    return maximum((_maximum_derivative_depth(argument, field, next_depth)
                    for argument in arguments); init=-1)
end

"""
Return every symbolic nesting that represents one commuting mixed partial.

Symbolics distinguishes `Dx(Dy(u))` from `Dy(Dx(u))` structurally, although
the smooth fields used by OPT satisfy Clairaut's theorem.  Searching only the
nesting produced by `makeMixPartials` silently loses PDE terms written in the
opposite order.  This was visible in 2-D elasticity as a missing λ or μ part of
the `(λ + μ) ∂xy` coupling.
"""
function commutingMixedPartialTerms(index, coordinates; field=identity)
    key = (_cartesian_orders(index), Tuple(coordinates), field)
    cached = lock(_commuting_cache_lock) do
        get(_commuting_terms_cache, key, nothing)
    end
    cached === nothing || return cached
    derivative_indices = Int[]
    for dimension in eachindex(coordinates)
        append!(derivative_indices, fill(dimension, index[dimension] - 1))
    end
    partials = Differential.(coordinates)
    terms = unique(map(_unique_derivative_orders(derivative_indices)) do order
        term = field
        for dimension in order
            term = partials[dimension](term)
        end
        term
    end)
    lock(_commuting_cache_lock) do
        _commuting_terms_cache[key] = terms
    end
    return terms
end

function _commuting_coefficient(expression, index, coordinates; field=identity)
    key = (expression, _cartesian_orders(index), Tuple(coordinates), field)
    cached = lock(_commuting_cache_lock) do
        get(_commuting_coefficient_cache, key, nothing)
    end
    cached === nothing || return cached
    coefficient = sum(myCoeff(expression, term) for term in
                      commutingMixedPartialTerms(index, coordinates; field=field))
    lock(_commuting_cache_lock) do
        _commuting_coefficient_cache[key] = coefficient
    end
    return coefficient
end

function varMmaker(maxPointsUsed,coordinates,vars,∂) 
    # this will make an array of material coeffs for with a local Cartesian grid (max points used for a node)
    vars = _as_variable_tuple(vars)
    Ndimension = length(coordinates)

    R = CartesianIndices(Tuple(maxPointsUsed))

    varM=Array{Any,2}(undef,length(vars),length(R))
    CartesianDependencies=Array{Int,2}(undef,Ndimension,length(vars))

    for iVar in eachindex(vars)


        newstring=split(string(vars[iVar]),"(")[1]
     
        
        CartesianDependency=findCartesianDependency(vars[iVar],Ndimension,∂)
        
        smallVarM=Symbolics.variables(Symbol(newstring),1:length(R))
        for j in R
            linearJ=LinearIndices(R)[j]
            realJ=(car2svec(j).-1).*CartesianDependency .+1 # if there is no dependence on a direction, it should get the same name
            linearRealJ=LinearIndices(R)[CartesianIndex(realJ...)]
            smallVarM[linearJ]=smallVarM[linearRealJ]
        end
        varM[iVar,:]=smallVarM
        CartesianDependencies[:,iVar]=CartesianDependency
    end
    
    return varM,CartesianDependencies
end

function PDECoefFinder(orders,coordinates,expr,field,vars)
    # PDECoefFinder cannot detect the material partials × material partials for the moment!! 
    # I know how to do it, but eq. 40 should be then more generalised (kind of the product of partials of different materials)

    # maxPolynomialOrderMaterial is also a chelou thing, that I need to work on more systematically
    # like the powers of partials should also be included but here search for Rm[1], yeah, that's what I am doing


    Ndimension = length(coordinates)
    alpha=[]
    
    maxPolynomialOrderMaterial = 2*(maximum(orders)-1)
    R=CartesianIndices(Tuple(orders))
    expr=mySimplify(expr)
    maximumFieldDerivative = _maximum_derivative_depth(expr, field)


    for i in R
        sum(_cartesian_orders(i)) > maximumFieldDerivative && continue
        # Mixed derivatives commute for the smooth fields represented here.
        # Sum every symbolic nesting because Symbolics does not canonicalise
        # Dx(Dy(field)) and Dy(Dx(field)) to the same expression.
        tmpCoeff = _commuting_coefficient(expr, i, coordinates; field=field)
        if tmpCoeff !== 0
            isTmpCoeffAConstant=true
            for iVar in eachindex(vars)
                
                Rm=CartesianIndices(Tuple(orders))
                maximumMaterialDerivative =
                    _maximum_derivative_depth(tmpCoeff, vars[iVar])
                for j in Rm
                    sum(_cartesian_orders(j)) > maximumMaterialDerivative && continue
                    tmpCoeffMaterial = _commuting_coefficient(
                        tmpCoeff, j, coordinates; field=vars[iVar])
                    
                    if tmpCoeffMaterial !==0
                        isTmpCoeffAConstant=false
                        isOKtoinclude =true
                        # This is to avoid partials of other material 
                        for jVar in eachindex(vars)
                            if jVar !== iVar
                                ∇n = makeMixPartials(orders,coordinates;field=vars[jVar]) 
                                for jj in Rm
                                    if jj !== Rm[1]
                                        sum(_cartesian_orders(jj)) >
                                            _maximum_derivative_depth(
                                                tmpCoeffMaterial, vars[jVar]) && continue
                                        differentialCoeff = _commuting_coefficient(
                                            tmpCoeffMaterial, jj, coordinates;
                                            field=vars[jVar])
                                        if differentialCoeff !==0
                                            isOKtoinclude = false
                                        end
                                    end
                                end
                            end
                        end
                        if isOKtoinclude
                            #tmpCoeffMaterial=substitute(tmpCoeffMaterial,mapping)
                            specificMaterialTerm=tmpCoeffMaterial*vars[iVar]
                            tmpAlphaIJ = (node=specificMaterialTerm,nᶜ=j,n=i) # the famous n prime and n in the equation 56 or 40
                            alpha = push!(alpha,tmpAlphaIJ)
                        end
                    end
                end
                for matPower in 2:maxPolynomialOrderMaterial
                    tmpCoeffMaterial = myCoeff(tmpCoeff,vars[iVar]^matPower)
                    #tmpCoeffMaterial=substitute(tmpCoeffMaterial,mapping)
                    if tmpCoeffMaterial !==0
                        specificMaterialTerm=tmpCoeffMaterial*vars[iVar]^matPower
                        tmpAlphaIJ = (node=specificMaterialTerm,nᶜ=Rm[1],n=i) # the famous n prime and n in the equation 56 or 40
                        alpha = push!(alpha,tmpAlphaIJ)
                    end
                end
            end
            if isTmpCoeffAConstant
                specificMaterialTerm = tmpCoeff
                tmpAlphaIJ=(node=specificMaterialTerm,nᶜ=R[1],n=i)
                alpha = push!(alpha,tmpAlphaIJ)
            end
        end
    end
    alpha=unique(alpha)

    return alpha # varM: iVar and linearised cartesian indices
end 

function TaylorCoefInversion(coefInversionDict::Dict)

    # this function will give Float64 array (not Any)

    # the user might want to have a look on illposedTaylorCoefficientsInversion_legend, which is deprecated as of 10/06/2025.

    # as of 26/04/2026, it can compute for an arbitrary μPoints

    # based on the equation 27 (of the version 10/06/2025 FD2025 : \psi_{;\mu,\nu}^{(l)}[\nu+\mu]=\sum_\eta C_{\mu+\eta;\mu,\nu}^{(l)} \psi[\nu+\mu+\eta]), we need to perform this inversion anyways for all the point \mu inside L(\nu) (the concerned points for \nu)


    # be careful that pointsIndices is now a 1D array of integer vectors!!

    @unpack multiOrdersIndices, pointsIndices, μpointsIndices, Δ  = coefInversionDict
    taylorInverseMode = Symbol(get(coefInversionDict, "taylor_inverse_mode", :hierarchical_constrained))
    exactTaylorTotalDegree = get(coefInversionDict, "exact_taylor_total_degree", nothing)


    @show typeof(pointsIndices)
    
    @show multiOrdersIndices, pointsIndices, μpointsIndices 
    
    numberOfEtas = length(pointsIndices)
    numberOfLs   = length(multiOrdersIndices)
    numberOfMus = length(μpointsIndices) 

    CˡηGlobal = Array{Float64,3}(undef,numberOfEtas,numberOfLs,numberOfMus)

    # this is the C^{(l)}_{\mu+\eta; μ, \nu}

    for μ_oneD in axes(CˡηGlobal,3)
        #@show typeof(multiOrdersIndices),typeof(pointsIndices)
        CˡηGlobal[:,:,μ_oneD]=TaylorCoefInversion(numberOfLs,numberOfEtas,multiOrdersIndices,pointsIndices,Δ,μpointsIndices[μ_oneD]; taylor_inverse_mode=taylorInverseMode, exact_total_degree=exactTaylorTotalDegree)
    end 

    return @strdict(CˡηGlobal)

end

function TaylorCoefInversion(numberOfLs,numberOfEtas,multiOrdersIndices,pointsIndices,Δ,μPoint; taylor_inverse_mode=:hierarchical_constrained, exact_total_degree=nothing)


    # the old version is : illposedTaylorCoefficientsInversionSingleCentre

    # in fact, available points depend on the position of μ (=k here), we need to 'mute' some points
    # with Y_μ
    
    # for this pointsIndices are filtered for every μ

    tmpTaylorExpansionCoeffs = Array{Any,2}(undef,numberOfLs,numberOfEtas)

    for iAvailablePoint in eachindex(pointsIndices)
        η_μ = pointsIndices[iAvailablePoint]
        η = Float64.(η_μ) .- μPoint
        distances= η .* Δ
        for j in multiOrdersIndices
            linearJ = LinearIndices(multiOrdersIndices)[j]
            orders = car2svec(j).-1
            numerator = prod(distances .^orders)
            denominator=prod(factorial.(orders))
            tmpTaylorCoeffs = numerator/denominator
            tmpTaylorExpansionCoeffs[linearJ,iAvailablePoint]=tmpTaylorCoeffs 
        end
    end

    # Use a scaled SVD pseudo-inverse of the rectangular Taylor matrix directly.
    # This avoids forming normal equations A'A, which square the condition number
    # and can break x/y symmetry for high supplementaryOrder.
    A = Num2Float64.(tmpTaylorExpansionCoeffs)
    tmpCˡηlocal = if taylor_inverse_mode == :hierarchical_constrained
        hierarchicalTaylorPseudoInverse(
            A, multiOrdersIndices, pointsIndices;
            exact_total_degree=exact_total_degree,
        )
    elseif taylor_inverse_mode == :weak_operator_optimized
        hierarchicalTaylorPseudoInverse(
            A, multiOrdersIndices, pointsIndices;
            exact_total_degree=exact_total_degree,
            supplementary_objective=:weak_moments,
        )
    elseif taylor_inverse_mode == :scaled_svd
        stableTaylorPseudoInverse(A)
    elseif taylor_inverse_mode == :moore_penrose_svd
        directTaylorPseudoInverse(A)
    else
        throw(ArgumentError("unknown taylor inverse mode: $taylor_inverse_mode"))
    end

    return tmpCˡηlocal
end

"""
    hierarchicalTaylorPseudoInverse(A, multiOrdersIndices, pointsIndices)

Construct a Taylor inverse that reproduces a downward-closed low-order
polynomial space exactly. Higher Taylor coefficients are fitted only inside
the null space of those exact constraints, so supplementary orders cannot
alter constants, linear terms, pure second derivatives, or mixed second
derivatives.

For a three-point tensor stencil the exact space is the set of monomials of
total degree at most two. This includes `x*z`, while retaining null-space
degrees of freedom with which supplementary terms can improve the weak
operator. The former unconstrained Moore-Penrose inverse mixed aliased
monomials (`x^2` with `x^4`, and `x*z` with three higher-order aliases).
"""
function hierarchicalTaylorPseudoInverse(
    A,
    multiOrdersIndices,
    pointsIndices;
    rtol=sqrt(eps(Float64)),
    exact_total_degree=nothing,
    supplementary_objective=:uniform,
)
    A = Matrix{Float64}(A)
    order_indices = vec(collect(multiOrdersIndices))
    length(order_indices) == size(A, 1) ||
        throw(DimensionMismatch("Taylor order count does not match the rows of A"))

    points = vec(collect(pointsIndices))
    coordinate_dimension = length(first(points))
    coordinate_counts = [
        length(unique(point[d] for point in points))
        for d in 1:coordinate_dimension
    ]
    active_dimensions = findall(>(1), coordinate_counts)
    requested_degree = isnothing(exact_total_degree) ?
        (isempty(active_dimensions) ? 0 :
            minimum(coordinate_counts[d] - 1 for d in active_dimensions)) :
        Int(exact_total_degree)
    requested_degree >= 0 ||
        throw(ArgumentError("exact Taylor total degree must be nonnegative"))
    total_degree(index) = sum(index[d] - 1 for d in active_dimensions)

    # Add low-order monomials in degree order, retaining only independent
    # rows whose direct predecessors are already present.  A clipped 3×2
    # stencil can therefore preserve x² and xz while correctly rejecting z².
    candidates = sort(
        findall(index -> total_degree(index) <= requested_degree, order_indices);
        by=row -> (total_degree(order_indices[row]), Tuple(order_indices[row])),
    )
    exact_rows = Int[]
    selected_orders = Set{Tuple}()
    current_rank = 0
    for row in candidates
        exponent = Tuple(order_indices[row][d] - 1 for d in 1:coordinate_dimension)
        predecessors_present = all(1:coordinate_dimension) do d
            exponent[d] == 0 || begin
                predecessor = ntuple(k -> exponent[k] - (k == d), coordinate_dimension)
                predecessor in selected_orders
            end
        end
        predecessors_present || continue
        trial_rows = [exact_rows; row]
        trial_rank = rank(A[trial_rows, :]; rtol=rtol)
        trial_rank > current_rank || continue
        push!(exact_rows, row)
        push!(selected_orders, exponent)
        current_rank = trial_rank
    end
    supplementary_rows = setdiff(collect(axes(A, 1)), exact_rows)
    Aexact = A[exact_rows, :]

    rank_exact = rank(Aexact; rtol=rtol)
    rank_exact == length(exact_rows) ||
        throw(ArgumentError(
            "the exact Taylor constraint space is rank deficient: " *
            "rank=$rank_exact for $(length(exact_rows)) constraints",
        ))

    # C has dimensions (nodes, Taylor terms), so A*C is the Taylor-order
    # reproduction matrix. First impose exact identity rows.
    C = zeros(Float64, size(A, 2), size(A, 1))
    C[:, exact_rows] .= directTaylorPseudoInverse(Aexact; rtol=rtol)

    null_exact = nullspace(Aexact; rtol=rtol)
    if !isempty(supplementary_rows) && size(null_exact, 2) > 0
        Ahigh = A[supplementary_rows, :]
        target_high = zeros(Float64, length(supplementary_rows), size(A, 1))
        for (row, global_row) in enumerate(supplementary_rows)
            target_high[row, global_row] = 1.0
        end

        # Corrections lie in null(Aexact); exact low-order reproduction is
        # therefore invariant under the supplementary least-squares fit.
        reduced = Ahigh * null_exact
        residual = target_high - Ahigh * C
        if supplementary_objective === :uniform
            C .+= null_exact * directTaylorPseudoInverse(reduced; rtol=rtol) * residual
        elseif supplementary_objective === :weak_moments
            # This objective is intrinsic to Taylor consistency: lower total
            # degree moments, which enter a compact weak operator first, carry
            # more weight. A small null-space penalty selects a stable member
            # of the constrained family. No external FD/Takeuchi coefficients
            # are used as a target.
            degrees = Float64[total_degree(order_indices[row]) for row in supplementary_rows]
            weights = 1.0 ./ (1.0 .+ degrees).^2
            weighted_reduced = weights .* reduced
            weighted_residual = weights .* residual
            ridge = sqrt(eps(Float64)) * max(opnorm(weighted_reduced), 1.0)
            augmented = [weighted_reduced; ridge .* Matrix{Float64}(I, size(reduced, 2), size(reduced, 2))]
            augmented_residual = [weighted_residual; zeros(size(reduced, 2), size(residual, 2))]
            correction = directTaylorPseudoInverse(augmented; rtol=rtol) * augmented_residual
            C .+= null_exact * correction
        else
            throw(ArgumentError("unknown supplementary objective: $supplementary_objective"))
        end
    end

    target_exact = zeros(Float64, length(exact_rows), size(A, 1))
    for (row, global_row) in enumerate(exact_rows)
        target_exact[row, global_row] = 1.0
    end
    exact_error = norm(Aexact * C - target_exact, Inf)
    exact_error <= 1e3 * rtol ||
        error("hierarchical Taylor constraints were not preserved; error=$exact_error")

    return C
end

function directTaylorPseudoInverse(A; rtol=sqrt(eps(Float64)))
    F = svd(Matrix{Float64}(A))
    cutoff = rtol * maximum(F.S)
    Sinv = map(s -> s <= cutoff ? 0.0 : inv(s), F.S)
    return F.Vt' * Diagonal(Sinv) * F.U'
end

function stableTaylorPseudoInverse(A; rtol=sqrt(eps(Float64)), ridge=0.0)
    A = Matrix{Float64}(A)

    rowScale = vec(maximum(abs.(A); dims=2))
    colScale = vec(maximum(abs.(A); dims=1))
    rowScale .= ifelse.(rowScale .== 0.0, 1.0, rowScale)
    colScale .= ifelse.(colScale .== 0.0, 1.0, colScale)

    As = A ./ rowScale
    As = As ./ colScale'

    F = svd(As)
    cutoff = rtol * maximum(F.S)
    Sinv = similar(F.S)
    for i in eachindex(F.S)
        if F.S[i] <= cutoff
            Sinv[i] = 0.0
        elseif ridge > 0
            Sinv[i] = F.S[i] / (F.S[i]^2 + ridge^2)
        else
            Sinv[i] = inv(F.S[i])
        end
    end

    pinv_scaled = F.Vt' * Diagonal(Sinv) * F.U'

    # If As = Dr^{-1} * A * Dc^{-1}, then pinv(A) = Dc^{-1} * pinv(As) * Dr^{-1}.
    return (pinv_scaled ./ colScale) ./ rowScale'
end

function numbersOfTheExpression(equationCharacteristics,
                               trialFunctionsCharacteristics,
                               TaylorOptions1, TaylorOptions2)

    numbersOfTheSystem1=numbersOfTheExpression(equationCharacteristics,trialFunctionsCharacteristics,TaylorOptions1)
    numbersOfTheSystem2=numbersOfTheExpression(equationCharacteristics,trialFunctionsCharacteristics,TaylorOptions2)
    @unpack timeMarching,NtypeofExpr,NtypeofMaterialVariables,NtypeofFields,nCoordinates,Ndimension,pointsUsed,pointsμUsed,offsetsμUsed = numbersOfTheSystem1

    return (
        timeMarching = timeMarching,
        NtypeofExpr = NtypeofExpr,
        NtypeofMaterialVariables = NtypeofMaterialVariables,
        NtypeofFields = NtypeofFields,
        nCoordinates = nCoordinates,
        Ndimension = Ndimension,   # 🔥 key change
        pointsUsed = pointsUsed,
        pointsμUsed = pointsμUsed,
        offsetsμUsed = offsetsμUsed,
        pointsμᶜUsed = numbersOfTheSystem2.pointsμUsed,
        offsetsμᶜUsed = numbersOfTheSystem2.offsetsμUsed,
    )
end

function numbersOfTheExpression(equationCharacteristics,
                               trialFunctionsCharacteristics,
                               TaylorOptions)

    @unpack exprs, fields, vars, coordinates = equationCharacteristics
    @unpack pointsInSpace, pointsInTime = trialFunctionsCharacteristics
    @unpack pointsμInSpace, pointsμInTime,
            offsetμInΔyInSpace, offsetμInΔyInTime = TaylorOptions

    timeMarching = any(a -> a === timeDimensionString, string.(coordinates))

    vars = _as_variable_tuple(vars)

    NtypeofExpr = length(exprs)
    NtypeofMaterialVariables = length(vars)
    NtypeofFields = length(fields)

    Ndimension = length(coordinates)

    # 🔥 cleaner + no broadcast
    nSpatial = timeMarching ? Ndimension - 1 : Ndimension
    spatialPoints = if pointsInSpace isa Integer
        fill(Int(pointsInSpace), nSpatial)
    else
        values = Int.(collect(pointsInSpace))
        length(values) == nSpatial || throw(DimensionMismatch(
            "pointsInSpace must contain one value per spatial coordinate",
        ))
        values
    end
    all(>(0), spatialPoints) ||
        throw(ArgumentError("all spatial point counts must be positive"))
    pointsUsed = timeMarching ? [spatialPoints; Int(pointsInTime)] : spatialPoints

    # Interpolation grids may follow the actual local geometry.  In
    # particular, a clipped surface stencil is 3×2 in space and must not be
    # represented by a fictitious scalar 3×3 interpolation grid.
    spatialInterpolationPoints = if pointsμInSpace isa Integer
        fill(Int(pointsμInSpace), nSpatial)
    else
        values = Int.(collect(pointsμInSpace))
        length(values) == nSpatial || throw(DimensionMismatch(
            "field/material ptsSpace must contain one value per spatial coordinate",
        ))
        values
    end
    all(>(0), spatialInterpolationPoints) || throw(ArgumentError(
        "all field/material spatial interpolation point counts must be positive",
    ))
    spatialInterpolationOffsets = if offsetμInΔyInSpace isa Real
        fill(Float64(offsetμInΔyInSpace), nSpatial)
    else
        values = Float64.(collect(offsetμInΔyInSpace))
        length(values) == nSpatial || throw(DimensionMismatch(
            "field/material offsetSpace must contain one value per spatial coordinate",
        ))
        values
    end
    pointsμUsed = timeMarching ?
        [spatialInterpolationPoints; Int(pointsμInTime)] :
        spatialInterpolationPoints
    offsetsμUsed = timeMarching ?
        [spatialInterpolationOffsets; Float64(offsetμInΔyInTime)] :
        spatialInterpolationOffsets

    return (
        timeMarching = timeMarching,
        NtypeofExpr = NtypeofExpr,
        NtypeofMaterialVariables = NtypeofMaterialVariables,
        NtypeofFields = NtypeofFields,
        nCoordinates = Ndimension,
        Ndimension = Val(Ndimension),   # 🔥 key change
        pointsUsed = pointsUsed,
        pointsμUsed = pointsμUsed,
        offsetsμUsed = offsetsμUsed
    )
end

function investigateDependencies(equationCharacteristics,
                                 numbersOfTheSystem,
                                 trialFunctionsCharacteristics,
                                 TaylorOptionsμ,TaylorOptionsμᶜ)
    dependencies,ordersForSplinesμ,configsTaylorμ=investigateDependencies(equationCharacteristics,numbersOfTheSystem,trialFunctionsCharacteristics,TaylorOptionsμ)
    #need to change this stuff but for the moment this is the least thing i can do ... (02/05/2026)
    numbersOfTheSystemTmp = merge(numbersOfTheSystem,(pointsμUsed = numbersOfTheSystem.pointsμᶜUsed,
        offsetsμUsed = numbersOfTheSystem.offsetsμᶜUsed,),)
    _,ordersForSplinesμᶜ,configsTaylorμᶜ= investigateDependencies(equationCharacteristics,numbersOfTheSystemTmp,trialFunctionsCharacteristics,TaylorOptionsμᶜ)
    return dependencies,ordersForSplinesμ,configsTaylorμ,ordersForSplinesμᶜ,configsTaylorμᶜ
end


function investigateDependencies(equationCharacteristics,
                                 numbersOfTheSystem,
                                 trialFunctionsCharacteristics,
                                 TaylorOptions)

        @unpack Ndimension = numbersOfTheSystem

    return _investigateDependencies(Ndimension,
                                   equationCharacteristics,
                                   numbersOfTheSystem,
                                   trialFunctionsCharacteristics,
                                   TaylorOptions)
end

function _investigateDependencies(::Val{N},
                                 equationCharacteristics,
                                 numbersOfTheSystem,
                                 trialFunctionsCharacteristics,
                                 TaylorOptions) where N

    @unpack exprs,fields,vars,∂ = equationCharacteristics
    @unpack timeMarching,NtypeofExpr,NtypeofMaterialVariables,NtypeofFields,
            pointsUsed,pointsμUsed,offsetsμUsed = numbersOfTheSystem
    @unpack orderBtime,orderBspace = trialFunctionsCharacteristics
    @unpack YorderBtime,YorderBspace,supplementaryOrder = TaylorOptions

    # ---------------- Dependencies ----------------

    variableDependency = ones(Int, N)
    fieldDependency    = ones(Int, N)

    eachVariableDependency = zeros(Int, N, NtypeofMaterialVariables)
    eachFieldDependency    = zeros(Int, N, NtypeofFields)

    for iFields in 1:NtypeofFields
        dep = findCartesianDependency(fields[iFields], N, ∂)
        eachFieldDependency[:,iFields] = dep
        fieldDependency .*= (1 .- dep)
    end

    for iVars in 1:NtypeofMaterialVariables
        @show vars, iVars
        dep = findCartesianDependency(vars[iVars], N, ∂)
        eachVariableDependency[:,iVars] = dep
        variableDependency .*= (1 .- dep)
    end

    fieldDependency .= 1 .- fieldDependency
    variableDependency .= (1 .- variableDependency) .* fieldDependency

    # ---------------- Orders ----------------

    orderBspline  = zeros(Int, N)
    YorderBspline = zeros(Int, N)

    if timeMarching
        orderBspline[end] = orderBtime * fieldDependency[end]
        orderBspline[1:end-1] .= orderBspace .* fieldDependency[1:end-1]

        YorderBspline[end] = YorderBtime * fieldDependency[end]
        YorderBspline[1:end-1] .= YorderBspace .* fieldDependency[1:end-1]
    else
        orderBspline  .= orderBspace .* fieldDependency
        YorderBspline .= YorderBspace .* fieldDependency
    end

    pointsUsedForFields = (pointsUsed .- 1) .* fieldDependency .+ 1
    orderExpressions = pointsUsedForFields
    orderU = (orderExpressions .- 1) .+ (supplementaryOrder .* fieldDependency) .+ 1

    # ---------------- Taylor grids ----------------

    orderTaylors = Array{Any,N}(undef, Tuple(orderU))
    pointsInSpaceTime = Array{Any,N}(undef, Tuple(pointsUsedForFields))

    multiOrdersIndices = CartesianIndices(orderTaylors)
    multiPointsIndices = CartesianIndices(pointsInSpaceTime)

    availablePointsConfigurations = Vector{Any}()
    availableμPoints = Vector{Any}()
    availableμaxes = Vector{Any}()
    centrePointConfigurations = Int[]


    # ---------------- Candidate positions of ν ----------------

    tmpVec = ((car2vec(multiPointsIndices[end]) .- 1) .÷ 2) .+ 1
    if timeMarching && car2vec(multiPointsIndices[end])[end] > 1
        tmpVec[end] = car2vec(multiPointsIndices[end])[end] - 1
    end
    requestedCentre = hasproperty(trialFunctionsCharacteristics, :nuCentre) ?
        trialFunctionsCharacteristics.nuCentre : nothing
    if !isnothing(requestedCentre)
        centre = Int.(collect(requestedCentre))
        if timeMarching && length(centre) == N - 1
            push!(centre, tmpVec[end])
        end
        length(centre) == N || throw(DimensionMismatch(
            "nuCentre must contain one value per spatial coordinate",
        ))
        all(d -> 1 <= centre[d] <= size(pointsInSpaceTime, d), eachindex(centre)) ||
            throw(ArgumentError("nuCentre lies outside the local stencil"))
        tmpVec .= centre
    end
    middleν = vec2car(tmpVec)
    nuGeometryMode = hasproperty(trialFunctionsCharacteristics, :nuGeometryMode) ?
        trialFunctionsCharacteristics.nuGeometryMode : :middle

    candidateCentres = if nuGeometryMode == :all
        if timeMarching
            # Move ν through every spatial position, while retaining the
            # established time level used to define one marching step.
            [I for I in multiPointsIndices if I[N] == middleν[N]]
        else
            collect(multiPointsIndices)
        end
    else
        [middleν]
    end
    numberGeometries = length(candidateCentres)

    for middleν in candidateCentres

        push!(availablePointsConfigurations, car2svec.(multiPointsIndices))
        push!(centrePointConfigurations,
            LinearIndices(multiPointsIndices)[middleν])

        # ---------------- μ coordinates ----------------

        tmpμCoordinates = Array{SVector{N,Float64}}(undef, pointsμUsed...)

        tmpDistances = Float64.(availablePointsConfigurations[1][end] .- 1)
        tmpΔμ = tmpDistances .- 2 .* offsetsμUsed

        for i in 1:N
            if pointsμUsed[i] > 1
                tmpΔμ[i] /= (pointsμUsed[i] - 1)
            end
        end

        μaxes = ntuple(d -> begin
            [1.0 + offsetsμUsed[d] + (j - 1) * tmpΔμ[d] for j in 1:size(tmpμCoordinates, d)]
        end, N)

        for I in CartesianIndices(tmpμCoordinates)
            tmpμCoordinates[I] = SVector{N}(ntuple(d -> μaxes[d][I[d]], N))
        end


        push!(availableμPoints, tmpμCoordinates)
        push!(availableμaxes,μaxes)
    end
    
    

    # ---------------- Outputs ----------------

    dependencies = (
        variableDependency = variableDependency,
        fieldDependency = fieldDependency,
        eachVariableDependency = eachVariableDependency,
        eachFieldDependency = eachFieldDependency
    )

    ordersForSplines = (
        orderBspline = orderBspline,
        YorderBspline = YorderBspline,
        orderExpressions = orderExpressions,
        orderU = orderU
    )

    configsTaylor = (
        numberGeometries = numberGeometries,
        multiOrdersIndices = multiOrdersIndices, # we still need this since available points can differ
        availablePointsConfigurations = availablePointsConfigurations,
        centrePointConfigurations = centrePointConfigurations,
        availableμPoints = availableμPoints,
        availableμaxes = availableμaxes,
        exactTaylorTotalDegree = hasproperty(
            trialFunctionsCharacteristics, :exactTaylorTotalDegree,
        ) ? trialFunctionsCharacteristics.exactTaylorTotalDegree : nothing,
    )

    return dependencies, ordersForSplines, configsTaylor
end

function bigαFinder(equationCharacteristics,numbersOfTheSystem,ordersForSplines)

    @unpack NtypeofExpr, NtypeofFields = numbersOfTheSystem
    @unpack exprs,fields,vars,coordinates,∂ = equationCharacteristics
    @unpack orderExpressions = ordersForSplines

    bigα=Array{Any,2}(missing,NtypeofExpr,NtypeofFields)
    varM=nothing
    pointsUsedForFields=orderExpressions

    for iExpr in eachindex(exprs)
        for iField in eachindex(fields)
            
            tmpNonZeroAlphas=PDECoefFinder(orderExpressions,coordinates,exprs[iExpr],fields[iField],vars) 
            # we assume that the pointsUsedForFields represent the highest order of partials
            bigα[iExpr,iField]=unique(tmpNonZeroAlphas)
        end
    end
    varM,CartesianDependencies=varMmaker(pointsUsedForFields,coordinates,vars,∂)
    return bigα,varM,CartesianDependencies
end
