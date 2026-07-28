using Symbolics
using SparseArrays
using LinearAlgebra
using JLD2
using DrWatson

# Optional plotting support inside timeMarchingSchemeLinear
# Uncomment if you want videoMode=true.
# using CairoMakie
# CairoMakie.activate!()


_numop_get(x, name::Symbol, default=nothing) = x isa AbstractDict ? get(x, String(name), get(x, name, default)) : (hasproperty(x, name) ? getproperty(x, name) : default)

function _linear_indices_from_operator(op)
    geometry = op.geometry
    NpointsSpace = length(geometry.νWhole)
    activeTimePoints = geometry.activeTimePoints
    spaceShape = Tuple(geometry.wholeRegionPointsSpace)
    return NpointsSpace, activeTimePoints, spaceShape
end

function _split_operator_columns(op, NField, NpointsSpace, activeTimePoints, spaceShape, target)
    nEq = op.size[1]
    if target == :unknown
        nCol = NpointsSpace * NField
    elseif target == :knownfield
        nCol = NpointsSpace * NField * max(activeTimePoints - 1, 0)
    elseif target == :knownforce
        nCol = NpointsSpace * NField * activeTimePoints
    else
        throw(ArgumentError("unknown target $target"))
    end

    rows = Int[]
    cols = Int[]
    vals = eltype(op.table.vals)[]
    fieldTimeSpace = CartesianIndices((NField, activeTimePoints, spaceShape...))
    pointLinear = LinearIndices(spaceShape)

    for k in eachindex(op.table.vals)
        ci = fieldTimeSpace[Int(op.table.cols[k])]
        iField = ci[1]
        iT = ci[2]
        jPoint = CartesianIndex(Tuple(ci)[3:end])
        lp = pointLinear[jPoint]

        if target == :unknown
            iT == activeTimePoints || continue
            col = lp + (iField - 1) * NpointsSpace
        elseif target == :knownfield
            iT < activeTimePoints || continue
            col = lp + (iField - 1) * NpointsSpace + (iT - 1) * NpointsSpace * NField
        else
            col = lp + (iField - 1) * NpointsSpace + (iT - 1) * NpointsSpace * NField
        end

        push!(rows, Int(op.table.rows[k]))
        push!(cols, col)
        push!(vals, op.table.vals[k])
    end

    return sparse(rows, cols, vals, nEq, nCol)
end

function _replace_sparse_rows(matrix, replacements::Dict{Int,<:Any})
    I, J, V = findnz(matrix)
    replaced = Set(keys(replacements))
    keep = [!(row in replaced) for row in I]
    rows, cols, vals = I[keep], J[keep], V[keep]
    for (row, entries) in replacements
        for (col, value) in entries
            iszero(value) && continue
            push!(rows, row); push!(cols, col); push!(vals, value)
        end
    end
    sparse(rows, cols, vals, size(matrix)...)
end

function _zero_sparse_rows(matrix, rows_to_zero)
    I, J, V = findnz(matrix)
    removed = Set(rows_to_zero)
    keep = [!(row in removed) for row in I]
    sparse(I[keep], J[keep], V[keep], size(matrix)...)
end

function _elastic_free_surface_rows_2d(left, NField, NpointsSpace,
                                       spaceShape, spacing)
    geometry = left.geometry
    boundary = geometry.freeSurfaceBoundary
    isnothing(boundary) && return Dict{Int,Vector{Tuple{Int,Float64}}}()
    NField == 2 ||
        throw(DimensionMismatch("the 2D elastic boundary needs two fields"))
    length(geometry.Models) >= 3 ||
        throw(ArgumentError("the elastic model must contain ρ, λ and μ"))
    left.size[1] == NField * NpointsSpace ||
        throw(ArgumentError(
            "free-surface row replacement currently requires one OPT test block",
        ))
    λmodel, μmodel = geometry.Models[2], geometry.Models[3]
    point_linear = LinearIndices(spaceShape)
    residual_linear = LinearIndices((NField, NpointsSpace))
    dx, dz = Float64.(spacing)
    replacements = Dict{Int,Vector{Tuple{Int,Float64}}}()

    for (point, normal) in zip(boundary.points, boundary.normals)
        model_point = geometry.conv.whole2model(point)
        λ = Float64(λmodel[model_point])
        μ = Float64(μmodel[model_point])
        nx, nz = Float64.(normal)
        lp = point_linear[point]
        rows = (residual_linear[1, lp], residual_linear[2, lp])
        coefficients = (Dict{Int,Float64}(), Dict{Int,Float64}())

        function add_derivative!(equation, field, axis, factor)
            iszero(factor) && return
            if axis == 1
                minus = point - CartesianIndex(1, 0)
                plus = point + CartesianIndex(1, 0)
                derivative = ((minus, -inv(2dx)), (plus, inv(2dx)))
            else
                # The elastic medium is below the extracted topographic surface.
                minus = point - CartesianIndex(0, 1)
                derivative = ((minus, -inv(dz)), (point, inv(dz)))
            end
            for (sample_point, weight) in derivative
                checkbounds(Bool, point_linear, sample_point) ||
                    throw(BoundsError(point_linear, sample_point))
                column = point_linear[sample_point] +
                         (field - 1) * NpointsSpace
                coefficients[equation][column] =
                    get(coefficients[equation], column, 0.0) + factor * weight
            end
        end

        # t_x = n_x σ_xx + n_z σ_xz
        add_derivative!(1, 1, 1, nx * (λ + 2μ))
        add_derivative!(1, 1, 2, nz * μ)
        add_derivative!(1, 2, 1, nz * μ)
        add_derivative!(1, 2, 2, nx * λ)
        # t_z = n_x σ_xz + n_z σ_zz
        add_derivative!(2, 1, 1, nz * λ)
        add_derivative!(2, 1, 2, nx * μ)
        add_derivative!(2, 2, 1, nx * μ)
        add_derivative!(2, 2, 2, nz * (λ + 2μ))
        replacements[rows[1]] = [(column, value) for (column, value) in coefficients[1]]
        replacements[rows[2]] = [(column, value) for (column, value) in coefficients[2]]
    end
    replacements
end

function prepareNumericalLinearSystem(numOperators; T=Float64,
                                      free_surface_spacing=(1.0, 1.0))
    numericalOperators = _numop_get(numOperators, :numericalOperators)
    numericalOperators === nothing && error("numOperators does not contain numericalOperators")

    left = _numop_get(numericalOperators, :left)
    right = _numop_get(numericalOperators, :right)
    left === nothing && error("numOperators.numericalOperators.left is missing")
    right === nothing && error("numOperators.numericalOperators.right is missing")

    NpointsSpace, activeTimePoints, spaceShape = _linear_indices_from_operator(left)
    NField = div(left.size[2], NpointsSpace * activeTimePoints)
    NForceField = div(right.size[2], NpointsSpace * activeTimePoints)
    timePointsUsedForOneStep = activeTimePoints
    NknownTime = max(timePointsUsedForOneStep - 1, 0)
    NforcePoints = NpointsSpace

    A_unknown = _split_operator_columns(left, NField, NpointsSpace, activeTimePoints, spaceShape, :unknown)
    L_known = _split_operator_columns(left, NField, NpointsSpace, activeTimePoints, spaceShape, :knownfield)
    R_force = _split_operator_columns(right, NForceField, NpointsSpace, activeTimePoints, spaceShape, :knownforce)

    boundary_conditions = _numop_get(numOperators, :boundaryConditions)
    free_surface_rows = Int[]
    if !isnothing(boundary_conditions) &&
       !isnothing(boundary_conditions.free_surface) &&
       length(spaceShape) == 2
        replacements = _elastic_free_surface_rows_2d(
            left, NField, NpointsSpace, spaceShape, free_surface_spacing,
        )
        free_surface_rows = sort!(collect(keys(replacements)))
        A_unknown = _replace_sparse_rows(A_unknown, replacements)
        L_known = _zero_sparse_rows(L_known, free_surface_rows)
        R_force = _zero_sparse_rows(R_force, free_surface_rows)
    end

    b_template = zeros(promote_type(T, eltype(left.table.vals), eltype(right.table.vals)), left.size[1])
    known_lhs_template = zeros(eltype(b_template), NpointsSpace, NField, NknownTime)
    known_rhs_template = zeros(eltype(b_template), NforcePoints, NForceField, timePointsUsedForOneStep)

    function b_fun!(b, knownInputs)
        nKnownField = length(known_lhs_template)
        knownFieldVec = @view knownInputs[1:nKnownField]
        knownForceVec = @view knownInputs[nKnownField+1:nKnownField+length(known_rhs_template)]
        b .= 0
        if nKnownField > 0
            mul!(b, L_known, knownFieldVec, -one(eltype(b)), one(eltype(b)))
        end
        mul!(b, R_force, knownForceVec, one(eltype(b)), one(eltype(b)))
        return b
    end

    function A_fun!(avec, knownInputs)
        avec .= vec(A_unknown)
        return avec
    end

    return (
        is_numerical = true,
        residual_operator = _numop_get(numericalOperators, :residual),
        left_operator = left,
        right_operator = right,
        A_unknown = A_unknown,
        L_known = L_known,
        R_force = R_force,
        A_fun! = A_fun!,
        b_fun! = b_fun!,
        A_template = copy(A_unknown),
        b_template = b_template,
        known_lhs_template = known_lhs_template,
        known_rhs_template = known_rhs_template,
        spaceShape = spaceShape,
        NpointsSpace = NpointsSpace,
        NforcePoints = NforcePoints,
        NField = NField,
        NForceField = NForceField,
        timePointsUsedForOneStep = timePointsUsedForOneStep,
        free_surface_rows = free_surface_rows,
    )
end

function prepareLinearSystem(numOperators; T=Float64, kwargs...)
    if _numop_get(numOperators, :numericalOperators) !== nothing
        return prepareNumericalLinearSystem(numOperators; T=T, kwargs...)
    end

    @unpack costfunctions, fieldLHS, fieldRHS, champsLimité = numOperators

    costvec = vec(costfunctions)

    timePointsUsedForOneStep = size(fieldLHS, 2)
    NField = size(fieldLHS, 1)

    spaceShape = size(fieldLHS[1, 1])
    Rcoord = CartesianIndices(spaceShape)
    linearR = LinearIndices(Rcoord)
    NpointsSpace = length(Rcoord)
    NknownTime = max(timePointsUsedForOneStep - 1, 0)

    symbUnknownField = Array{Num,2}(undef, NpointsSpace, NField)
    symbKnownField = Array{Num,3}(undef, NpointsSpace, NField, NknownTime)

    for iT in 1:NknownTime
        for iField in 1:NField
            A = fieldLHS[iField, iT]
            for j in Rcoord
                lj = linearR[j]
                symbKnownField[lj, iField, iT] = A[j]
            end
        end
    end

    for iField in 1:NField
        A = fieldLHS[iField, timePointsUsedForOneStep]
        for j in Rcoord
            lj = linearR[j]
            symbUnknownField[lj, iField] = A[j]
        end
    end

    if champsLimité === nothing
        NforcePoints = NpointsSpace
        symbKnownForce = Array{Num,3}(undef, NforcePoints, NField, timePointsUsedForOneStep)
        for iT in 1:timePointsUsedForOneStep
            for iField in 1:NField
                A = fieldRHS[iField, iT]
                for j in Rcoord
                    lj = linearR[j]
                    symbKnownForce[lj, iField, iT] = A[j]
                end
            end
        end
    else
        NforcePoints = length(champsLimité[1, 1])
        symbKnownForce = Array{Num,3}(undef, NforcePoints, NField, timePointsUsedForOneStep)
        for iT in 1:timePointsUsedForOneStep
            for iField in 1:NField
                for j in 1:NforcePoints
                    symbKnownForce[j, iField, iT] = champsLimité[iField, iT][j]
                end
            end
        end
    end

    unknownInputs = vec(symbUnknownField)
    knownFieldInputs = NknownTime == 0 ? Num[] : vec(symbKnownField)
    knownForceInputs = vec(symbKnownForce)
    knownInputs = vcat(knownFieldInputs, knownForceInputs)

    # residual = A*u + c  =>  A*u = -c
    A_sym = Symbolics.jacobian(costvec, unknownInputs)
    zero_map = Dict(u => 0 for u in unknownInputs)
    c_sym = Symbolics.substitute.(costvec, Ref(zero_map))
    b_sym = .-c_sym

    A_expr = build_function(vec(A_sym), knownInputs; expression=Val(false))
    b_expr = build_function(b_sym, knownInputs; expression=Val(false))

    A_tuple = eval(A_expr)
    b_tuple = eval(b_expr)

    A_fun = A_tuple[1]
    A_fun! = A_tuple[2]

    b_fun = b_tuple[1]
    b_fun! = b_tuple[2]

    nUnknown = length(unknownInputs)
    nEq = length(costvec)

    known_lhs_template = zeros(T, NpointsSpace, NField, NknownTime)
    known_rhs_template = zeros(T, NforcePoints, NField, timePointsUsedForOneStep)

    A_template = zeros(T, nEq, nUnknown)
    b_template = zeros(T, nEq)

    return (
        A_fun = A_fun,
        A_fun! = A_fun!,
        b_fun = b_fun,
        b_fun! = b_fun!,
        A_template = A_template,
        b_template = b_template,
        known_lhs_template = known_lhs_template,
        known_rhs_template = known_rhs_template,
        spaceShape = spaceShape,
        NpointsSpace = NpointsSpace,
        NforcePoints = NforcePoints,
        NField = NField,
        timePointsUsedForOneStep = timePointsUsedForOneStep,
    )
end

function makeLinearEvaluator(prepared)
    nKnownField = length(prepared.known_lhs_template)
    nKnownForce = length(prepared.known_rhs_template)
    knownInputs = zeros(eltype(prepared.b_template), nKnownField + nKnownForce)

    function eval!(A, b, knownField, knownForce)
        knownInputs[1:nKnownField] .= vec(knownField)
        knownInputs[nKnownField+1:nKnownField+nKnownForce] .= vec(knownForce)

        if hasproperty(prepared, :is_numerical) && prepared.is_numerical
            prepared.b_fun!(b, knownInputs)
        else
            prepared.A_fun!(vec(A), knownInputs)
            prepared.b_fun!(b, knownInputs)
        end
        return A, b
    end

    return eval!
end

function evaluateLinearSystem(prepared, knownField, knownForce; sparse_output=false)
    eval! = makeLinearEvaluator(prepared)
    A = copy(prepared.A_template)
    b = copy(prepared.b_template)
    eval!(A, b, knownField, knownForce)
    if sparse_output
        A = sparse(A)
    end
    return A, b
end

function applyBoundaryConditionForced!(A, b; leftValue=1.0, rightValue=1.0)
    A[1, :] .= 0
    A[1, 1] = 1
    b[1] = leftValue

    A[end, :] .= 0
    A[end, end] = 1
    b[end] = rightValue

    return A, b
end

function prepareConstantMatrix(preparedLin;
    boundaryConditionForced=false,
    sparse_output=true,
    leftValue=1.0,
    rightValue=1.0,
)
    eval! = makeLinearEvaluator(preparedLin)

    knownField0 = zero(preparedLin.known_lhs_template)
    knownForce0 = zero(preparedLin.known_rhs_template)

    A = copy(preparedLin.A_template)
    btmp = copy(preparedLin.b_template)

    eval!(A, btmp, knownField0, knownForce0)

    if sparse_output
        A = sparse(A)
    elseif hasproperty(preparedLin, :is_numerical) && preparedLin.is_numerical
        A = Matrix(A)
    end

    if boundaryConditionForced
        applyBoundaryConditionForced!(A, btmp; leftValue=leftValue, rightValue=rightValue)
    end

    factor = lu(A)
    return A, factor
end

"""
    propagateLinearSystem(prepared, Nt, dt; sourceFull, output_stride=1)

Run an already prepared OPT time-marching system in memory and return field
snapshots with layout `(space..., field, output_time)`. `sourceFull` uses
`(space_point, force_field, source_time)` and must include the extra time
levels required by the OPT stencil.
"""
function propagateLinearSystem(
    prepared,
    Nt::Integer,
    dt::Real;
    sourceFull,
    output_stride::Integer=1,
    initialCondition=0.0,
)
    Nt > 0 || throw(ArgumentError("Nt must be positive"))
    output_stride > 0 || throw(ArgumentError("output_stride must be positive"))
    minimum_source_times = Nt + prepared.timePointsUsedForOneStep - 1
    size(sourceFull) == (
        prepared.NforcePoints,
        prepared.NForceField,
        size(sourceFull, 3),
    ) || throw(DimensionMismatch("sourceFull has incompatible space/field axes"))
    size(sourceFull, 3) >= minimum_source_times ||
        throw(DimensionMismatch(
            "sourceFull needs at least $minimum_source_times time samples",
        ))

    knownField = fill(Float64(initialCondition), size(prepared.known_lhs_template))
    knownForce = fill(Float64(initialCondition), size(prepared.known_rhs_template))
    unknownField = fill(Float64(initialCondition), prepared.NpointsSpace, prepared.NField)
    _, factor = prepareConstantMatrix(prepared; sparse_output=true)
    b = copy(prepared.b_template)
    stored_steps = unique(vcat(0, collect(output_stride:output_stride:Nt), Nt))
    stored_lookup = Set(stored_steps)
    snapshots = Array{Float64}[]
    times = Float64[]
    push!(snapshots, reshape(copy(unknownField), prepared.spaceShape..., prepared.NField))
    push!(times, 0.0)
    nKnownTime = size(knownField, 3)

    for step in 1:Nt
        knownForce .= sourceFull[:, :, step:step+prepared.timePointsUsedForOneStep-1]
        prepared.b_fun!(b, vcat(vec(knownField), vec(knownForce)))
        unknownField .= reshape(factor \ b, prepared.NpointsSpace, prepared.NField)
        if nKnownTime > 0
            nKnownTime > 1 &&
                (knownField[:, :, 1:end-1] .= knownField[:, :, 2:end])
            knownField[:, :, end] .= unknownField
        end
        if step in stored_lookup
            push!(snapshots,
                  reshape(copy(unknownField), prepared.spaceShape..., prepared.NField))
            push!(times, step * Float64(dt))
        end
    end
    history = cat(snapshots...; dims=length(prepared.spaceShape) + 2)
    (; history, times, stored_steps, dt=Float64(dt))
end

function timeMarchingSchemeLinear(
    preparedLin,
    Nt,
    Δnum,
    modelName;
    videoMode=true,
    sourceType="Ricker",
    t₀=50,
    f₀=0.04,
    initialCondition=0.0,
    sourceFull=nothing,
    iExperiment=nothing,
    sparse_output=true,
    boundaryConditionForced=false,
    assume_constant_matrix=true,
)
    mkpath(datadir("fieldResults"))


    if isnothing(iExperiment)
        iExperiment = (iH=1, iCase=1, iPointsUsed=1)
    end

    @unpack iH, iCase, iPointsUsed = iExperiment
    experiment_name = "$(iH)_$(iCase)_$(iPointsUsed)"

    sequentialFileName = datadir("fieldResults", savename((Nt, Δnum..., modelName, sourceType, experiment_name, "linear"), "jld2"))
    compactFileName = datadir("fieldResults", savename("compact", (Nt, Δnum..., modelName, sourceType, experiment_name, "linear"), "jld2"))

    NField = preparedLin.NField
    NForceField = hasproperty(preparedLin, :NForceField) ? preparedLin.NForceField : preparedLin.NField
    NpointsSpace = preparedLin.NpointsSpace
    NforcePoints = preparedLin.NforcePoints
    timePointsUsedForOneStep = preparedLin.timePointsUsedForOneStep
    spaceShape = preparedLin.spaceShape
    NknownTime = max(timePointsUsedForOneStep - 1, 0)

    knownField = fill(Float64(initialCondition), size(preparedLin.known_lhs_template))
    knownForce = fill(Float64(initialCondition), size(preparedLin.known_rhs_template))
    unknownField = fill(Float64(initialCondition), NpointsSpace, NField)

    itVec = collect(1:Nt)
    t = (itVec .- 1) .* Δnum[end]

    sourceTime = nothing
    if sourceType == "Ricker"
        myRicker(x) = Ricker(x, t₀, f₀)
        sourceTime = myRicker.(t)
        if timePointsUsedForOneStep != 1
            prepend!(sourceTime, zeros(timePointsUsedForOneStep))
        end
    elseif sourceType == "Explicit"
        size(sourceFull, 1) == NforcePoints || error("sourceFull first dimension should be NforcePoints=$NforcePoints")
        size(sourceFull, 2) == NForceField || error("sourceFull second dimension should be NForceField=$NForceField")
        minimumTimeSamples = Nt + timePointsUsedForOneStep - 1
        size(sourceFull, 3) >= minimumTimeSamples ||
            error("sourceFull third dimension should be at least $minimumTimeSamples for Nt=$Nt and timePointsUsedForOneStep=$timePointsUsedForOneStep")
    end

    if !isfile(sequentialFileName)
        fieldFile = jldopen(sequentialFileName, "w")

        hm = nothing
        fig = nothing
        if videoMode
            fig = Figure()
            ax = Axis(fig[1, 1])
            firstField = reshape(unknownField[:, 1], spaceShape...)
            hm = heatmap!(ax, Float32.(firstField), colormap=:deep, colorrange=(-1e-5, 1e-5))
            display(fig)
        end

        eval! = makeLinearEvaluator(preparedLin)
        Awork = copy(preparedLin.A_template)
        b = copy(preparedLin.b_template)

        factor = nothing
        if assume_constant_matrix
            _, factor = prepareConstantMatrix(
                preparedLin;
                boundaryConditionForced=boundaryConditionForced,
                sparse_output=sparse_output,
            )
        end

        for it in itVec
            if sourceType == "Ricker"
                knownForce .= 0.0
                knownForce[1, 1, :] .= sourceTime[it:it+timePointsUsedForOneStep-1]
            else
                knownForce .= sourceFull[:, :, it:it+timePointsUsedForOneStep-1]
            end

            if assume_constant_matrix
                preparedLin.b_fun!(b, vcat(vec(knownField), vec(knownForce)))
                if boundaryConditionForced
                    b[1] = 1.0
                    b[end] = 1.0
                end
            else
                eval!(Awork, b, knownField, knownForce)
                A = sparse_output ? sparse(Awork) : (hasproperty(preparedLin, :is_numerical) && preparedLin.is_numerical ? Matrix(Awork) : copy(Awork))
                if boundaryConditionForced
                    applyBoundaryConditionForced!(A, b; leftValue=1.0, rightValue=1.0)
                end
                factor = lu(A)
            end

            u = factor \ b
            unknownField .= reshape(u, NpointsSpace, NField)

            if NknownTime > 0
                if NknownTime > 1
                    knownField[:, :, 1:end-1] .= knownField[:, :, 2:end]
                end
                knownField[:, :, end] .= unknownField
            end

            newField = reshape(unknownField, NpointsSpace, NField)
            fieldFile["timestep_$it"] = reshape(newField, spaceShape..., NField)

            if videoMode
                firstField = reshape(unknownField[:, 1], spaceShape...)
                hm[1] = Float32.(firstField)
            end
        end

        close(fieldFile)
    end

    file = jldopen(sequentialFileName, "r")
    a = [file["timestep_$it"] for it in itVec]
    close(file)

    acompact = reduce(hcat, a)
    jldsave(compactFileName; acompact=acompact)

    return acompact
end
