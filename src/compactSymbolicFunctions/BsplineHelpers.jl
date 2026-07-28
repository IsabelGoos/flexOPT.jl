using Symbolics

# this source file was truly a fruit of some weeks of optimisations with codex (chatGPT) March-April 2026
# Nobuaki Fuji, IPGP/UPC/IUF

# Public API still uses integer orders: -1, 0, 1, 2, ...
# Internally we immediately map that to a semantic mode so that -1 never
# becomes an accidental array slot through order + 1 arithmetic.
function bspline_order_spec(order::Integer)
    order < -1 && throw(ArgumentError("B-spline order must be -1 for :indicator or >= 0 for :bspline; got $order"))
    return (kind = order == -1 ? :indicator : :bspline, order = Int(order))
end

bspline_order_spec(params) = bspline_order_spec(params.maximumOrder)
bspline_order_slot(spec::NamedTuple) = spec.kind == :bspline ? spec.order + 1 :
    throw(ArgumentError(":indicator has no B-spline storage slot"))
is_indicator_family(family) = family !== nothing && hasproperty(family, :kind) && family.kind == :indicator
is_bspline_family(family) = family !== nothing && hasproperty(family, :kind) && family.kind == :bspline

"""
    constructBsplineFamily(params; simplify_expr=Symbolics.simplify, boundary_mode=:clamped)

Build a limited 1D B-spline family on the node grid described by
`(maximumOrder, allNodes, idx_nodesNum, idx_refPoints)`.

The construction:

- augments the requested reference points with numerical left/right closures
  when needed,
- pads the resulting knot list with `maximumOrder` auxiliary knots on each side,
- constructs the full auxiliary family,
- returns only the centre functions corresponding to the original
  `idx_refPoints`.

The returned `b` array stores piecewise symbolic expressions with layout
`(segment, function, order_slot)`, where `order_slot = p + 1` corresponds to
degree `p`.
"""
function constructBsplineFamily(params;simplify_expr=mySimplify,boundary_mode=:ghost,correction_truncation=true,)
    
    orderSpec = bspline_order_spec(params)
    maximumOrder = orderSpec.order
    allNodes = collect(params.allNodes)
    idx_nodesNum = collect(params.idx_nodesNum)
    idx_refPoints_original = collect(params.idx_refPoints)
    idx_selectedPoints = collect(params.idx_selectedPoints)

    issubset(Set(idx_selectedPoints), Set(idx_refPoints_original)) ||
        error("idx_selectedPoints should be a subset of idx_refPoints_original")

    selected_positions = Int.(indexin(idx_selectedPoints, idx_refPoints_original))
    selected_positions = collect(selected_positions)

    @variables x Δx
    variables = Num[x, Δx]

    rangeSegments = idx_nodesNum[1]:idx_nodesNum[end]
    segmentNodes = allNodes[rangeSegments]
    numberFunctionsOriginal = length(idx_refPoints_original)


    if orderSpec.kind == :indicator
        return (
            kind = :indicator,
            order = -1,
            b = nothing,
            b_full = nothing,
            b_support = nothing,
        )
    end

    NorderShiftedByOne = maximumOrder + 1

    core_idx_refPoints = copy(idx_refPoints_original)
    leftPlusNumberFunctions = 0
    rightPlusNumberFunctions = 0

    if allNodes[core_idx_refPoints[1]] > allNodes[idx_nodesNum[1]]
        core_idx_refPoints = vcat(idx_nodesNum[1], core_idx_refPoints)
        leftPlusNumberFunctions = 1
    end

    if allNodes[core_idx_refPoints[end]] < allNodes[idx_nodesNum[end]]
        core_idx_refPoints = vcat(core_idx_refPoints, idx_nodesNum[end])
        rightPlusNumberFunctions = 1
    end

    coreKnots = allNodes[core_idx_refPoints]
    refKnots = extend_knots(coreKnots, maximumOrder; mode=boundary_mode)
    refKnotsSymbolic = Δx .* refKnots

    numberFunctions = length(refKnots) - 1

    denominators = zeros(Num, 2, numberFunctions, maximumOrder)
    for p in 1:maximumOrder
        for idx in 1:numberFunctions-p
            denominators[1, idx, p] = refKnotsSymbolic[idx + p] - refKnotsSymbolic[idx]
            denominators[2, idx, p] = refKnotsSymbolic[idx + p + 1] - refKnotsSymbolic[idx + 1]
        end
    end

    b_basis_full = CompactSymbolicFunctions(
        allNodes,
        numberFunctions;
        auxDims=(NorderShiftedByOne,),
        variables=variables,
    )

    b_basis_support = CompactSymbolicFunctions(
        segmentNodes,
        numberFunctionsOriginal;
        auxDims=(NorderShiftedByOne,),
        variables=variables,
    )

    for p in 0:maximumOrder
        slot = p + 1

        if p == 0
            for idx in 1:numberFunctions
                leftKnot = refKnots[idx]
                rightKnot = refKnots[idx + 1]
                support_segments = support_segment_indices(leftKnot, rightKnot, allNodes)

                for seg in support_segments
                    b_basis_full.data[seg, idx, slot] = 1.0
                end
            end
        else
            b_basis_full.data[:, :, slot] .= 0.0

            for idx in 1:numberFunctions-p
                up = denominators[1, idx, p]
                down = denominators[2, idx, p]
                leftKnot = refKnotsSymbolic[idx]
                rightKnot = refKnotsSymbolic[idx + p + 1]

                active_segments_left = support_segment_indices(refKnots[idx], refKnots[idx + p], allNodes)
                active_segments_right = support_segment_indices(refKnots[idx + 1], refKnots[idx + p + 1], allNodes)

                for seg in active_segments_left
                    if !iszero_denominator(up)
                        b_basis_full.data[seg, idx, slot] += simplify_expr(
                            (x - leftKnot) / up * b_basis_full.data[seg, idx, slot - 1]
                        )
                    end
                end

                for seg in active_segments_right
                    if !iszero_denominator(down)
                        b_basis_full.data[seg, idx, slot] += simplify_expr(
                            (rightKnot - x) / down * b_basis_full.data[seg, idx + 1, slot - 1]
                        )
                    end
                end
            end

            if boundary_mode == :clamped
                rebuild_clamped_residuals!(
                    b_basis_full.data,
                    allNodes,
                    refKnots,
                    numberFunctions,
                    slot,
                    p,
                    simplify_expr,
                )
            end
        end

        firstCentral = leftPlusNumberFunctions + maximumOrder + 1 - div(p + 1, 2)
        lastCentral = firstCentral + numberFunctionsOriginal - 1
        centre_range = firstCentral:lastCentral

        b_basis_support.data[:, :, slot] .= b_basis_full.data[rangeSegments, centre_range, slot]

        if correction_truncation && p > 0
            left_support = support_segment_indices(refKnots[firstCentral], refKnots[firstCentral + p + 1], allNodes)
            right_support = support_segment_indices(refKnots[lastCentral], refKnots[lastCentral + p + 1], allNodes)

            b_basis_support.data[:, 1, slot] .= 0.0
            b_basis_support.data[:, end, slot] .= 0.0

            for seg in left_support
                seg_local = local_segment_index(seg, rangeSegments)
                isnothing(seg_local) && continue

                s = zero(Num)
                for j in 2:numberFunctionsOriginal
                    s += b_basis_support.data[seg_local, j, slot]
                end
                b_basis_support.data[seg_local, 1, slot] = simplify_expr(1 - s)
            end

            for seg in right_support
                seg_local = local_segment_index(seg, rangeSegments)
                isnothing(seg_local) && continue

                s = zero(Num)
                for j in 1:numberFunctionsOriginal-1
                    s += b_basis_support.data[seg_local, j, slot]
                end
                b_basis_support.data[seg_local, end, slot] = simplify_expr(1 - s)
            end
        end
    end

    b_full = differentiate(b_basis_full, x, maximumOrder; simplify_expr=simplify_expr)
    b_support = differentiate(b_basis_support, x, maximumOrder; simplify_expr=simplify_expr)

    b = CompactSymbolicFunctions(
        segmentNodes,
        length(selected_positions);
        auxDims=b_support.auxDims,
        variables=variables,
        init=b_support.data[:, selected_positions, :, :],
    )

    return (
        kind = :bspline,
        order = maximumOrder,
        b = b,
        b_full = b_full,
        b_support = b_support,
    )
end

function local_segment_index(seg, rangeSegments)
    seg < first(rangeSegments) && return nothing
    seg > last(rangeSegments) && return nothing
    return seg - first(rangeSegments) + 1
end

function rebuild_clamped_residuals!(b, allNodes, refKnots, numberFunctions, slot, p, simplify_expr)
    left_support = support_segment_indices(refKnots[1], refKnots[p + 2], allNodes)
    for seg in left_support
        s = zero(Num)
        for idx in 2:numberFunctions
            s += b[seg, idx, slot]
        end
        b[seg, 1, slot] = simplify_expr(1 - s)
    end

    right_support = support_segment_indices(refKnots[numberFunctions], refKnots[numberFunctions + p + 1], allNodes)
    for seg in right_support
        s = zero(Num)
        for idx in 1:numberFunctions-1
            s += b[seg, idx, slot]
        end
        b[seg, numberFunctions, slot] = simplify_expr(1 - s)
    end

    return b
end




function extend_knots(coreKnots, maximumOrder; mode=:clamped)
    if mode == :clamped
        return vcat(
            fill(coreKnots[1], maximumOrder),
            coreKnots,
            fill(coreKnots[end], maximumOrder),
        )
    elseif mode == :ghost
        left_step, right_step = edge_steps(coreKnots)
        left_ghosts = [coreKnots[1] - k * left_step for k in maximumOrder:-1:1]
        right_ghosts = [coreKnots[end] + k * right_step for k in 1:maximumOrder]
        return vcat(left_ghosts, coreKnots, right_ghosts)
    else
        throw(ArgumentError("unknown boundary_mode=$mode; use :clamped or :ghost"))
    end
end

function edge_steps(coreKnots)
    if length(coreKnots) == 1
        return (1.0, 1.0)
    end

    return (coreKnots[2] - coreKnots[1], coreKnots[end] - coreKnots[end - 1])
end

function support_segment_indices(leftKnot, rightKnot, allNodes)
    if rightKnot <= leftKnot
        return Int[]
    end

    segments = Int[]
    for seg in 1:length(allNodes)-1
        leftNode = allNodes[seg]
        rightNode = allNodes[seg + 1]
        if max(leftNode, leftKnot) < min(rightNode, rightKnot)
            push!(segments, seg)
        end
    end
    return segments
end


function iszero_denominator(expr)
    try
        return Bool(Symbolics.value(iszero(expr)))
    catch
        return false
    end
end

function rebuild_clamped_residuals!(b, allNodes, refKnots, numberFunctions, slot, p, simplify_expr)
    left_support = support_segment_indices(1, refKnots[p + 2], allNodes)
    for seg in left_support
        s = zero(Num)
        for idx in 2:numberFunctions
            s += b[seg, idx, slot]
        end
        b[seg, 1, slot] = simplify_expr(1 - s)
    end

    rightExtreme = length(allNodes)

    right_support = support_segment_indices(refKnots[numberFunctions-p-1], rightExtreme, allNodes)
    for seg in right_support
        s = zero(Num)
        for idx in 1:numberFunctions-1
            s += b[seg, idx, slot]
        end
        b[seg, numberFunctions, slot] = simplify_expr(1 - s)
    end

    return b
end

"""
    evaluate_bspline_piecewise(b, idx, orderSlot, ξ, allNodes, x, Δx; Δxval=1.0)

Evaluate the piecewise symbolic basis function stored in `b[:, idx, orderSlot]`
at the normalized coordinate `ξ = x / Δx`.
"""
function evaluate_bspline_piecewise(b, idx, orderSlot, ξ, allNodes, x, Δx; Δxval=1.0)
    segment = find_segment(ξ, allNodes)
    expr = b[segment, idx, orderSlot]
    value = Symbolics.substitute(expr, Dict(x => ξ * Δxval, Δx => Δxval))
    return Symbolics.value(value)
end

function find_segment(ξ, allNodes)
    for k in 1:length(allNodes)-1
        if allNodes[k] <= ξ < allNodes[k + 1]
            return k
        end
    end
    return length(allNodes)
end

"""
    plot_bspline_family(result; order=result_order(result), N=10, Δxval=1.0, show_full=false)

Sample and plot the symbolic B-spline family returned by `construct_bspline_family`.
The x-axis is the normalized coordinate `ξ = x / Δx`, with `N` samples per unit span.
"""
function segment_grid(allNodes; N=10)
    ξgrid = Float64[]
    for k in 1:length(allNodes)-1
        a = allNodes[k]
        b = allNodes[k+1]
        localgrid = collect(range(a, b; length=N+1))
        if k > 1
            localgrid = localgrid[2:end]  # avoid duplicating left endpoint
        end
        append!(ξgrid, localgrid)
    end
    return ξgrid
end

function evaluate_bspline_piecewise_deriv(b_deriv, idx, derivSlot, orderSlot, ξ, allNodes, x, Δx; Δxval=1.0)
    segment = find_segment(ξ, allNodes)
    expr = b_deriv[segment, idx, derivSlot, orderSlot]
    value = Symbolics.substitute(expr, Dict(x => ξ * Δxval, Δx => Δxval))
    return Symbolics.value(value)
end


function plotBSpline(csf::CompactSymbolicFunctions;derivOrder=0,order=0,N=10, Δxval=1.0) 
    slots=(derivOrder+1,order+1,)
    ylabel="Derivative order $derivOrder"
    return plotCompactSymbolicFunctions(csf, slots::NTuple; N=N, Δxval=Δxval,ylabel=ylabel)
end
