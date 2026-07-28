module fieldOutput

export FieldDerivative, FieldExpression, FieldOutputRequest, sample_field_output
export spatial_derivative, time_derivative, differentiate_field

"""
    FieldDerivative(field, coordinates...)

Describe a derivative without assuming a particular number or names of
coordinates. Examples: `FieldDerivative(:ux, :x)` and
`FieldDerivative(:uz, :x, :t)`.
"""
struct FieldDerivative
    field
    coordinates::Tuple{Vararg{Symbol}}
end
FieldDerivative(field, coordinates::Symbol...) =
    FieldDerivative(field, coordinates)

"""
    FieldExpression(name, evaluator)

Define a derived result such as strain, stress, velocity, or acceleration.
`evaluator(fields, coordinates)` must return an array whose axes follow the
coordinate `NamedTuple`. This keeps the output layer independent from the
field symbols and dimensionality selected by `famousEquations`.
"""
struct FieldExpression{N,F}
    name::N
    evaluator::F
end

"""
    FieldOutputRequest(expression; selection=NamedTuple(), label=nothing)

`expression` is a field symbol, `FieldDerivative`, or a function accepting
`(fields, coordinates)`. Each selection value may be `:all`, a scalar nearest
coordinate, a range, or `(minimum, maximum, sampling)`.
"""
Base.@kwdef struct FieldOutputRequest
    expression
    selection::NamedTuple = NamedTuple()
    label::Union{Nothing,Symbol,String} = nothing
    interpolation::Symbol = :linear
end
FieldOutputRequest(expression; selection=NamedTuple(), label=nothing,
                   interpolation=:linear) =
    FieldOutputRequest(expression, selection, label, interpolation)

function _differentiate(values, axis::Int, coordinates)
    result = similar(values)
    n = size(values, axis)
    n >= 2 || throw(ArgumentError("cannot differentiate a singleton axis"))
    select(index) = ntuple(d -> d == axis ? index : Colon(), ndims(values))
    result[select(1)...] .=
        (values[select(2)...] .- values[select(1)...]) /
        (coordinates[2] - coordinates[1])
    for i in 2:n-1
        result[select(i)...] .=
            (values[select(i+1)...] .- values[select(i-1)...]) /
            (coordinates[i+1] - coordinates[i-1])
    end
    result[select(n)...] .=
        (values[select(n)...] .- values[select(n-1)...]) /
        (coordinates[n] - coordinates[n-1])
    result
end

_fielddata(fields::NamedTuple, key::Symbol) = getproperty(fields, key)
_fielddata(fields::AbstractDict, key) = fields[key]
_fielddata(fields, key) = getindex(fields, key)

function differentiate_field(values, coordinates::NamedTuple, name::Symbol)
    axis = findfirst(==(name), propertynames(coordinates))
    isnothing(axis) && throw(ArgumentError("coordinate $name is unavailable"))
    _differentiate(values, axis, getproperty(coordinates, name))
end

function _evaluate(expression, fields, coordinates)
    expression isa FieldExpression &&
        return expression.evaluator(fields, coordinates)
    expression isa Function && return expression(fields, coordinates)
    !(expression isa FieldDerivative) && return _fielddata(fields, expression)
    expression isa FieldDerivative ||
        throw(ArgumentError("expression must be a field, FieldDerivative or function"))
    values = _fielddata(fields, expression.field)
    for name in expression.coordinates
        values = differentiate_field(values, coordinates, name)
    end
    values
end

function _targets(axis, specification)
    source = Float64.(collect(axis))
    specification === :all && return source
    specification isa Real && return [Float64(specification)]
    if specification isa Tuple && length(specification) == 3
        low, high, sampling = specification
        sampling > 0 || throw(ArgumentError("sampling must be positive"))
        return Float64.(collect(low:sampling:high))
    end
    specification isa AbstractRange || specification isa AbstractVector ||
        throw(ArgumentError("invalid coordinate selection $specification"))
    Float64.(collect(specification))
end

function _interpolate_axis(values, source, targets, dimension)
    source_values = Float64.(collect(source))
    all(diff(source_values) .> 0) ||
        throw(ArgumentError("interpolation coordinates must be increasing"))
    all(target -> first(source_values) <= target <= last(source_values), targets) ||
        throw(BoundsError(source_values, targets))
    output_shape = collect(size(values))
    output_shape[dimension] = length(targets)
    output = similar(values, Tuple(output_shape))
    source_selector = Any[Colon() for _ in 1:ndims(values)]
    target_selector = Any[Colon() for _ in 1:ndims(values)]
    for (target_index, target) in pairs(targets)
        upper = searchsortedfirst(source_values, target)
        if upper == 1
            lower = upper = 1
            weight = 0.0
        elseif upper > length(source_values)
            lower = upper = length(source_values)
            weight = 0.0
        elseif source_values[upper] == target
            lower = upper
            weight = 0.0
        else
            lower = upper - 1
            weight = (target - source_values[lower]) /
                     (source_values[upper] - source_values[lower])
        end
        target_selector[dimension] = target_index
        source_selector[dimension] = lower
        if lower == upper
            output[target_selector...] .= values[source_selector...]
        else
            lower_values = values[source_selector...]
            source_selector[dimension] = upper
            output[target_selector...] .=
                (1 - weight) .* lower_values .+ weight .* values[source_selector...]
        end
    end
    output
end

function _nearest_axis(values, source, targets, dimension)
    source_values = Float64.(collect(source))
    indices = [argmin(abs.(source_values .- target)) for target in targets]
    selectors = ntuple(d -> d == dimension ? indices : Colon(), ndims(values))
    values[selectors...]
end

"""
    sample_field_output(fields, coordinates, requests)

Evaluate and subset any number of field/derivative requests. Field arrays must
have axes in exactly `propertynames(coordinates)` order; no spatial dimension
is hard-coded.
"""
function sample_field_output(fields, coordinates::NamedTuple, requests)
    names = propertynames(coordinates)
    map(requests) do request
        values = _evaluate(request.expression, fields, coordinates)
        ndims(values) == length(names) ||
            throw(DimensionMismatch("field dimensions and coordinate axes differ"))
        target_axes = ntuple(length(names)) do dimension
            name = names[dimension]
            specification = hasproperty(request.selection, name) ?
                            getproperty(request.selection, name) : :all
            _targets(getproperty(coordinates, name), specification)
        end
        sampled = values
        current_axes = coordinates
        request.interpolation in (:linear, :nearest) ||
            throw(ArgumentError("interpolation must be :linear or :nearest"))
        for dimension in eachindex(names)
            sampler = request.interpolation === :linear ?
                      _interpolate_axis : _nearest_axis
            sampled = sampler(sampled,
                              getproperty(current_axes, names[dimension]),
                              target_axes[dimension], dimension)
            current_axes = merge(
                current_axes,
                NamedTuple{(names[dimension],)}((target_axes[dimension],)),
            )
        end
        sampled_coordinates = NamedTuple{names}(target_axes)
        default_label = request.expression isa FieldExpression ?
                        request.expression.name : string(request.expression)
        label = isnothing(request.label) ? default_label : request.label
        (; label, values=sampled, coordinates=sampled_coordinates,
           expression=request.expression, selection=request.selection)
    end
end

spatial_derivative(field, coordinate::Symbol) =
    FieldDerivative(field, coordinate)
time_derivative(field, time_coordinate::Symbol=:t) =
    FieldDerivative(field, time_coordinate)

end
