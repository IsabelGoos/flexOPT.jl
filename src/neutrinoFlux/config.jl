# ============================================================================
# Neutrino flux model parameter structs
# ============================================================================

"""
    NeutrinoFluxParameters

Abstract supertype of the parameter sets of the supported neutrino flux
models. Each model has its own concrete subtype, since each depends on a
different set of parameters.
"""
abstract type NeutrinoFluxParameters end

"""
    HondaParameters

Parameters specific to the Honda atmospheric neutrino flux model tables.
"""
Base.@kwdef struct HondaParameters <: NeutrinoFluxParameters
    location::Symbol            = :gran_sasso
    season::Symbol              = :all_year
    angles_separation::Symbol   = :averaged_ϕ
    solar_activity::Symbol      = :minimum
    mountain_overburden::Symbol = :without
    return_uncertainties::Bool  = false
end

"""
    DaemonfluxParameters

Parameters specific to the Daemonflux model.
"""
Base.@kwdef struct DaemonfluxParameters <: NeutrinoFluxParameters
    flux_mode::Symbol          = :total
    return_uncertainties::Bool = false
    params_df::Union{Nothing, Dict{String, Any}} = nothing
end

"""
    ChengParameters

Parameters specific to the Cheng atmospheric neutrino flux model tables.
"""
Base.@kwdef struct ChengParameters <: NeutrinoFluxParameters
    location::Symbol           = :orca
    season::Symbol             = :all_year
    angles_separation::Symbol  = :averaged_ϕ
    solar_activity::Symbol     = :minimum
    muon_propa_earth::Symbol   = :with 
    return_uncertainties::Bool = false
end

# Associate to each model the corresponding parameter set subtype.
const MODEL_2_PARAMSTYPE = Dict{Symbol, Type{<:NeutrinoFluxParameters}}(
    :Honda      => HondaParameters,
    :Daemonflux => DaemonfluxParameters,
    :Cheng      => ChengParameters,
)

# ============================================================================
# Defaults
# ============================================================================

const DEFAULT_NUFLUX_MODEL  = Ref{Symbol}(:Honda)
const DEFAULT_NUFLUX_PARAMS = Ref{NeutrinoFluxParameters}(HondaParameters())

"""
    set_default_neutrino_flux!(model::Symbol; kwargs...)

Set the default neutrino flux model and its parameters (default=Honda),
analogous to `set_default_planet!`.

The parameter set is built from the defaults of `model`, and the keyword
arguments override the fields you want to change. Fields not given as keywords
take the model's default values, they do **not** keep settings from earlier
calls. Optional keywords update the typed default parameters.

No network access or flux computation occurs here.

# Returns
A `NamedTuple` `(; model, params)` with the new default model and its
validated parameter set.
"""
function set_default_neutrino_flux!(model::Symbol=:Honda; kwargs...)

    # Check if model is legit.
    model_validated = _validate_neutrino_flux_model(model)

    # Define default parameters corresponding to the model.
    params_type = MODEL_2_PARAMSTYPE[model_validated] 
    params      = params_type()

    # Update any parameters and validate the new parameter set.
    params_updated   = _update_parameters(params; kwargs...)
    params_validated = _validate_parameters(params_updated)

    # Commit only after validation succeeded.
    DEFAULT_NUFLUX_MODEL[]  = model_validated
    DEFAULT_NUFLUX_PARAMS[] = params_validated
    return (; model=DEFAULT_NUFLUX_MODEL[], params=DEFAULT_NUFLUX_PARAMS[])
end

"""
    _update_parameters(params::NeutrinoFluxParameters; kwargs...)

Return a copy of the parameter set `params` in which (only) the fields given as
keyword arguments are replaced by their new values.
"""
function _update_parameters(params::ParamsType; kwargs...) where {ParamsType<:NeutrinoFluxParameters}

    # Names of the fields of this parameter type.
    names = fieldnames(ParamsType) # Note that fieldnames gives Symbols
                                   # (except in the case of tuples).

    # Catch typos and keywords that belong to another model.
    bad_keys = [k for k in keys(kwargs) if k ∉ names]
    isempty(bad_keys) || throw(ArgumentError(
        "Unknown $(nameof(ParamsType)): $(repr(bad_keys)). Valid: $(collect(names))."))

    # Convert the struct into a NamedTuple...
    current_params = NamedTuple{names}(ntuple(i -> getfield(params, i), length(names)))

    # Convert the given values to the types the fields expect (e.g. "juno" → :juno).
    normalized_kwargs = (; (k => _normalize_choice(fieldtype(ParamsType, k), v) for (k, v) in pairs(kwargs))...)

    # ... so that it can be merged with the given normalized keyword arguments.
    return ParamsType(; merge(current_params, normalized_kwargs)...)
end

















const LOCATION_MAPPING = Dict{Symbol,String}(
    :kamioka => "kam",
    :gran_sasso => "grn",
    :sudbury => "sno",
    :frejus => "frj",
    :ino => "ino",
    :south_pole => "spl",
    :pythasalmi => "pyh",
    :homestake => "hms",
    :juno => "juno",
)

const SEASON_MAPPING = Dict{Symbol,String}(
    :all_year => "ally",
    :march_may => "0305",
    :june_aug => "0608",
    :sept_nov => "0911",
    :dec_feb => "1202",
)

const ANGLES_MAPPING = Dict{Symbol,String}(
    :variable_ϕ => "20-12",
    :averaged_ϕ => "20-01",
    :averaged_ϕθ => "01-01",
)

const MOUNTAIN_MAPPING = Dict{Symbol,String}(
    :with => "-mtn",
    :without => "",
)

const SOLAR_MAPPING = Dict{Symbol,String}(
    :minimum => "solmin",
    :maximum => "solmax",
)

function _validate_neutrino_flux_model(model::Symbol)
    model in (:Honda, :Daemonflux) ||
        throw(ArgumentError("model must be :Honda or :Daemonflux; got $model"))
    return model
end

function _validate_neutrino_flux_parameters(params::NeutrinoFluxParameters)
    params.flux_mode in (:total, :conventional) ||
        throw(ArgumentError("flux_mode must be :total or :conventional"))
    haskey(LOCATION_MAPPING, params.location_hf) ||
        throw(ArgumentError("unsupported Honda location: $(params.location_hf)"))
    haskey(SEASON_MAPPING, params.season_hf) ||
        throw(ArgumentError("unsupported Honda season: $(params.season_hf)"))
    haskey(ANGLES_MAPPING, params.angles_separation_hf) ||
        throw(ArgumentError(
            "unsupported Honda angular separation: $(params.angles_separation_hf)",
        ))
    haskey(MOUNTAIN_MAPPING, params.mountain_overburden_hf) ||
        throw(ArgumentError(
            "mountain_overburden_hf must be :with or :without",
        ))
    haskey(SOLAR_MAPPING, params.solar_activity_hf) ||
        throw(ArgumentError("solar_activity_hf must be :minimum or :maximum"))
    return params
end

_parameters_namedtuple(params::NeutrinoFluxParameters) = (;
    flux_mode=params.flux_mode,
    return_uncertainties=params.return_uncertainties,
    params_df=params.params_df,
    location_hf=params.location_hf,
    season_hf=params.season_hf,
    angles_separation_hf=params.angles_separation_hf,
    mountain_overburden_hf=params.mountain_overburden_hf,
    solar_activity_hf=params.solar_activity_hf,
)

function configure_neutrino_flux(
    params::NeutrinoFluxParameters=DEFAULT_NUFLUX_PARAMS[];
    kwargs...,
)
    updated = NeutrinoFluxParameters(;
        merge(_parameters_namedtuple(params), NamedTuple(kwargs))...,
    )
    return _validate_neutrino_flux_parameters(updated)
end



const _NUFLUX_PARAMETER_KEYS = Set(fieldnames(NeutrinoFluxParameters))
const _NUFLUX_CONFIGURATION_KEYS = union(_NUFLUX_PARAMETER_KEYS, Set((:model,)))

function _normalize_choice(value, choices, key::Symbol)
    candidate = lowercase(String(Symbol(value)))
    match = findfirst(choice -> lowercase(String(choice)) == candidate, choices)
    isnothing(match) &&
        throw(ArgumentError(
            "invalid $key=$(repr(value)); supported values are $(collect(choices))",
        ))
    return choices[match]
end

function _normalize_boolean(value, key::Symbol)
    value isa Bool && return value
    if value isa AbstractString || value isa Symbol
        normalized = lowercase(String(value))
        normalized == "true" && return true
        normalized == "false" && return false
    end
    throw(ArgumentError("$key must be true or false; got $(repr(value))"))
end

function _normalize_daemonflux_parameters(value)
    isnothing(value) && return nothing
    value isa AbstractDict ||
        throw(ArgumentError("params_df must be nothing or a dictionary"))
    return Dict{String,Any}(string(key) => item for (key, item) in value)
end

function _normalized_configuration_dict(raw::AbstractDict)
    normalized = Dict{Symbol,Any}()
    for (raw_key, value) in raw
        key = Symbol(raw_key)
        key in _NUFLUX_CONFIGURATION_KEYS ||
            throw(ArgumentError(
                "unknown neutrino-flux parameter $(repr(raw_key)); " *
                "supported keys are $(sort!(collect(_NUFLUX_CONFIGURATION_KEYS)))",
            ))
        haskey(normalized, key) &&
            throw(ArgumentError(
                "parameter $key was provided more than once using String/Symbol keys",
            ))
        normalized[key] = value
    end
    return normalized
end

function _parameters_from_dict(
    raw::AbstractDict,
    base::NeutrinoFluxParameters=DEFAULT_NUFLUX_PARAMS[],
)
    updates = _normalized_configuration_dict(raw)
    current = Dict{Symbol,Any}(pairs(_parameters_namedtuple(base)))

    if haskey(updates, :flux_mode)
        current[:flux_mode] = _normalize_choice(
            updates[:flux_mode],
            (:total, :conventional),
            :flux_mode,
        )
    end
    if haskey(updates, :return_uncertainties)
        current[:return_uncertainties] = _normalize_boolean(
            updates[:return_uncertainties],
            :return_uncertainties,
        )
    end
    if haskey(updates, :params_df)
        current[:params_df] = _normalize_daemonflux_parameters(updates[:params_df])
    end

    choice_fields = (
        location_hf=Tuple(keys(LOCATION_MAPPING)),
        season_hf=Tuple(keys(SEASON_MAPPING)),
        angles_separation_hf=Tuple(keys(ANGLES_MAPPING)),
        mountain_overburden_hf=Tuple(keys(MOUNTAIN_MAPPING)),
        solar_activity_hf=Tuple(keys(SOLAR_MAPPING)),
    )
    for (key, choices) in pairs(choice_fields)
        if haskey(updates, key)
            current[key] = _normalize_choice(updates[key], choices, key)
        end
    end

    params = NeutrinoFluxParameters(;
        (key => current[key] for key in fieldnames(NeutrinoFluxParameters))...,
    )
    return _validate_neutrino_flux_parameters(params)
end




function _parameters_dict(
    params::NeutrinoFluxParameters;
    model::Symbol,
    bin_centers_arrays=nothing,
)
    result = Dict{String,Any}(
        "model" => _validate_neutrino_flux_model(model),
        "flux_mode" => params.flux_mode,
        "return_uncertainties" => params.return_uncertainties,
        "params_df" => params.params_df,
        "location_hf" => params.location_hf,
        "season_hf" => params.season_hf,
        "angles_separation_hf" => params.angles_separation_hf,
        "mountain_overburden_hf" => params.mountain_overburden_hf,
        "solar_activity_hf" => params.solar_activity_hf,
    )
    if !isnothing(bin_centers_arrays)
        result["bin_centers_arrays"] = bin_centers_arrays
    end
    return result
end
