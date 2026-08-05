# ===================================================
# Abstract type and model-specific parameter structs
# ===================================================

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
end

"""
    DaemonfluxParameters

Parameters specific to the Daemonflux model with continuous parameters.
"""
Base.@kwdef struct DaemonfluxParameters <: NeutrinoFluxParameters
    flux_mode::Symbol          = :total
    return_uncertainties::Bool = false
    params::Union{Nothing,Dict{String,Any}} = nothing
end

"""
    ChengParameters

Parameters specific to the Cheng atmospheric neutrino flux model tables.
"""
Base.@kwdef struct ChengParameters <: NeutrinoFluxParameters
    location::Symbol            = :gran_sasso
    season::Symbol              = :all_year
    angles_separation::Symbol   = :averaged_ϕ
    solar_activity::Symbol      = :minimum
    muon_propa_earth::Symbol    = :with 

end

# =============================================================
# Default values container (model + model-specific parameters)
# =============================================================

const DEFAULT_NUFLUX_MODEL  = Ref{Symbol}(:Honda)
const DEFAULT_NUFLUX_PARAMS = Ref{NeutrinoFluxParameters}(HondaParameters())

# =========
# Mappings
# =========

const LOCATION_MAPPING_H = Dict{Symbol, String}(
    :kamioka    => "kam",
    :gran_sasso => "grn",
    :sudbury    => "sno",
    :frejus     => "frj",
    :ino        => "ino",
    :south_pole => "spl",
    :pythasalmi => "pyh",
    :homestake  => "hms",
    :juno       => "juno",

)

const LOCATION_MAPPING_C = Dict{Symbol, String}(
    :juno       => ["juno", "JUNO"],
    :superk     => ["kam", "SK"],
    :orca       => ["orca", "KM3NeT-ORCA"],
    :icecube    => ["icecube", "IceCube"],
    :dune       => ["dune", "DUNE"],
    :trident    => ["trident", "TRIDENT"],
    :cjpl       => ["jp", "CJPL"],
)

const SEASON_MAPPING = Dict{Symbol, String}(
    :all_year  => "ally",
    :march_may => "0305",
    :june_aug  => "0608",
    :sept_nov  => "0911",
    :dec_feb   => "1202",
)

const ANGLES_MAPPING = Dict{Symbol, String}(
    :variable_ϕ  => "20-12",
    :averaged_ϕ  => "20-01",
    :averaged_ϕθ => "01-01",
)

const SOLAR_MAPPING = Dict{Symbol, String}(
    :minimum => ["solmin", "solar-min"],
    :maximum => ["solmax", "solar-max"],
)

const MOUNTAIN_MAPPING = Dict{Symbol, String}(
    :with    => "-mtn",
    :without => "",
)

const MUON_PROPA_MAPPING = Dict{Symbol, String}(
    :with    => "with-muon-in-earth",
    :without => "without-muon-in-earth",
)

# =====================================
# Allowed parameter keys in Daemonflux
# =====================================

const DAEMONFLUX_ALLOWED_KEYS = Set{String}([
    "K+_158G", "K+_2P", "K+_31G", "K-_158G", "K-_2P", "K-_31G",
    "n_158G", "n_2P", "p_158G", "p_2P",
    "pi+_158G", "pi+_20T", "pi+_2P", "pi+_31G",
    "pi-_158G", "pi-_20T", "pi-_2P", "pi-_31G",
    "GSF_1", "GSF_2", "GSF_3", "GSF_4", "GSF_5", "GSF_6"
])

# ===================
# Validation methods
# ===================

function _validate_neutrino_flux_model(model::Symbol)
    model in (:Honda, :Daemonflux, :Cheng) ||
        throw(ArgumentError("Model must be :Honda, :Daemonflux, or :Cheng; got $model ."))
    return model
end

function _validate_parameters(p::HondaParameters)
    haskey(LOCATION_MAPPING_H, p.location) ||
        throw(ArgumentError("Unsupported Honda location: $(p.location) ."))
    haskey(SEASON_MAPPING, p.season) ||
        throw(ArgumentError("Unsupported Honda season: $(p.season) ."))
    haskey(ANGLES_MAPPING, p.angles_separation) ||
        throw(ArgumentError("Unsupported Honda angular separation: $(p.angles_separation) ."))
    haskey(SOLAR_MAPPING, p.solar_activity) ||
        throw(ArgumentError("Solar_activity must be :minimum or :maximum ."))
    haskey(MOUNTAIN_MAPPING, p.mountain_overburden) ||
        throw(ArgumentError("Mountain_overburden must be :with or :without ."))
    return p
end

function _validate_parameters(p::DaemonfluxParameters)
    p.flux_mode in (:total, :conventional) ||
        throw(ArgumentError("Flux_mode must be :total or :conventional ."))
    p.return_uncertainties in (false, true) ||
        throw(ArgumentError("Return_uncertainties must be true or false ."))
    if !isnothing(p.params)
        for key in keys(p.params)
            key_str = string(key)
            if !(key_str in DAEMONFLUX_ALLOWED_KEYS)
                throw(ArgumentError("Invalid Daemonflux parameter: $(repr(key)) ."))
            end
        end
    end
    return p
end

function _validate_parameters(p::ChengParameters)
    haskey(LOCATION_MAPPING_C, p.location) ||
        throw(ArgumentError("Unsupported Honda location: $(p.location) ."))
    haskey(SEASON_MAPPING, p.season) ||
        throw(ArgumentError("Unsupported Honda season: $(p.season) ."))
    haskey(ANGLES_MAPPING, p.angles_separation) ||
        throw(ArgumentError("Unsupported Honda angular separation: $(p.angles_separation) ."))
    haskey(SOLAR_MAPPING, p.solar_activity) ||
        throw(ArgumentError("Solar_activity must be :minimum or :maximum ."))
    haskey(MUON_PROPA_MAPPING, p.muon_propa_earth) ||
        throw(ArgumentError("Muon_propa_earth must be :with or :without ."))
    return p
end

# ================================
# Configuration and state setters
# ================================

"""
    set_default_neutrino_flux!(params::NeutrinoFluxParameters, model::Symbol)

Set default parameters to update the default model.
No network access or flux computation occurs here.
"""
function set_default_neutrino_flux!(
    model::Symbol = DEFAULT_NUFLUX_MODEL[],
    params::NeutrinoFluxParameters = DEFAULT_NUFLUX_PARAMS[] 
)
    validated_m = _validate_neutrino_flux_model(model)
    validated_p = _validate_parameters(params)

    @assert validated_m == _model_symbol_from_type(params) "The model :$validated_m and the parameters type $(typeof(params)) do not go together."

    DEFAULT_NUFLUX_MODEL[]  = validated_m
    DEFAULT_NUFLUX_PARAMS[] = validated_p
    return (; model=DEFAULT_NUFLUX_MODEL[], params=DEFAULT_NUFLUX_PARAMS[])
end

_model_symbol_from_type(::HondaParameters)      = :Honda
_model_symbol_from_type(::DaemonfluxParameters) = :Daemonflux
_model_symbol_from_type(::ChengParameters)      = :Cheng

"""
    set_default_neutrino_flux!()

Return the current default model symbol and parameters without making changes.
"""
function set_default_neutrino_flux!()
    @info "No parameters provided. Returning current default configuration."
    return (; model=DEFAULT_NUFLUX_MODEL[], params=DEFAULT_NUFLUX_PARAMS[])
end

"""
    set_default_neutrino_flux!(raw::AbstractDict)

Return the current default model symbol and parameters without making changes.
"""
function set_default_neutrino_flux!(raw::AbstractDict)
    if isempty(raw)
        @info "Received empty dictionary. Returning current default configuration."
        return set_default_neutrino_flux!()
    end
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

"""
    set_default_neutrino_flux!(updates::AbstractDict)

Partially update the current default model and/or any flux parameters. Keys may
be `String` or `Symbol`; symbolic choices may likewise be strings or symbols.
Unspecified values are preserved.

# Examples

    set_default_neutrino_flux!(Dict(:model => :Honda))

    set_default_neutrino_flux!(Dict(
        "location_hf" => "juno",
        "solar_activity_hf" => "maximum",
    ))

    set_default_neutrino_flux!(Dict(
        :model => :Daemonflux,
        :flux_mode => :conventional,
        :return_uncertainties => true,
        :params_df => Dict("K+_158G" => 1.0),
    ))
"""
function set_default_neutrino_flux!(raw::AbstractDict)
    updates = _normalized_configuration_dict(raw)
    model = haskey(updates, :model) ?
            _normalize_choice(
                updates[:model],
                (:Honda, :Daemonflux),
                :model,
            ) :
            DEFAULT_NUFLUX_MODEL[]
    return set_default_neutrino_flux!(
        _parameters_from_dict(raw, DEFAULT_NUFLUX_PARAMS[]);
        model=model,
    )
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
