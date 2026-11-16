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

set_default_neutrino_flux!(model::AbstractString; kwargs...) =
    set_default_neutrino_flux!(Symbol(model); kwargs...)

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

# ============================================================================
# Mappings
# ============================================================================
# These mappings are used to build the file names that are then used to
# retrieve the corresponding models from the web.

# Honda location mapping
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

# Cheng Location mapping
const LOCATION_MAPPING_C = Dict{Symbol, Vector{String}}(
    :cjpl       => ["CJPL", "jp"],
    :dune       => ["DUNE", "dune"],
    :icecube    => ["IceCube", "icecube"],
    :juno       => ["JUNO", "juno"],
    :orca       => ["KM3NeT-ORCA", "orca"],
    :superk     => ["SK", "kam"],
    :trident    => ["TRIDENT", "trident"],
)

# Seasons, used by Honda and Cheng.
const SEASON_MAPPING = Dict{Symbol, String}(
    :all_year  => "ally",
    :march_may => "0305",
    :june_aug  => "0608",
    :sept_nov  => "0911",
    :dec_feb   => "1202",
)

# Angles separation, used by Honda and Cheng.
const ANGLES_MAPPING = Dict{Symbol, String}(
    :variable_ϕ  => "20-12",
    :averaged_ϕ  => "20-01",
    :averaged_ϕθ => "01-01",
)

# The first/second option is used by Honda/Cheng.
const SOLAR_MAPPING = Dict{Symbol, @NamedTuple{honda::String, cheng::String}}(
    :minimum => (honda = "solmin", cheng = "solar-min"),
    :maximum => (honda = "solmax", cheng = "solar-max"),
)

# Mountain overburden, used by Honda only.
const MOUNTAIN_MAPPING = Dict{Symbol, String}(
    :with    => "-mtn",
    :without => "",
)

# Muon propagation, used by Cheng only.
const MUON_PROPA_MAPPING = Dict{Symbol, String}(
    :with    => "with-muon-in-earth",
    :without => "without-muon-in-earth",
)

# ============================================================================
# Allowed parameters
# ============================================================================

const ALLOWED_MODELS = Set{Symbol}([
    :Honda, :Daemonflux, :Cheng,
])

const DAEMONFLUX_ALLOWED_FLUX_MODES = Set{Symbol}([
    :total, :conventional,
])

const DAEMONFLUX_ALLOWED_KEYS = Set{String}([
    "K+_158G", "K+_2P", "K+_31G", "K-_158G", "K-_2P", "K-_31G",
    "n_158G", "n_2P", "p_158G", "p_2P",
    "pi+_158G", "pi+_20T", "pi+_2P", "pi+_31G",
    "pi-_158G", "pi-_20T", "pi-_2P", "pi-_31G",
    "GSF_1", "GSF_2", "GSF_3", "GSF_4", "GSF_5", "GSF_6",
])

const FLUX_KEYS      = ("NuMu_flux",     "Nue_flux",     "AntiNuMu_flux",     "AntiNue_flux")
const FLUX_ERRS_KEYS = ("NuMu_flux_err", "Nue_flux_err", "AntiNuMu_flux_err", "AntiNue_flux_err")

# ============================================================================
# Validation and normalization methods
# ============================================================================

function _validate_neutrino_flux_model(model::Symbol)
    model in ALLOWED_MODELS ||
        throw(ArgumentError(
            "Unsupported neutrino flux model $(repr(model)). " *
            "Valid: $(collect(ALLOWED_MODELS))."
        ))
    return model
end

function _validate_parameters(p::HondaParameters)
    haskey(LOCATION_MAPPING_H, p.location) ||
        throw(ArgumentError(
            "Unsupported Honda location $(repr(p.location)). " *
            "Valid: $(collect(keys(LOCATION_MAPPING_H)))."
        ))
    haskey(SEASON_MAPPING, p.season) ||
        throw(ArgumentError(
            "Unsupported Honda season $(repr(p.season)). " *
            "Valid: $(collect(keys(SEASON_MAPPING)))."
        ))
    haskey(ANGLES_MAPPING, p.angles_separation) ||
        throw(ArgumentError(
            "Unsupported Honda angular separation $(repr(p.angles_separation)). " * 
            "Valid: $(collect(keys(ANGLES_MAPPING)))."
        ))
    haskey(SOLAR_MAPPING, p.solar_activity) ||
        throw(ArgumentError(
            "Unsupported Honda solar activity $(repr(p.solar_activity)). " *
            "Valid: $(collect(keys(SOLAR_MAPPING)))."
        ))
    haskey(MOUNTAIN_MAPPING, p.mountain_overburden) ||
        throw(ArgumentError(
            "Unsupported Honda mountain overburden $(repr(p.mountain_overburden)). " *
            "Valid: $(collect(keys(MOUNTAIN_MAPPING)))."
        ))
    return p
end

function _validate_parameters(p::DaemonfluxParameters)
    p.flux_mode in DAEMONFLUX_ALLOWED_FLUX_MODES ||
        throw(ArgumentError(
            "Unsupported Daemonflux mode $(repr(p.flux_mode)). " *
            "Valid: $(collect(DAEMONFLUX_ALLOWED_FLUX_MODES))."
        ))
    if !isnothing(p.params_df)
       bad_keys = [k for k in keys(p.params_df) if k ∉ DAEMONFLUX_ALLOWED_KEYS]
       isempty(bad_keys) || throw(ArgumentError(
            "Invalid Daemonflux parameter(s): $(repr.(bad_keys)). " *
            "Valid: $(collect(DAEMONFLUX_ALLOWED_KEYS))."
        ))
    end
    return p
end

function _validate_parameters(p::ChengParameters)
    haskey(LOCATION_MAPPING_C, p.location) ||
        throw(ArgumentError(
            "Unsupported Cheng location $(repr(p.location)). " *
            "Valid: $(collect(keys(LOCATION_MAPPING_C)))."
        ))
    haskey(SEASON_MAPPING, p.season) ||
        throw(ArgumentError(
            "Unsupported Cheng season $(repr(p.season)). " *
            "Valid: $(collect(keys(SEASON_MAPPING)))."
        ))
    haskey(ANGLES_MAPPING, p.angles_separation) ||
        throw(ArgumentError(
            "Unsupported Cheng angular separation $(repr(p.angles_separation)). " * 
            "Valid: $(collect(keys(ANGLES_MAPPING)))."
        ))
    haskey(SOLAR_MAPPING, p.solar_activity) ||
        throw(ArgumentError(
            "Unsupported Cheng solar activity $(repr(p.solar_activity)). " *
            "Valid: $(collect(keys(SOLAR_MAPPING)))."
        ))
    haskey(MUON_PROPA_MAPPING, p.muon_propa_earth) ||
        throw(ArgumentError(
            "Unsupported Cheng muon propagation $(repr(p.mountain_overburden)). " * 
            "Valid: $(collect(keys(MUON_PROPA_MAPPING)))."
        ))
    return p
end

"""
    _normalize_choice(T, value)

Convert a user-supplied `value` to the type `T` expected by a parameter field,
so that strings can be used in place of symbols, booleans, or dictionary keys.
"""
function _normalize_choice end

# Symbol fields: accept String or Symbol.
_normalize_choice(::Type{Symbol}, v::Union{AbstractString, Symbol}) =
    Symbol(lowercase(String(v)))

# Bool fields: accept "true" / "false" and similar, or Symbol.
_normalize_choice(::Type{Bool}, v::Union{AbstractString, Symbol}) = parse(Bool, lowercase(String(v)))

# Dictionary values (e.g. params_df): make the keys Strings.
_normalize_choice(::Any, v::AbstractDict) =
    Dict{String, Any}(string(k) => x for (k, x) in v)

# Everything else is left as it is.
_normalize_choice(::Any, v) = v

# ============================================================================
# Base urls of the different neutrino flux models
# ============================================================================

const HONDA_BASE_URL = "http://www-rccn.icrr.u-tokyo.ac.jp/mhonda/public/nflx2014"
const CHENG_BASE_URL = "https://raw.githubusercontent.com/JIECheng2021/atm_nu_flux_data/main"
         
# ============================================================================
# Valid energy and angular ranges of each model's tables
# ============================================================================

const HONDA_E_RANGE    = (0.1, 1.0e4)
const HONDA_COSθ_RANGE = (-1.0, 1.0)

const CHENG_E_RANGE    = (0.01, 1.0e4)
const CHENG_COSθ_RANGE = (-1.0, 1.0)

# Associate to each model the corresponding energy and angular ranges.
const MODEL_2_RANGETYPE = Dict{Symbol, NTuple{2, NTuple{2, Float64}}}(
    :Honda      => (HONDA_E_RANGE, HONDA_COSθ_RANGE),
    :Cheng      => (CHENG_E_RANGE, CHENG_COSθ_RANGE),
)







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




A
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
