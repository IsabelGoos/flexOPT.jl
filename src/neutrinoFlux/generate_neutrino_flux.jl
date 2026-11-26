# ============================================================================
# Main function
# ============================================================================

"""
    produce_neutrino_flux(bin_centers_arrays::Tuple)
    produce_neutrino_flux(user_params::AbstractDict{String})

Produce the atmospheric neutrino flux on a given binning, using either the
default neutrino flux model and parameters or substituting with `user_params`.

# Arguments
- `bin_centers_arrays::Tuple`: binning `(energies, cosθgrid)`. The flux is
  computed with the current default model and parameters. The bin centers for 
  cosθ have to be given in increasing order.
- `user_params::AbstractDict{String}`: dictionary with the following keys:
  - `bin_centers_arrays` (required): binning `(energies, cosθgrid)`.
  - `model` (optional): replaces the default model.
  - any other key (optional): a parameter of the chosen model that
    replaces its default value.

# Returns 
a 6-tuple (or 10-tuple, if `return_uncertainties::Bool` is true) with the following matrices in the shape (nθ, nE):
- `NuMu_flux` 
- `Nue_flux`
- `AntiNuMu_flux` 
- `AntiNue_flux`
- `Energies`
- `cosθs`
(followed by the corresponding uncertainties if `return_uncertainties::Bool` is true, in the same order).
"""
produce_neutrino_flux(bin_centers_arrays::NamedTuple{names, <:NTuple{2, AbstractVector{<:Real}}}) where {names} =
    produce_neutrino_flux(Dict{String, Any}("bin_centers_arrays" => bin_centers_arrays))

function produce_neutrino_flux(user_params::Dict{String, Any})

    # Make sure that user_params contains the bin centers.
    haskey(user_params, "bin_centers_arrays") ||
        throw(ArgumentError("The key \"bin_centers_arrays\" is required: (energies, cosθgrid)."))

    # Pick the model and parameter set.
    if haskey(user_params, "model")
        model  = _validate_neutrino_flux_model(Symbol(uppercasefirst(lowercase(string(user_params["model"])))))
        params = MODEL_2_PARAMSTYPE[model]()      
    else
        model  = DEFAULT_NUFLUX_MODEL[]
        params = DEFAULT_NUFLUX_PARAMS[]          
    end

    # Update parameters. Any problem will be raised by _update_parameters.
    updates = Dict{Symbol, Any}(
        Symbol(k) => v for (k, v) in user_params if k ∉ ("model", "bin_centers_arrays")
    )  
    params_validated = _validate_parameters(_update_parameters(params; updates...))

    # Dictionary expected by the _produce_*_neutrino_flux functions.
    backend_params = Dict{String, Any}(
        String(f) => getfield(params_validated, f) for f in fieldnames(typeof(params_validated))
    )
    backend_params["model"]              = model
    backend_params["bin_centers_arrays"] = _validate_flux_binning(user_params["bin_centers_arrays"])

    return_uncertainties                = backend_params["return_uncertainties"]
    angles_separation                   = backend_params["angles_separation"]

    for (key, value) in backend_params
        println(key, " => ", value)
    end

    # Produce and return atmospheric neutrino flux.
    if model === :Honda
        nuflux = _produce_honda_neutrino_flux(backend_params)
    elseif model === :Daemonflux
        nuflux = _produce_daemonflux_neutrino_flux(backend_params)
    else  # :Cheng
        nuflux = _produce_cheng_neutrino_flux(backend_params)
    end

    # Interpolate if the models are Honda or Cheng.
    if model in (:Honda, :Cheng)  
        E_orig    = nuflux["Energies"]
        cosθ_orig = nuflux["cosθgrid"]
        bins_new  = user_params["bin_centers_arrays"]
        if angles_separation === :averaged_ϕ
            nuflux_new = Dict{String, Any}( 
                key => _get_interpolated_flux(nuflux[key], E_orig, cosθ_orig, model, bins_new)
                for key in FLUX_KEYS
            )
        end
        if angles_separation === :variable_ϕ
            ϕ_orig = nuflux["ϕs"]
            nuflux_new = Dict{String, Any}( 
                key => _get_interpolated_flux(nuflux[key], E_orig, cosθ_orig, ϕ_orig, model, bins_new)
                for key in FLUX_KEYS
            )
            nuflux_new["ϕs"]    = collect(bins_new[3])
        end
        nuflux_new["Energies"] = collect(bins_new[1])
        nuflux_new["cosθgrid"] = collect(bins_new[2])

        # Now that interpolation has been done, uncertainties can be included.
        if return_uncertainties

            if model === :Honda
                rel_err = map(E -> E < 10.0 ? 0.07 : (E <= 100.0 ? 0.14 : 0.25), nuflux_new["Energies"])
                for (key, key_err) in zip(FLUX_KEYS, FLUX_ERRS_KEYS)
                    nuflux_new[key_err] = nuflux_new[key] .* rel_err
                end
            end

            if model === :Cheng
                E_new    = nuflux_new["Energies"]
                cosθ_new = nuflux_new["cosθgrid"]

                E_pts = [0.35, 0.6, 1.0, 2.0, 5.0, 10.0, 30.0, 50.0, 100.0]
                #                           0.35   0.6    1      2      5     10     30     50    100
                errs_vert = Dict(
                    "NuMu_flux_err"     => [0.175, 0.095, 0.058, 0.052, 0.072, 0.090, 0.115, 0.140, 0.180],
                    "Nue_flux_err"      => [0.225, 0.120, 0.075, 0.050, 0.035, 0.042, 0.075, 0.115, 0.225],
                    "AntiNuMu_flux_err" => [0.170, 0.090, 0.060, 0.043, 0.038, 0.047, 0.065, 0.075, 0.100],
                    "AntiNue_flux_err"  => [0.225, 0.120, 0.080, 0.055, 0.042, 0.050, 0.085, 0.125, 0.275],
                )
                errs_hor = Dict(
                    "Nue_flux_err"      => [0.280, 0.130, 0.072, 0.045, 0.030, 0.022, 0.020, 0.020, 0.020],
                    "NuMu_flux_err"     => [0.200, 0.100, 0.057, 0.048, 0.047, 0.048, 0.054, 0.058, 0.066],
                    "AntiNuMu_flux_err" => [0.200, 0.100, 0.058, 0.038, 0.030, 0.027, 0.030, 0.033, 0.037],
                    "AntiNue_flux_err"  => [0.270, 0.130, 0.072, 0.048, 0.032, 0.022, 0.020, 0.020, 0.020],
                )
                w = abs.(cosθ_new)'

                for key in FLUX_KEYS
                    itp_v = extrapolate(interpolate(log10.(E_pts), errs_vert[key * "_err"], FritschCarlsonMonotonicInterpolation()), Flat())
                    itp_h = extrapolate(interpolate(log10.(E_pts), errs_hor[key * "_err"],  FritschCarlsonMonotonicInterpolation()), Flat())

                    errs_v = itp_v.(log10.(E_new))    
                    errs_h = itp_h.(log10.(E_new))

                    # Linear interpolation wrt angle.
                    rel_err = errs_h .+ (errs_v .- errs_h) .* w    

                    nuflux_new[key * "_err"] = nuflux_new[key] .* rel_err
                end

            end

        end

        nuflux = nuflux_new
    end

    return nuflux
end

"""
    produce_or_load_neutrino_flux(bin_centers_arrays;
        directory="neutrinoFlux", prefix="flux")

Load a previously generated flux table, or generate it and cache it with
`myProduceOrLoad`.

# Arguments
- `bin_centers_arrays::NamedTuple`: two real-valued vectors named `energies`
  and `cosθgrid`.

# Keyword arguments
- `directory::AbstractString`: folder holding the cache files.
- `prefix::AbstractString`: prefix of the cache file names.

The default model and its parameters are used. The cache key contains the complete binning, the
model and the typed parameter set, so changing any of them triggers a new
computation. The returned value is the same `NeutrinoFluxTable` named tuple
as [`produce_neutrino_flux`](@ref).
"""
function produce_or_load_neutrino_flux(
    bin_centers_arrays::NamedTuple{names, <:NTuple{2, AbstractVector{<:Real}}};
    directory::String="neutrinoFlux",
    prefix::String="flux",
) where {names}

    # Read the globally selected model and its matching parameter type.
    model       = DEFAULT_NUFLUX_MODEL[]
    params_type = MODEL_2_PARAMSTYPE[model] 
    params      = params_type()   
    # Validate the grids.
    validated_binning = _validate_flux_binning(bin_centers_arrays)

    # Everything that determines the flux goes into the cache key.
    config_dict = Dict{String, Any}(
        "model"    => model,
        "params"   => params,
        "energies" => validated_binning.energies,
        "cosθgrid" => validated_binning.cosθgrid,
    )
 
    # Load from disk if an entry with this config_dict exists; 
    # otherwise compute it with `_cached_neutrino_flux` and store it.
    cached = myProduceOrLoad(
        _cached_neutrino_flux,
        config_dict,
        directory,
        prefix,
    )

    # Guard against stale or corrupted cache files.
    haskey(cached, "table") ||
        error("invalid neutrino-flux cache entry: missing \"table\"")
    return cached["table"]
end

"""
    _cached_neutrino_flux(config::AbstractDict) 

Compute a neutrino flux table from a cache-key dictionary. This is the function
that `myProduceOrLoad` calls on a cache miss, so its input is exactly the cache
key built in [`produce_or_load_neutrino_flux`](@ref).

# Arguments
- `config::AbstractDict`: must contain the keys
  - `"energies"`: energy bin centers,
  - `"cosθgrid"`: cosθ bin centers,
  - `"model"`: flux model (e.g. `:Honda`),
  - `"params"`: parameter set of that model.

# Returns
A `Dict{String, Any}` with a single entry `"table"` holding the
`NeutrinoFluxTable` named tuple from [`produce_neutrino_flux`](@ref).
"""
function _cached_neutrino_flux(config_dict::AbstractDict)

    # Rebuild the binning named tuple from the cache key.
    binning = (
        energies = config_dict["energies"],
        cosθgrid = config_dict["cosθgrid"],
    )

    # Generate the flux table with the model and parameters stored in the key,
    # not with the current global defaults, so the result always matches the key.
    table = produce_neutrino_flux(binning)
    
    # DrWatson's metadata/tagging layer requires an AbstractDict. Keep that
    # implementation detail inside the cache and expose only the named tuple.
    return Dict{String, Any}("table" => table)
end

# ============================================================================
# Validation methods
# ============================================================================

"""
    _validate_flux_binning(bin_centers_arrays)

Validate the binning of a neutrino flux table and return clean grids.

# Arguments
- `bin_centers_arrays::NamedTuple`: must contain exactly two real-valued vectors,
  named `energies` and `cosθgrid`.
"""
function _validate_flux_binning(
    bin_centers_arrays::NamedTuple{names, <:NTuple{2, AbstractVector{<:Real}}},
    model::Symbol) where {names}

    # Check the names.
    hasproperty(bin_centers_arrays, :energies) ||
        throw(ArgumentError("`bin_centers_arrays` must contain `energies`."))
    hasproperty(bin_centers_arrays, :cosθgrid) ||
        throw(ArgumentError("`bin_centers_arrays` must contain `cosθgrid`."))

    # Check that the arrays are not empty.
    energies = Float64.(collect(bin_centers_arrays.energies))
    cosθgrid = Float64.(collect(bin_centers_arrays.cosθgrid))
    isempty(energies) && throw(ArgumentError("`energies` cannot be empty."))
    isempty(cosθgrid) && throw(ArgumentError("`cosθgrid` cannot be empty."))
    
    # Check that energies and cosθgrid are strictly increasing.
    all(diff(energies) .> 0) || throw(ArgumentError("`energies` must be strictly increasing."))
    all(diff(cosθgrid) .> 0) || throw(ArgumentError("`cosθgrid` must be strictly increasing."))

    # Check if the new bin centers are within the allowed range.
    E_range, cosθ_range = MODEL_2_RANGETYPE[model]
    _check_bins_in_range(energies, cosθgrid, E_range, cosθ_range)
    
    return (; energies, cosθgrid)
end

"""
    _check_bins_in_range(E_new, cosθ_new, E_range, cosθ_range; rtol=1e-8)

Check that a requested binning lies inside the allowed range, 
so that interpolation never needs to extrapolate.

# Arguments
- `E_new`: requested energies, in GeV.
- `cosθ_new`: requested cosines of the zenith angle.
- `E_range`: allowed energy range.
- `cosθ_range`: allowed angular range.

# Keywords
- `rtol::Real = 1e-8`: tolerance.
"""
function _check_bins_in_range(E_new, cosθ_new, E_range, cosθ_range; rtol=1e-8)

    # Accepted energy and angular ranges.
    E_lo,  E_hi  = E_range
    ct_lo, ct_hi = cosθ_range

    # Proposed energy and angular ranges.
    Emin,  Emax  = extrema(E_new)
    ctmin, ctmax = extrema(cosθ_new)

    # Check, with some tolerance, if the proposed ranges are within the accepted ones.
    (Emin >= E_lo * (1 - rtol) && Emax <= E_hi * (1 + rtol)) ||
        throw(ArgumentError("Energies [$Emin, $Emax] GeV are outside the allowed range [$E_lo, $E_hi] GeV."))
    (ctmin >= ct_lo - rtol && ctmax <= ct_hi + rtol) ||
        throw(ArgumentError("cosθ values [$ctmin, $ctmax] are outside the allowed range [$ct_lo, $ct_hi]."))
    return nothing
end

# ============================================================================
# Individual functions for each neutrino flux model
# ============================================================================

# ============================================================================
# Honda model
# ============================================================================

"""
    _produce_honda_neutrino_flux(nuflux_params::Dict{String, Any})

Produce the atmospheric neutrino flux for all relevant neutrino types using the Honda model.
"""
function _produce_honda_neutrino_flux(nuflux_params::Dict{String, Any})

    # Read parameters and build url from where to fetch the data.
    url_honda  = _build_flux_url(nuflux_params)
    angles_sep = nuflux_params["angles_separation"]
    return_uncertainties = nuflux_params["return_uncertainties"]
    
    # Download.
    decompressed_stream = download_nuflux(url_honda)
    
    # Read the text grid into a matrix, omitting any lines with text.
    lines = readlines(decompressed_stream)
    close(decompressed_stream)
    numeric_lines = filter(line -> occursin(r"^\s*[0-9.-]", line), lines)
    clean_text_block = join(numeric_lines, "\n")
    matrix_2d = readdlm(IOBuffer(clean_text_block), Float64)
    
    # Reshape, based on angular separation mode.

    # :averaged_ϕ-grid layout: 5 columns x 101 energies x 20 cosθ steps
    if angles_sep === :averaged_ϕ

        n_E    = 101
        n_cosθ = 20
        flux_reshaped = reshape(matrix_2d, n_E, n_cosθ, 5)

        # Extract matrices.
        results = Dict{String, Any}(
            "NuMu_flux"     => collect(flux_reshaped[:, :, 2]),
            "Nue_flux"      => collect(flux_reshaped[:, :, 4]),
            "AntiNuMu_flux" => collect(flux_reshaped[:, :, 3]),
            "AntiNue_flux"  => collect(flux_reshaped[:, :, 5]),
            "Energies"      => flux_reshaped[:, 1, 1],
            "cosθs"         => range(-1.0, .9, length=n_cosθ).+0.05
        )

    elseif angles_sep === :variable_ϕ
        # TODO: Implement the reshaping logic for the azimuth-dependent grid
        error("angles_separation_hf mode :variable_ϕ is not yet implemented.")

    elseif angles_sep === :averaged_ϕθ
        # TODO: Implement the reshaping logic for the fully integrated/averaged grid
        error("angles_separation_hf mode :averaged_ϕθ is not yet implemented.")

    else
        error("Invalid angles_separation_hf option: $angles_sep. Use :variable_ϕ, :averaged_ϕ, or :averaged_ϕθ.")
    end

    return results
end

# ============================================================================
# Daemonflux model
# ============================================================================

"""
    _produce_daemonflux_neutrino_flux(nuflux_params::Dict{String, Any})

Produce the atmospheric neutrino flux for all relevant neutrino types using the Daemonflux model.
"""
function _produce_daemonflux_neutrino_flux(nuflux_params::Dict{String, Any})
    pythoncall = _require_pythoncall()
    return Base.invokelatest(
        _produce_daemonflux_neutrino_flux,
        nuflux_params,
        pythoncall
    )
end

function _produce_daemonflux_neutrino_flux(nuflux_params::Dict{String, Any}, pythoncall::Module)

    # Read some inputs.
    E_new, cosθ_new      = nuflux_params["bin_centers_arrays"]
    flux_mode            = nuflux_params["flux_mode"]
    return_uncertainties = nuflux_params["return_uncertainties"]

    # Import and instantiate Python daemonflux class.
    nuflux = _get_daemonflux_instance(pythoncall)

    # Convert to degrees, as needed by daemonflux.
    θ_new = string.(rad2deg.(acos.(cosθ_new)))

    # Pre-allocate output arrays.
    nE = length(E_new)
    nθ = length(θ_new)
    NuMu_flux     = zeros(Float64, nθ, nE)
    Nue_flux      = zeros(Float64, nθ, nE)
    AntiNuMu_flux = zeros(Float64, nθ, nE)
    AntiNue_flux  = zeros(Float64, nθ, nE)
    if return_uncertainties
        NuMu_flux_err     = zeros(Float64, nθ, nE)
        Nue_flux_err      = zeros(Float64, nθ, nE)
        AntiNuMu_flux_err = zeros(Float64, nθ, nE)
        AntiNue_flux_err  = zeros(Float64, nθ, nE)
    end

    # Define key names depending on the chosen flux_mode. 
    prefix = flux_mode === :total ? "total_" : ""
    keys = (
        numu    = prefix * "numu",
        nue     = prefix * "nue",
        anumu   = prefix * "antinumu",
        anue    = prefix * "antinue" 
    )

    # Set daemonflux parameters.
    raw_params_df = get(nuflux_params, "params_df", nothing)
    params = (isnothing(raw_params_df) || isempty(raw_params_df)) ? 
             _set_dflux_params(nuflux) : 
             _set_dflux_params(nuflux, raw_params_df)

    # Energy scaling.
    E_scaling = 1.0e4 ./ (E_new .^ 3)

    # Fetch raw daemonflux values.
    @views for (i, θ) in enumerate(θ_new)
        # flux values
        NuMu_flux[i, :]     .= pyconvert(Vector{Float64}, nuflux.flux(E_new, θ, keys.numu,  params=params)) .* E_scaling
        Nue_flux[i, :]      .= pyconvert(Vector{Float64}, nuflux.flux(E_new, θ, keys.nue,   params=params)) .* E_scaling
        AntiNuMu_flux[i, :] .= pyconvert(Vector{Float64}, nuflux.flux(E_new, θ, keys.anumu, params=params)) .* E_scaling
        AntiNue_flux[i, :]  .= pyconvert(Vector{Float64}, nuflux.flux(E_new, θ, keys.anue,  params=params)) .* E_scaling
        # flux uncertainties
        # E_scaling relevant here too? TODO: check
        if return_uncertainties
            NuMu_flux_err[i, :]     .= pyconvert(Vector{Float64}, nuflux.error(E_new, θ, keys.numu))
            Nue_flux_err[i, :]      .= pyconvert(Vector{Float64}, nuflux.error(E_new, θ, keys.nue))
            AntiNuMu_flux_err[i, :] .= pyconvert(Vector{Float64}, nuflux.error(E_new, θ, keys.anumu))
            AntiNue_flux_err[i, :]  .= pyconvert(Vector{Float64}, nuflux.error(E_new, θ, keys.anue))
        end
    end

    # Prepare dictionary for output.
    results = Dict{String, Any}(
        "NuMu_flux"     => NuMu_flux,
        "Nue_flux"      => Nue_flux,
        "AntiNuMu_flux" => AntiNuMu_flux,
        "AntiNue_flux"  => AntiNue_flux
    )
    if return_uncertainties
        results["NuMu_flux_err"]     = NuMu_flux_err
        results["Nue_flux_err"]      = Nue_flux_err
        results["AntiNuMu_flux_err"] = AntiNuMu_flux_err
        results["AntiNue_flux_err"]  = AntiNue_flux_err
    end

    results["Energies"] = E_new
    results["cosθs"]    = cosθ_new

    return results

end

"""
    _get_daemonflux_instance(pythoncall::Module; retries=3, retry_delay=2)

Create and return the Python `daemonflux.Flux` object for IceCube, using the `pythoncall` module. 
"""
function _get_daemonflux_instance(
    pythoncall::Module;
    retries::Integer=3,
    retry_delay::Real=2,
)
    retries >= 1 || throw(ArgumentError("`retries` must be at least 1."))
    daemonflux = pythoncall.pyimport("daemonflux")
    for attempt in 1:retries
        try
            return daemonflux.Flux(
                location="IceCube",
                use_calibration=true,
                debug=1,
            )
        catch error_pc
            if attempt == retries
                message = sprint(showerror, error_pc)
                throw(ErrorException(
                    "Daemonflux could not initialize after $retries attempts. " *
                    "Its official GitHub release may be temporarily unreachable. " *
                    "Check access to github repository, " *
                    "and retry the computation.\n" *
                    "Last Python error:\n$message",
                ))
            end
            sleep(retry_delay * 2.0^(attempt - 1))
        end
    end
    error("Daemonflux unreachable.")
end

"""
    _set_dflux_params(flux_obj)

Initialize a dictionary of Daemonflux parameters set to their default values.

# Arguments
- `flux_obj`: A Python object representing a Daemonflux instance.
"""
function _set_dflux_params(flux_obj)
    pythoncall = _require_pythoncall()
    return Base.invokelatest(_set_dflux_params, flux_obj, pythoncall)
end

function _set_dflux_params(flux_obj, pythoncall::Module)

    # Fetch all the Daemonflux parameter names. 
    params = pythoncall.pyconvert(
        Vector{String},
        flux_obj.params.known_parameters,
    )

    # Set all parameters to 0 
    # (meaning that their difference with respect to the default parameters is zero).
    params_dict = Dict(p => 0.0 for p in params) 
    
    return params_dict

end

"""
    _set_dflux_params(flux_obj, params_dict)

Initialize a dictionary of Daemonflux parameters to those values given in `params_dict`.
Those parameters not mentioned there are set to default values.

# Arguments
- `flux_obj`: A Python object representing a Daemonflux instance.
- `params_dict`: A dictionary with those pairs of parameter names and values which should be modified.
The possible parameter names are K+_158G, K+_2P, K+_31G, K-_158G, K-_2P, K-_31G, n_158G, n_2P, p_158G, p_2P, 
pi+_158G, pi+_20T, pi+_2P, pi+_31G, pi-_158G, pi-_20T, pi-_2P, pi-_31G, GSF_1, GSF_2, GSF_3, GSF_4, GSF_5, GSF_6.
The get_dflux_param_value function helps find specific values.

# Returns
- `Dict{String, Float64}` 
"""
function _set_dflux_params(flux_obj, params_dict::Union{Nothing, Dict{String, Any}})
    pythoncall = _require_pythoncall()
    return Base.invokelatest(
        _set_dflux_params,
        flux_obj,
        params_dict,
        pythoncall
    )
end

function _set_dflux_params(flux_obj, params_dict::Union{Nothing, Dict{String, Any}}, pythoncall::Module)

    # Set default Daemonflux parameters.
    full_params = _set_dflux_params(flux_obj, pythoncall)

    # If params_dict is nothing, just return the defaults without trying to merge.
    if isnothing(params_dict)
        return full_params
    end

    # Otherwise, update with pairs from params_dict.
    merge!(full_params, params_dict)

    return full_params

end

"""
    get_dflux_param_value(flux_obj, param_name, signed_sigma; max_iterations=10000, max_param=10.0, tol=0.01)

Find the value of the parameter `param_name` which is `signed_sigma` away from the default value.

# Arguments
- `flux_obj`: A Python object representing a Daemonflux instance.
- `param_name`: The name of the parameter.
- `signed_sigma`: Amount of sigmas you want to be away from the default value, up or down.
- `max_iterations`: The maximum number of iterations in the bisection method.
- `max_param`: The maximum search limit for the absolute value of the parameter value.
- `tol`: The relative tolerance on the target distance from the default value.
"""
function get_dflux_param_value(
    flux_obj,
    param_name::String,
    signed_sigma::Real;
    max_iterations::Int=10000,
    max_param::Real=10.0,
    tol::Real=0.01,
)
    pythoncall = _require_pythoncall()
    return Base.invokelatest(
        _get_dflux_param_value,
        flux_obj,
        param_name,
        signed_sigma;
        max_iterations=max_iterations,
        max_param=max_param,
        tol=tol,
        pythoncall
    )
end

function get_dflux_param_value(
    flux_obj,
    param_name::String,
    signed_sigma::Real;
    max_iterations::Int=10000,
    max_param::Real=10.0,
    tol::Real=0.01,
    pythoncall::Module
)

    # Determine the sign of the parameter based on if it should be above or below the default value.
    p_sign = sign(signed_sigma)
    # flux_obj.chi2 is symmetric and positive, the search is performed on the positive axis.
    sigma = abs(signed_sigma)

    # Define the starting search bounds.
    p_low  = 0.
    p_high = max_param

    # Evaluate the function at the upper boundary to ensure a root exists.
    chi2_high = pythoncall.pyconvert(
        Float64,
        flux_obj.chi2(params=Dict(param_name => p_high)),
    )
    diff_high = sqrt(chi2_high) - sigma
    # If the target is out of bounds, fail immediately.
    if diff_high < 0.0
        error("Target sigma ($signed_sigma) is outside the search bounds. At p = $max_param, chi2 is only $chi2_high. Increase `max_param`.")
    end

    # Run bisection search.
    for i in 1:max_iterations

        p_mid = (p_low + p_high) / 2
        chi2_mid = pythoncall.pyconvert(
            Float64,
            flux_obj.chi2(params=Dict(param_name => p_mid)),
        )
        sigma_diff = sqrt(chi2_mid) - sigma

        # Check if we are within the relative tolerance threshold
        if abs(sigma_diff) <= tol * sigma
            return p_sign * p_mid  
        end

        # Adjust search bounds
        if sigma_diff < 0.0 
            p_low = p_mid
        else
            p_high = p_mid
        end

    end

    error("Bisection method failed to converge within $max_iterations iterations. Check your parameters.")

end

# ============================================================================
# Cheng model
# ============================================================================

"""
    _produce_cheng_neutrino_flux(nuflux_params::Dict{String, Any})

Produce the atmospheric neutrino flux for all relevant neutrino types using the Cheng model.
"""
function _produce_cheng_neutrino_flux(nuflux_params::Dict{String, Any})

    # Read parameters and build url from where to fetch the data.
    url_cheng  = _build_flux_url(nuflux_params)
    angles_sep = nuflux_params["angles_separation"]

    # Download.
    decompressed_stream = download_nuflux(url_cheng)

    # Read the text grid into a matrix, omitting any lines with text.
    lines = readlines(decompressed_stream)
    close(decompressed_stream)
    numeric_lines = filter(line -> occursin(r"^\s*[0-9.-]", line), lines)
    clean_text_block = join(numeric_lines, "\n")
    matrix_2d = readdlm(IOBuffer(clean_text_block), Float64)

    # Reshape, based on angular separation mode.

    # :averaged_ϕ-grid layout: 5 columns x 101 energies x 20 cosθ steps
    if angles_sep === :averaged_ϕ

        n_E    = 121
        n_cosθ = 20
        flux_reshaped = reshape(matrix_2d, n_E, n_cosθ, 5)

        # Extract matrices.
        results = Dict{String, Any}(
            "NuMu_flux"     => collect(flux_reshaped[:, :, 2]),
            "Nue_flux"      => collect(flux_reshaped[:, :, 4]),
            "AntiNuMu_flux" => collect(flux_reshaped[:, :, 3]),
            "AntiNue_flux"  => collect(flux_reshaped[:, :, 5]),
            "Energies"      => flux_reshaped[:, 1, 1],
            "cosθs"         => range(-1.0, .9, length=n_cosθ).+0.05
        )

    elseif angles_sep === :variable_ϕ

        n_E    = 121
        n_cosθ = 20
        n_ϕ    = 12

        A = reshape(matrix_2d, n_E, n_ϕ, n_cosθ, 5)   
        A = A[:, :, end:-1:1, :]                     
        A = permutedims(A, (1, 3, 2, 4))              

        results = Dict{String, Any}(
            "NuMu_flux"     => collect(A[:, :, :, 2]),
            "Nue_flux"      => collect(A[:, :, :, 4]),
            "AntiNuMu_flux" => collect(A[:, :, :, 3]),
            "AntiNue_flux"  => collect(A[:, :, :, 5]),
            "Energies"      => A[:, 1, 1, 1],
            "cosθs"         => range(-1.0, .9, length=n_cosθ) .+0.05,
            "ϕs"            => range(0, 330.0, length=n_ϕ) .+15.
        )

    elseif angles_sep === :averaged_ϕθ
        # TODO: Implement the reshaping logic for the fully integrated/averaged grid
        error("angles_separation_hf mode :averaged_ϕθ is not yet implemented.")

    else
        error("Invalid angles_separation_hf option: $angles_sep. Use :variable_ϕ, :averaged_ϕ, or :averaged_ϕθ.")
    end

    return results
end

# ============================================================================
# Function that builds the urls of the different neutrino flux models
# ============================================================================

"""
    _params_from_dict(nuflux_params::Dict{String, Any})

Rebuild the parameter struct from the dictionary `nuflux_params`.
"""
function _params_from_dict(nuflux_params::Dict{String, Any})
    ParamsType = MODEL_2_PARAMSTYPE[nuflux_params["model"]]
    fields = (Symbol(k) => v for (k, v) in nuflux_params
              if k ∉ ("model", "bin_centers_arrays"))
    return ParamsType(; fields...)
end

_build_flux_url(nuflux_params::Dict{String, Any}) =
    _build_flux_url(_params_from_dict(nuflux_params))

function _build_flux_url(p::HondaParameters)

    loc   = LOCATION_MAPPING_H[p.location]
    seas  = SEASON_MAPPING[p.season]
    ang   = ANGLES_MAPPING[p.angles_separation]
    sol   = SOLAR_MAPPING[p.solar_activity].honda 
    mtn   = MOUNTAIN_MAPPING[p.mountain_overburden]

    return "$HONDA_BASE_URL/$loc-$seas-$ang$mtn-$sol.d.gz"
end

function _build_flux_url(p::ChengParameters)

    loc1  = LOCATION_MAPPING_C[p.location][1]
    loc2  = LOCATION_MAPPING_C[p.location][2]
    seas  = SEASON_MAPPING[p.season]
    ang   = ANGLES_MAPPING[p.angles_separation]
    sol1  = SOLAR_MAPPING[p.solar_activity].cheng 
    sol2  = SOLAR_MAPPING[p.solar_activity].honda 
    mu    = MUON_PROPA_MAPPING[p.muon_propa_earth]

    return "$CHENG_BASE_URL/$loc1/$(sol1)_$mu/$loc2-$seas-$ang-$sol2.d"
end

"""
    download_nuflux(url::String)

Download a neutrino flux and return a stream of its decompressed contents.
"""
function download_nuflux(url::String)
    
    # Download into an in-memory buffer and get the raw bytes.
    io_buffer = IOBuffer()
    try
        Downloads.download(url, io_buffer; timeout = 60)
    catch e 
        if e isa Downloads.RequestError
            error("Failed to download the flux table. The URL does not exist or the server is down.\n" *
                  "URL attempted: $url\nDetails: $(sprint(showerror, e))")        
        else
            rethrow(e) 
        end
    end
    data = take!(io_buffer)
    
    # Gzip files start with the bytes 0x1f 0x8b.
    is_gzip = length(data) >= 2 && data[1] == 0x1f && data[2] == 0x8b

    # Decompress, if possible.
    decompressed_stream = is_gzip ? GzipDecompressorStream(IOBuffer(data)) : IOBuffer(data)

    return decompressed_stream
end

# ============================================================================
# Functions to interpolate neutrino fluxes using the provided bins
# ============================================================================

"""
    _get_interpolated_flux(F_orig, E_orig, cosθ_orig, model, bin_centers_arrays)

Interpolate a flux table onto a new binning.

`F_orig` is a matrix defined on the original grid `(E_orig, cosθ_orig)`.
It is interpolated onto the grid `(E_new, cosθ_new)` given by `bin_centers_arrays`.

The interpolation is linear in `log(E)` and `log(Phi)` and linear in `cosθ`.

# Arguments
- `F_orig::AbstractMatrix{<:Real}`: flux table.
- `E_orig::AbstractVector{<:Real}`: energies of the original table, in GeV.
- `cosθ_orig::AbstractVector{<:Real}`: cosθ values of the original table
  (ascending or descending, in the same order as the columns of `F_orig`).
- `model::Symbol`: atmospheric neutrino flux model used.
- `bin_centers_arrays::NTuple{2, AbstractVector{<:Real}}`: new binning `(E_new, cosθ_new)`.

# Returns
A matrix of size `length(E_new) × length(cosθ_new)`.
"""
function _get_interpolated_flux(F_orig::AbstractMatrix{<:Real}, 
                                E_orig::AbstractVector{<:Real}, 
                                cosθ_orig::AbstractVector{<:Real}, 
                                model::Symbol, 
                                bin_centers_arrays::NTuple{2, AbstractVector{<:Real}})

    # Check if the new bin centers are within the allowed range.
    E_range, cosθ_range = MODEL_2_RANGETYPE[model]
    E_new,   cosθ_new   = bin_centers_arrays
    _check_bins_in_range(E_new, cosθ_new, E_range, cosθ_range)

    # The interpolation nodes must be in ascending order.
    order     = sortperm(cosθ_orig)
    cosθ_orig = cosθ_orig[order]
    F_orig    = F_orig[:, order]

    # Interpolate in log-log space for energy and flux. Atmospheric neutrino
    # spectra fall roughly as a power law in terms of energy. 
    logE_orig = log.(E_orig)
    logF_orig = log.(F_orig)
    itp = linear_interpolation((logE_orig, cosθ_orig), logF_orig;
                               extrapolation_bc = Flat())

    # Evaluate on the new grid and go back from log(Phi) to Phi.
    Phi_new = [exp(itp(log(E), cosθ)) for E in E_new, cosθ in cosθ_new]
    return Phi_new
end

function _get_interpolated_flux(F_orig::AbstractArray{<:Real, 3},
                                E_orig::AbstractVector{<:Real},
                                cosθ_orig::AbstractVector{<:Real},
                                ϕ_orig::AbstractVector{<:Real},        
                                model::Symbol,
                                bin_centers_arrays::NTuple{3, AbstractVector{<:Real}})

    # Check if the new bin centers are within the allowed range.
    E_range, cosθ_range    = MODEL_2_RANGETYPE[model]
    E_new, cosθ_new, ϕ_new = bin_centers_arrays
    _check_bins_in_range(E_new, cosθ_new, E_range, cosθ_range)

    # The interpolation nodes must be in ascending order.
    ocosθ     = sortperm(cosθ_orig)
    oϕ        = sortperm(ϕ_orig)
    cosθ_orig = cosθ_orig[ocosθ]
    ϕ_orig    = ϕ_orig[oϕ]
    F_orig    = F_orig[:, ocosθ, oϕ]

    # ϕ is periodic (360°): pad one bin on each side so that interpolation
    # between the last and first bin wraps around instead of being flat.
    ϕ_ext   = vcat(ϕ_orig[end] - 360, ϕ_orig, ϕ_orig[1] + 360)
    Phi_ext = cat(F_orig[:, :, end:end], F_orig, F_orig[:, :, 1:1]; dims = 3)

    # Interpolate in log-log space for energy and flux, linearly in cosθ and ϕ.
    logE_orig   = log.(E_orig)
    logPhi_ext  = log.(Phi_ext)
    itp = linear_interpolation((logE_orig, cosθ_orig, ϕ_ext), logPhi_ext;
                               extrapolation_bc = Flat())

    # Evaluate on the new grid (ϕ wrapped into [0, 360)) and go back to Phi.
    Phi_new = [exp(itp(log(E), c, mod(p, 360.0)))
               for E in E_new, c in cosθ_new, p in ϕ_new]
    return Phi_new
end

# ============================================================================
# Python call
# ============================================================================

# Identifier of the PythonCall package.
const _PYTHONCALL_PKGID = Base.PkgId(
    Base.UUID("6099a3de-0909-46bc-b1f4-468b9a2dfc0d"),
    "PythonCall",
)

# Load PythonCall on first call and return the module (cheap afterwards).
_require_pythoncall() = Base.require(_PYTHONCALL_PKGID)
