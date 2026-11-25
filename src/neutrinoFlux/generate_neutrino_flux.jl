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
a 4-tuple (or 8-tuple, if `return_uncertainties::Bool` is true) with the following matrices in the shape (nθ, nE):
- `NuMu_flux` 
- `Nue_flux`
- `AntiNuMu_flux` 
- `AntiNue_flux`
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
function _set_dflux_params(flux_obj, params_df_dict::Union{Nothing, Dict{String, Any}})
    pythoncall = _require_pythoncall()
    return Base.invokelatest(
        _set_dflux_params,
        flux_obj,
        params_dict,
        pythoncall
    )
end

function _set_dflux_params(flux_obj, params_df_dict::Union{Nothing, Dict{String, Any}}, pythoncall::Module)

    # Set default Daemonflux parameters.
    full_params = _set_dflux_params(flux_obj, pythoncall)

    # If params_df_dict is nothing, just return the defaults without trying to merge.
    if isnothing(params_df_dict)
        return full_params
    end

    # Otherwise, update with pairs from params_dict.
    merge!(full_params, params_df_dict)

    return full_params

end





















"""
    interpolate_flux_at_bin_centers(energies::Vector{Float64}, flux::Matrix{Float64}, energies_in_log::Bool)

Calculates energy bin centers, interpolates the log(flux) at log(energies).

# Arguments
- `energies`: The bin edges for the energy axis.
- `flux`: The flux data matrix.
- `energies_in_log`: If `true`, it assumes the input `energies` are evenly spaced in log-scale 
  andit calculates geometric centers. Otherwise, it calculates arithmetic centers.

# Returns
A tuple containing:
1. `bin_centers`: The calculated centers of the energy bins.
2. `interpolated_flux`: The flux evaluated at the bin centers.
"""
function interpolate_flux_at_bin_centers(
    energies::AbstractVector{<:Real},
    flux::AbstractMatrix{<:Real},
    energies_in_log::Bool,
)
    
    # calculate the bin centers    
    if energies_in_log == true
        bin_centers = sqrt.(energies[1:end-1] .* energies[2:end]) 
    else
        bin_centers = 0.5 .* (energies[1:end-1] .+ energies[2:end])
    end

    # we interpolate log(flux) vs log(energies)
    interpolators = [interpolate(log.(energies), log.(col), FritschCarlsonMonotonicInterpolation()) for col in eachcol(flux)]
    # evaluate the interpolators at the log of the bin centers, then exponentiate back
    interpolated_flux = exp.(hcat([interp.(log.(bin_centers)) for interp in interpolators]...))   
    return bin_centers, interpolated_flux

end

"""
    set_dflux_params(flux_obj)

Initialize a dictionary of Daemonflux parameters set to their default values.

# Arguments
- `flux_obj`: A Python object representing a Daemonflux instance.

# Returns
- `Dict{String, Float64}`: A dictionary mapping parameter names to `0.0`.
"""
const _PYTHONCALL_PKGID = Base.PkgId(
    Base.UUID("6099a3de-0909-46bc-b1f4-468b9a2dfc0d"),
    "PythonCall",
)

_require_pythoncall() = Base.require(_PYTHONCALL_PKGID)

function _daemonflux_instance(
    pythoncall::Module;
    retries::Integer=3,
    retry_delay::Real=2,
)
    retries >= 1 || throw(ArgumentError("retries must be at least 1"))
    daemonflux = pythoncall.pyimport("daemonflux")
    for attempt in 1:retries
        try
            return daemonflux.Flux(
                location="IceCube",
                use_calibration=true,
                debug=1,
            )
        catch error
            if attempt == retries
                message = sprint(showerror, error)
                throw(ErrorException(
                    "Daemonflux could not initialize its IceCube data after " *
                    "$retries attempts. Its official GitHub release assets may " *
                    "be temporarily unreachable. Check access to github.com and " *
                    "release-assets.githubusercontent.com, then call " *
                    "`prefetch_daemonflux_data!()` or retry the computation.\n" *
                    "Last Python error:\n$message",
                ))
            end
            sleep(retry_delay * 2.0^(attempt - 1))
        end
    end
    error("unreachable")
end

"""
    prefetch_daemonflux_data!(; retries=3, retry_delay=2)

Download, cache and validate the official Daemonflux IceCube spline and
calibration assets before a long computation. Existing cached files are reused.
This is optional because Daemonflux also initializes lazily on first use.
"""
function prefetch_daemonflux_data!(;
    retries::Integer=3,
    retry_delay::Real=2,
)
    pythoncall = _require_pythoncall()
    return Base.invokelatest(
        _prefetch_daemonflux_data!,
        pythoncall;
        retries=retries,
        retry_delay=retry_delay,
    )
end

function _prefetch_daemonflux_data!(
    pythoncall::Module;
    retries::Integer=3,
    retry_delay::Real=2,
)
    flux = _daemonflux_instance(
        pythoncall;
        retries=retries,
        retry_delay=retry_delay,
    )
    parameter_count = length(
        pythoncall.pyconvert(Vector{String}, flux.params.known_parameters),
    )
    return (
        model=:Daemonflux,
        location=:IceCube,
        ready=true,
        parameter_count=parameter_count,
    )
end




"""
    get_dflux_param_value(flux_obj, param_name, signed_sigma; max_iterations=10000, max_param=10., tol=0.01)

Find the value of the parameter `param_name` which is `signed_sigma` away from the default value.

# Arguments
- `flux_obj`: A Python object representing a Daemonflux instance.
- `param_name`: The name of the parameter (see the set_dflux_params function for possible parameter names).
- `signed_sigma`
- `max_iterations`: The maximum number of iterations in the Bisection method.
- `max_param`: The maximum search limit for the absolute value of the parameter value.
- `tol`: The relative tolerance on the target distance from the default value.
"""
function get_dflux_param_value(
    flux_obj,
    param_name,
    signed_sigma;
    max_iterations=10000,
    max_param=10.0,
    tol=0.01,
)
    pythoncall = _require_pythoncall()
    return Base.invokelatest(
        _get_dflux_param_value,
        pythoncall,
        flux_obj,
        param_name,
        signed_sigma;
        max_iterations=max_iterations,
        max_param=max_param,
        tol=tol,
    )
end

function _get_dflux_param_value(
    pythoncall::Module,
    flux_obj,
    param_name,
    signed_sigma;
    max_iterations=10000,
    max_param=10.0,
    tol=0.01,
)
    
    # Determine the sign of the parameter based on if it should be above or below the default value
    p_sign = sign(signed_sigma)
    # flux_obj.chi2 is symmetric and positive, the search is performed on the positive axis
    sigma = abs(signed_sigma)

    # Define the starting search bounds
    p_low  = 0.
    p_high = max_param

    # Evaluate the function at the upper boundary to ensure a root exists
    chi2_high = pythoncall.pyconvert(
        Float64,
        flux_obj.chi2(params=Dict(param_name => p_high)),
    )
    diff_high = sqrt(chi2_high) - sigma
    # If the target is out of bounds, fail immediately 
    if diff_high < 0.0
        error("Target sigma ($signed_sigma) is outside the search bounds. At p = $max_param, chi2 is only $chi2_high. Increase max_p.")
    end

    # Run Bisection search
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

    error("Bisection failed to converge within $max_iterations iterations. Check your parameters.")

end

function set_hflux_params(nuflux_params::Dict{String, Any})

    loc        = Symbol(nuflux_params["location_hf"])
    season     = Symbol(nuflux_params["season_hf"])
    angles_sep = Symbol(nuflux_params["angles_separation_hf"])
    mountain   = Symbol(nuflux_params["mountain_overburden_hf"])
    solar      = Symbol(nuflux_params["solar_activity_hf"])

    loc_str    = LOCATION_MAPPING[loc]        
    season_str = SEASON_MAPPING[season]
    angles_str = ANGLES_MAPPING[angles_sep]
    mtn_str    = MOUNTAIN_MAPPING[mountain]
    solar_str  = SOLAR_MAPPING[solar]

    filename = "http://www-rccn.icrr.u-tokyo.ac.jp/mhonda/public/nflx2014/$(loc_str)-$(season_str)-$(angles_str)$(mtn_str)-$(solar_str).d.gz"
    return filename
    
end

"""
    NeutrinoFluxTable

Named-tuple flux table evaluated on the exact flexOPT binning. Every matrix has
shape `(length(energies), length(cosθgrid))`, matching the first two dimensions
of the oscillation-probability arrays.
"""
const NeutrinoFluxTable = NamedTuple{
    (:model, :energies, :cosθgrid, :νμ, :νe, :antiνμ, :antiνe, :uncertainties),
}


function _resample_honda_flux(
    raw_flux,
    native_energies,
    native_cosθ,
    target_energies,
    target_cosθ,
)
    flux = Float64.(raw_flux)
    energies = Float64.(collect(native_energies))
    cosines = Float64.(collect(native_cosθ))
    size(flux) == (length(cosines), length(energies)) ||
        throw(DimensionMismatch(
            "Honda flux size must be (cosθ, energy); got $(size(flux))",
        ))
    all(>(0), energies) ||
        throw(ArgumentError("Honda energy coordinates must be positive"))

    # Interpolate log(flux) over (log(E), cosθ) to preserve positivity and the
    # many-decade energy dependence.
    log_flux = log.(max.(permutedims(flux), floatmin(Float64)))
    interpolator = LinearInterpolation(
        (log.(energies), cosines),
        log_flux;
        extrapolation_bc=Flat(),
    )
    return [
        exp(interpolator(log(energy), cosine))
        for energy in target_energies, cosine in target_cosθ
    ]
end

function _flux_table_from_raw(
    raw::AbstractDict,
    binning::NamedTuple,
    model::Symbol,
)
    (; energies, cosθgrid) = _validate_flux_binning(binning)
    flux_names = ("NuMu_flux", "Nue_flux", "AntiNuMu_flux", "AntiNue_flux")
    all(name -> haskey(raw, name), flux_names) ||
        throw(ArgumentError("flux backend returned an incomplete flux table"))

    matrices = if model === :Honda
        native_energies = raw["Energies"]
        native_cosθ = raw["cosθgrid"]
        map(flux_names) do name
            _resample_honda_flux(
                raw[name],
                native_energies,
                native_cosθ,
                energies,
                cosθgrid,
            )
        end
    else
        expected_size = (length(cosθgrid), length(energies))
        map(flux_names) do name
            matrix = Float64.(raw[name])
            size(matrix) == expected_size ||
                throw(DimensionMismatch(
                    "$name has size $(size(matrix)); expected $expected_size",
                ))
            permutedims(matrix)
        end
    end

    uncertainty_names = (
        "NuMu_flux_err",
        "Nue_flux_err",
        "AntiNuMu_flux_err",
        "AntiNue_flux_err",
    )
    uncertainties = if all(name -> haskey(raw, name), uncertainty_names)
        converted = map(uncertainty_names) do name
            matrix = Float64.(raw[name])
            size(matrix) == (length(cosθgrid), length(energies)) ||
                throw(DimensionMismatch("invalid uncertainty size for $name"))
            permutedims(matrix)
        end
        (νμ=converted[1], νe=converted[2],
         antiνμ=converted[3], antiνe=converted[4])
    else
        nothing
    end

    return (
        model=model,
        energies=energies,
        cosθgrid=cosθgrid,
        νμ=matrices[1],
        νe=matrices[2],
        antiνμ=matrices[3],
        antiνe=matrices[4],
        uncertainties=uncertainties,
    )
end



function _cached_neutrino_flux(config::AbstractDict)
    binning = (
        energies=config["energies"],
        cosθgrid=config["cosθgrid"],
    )
    table = produce_neutrino_flux(
        binning;
        model=config["model"],
        params=config["params"],
    )
    # DrWatson's metadata/tagging layer requires an AbstractDict. Keep that
    # implementation detail inside the cache and expose only the named tuple.
    return Dict{String,Any}("table" => table)
end

"""
    produce_or_load_neutrino_flux(binning;
        model=DEFAULT_NUFLUX_MODEL[], params=DEFAULT_NUFLUX_PARAMS[],
        directory="neutrinoFlux", prefix="flux")

Load a previously generated flux table or generate and cache it with
`commonBatchs.myProduceOrLoad`. The cache key contains the complete binning,
backend and typed parameter set. The returned value is the same
`NeutrinoFluxTable` named tuple as [`produce_neutrino_flux`](@ref).
"""
function produce_or_load_neutrino_flux(
    binning::NamedTuple;
    model::Symbol=DEFAULT_NUFLUX_MODEL[],
    params::NeutrinoFluxParameters=DEFAULT_NUFLUX_PARAMS[],
    directory::AbstractString="neutrinoFlux",
    prefix::AbstractString="flux",
)
    validated_binning = _validate_flux_binning(binning)
    selected_model = _validate_neutrino_flux_model(model)
    selected_params = _validate_neutrino_flux_parameters(params)
    config = Dict{String,Any}(
        "model" => selected_model,
        "params" => selected_params,
        "energies" => validated_binning.energies,
        "cosθgrid" => validated_binning.cosθgrid,
    )
    cached = myProduceOrLoad(
        _cached_neutrino_flux,
        config,
        String(directory),
        String(prefix),
    )
    haskey(cached, "table") ||
        error("invalid neutrino-flux cache entry: missing \"table\"")
    return cached["table"]
end

function _flux_component(table::NeutrinoFluxTable, flavor::Symbol)
    aliases = (
        νe=:νe,
        nue=:νe,
        νμ=:νμ,
        numu=:νμ,
        antiνe=:antiνe,
        antinue=:antiνe,
        antiνμ=:antiνμ,
        antinumu=:antiνμ,
    )
    hasproperty(aliases, flavor) ||
        throw(ArgumentError(
            "unknown flavor $flavor; use :νe, :νμ, :antiνe or :antiνμ",
        ))
    canonical = getproperty(aliases, flavor)
    return getproperty(table, canonical)
end

"""
    plot_neutrino_flux_spectra(table; energy_power=3, angular=:mean)

Create a Makie figure of the four atmospheric-flux spectra. `angular` may be
`:mean` or `:sum`; the former is independent of the number of angular bins.
"""
function plot_neutrino_flux_spectra(
    table::NeutrinoFluxTable;
    energy_power::Real=3,
    angular::Symbol=:mean,
)
    angular in (:mean, :sum) ||
        throw(ArgumentError("angular must be :mean or :sum"))
    reduce_angles = angular === :mean ?
                    matrix -> vec(sum(matrix; dims=2)) ./ size(matrix, 2) :
                    matrix -> vec(sum(matrix; dims=2))
    figure = Figure()
    axis = Axis(
        figure[1, 1];
        xscale=log10,
        yscale=log10,
        xlabel="E / GeV",
        ylabel="ϕ E^$(energy_power)",
    )
    scale = table.energies .^ energy_power
    for (field, label) in (
        (:νe, "νₑ"),
        (:νμ, "νμ"),
        (:antiνe, "ν̄ₑ"),
        (:antiνμ, "ν̄μ"),
    )
        lines!(
            axis,
            table.energies,
            scale .* reduce_angles(getproperty(table, field));
            label=label,
        )
    end
    axislegend(axis; position=:lb)
    return (; figure, axis)
end

"""
    plot_neutrino_flux_heatmap(table; flavor=:νμ, logscale=true)

Create a Makie `(energy, cosθ)` heatmap for one flavor. The returned named
tuple contains `figure`, `axis`, `plot` and `colorbar`.
"""
function plot_neutrino_flux_heatmap(
    table::NeutrinoFluxTable;
    flavor::Symbol=:νμ,
    logscale::Bool=true,
)
    values = _flux_component(table, flavor)
    plotted = logscale ? log10.(max.(values, floatmin(Float64))) : values
    figure = Figure()
    axis = Axis(
        figure[1, 1];
        xlabel="log₁₀(E / GeV)",
        ylabel="cosθ",
        title=string(flavor),
    )
    heatmap_plot = heatmap!(
        axis,
        log10.(table.energies),
        table.cosθgrid,
        plotted,
    )
    colorbar = Colorbar(
        figure[1, 2],
        heatmap_plot;
        label=logscale ? "log₁₀(ϕ)" : "ϕ",
    )
    return (; figure, axis, plot=heatmap_plot, colorbar)
end







