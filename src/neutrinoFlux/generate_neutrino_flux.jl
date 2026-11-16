function _resolve_neutrino_flux_csv(path::AbstractString)
    candidates = String[path]
    endswith(lowercase(path), ".csv") || push!(candidates, path * ".csv")
    data_dir = normpath(joinpath(@__DIR__, "..", "..", "dataInput", "neutrinoFluxData"))
    append!(candidates, joinpath.(Ref(data_dir), basename.(copy(candidates))))
    match = findfirst(isfile, candidates)
    isnothing(match) &&
        throw(ArgumentError(
            "neutrino-flux CSV not found; tried $(join(candidates, ", "))",
        ))
    return abspath(candidates[match])
end

# Positional compatibility without the former ten-element return value.
function read_neutrino_flux_table(
    path::AbstractString,
    nEbins::Integer,
    nθbins::Integer,
    has_header::Bool,
)
    return read_neutrino_flux_table(
        path;
        nEbins=nEbins,
        nθbins=nθbins,
        has_header=has_header,
    )
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

function set_dflux_params(flux_obj)
    pythoncall = _require_pythoncall()
    return Base.invokelatest(_set_dflux_params, pythoncall, flux_obj)
end

function _set_dflux_params(pythoncall::Module, flux_obj)

    # Fetch all the Daemonflux parameter names 
    params = pythoncall.pyconvert(
        Vector{String},
        flux_obj.params.known_parameters,
    )

    # Set all parameters to 0 (meaning that their difference with respect to the
    # default parameters is zero)
    params_dict = Dict(p => 0.0 for p in params) 
    
    return params_dict

end

"""
    set_dflux_params(flux_obj, params_dict)

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
function set_dflux_params(
    flux_obj,
    params_dict::Union{Nothing,Dict{String,Any}},
)
    pythoncall = _require_pythoncall()
    return Base.invokelatest(
        _set_dflux_params,
        pythoncall,
        flux_obj,
        params_dict,
    )
end

function _set_dflux_params(
    pythoncall::Module,
    flux_obj,
    params_dict::Union{Nothing,Dict{String,Any}},
)

    # Set default Daemonflux parameters
    full_params = _set_dflux_params(pythoncall, flux_obj)

    # If params_dict is nothing, just return the defaults without trying to merge
    if isnothing(params_dict)
        return full_params
    end

    # Otherwise, update with pairs from params_dict
    merge!(full_params, params_dict)

    return full_params

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

function _validate_flux_binning(binning::NamedTuple)
    hasproperty(binning, :energies) ||
        throw(ArgumentError("binning must contain energies"))
    hasproperty(binning, :cosθgrid) ||
        throw(ArgumentError("binning must contain cosθgrid"))
    energies = Float64.(collect(binning.energies))
    cosθgrid = Float64.(collect(binning.cosθgrid))
    isempty(energies) && throw(ArgumentError("energies cannot be empty"))
    isempty(cosθgrid) && throw(ArgumentError("cosθgrid cannot be empty"))
    all(>(0), energies) || throw(ArgumentError("energies must be positive"))
    all(diff(energies) .> 0) ||
        throw(ArgumentError("energies must be strictly increasing"))
    all(-1 .<= cosθgrid .<= 1) ||
        throw(ArgumentError("cosθgrid values must lie in [-1, 1]"))
    all(diff(cosθgrid) .> 0) ||
        throw(ArgumentError("cosθgrid must be strictly increasing"))
    return (; energies, cosθgrid)
end

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
        native_cosθ = raw["cosθs"]
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

"""
    produce_neutrino_flux(binning;
        model=DEFAULT_NUFLUX_MODEL[], params=DEFAULT_NUFLUX_PARAMS[])

Produce the selected default flux directly on a flexOPT binning NamedTuple:

    (energies=..., cosθgrid=...)

Unlike the legacy dictionary method, this returns a `NeutrinoFluxTable` whose
matrices use `(energy, cosθ)` ordering.
"""
function produce_neutrino_flux(
    binning::NamedTuple;
    model::Symbol=DEFAULT_NUFLUX_MODEL[],
    params::NeutrinoFluxParameters=DEFAULT_NUFLUX_PARAMS[],
)
    selected_model = _validate_neutrino_flux_model(model)
    selected_params = _validate_neutrino_flux_parameters(params)
    validated_binning = _validate_flux_binning(binning)
    backend_params = _parameters_dict(
        selected_params;
        model=selected_model,
        bin_centers_arrays=(
            validated_binning.energies,
            validated_binning.cosθgrid,
        ),
    )
    raw = produce_neutrino_flux(backend_params)
    return _flux_table_from_raw(raw, validated_binning, selected_model)
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

"""
    produce_neutrino_flux(bin_centers_arrays; model, flux_mode)

Produce the atmospheric neutrino flux for all relevant neutrino types following a given `model`.

If the `model` `Daemonflux` is chosen, compute the atmospheric neutrino fluxes at `IceCube` 
using the Daemonflux model (arXiv:2303.00022).
`Icecube` is, for the moment, the only location for which upgoing neutrinos are computed (July 2026).
`use_calibration` is set to true so that the Daemonflux calibration is used.
If the `model` `Honda` is chosen, fetch the atmospheric neutrino fluxes from 
<http://www-rccn.icrr.u-tokyo.ac.jp/mhonda/public/nflx2014/index.html> (arXiv:1502.03916).
Neither model provides tau-neutrino or -antineutrino fluxes since they are negligible.

# Arguments
- `model::String`: `Daemonflux` (default) or `Honda`.
- `bin_centers_arrays::Tuple{AbstractArray, AbstractArray}`: a 2-tuple containing the bin centers 
for energy and cosθ. The bin centers for cosθ have to be given in increasing order.
- `params_df::Union{Nothing, Dict{String, Any}}`: Daemonflux parameters. By default all Daemonflux parameters
are set to their default values. Only those provided in the `params_df` dictionary are updated to the desired values.
- `flux_mode::String`: `conventional` or `total` (default). `total` means conventional + prompt (only the 
conventional componennt is calibrated).
- `return_uncertainties::Bool`: if true, return neutrino flux uncertainties. Default is false. 
- `location_hf::String`: location for honda flux, options are `kamioka`, `gran_sasso`, `sudbury`, `frejus`, `ino`, 
`south_pole`, `pythasalmi`, `homestake`, `juno`
- `season_hf::String`: season for honda flux, options are `all_year`, `march_may`, `june_aug`, `sept_nov`, `dec_feb`
- `angles_separation_hf::String`: the way the angles are organized in the honda flux, options are 
`variable_ϕ`, `averaged_ϕ`, `averaged_ϕθ`
- `mountain_overburden_hf::String`: `with` or `without`
- `solar_activity_hf::String`: `minimum` or `maximum`

# Returns 
a 4-tuple (or 8-tuple, if `return_uncertainties::Bool` is true) with the following matrices 
in the shape (nθ, nE):
- `NuMu_flux` 
- `Nue_flux`
- `AntiNuMu_flux` 
- `AntiNue_flux`
(followed by the corresponding uncertanties if `return_uncertainties::Bool` is true, in the same order).
"""
function produce_neutrino_flux(nuflux_params::Dict{String, Any})
  
    model = Symbol(nuflux_params["model"])

    if (model === :Daemonflux)
        results = produce_daemonflux_neutrino_flux(nuflux_params)
    elseif (model === :Honda)
        results = produce_honda_neutrino_flux(nuflux_params)
    else 
        error("model $model is not supported. Use :Daemonflux or :Honda.")
    end

    return results

end

"""
    produce_daemonflux_neutrino_flux(nuflux_params::Dict{String, Any})

Produce the atmospheric neutrino flux for all relevant neutrino types using the Daemonflux model.

See produce_neutrino_flux for more information.
"""
function produce_daemonflux_neutrino_flux(nuflux_params::Dict{String, Any})
    pythoncall = _require_pythoncall()
    return Base.invokelatest(
        _produce_daemonflux_neutrino_flux,
        pythoncall,
        nuflux_params,
    )
end

function _produce_daemonflux_neutrino_flux(
    pythoncall::Module,
    nuflux_params::Dict{String,Any},
)
    
    # Read some inputs
    Ebin_centers, cosθbin_centers = nuflux_params["bin_centers_arrays"]
    flux_mode            = Symbol(nuflux_params["flux_mode"])
    return_uncertainties = Bool(nuflux_params["return_uncertainties"])

    # Instantiate Daemonflux and retry transient GitHub release-asset failures.
    # Once downloaded, Daemonflux reuses its package-local spline cache.
    flux = _daemonflux_instance(pythoncall)

    # Convert to degrees, as needed by daemonflux
    θbin_centers = string.(rad2deg.(acos.(cosθbin_centers)))

    # Pre-allocate output arrays
    nE = length(Ebin_centers)
    nθ = length(θbin_centers)
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

    # Check mode and define key names depending on the chosen flux_mode 
    if flux_mode !== :total && flux_mode !== :conventional
        error("flux_mode $flux_mode is not a valid option. Use :total or :conventional.")
    end
    prefix = flux_mode === :total ? "total_" : ""
    keys = (
        numu    = prefix * "numu",
        nue     = prefix * "nue",
        anumu   = prefix * "antinumu",
        anue    = prefix * "antinue" 
    )

    # Set daemonflux parameters
    raw_params_df = get(nuflux_params, "params_df", nothing)
    params = (isnothing(raw_params_df) || isempty(raw_params_df)) ? 
             _set_dflux_params(pythoncall, flux, Dict{String, Any}()) :
             _set_dflux_params(pythoncall, flux, raw_params_df)
            
    # Energy scaling
    E_scaling = 1.e4 ./ (Ebin_centers .^ 3)
    # Fetch raw daemonflux values
    @views for (i, θ) in enumerate(θbin_centers)
        # flux values
        NuMu_flux[i, :]     .= pythoncall.pyconvert(Vector{Float64}, flux.flux(Ebin_centers, θ, keys.numu,  params=params)) .* E_scaling
        Nue_flux[i, :]      .= pythoncall.pyconvert(Vector{Float64}, flux.flux(Ebin_centers, θ, keys.nue,   params=params)) .* E_scaling
        AntiNuMu_flux[i, :] .= pythoncall.pyconvert(Vector{Float64}, flux.flux(Ebin_centers, θ, keys.anumu, params=params)) .* E_scaling
        AntiNue_flux[i, :]  .= pythoncall.pyconvert(Vector{Float64}, flux.flux(Ebin_centers, θ, keys.anue,  params=params)) .* E_scaling
        if return_uncertainties
            # flux uncertainties
            NuMu_flux_err[i, :]     .= pythoncall.pyconvert(Vector{Float64}, flux.error(Ebin_centers, θ, keys.numu))
            Nue_flux_err[i, :]      .= pythoncall.pyconvert(Vector{Float64}, flux.error(Ebin_centers, θ, keys.nue))
            AntiNuMu_flux_err[i, :] .= pythoncall.pyconvert(Vector{Float64}, flux.error(Ebin_centers, θ, keys.anumu))
            AntiNue_flux_err[i, :]  .= pythoncall.pyconvert(Vector{Float64}, flux.error(Ebin_centers, θ, keys.anue))
        end
    end

    # Prepare dictionary for output
    results = Dict{String, Matrix{Float64}}(
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

    return results

end

"""
    produce_honda_neutrino_flux(nuflux_params::Dict{String, Any})

Produce the atmospheric neutrino flux for all relevant neutrino types using the honda model.

See produce_neutrino_flux for more information.
"""
function produce_honda_neutrino_flux(nuflux_params::Dict{String, Any})

    # Read parameters and define url from where to fetch the data
    url_honda  = set_hflux_params(nuflux_params)
    angles_sep = Symbol(nuflux_params["angles_separation_hf"])
    
    # Download and decompression
    io_buffer = IOBuffer()
    try
        Downloads.download(url_honda, io_buffer)
    catch e 
        if e isa Downloads.RequestError
            error("Failed to download Honda flux. The URL does not exist or the server is down.\nURL attempted: $url_honda\nDetails: $(e.message)")
        else
            rethrow(e) 
        end
    end
    seekstart(io_buffer)
    
    # Decompress
    decompressed_stream = GzipDecompressorStream(io_buffer)
    
    # Read the text grid into a Matrix{Float64}, omit any lines with text 
    lines = readlines(decompressed_stream)
    close(decompressed_stream)
    numeric_lines = filter(line -> occursin(r"^\s*[0-9.-]", line), lines)
    clean_text_block = join(numeric_lines, "\n")
    matrix_2d = readdlm(IOBuffer(clean_text_block), Float64)
    
    # Reshaping based on angular separation mode
    if angles_sep === :averaged_ϕ

        # Grid layout: 5 columns x 101 energies x 20 cosθ steps
        n_E    = 101
        n_cosθ = 20
        flux_reshaped = reshape(matrix_2d', 5, n_E, n_cosθ)

        # Extract matrices and transpose using optimization views ('.') 
        # to match your required shape without unnecessary memory allocation.
        results = Dict{String, Any}(
            "NuMu_flux"     => collect(flux_reshaped[2, :, :]'),
            "Nue_flux"      => collect(flux_reshaped[4, :, :]'),
            "AntiNuMu_flux" => collect(flux_reshaped[3, :, :]'),
            "AntiNue_flux"  => collect(flux_reshaped[5, :, :]'),
            "Energies"      => flux_reshaped[1, :, 1],
            "cosθs"         => range(-1.0, 1.0, length=20)
        )
        return results

    elseif angles_sep === :variable_ϕ
        # TODO: Implement the reshaping logic for the azimuth-dependent grid
        error("angles_separation_hf mode :variable_ϕ is not yet implemented.")

    elseif angles_sep === :averaged_ϕθ
        # TODO: Implement the reshaping logic for the fully integrated/averaged grid
        error("angles_separation_hf mode :averaged_ϕθ is not yet implemented.")

    else
        error("Invalid angles_separation_hf option: $angles_sep. Use :variable_ϕ, :averaged_ϕ, or :averaged_ϕθ.")
    end

end
