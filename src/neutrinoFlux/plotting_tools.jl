

"""
    NeutrinoFluxTable

Named-tuple flux table evaluated on the exact flexOPT binning. Every matrix has
shape `(length(energies), length(cosθgrid))`, matching the first two dimensions
of the oscillation-probability arrays.
"""
const NeutrinoFluxTable = NamedTuple{
    (:model, :energies, :cosθgrid, :νμ, :νe, :antiνμ, :antiνe, :uncertainties),
}



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







