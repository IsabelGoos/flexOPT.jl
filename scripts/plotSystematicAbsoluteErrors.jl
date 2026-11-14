#!/usr/bin/env julia

import Pkg
const FLEXOPT_ROOT = normpath(joinpath(@__DIR__, ".."))
Pkg.activate(FLEXOPT_ROOT)

using CairoMakie
using JLD2
using Statistics

const INPUT = get(ENV, "FLEXOPT_SYSTEMATIC_INPUT", joinpath(
    FLEXOPT_ROOT, "scripts", "tmp", "famous_equation_periodic_benchmarks",
    "wave_fd_vs_opt.jld2"))
const OUTPUT_DIR = get(ENV, "FLEXOPT_SYSTEMATIC_FIGURES", joinpath(
    FLEXOPT_ROOT, "scripts", "tmp", "famous_equation_periodic_benchmarks",
    "systematic_figures"))

const SCHEMES = ["explicitFD", "convFD", "OPT3", "OPT4", "OPT5"]
const SCENARIOS = [
    "homogeneous",
    "same_wave_phase0",
    "long_material_phase0",
    "short_material_phase0",
]
const SCENARIO_TITLES = Dict(
    "homogeneous" => "homogeneous",
    "same_wave_phase0" => "material wavelength = field wavelength",
    "long_material_phase0" => "longer material wavelength",
    "short_material_phase0" => "shorter material wavelength",
)
const COLORS = Dict(
    "explicitFD" => :black,
    "convFD" => :gray45,
    "OPT3" => :dodgerblue3,
    "OPT4" => :darkorange2,
    "OPT5" => :seagreen3,
)
const MARKERS = Dict(
    "explicitFD" => :circle,
    "convFD" => :rect,
    "OPT3" => :utriangle,
    "OPT4" => :diamond,
    "OPT5" => :star5,
)

function rows_for(rows, equation_id, scenario, scheme)
    selected = filter(rows) do row
        row.equation_id == equation_id && row.scenario == scenario &&
            row.scheme == scheme && isfinite(row.absolute_error) &&
            row.absolute_error > 0
    end
    sort!(selected; by=row -> row.dx)
    return selected
end

function slope_guide!(axis, panel_rows, order, fraction, color)
    isempty(panel_rows) && return
    xs = sort(unique(row.dx for row in panel_rows))
    # With only two resolutions this is merely a segment, and its label tends
    # to collide with the measured curves. Keep guides for proper series.
    length(xs) < 3 && return
    x1, x2 = xs[1], xs[min(end, 3)]
    finite_errors = [row.absolute_error for row in panel_rows if row.absolute_error > 0]
    isempty(finite_errors) && return
    y1 = exp(mean(log.(finite_errors))) * fraction
    ys = y1 .* (Float64[x1, x2] ./ x1) .^ order
    lines!(axis, [x1, x2], ys; color, linestyle=:dot, linewidth=1.5)
    text!(axis, x2, ys[2]; text="O($(order))", color, fontsize=13,
        align=(:left, :center))
end

function make_figure(rows, equations, filename; row_height=390)
    figure = Figure(size=(2300, row_height * length(equations) + 170),
        fontsize=20)
    axis_grid = Matrix{Axis}(undef, length(equations), length(SCENARIOS))
    for (row_index, equation) in pairs(equations)
        for (column_index, scenario) in pairs(SCENARIOS)
            axis = Axis(figure[row_index, column_index];
                xscale=log10, yscale=log10,
                title=SCENARIO_TITLES[scenario],
                xlabel=row_index == length(equations) ? "Δx" : "",
                ylabel=column_index == 1 ?
                    "$(equation.label)\nabsolute RMS residual" : "",
            )
            axis_grid[row_index, column_index] = axis
            panel_rows = NamedTuple[]
            for scheme in SCHEMES
                selected = rows_for(rows, equation.id, scenario, scheme)
                append!(panel_rows, selected)
                isempty(selected) && continue
                scatterlines!(axis,
                    [row.dx for row in selected],
                    [row.absolute_error for row in selected];
                    label=scheme, color=COLORS[scheme], marker=MARKERS[scheme],
                    markersize=10, linewidth=2.4)
            end
            slope_guide!(axis, panel_rows, 2, 0.18, :gray55)
            # Separate the fourth-order reference vertically from O(2), even
            # on panels whose measured curves span only one decade.
            slope_guide!(axis, panel_rows, 4, 0.008, :gray70)
            axis.xreversed = true
        end
    end
    Legend(figure[0, :], axis_grid[1, 1]; orientation=:horizontal,
        tellheight=true, framevisible=false, nbanks=1)
    for row in 1:size(axis_grid, 1)
        linkyaxes!(axis_grid[row, :]...)
    end
    mkpath(OUTPUT_DIR)
    save(joinpath(OUTPUT_DIR, filename * ".png"), figure; px_per_unit=1.5)
    save(joinpath(OUTPUT_DIR, filename * ".pdf"), figure)
    return figure
end

data = load(INPUT)
rows = data["convergence_rows"]

poisson = [
    (label="Poisson 1-D", id="1DpoissonHetero"),
    (label="Poisson 2-D", id="2DpoissonHetero"),
]
seismic = [
    (label="SH frequency 1-D", id="1DsismoFreqHetero"),
    (label="SH time 1-D", id="1DsismoTime"),
    (label="Elastic time 2-D (P–S RMS)", id="2DsismoTimeIsoHeteroForce"),
    (label="Elastic time 3-D (P–S RMS)", id="3DsismoTimeIsoHeteroForce"),
]

available_equations = Set(String(row.equation_id) for row in rows)
filter!(equation -> equation.id in available_equations, poisson)
filter!(equation -> equation.id in available_equations, seismic)

isempty(poisson) || make_figure(rows, poisson, "systematic_poisson_absolute_errors")
isempty(seismic) || make_figure(rows, seismic, "systematic_elastic_absolute_errors"; row_height=410)
println("Saved systematic figures in ", OUTPUT_DIR)
