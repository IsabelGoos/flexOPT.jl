#!/usr/bin/env julia

import Pkg
const FLEXOPT_ROOT = get(ENV, "FLEXOPT_ROOT", normpath(joinpath(@__DIR__, "..")))
Pkg.activate(FLEXOPT_ROOT)

using CairoMakie
using JLD2
using LinearAlgebra
using Statistics

const DATA_DIR = joinpath(FLEXOPT_ROOT, "scripts", "tmp",
    "famous_equation_periodic_benchmarks")
const REFERENCE_FILE = get(ENV, "FLEXOPT_POISSON_REFERENCE",
    joinpath(DATA_DIR, "wave_fd_vs_opt_weak.jld2"))
const OUTPUT_FILE = joinpath(DATA_DIR, "poisson1d_abide_spotz.jld2")
const FIGURE_DIR = joinpath(DATA_DIR, "systematic_figures")
const ALPHA_ABIDE = 1 / 22
const A_ABIDE = 12 / 11

const FIELD_WAVE = 2.0
const KAPPA0 = 2.5
const SCENARIOS = (
    (name="homogeneous", q=0.0, amplitude=0.0),
    (name="same_wave_phase0", q=2.0, amplitude=0.15),
    (name="long_material_phase0", q=1.0, amplitude=0.15),
    (name="short_material_phase0", q=6.0, amplitude=0.15),
)
const SCENARIO_TITLES = Dict(
    "homogeneous" => "homogeneous",
    "same_wave_phase0" => "material wavelength = field wavelength",
    "long_material_phase0" => "longer material wavelength",
    "short_material_phase0" => "shorter material wavelength",
)
const COLORS = Dict(
    "explicitFD" => :black, "convFD" => :gray50,
    "OPT3" => :dodgerblue3, "OPT4" => :darkorange2,
    "OPT5" => :mediumseagreen, "Spotz" => :purple3,
    "Abide" => :firebrick2,
)

periodic_index(i, n) = mod1(i, n)

function cyclic_compact_derivative(rhs; alpha=ALPHA_ABIDE)
    n = length(rhs)
    matrix = Matrix(Tridiagonal(fill(alpha, n - 1), ones(n),
        fill(alpha, n - 1)))
    matrix[1, end] = alpha
    matrix[end, 1] = alpha
    return matrix \ rhs
end

function periodic_errors(n, scenario)
    dx = 2pi / n
    x = (0:n-1) .* dx
    xf = ((0:n-1) .+ 0.5) .* dx
    u = cos.(FIELD_WAVE .* x)
    kface = KAPPA0 .* (1 .+ scenario.amplitude .* cos.(scenario.q .* xf))
    gradk = -KAPPA0 .* scenario.amplitude .* scenario.q .* sin.(scenario.q .* x)
    kcenter = KAPPA0 .* (1 .+ scenario.amplitude .* cos.(scenario.q .* x))
    exact = -gradk .* FIELD_WAVE .* sin.(FIELD_WAVE .* x) .-
        kcenter .* FIELD_WAVE^2 .* cos.(FIELD_WAVE .* x)

    # q_{i+1/2} = -k_{i+1/2}(u_{i+1}-u_i)/dx; -div(q) approximates
    # (k u')', consistently with ExamplePoisson1D_VarCoeff.jl.
    qspotz = similar(u)
    for i in 1:n
        qspotz[i] = -kface[i] * (u[periodic_index(i + 1, n)] - u[i]) / dx
    end
    raw_spotz = [-(qspotz[i] - qspotz[periodic_index(i - 1, n)]) / dx
                  for i in 1:n]
    lap_exact = [(exact[periodic_index(i - 1, n)] - 2exact[i] +
                  exact[periodic_index(i + 1, n)]) / dx^2 for i in 1:n]
    spotz = raw_spotz .- dx^2 / 12 .* lap_exact

    # Abide's first compact derivative lives at faces and its second compact
    # derivative at cell centres. Both systems are cyclic in this benchmark.
    rhs_du = [A_ABIDE * (u[periodic_index(i + 1, n)] - u[i]) / dx
              for i in 1:n]
    du_face = cyclic_compact_derivative(rhs_du)
    qabide = -kface .* du_face
    rhs_dq = [A_ABIDE * (qabide[i] - qabide[periodic_index(i - 1, n)]) / dx
              for i in 1:n]
    abide = -cyclic_compact_derivative(rhs_dq)

    rms(v) = sqrt(mean(abs2, v))
    return (
        dx=dx,
        spotz_absolute=rms(spotz .- exact),
        abide_absolute=rms(abide .- exact),
        force_rms=rms(exact),
    )
end

function original_example_errors(n)
    # Same manufactured fields as tmp/ExamplePoisson1D_VarCoeff.jl.
    dx = 1 / n
    xc = collect(range(-0.5 - dx / 2, 0.5 + dx / 2; length=n + 2))
    xf = collect(range(-0.5 - dx, 0.5 + dx; length=n + 3))
    u = cos.(xc)
    k = cos.(xf)
    # q = cos(x)sin(x), hence b = -q' = sin²(x)-cos²(x).
    f(x) = sin(x)^2 - cos(x)^2
    b = f.(xc)

    q = -k .* [(-sin(x)) for x in xf]
    # Replace interior face fluxes by the discrete Spotz flux.
    q[2:end-1] .= -k[2:end-1] .* diff(u) ./ dx
    raw = .-diff(q[2:end-1]) ./ dx
    lap_b = (b[1:end-2] .- 2b[2:end-1] .+ b[3:end]) ./ dx^2
    spotz = raw .- dx^2 / 12 .* lap_b

    # Compact derivatives with exact ghost derivative values, matching the
    # assumption made explicitly in the original script.
    m = n + 1
    matrix_face = Matrix(Tridiagonal(fill(ALPHA_ABIDE, m - 1), ones(m),
        fill(ALPHA_ABIDE, m - 1)))
    rhs_face = A_ABIDE .* diff(u) ./ dx
    rhs_face[1] -= ALPHA_ABIDE * (-sin(xf[1]))
    rhs_face[end] -= ALPHA_ABIDE * (-sin(xf[end]))
    du_inner = matrix_face \ rhs_face
    q_inner = -k[2:end-1] .* du_inner

    matrix_center = Matrix(Tridiagonal(fill(ALPHA_ABIDE, n - 1), ones(n),
        fill(ALPHA_ABIDE, n - 1)))
    qleft, qright = -cos(xf[1]) * (-sin(xf[1])), -cos(xf[end]) * (-sin(xf[end]))
    qext = vcat(qleft, q_inner, qright)
    rhs_center = A_ABIDE .* diff(qext[2:end-1]) ./ dx
    # Exact dq/dx at the two centre ghost nodes.
    dq(x) = sin(x)^2 - cos(x)^2
    rhs_center[1] -= ALPHA_ABIDE * dq(xc[1])
    rhs_center[end] -= ALPHA_ABIDE * dq(xc[end])
    abide = -(matrix_center \ rhs_center)

    exact = b[2:end-1]
    rms(v) = sqrt(mean(abs2, v))
    return (dx=dx, spotz_absolute=rms(spotz .- exact),
        abide_absolute=rms(abide .- exact))
end

function plot_periodic(rows, reference_rows)
    figure = Figure(size=(2100, 560), fontsize=19)
    schemes = ("explicitFD", "convFD", "OPT3", "OPT4", "OPT5")
    for (column, scenario) in pairs(SCENARIOS)
        axis = Axis(figure[1, column], xscale=log10, yscale=log10,
            title=SCENARIO_TITLES[scenario.name], xlabel="Δx",
            ylabel=column == 1 ? "absolute RMS residual" : "")
        for scheme in schemes
            selected = sort(filter(r -> r.equation_id == "1DpoissonHetero" &&
                r.scenario == scenario.name && r.scheme == scheme,
                reference_rows); by=r -> r.dx)
            isempty(selected) && continue
            scatterlines!(axis, getproperty.(selected, :dx),
                getproperty.(selected, :absolute_error); color=COLORS[scheme],
                label=scheme, linewidth=2, markersize=9)
        end
        selected = filter(r -> r.scenario == scenario.name, rows)
        for (scheme, field, marker) in (("Spotz", :spotz_absolute, :diamond),
                                        ("Abide", :abide_absolute, :rect))
            scatterlines!(axis, getproperty.(selected, :dx),
                getproperty.(selected, field); color=COLORS[scheme], marker,
                label=scheme, linewidth=2.5, markersize=11)
        end
        axis.xreversed = true
        column == 1 && axislegend(axis; position=:lb, nbanks=2, framevisible=false)
    end
    save(joinpath(FIGURE_DIR, "poisson1d_abide_spotz_periodic.png"), figure;
        px_per_unit=1.5)
    save(joinpath(FIGURE_DIR, "poisson1d_abide_spotz_periodic.pdf"), figure)
end

function plot_original(rows)
    figure = Figure(size=(850, 600), fontsize=20)
    axis = Axis(figure[1, 1], xscale=log10, yscale=log10,
        title="Original ExamplePoisson1D_VarCoeff.jl model",
        xlabel="Δx", ylabel="absolute RMS residual")
    for (scheme, field, marker) in (("Spotz", :spotz_absolute, :diamond),
                                    ("Abide", :abide_absolute, :rect))
        scatterlines!(axis, getproperty.(rows, :dx), getproperty.(rows, field);
            color=COLORS[scheme], marker, label=scheme, linewidth=2.5,
            markersize=11)
    end
    axis.xreversed = true
    axislegend(axis; position=:lb, framevisible=false)
    save(joinpath(FIGURE_DIR, "poisson1d_abide_spotz_original.png"), figure;
        px_per_unit=1.5)
    save(joinpath(FIGURE_DIR, "poisson1d_abide_spotz_original.pdf"), figure)
end

function main()
    mkpath(FIGURE_DIR)
    periodic_rows = NamedTuple[]
    for scenario in SCENARIOS, n in (16, 24, 32)
        push!(periodic_rows, merge((scenario=scenario.name, n=n),
            periodic_errors(n, scenario)))
    end
    original_rows = [merge((n=n,), original_example_errors(n))
                     for n in (10, 20, 40, 80, 160)]
    jldsave(OUTPUT_FILE; periodic_rows, original_rows,
        abide_coefficients=(alpha=ALPHA_ABIDE, a=A_ABIDE),
        reference_file=REFERENCE_FILE)
    reference_rows = load(REFERENCE_FILE)["convergence_rows"]
    plot_periodic(periodic_rows, reference_rows)
    plot_original(original_rows)
    println("Saved: ", OUTPUT_FILE)
end

abspath(PROGRAM_FILE) == abspath(@__FILE__) && main()
