# dev/spherical1Delastic/stages/plots.jl
# Plot/print functions of notebooks/SphericalTI_SH_oct2026.ipynb (cells
# sh-step2-10, sh-step2-ratio, sh-step3-03, sh-step3b-04), unchanged except
# display(fig) -> devSave(fig, ...) (PNG under outputs/<step>/).

function plotSHConvergence(rows; title="")
    fig = Figure(size=(1000, 420))
    ax1 = Axis(fig[1, 1]; xscale=log10, yscale=log10, xlabel="h (km)",
               ylabel="max|A W| / max|A_Ω² W|", title="residual of the exact solution " * title)
    ax2 = Axis(fig[1, 2]; xscale=log10, yscale=log10, xlabel="h (km)",
               ylabel="‖A - sK‖/‖sK‖", title="distance to conventional DSM " * title)
    names = unique(getfield.(rows, :config))
    for nm in names
        rr = filter(r -> r.config == nm, rows)
        h = getfield.(rr, :h_km)
        scatterlines!(ax1, h, getfield.(rr, :resFlexOPT); label=nm)
        scatterlines!(ax2, h, getfield.(rr, :relDiffDSM); label=nm)
    end
    rr = filter(r -> r.config == first(names), rows)
    h = getfield.(rr, :h_km)
    lines!(ax1, h, getfield.(rr, :resDSM); color=:black, linewidth=3, label="DSM conventional")
    h0, r0 = h[end], rr[end].resDSM
    lines!(ax1, h, r0 .* (h ./ h0) .^ 2; color=:gray, linestyle=:dash, label="h²")
    lines!(ax1, h, r0 .* (h ./ h0) .^ 4; color=:gray, linestyle=:dot, label="h⁴")
    axislegend(ax1; position=:rb)
    axislegend(ax2; position=:rb)
    devSave(fig, "convergence " * title)
    return fig
end

# Error constant relative to the conventional FD scheme, at the same h
function shErrorRatios(rows; reference="convFD3")
    for gp in unique(getfield.(rows, :gridParam))
        rr = filter(r -> r.gridParam == gp, rows)
        ref = only(filter(r -> r.config == reference, rr)).resFlexOPT
        @printf("h = %8.3f km   ", rr[1].h_km)
        for r in rr
            r.config == reference && continue
            @printf("%s/%s = %.3e   ", r.config, reference, r.resFlexOPT / ref)
        end
        @printf("DSMconv/%s = %.3e\n", reference, rr[1].resDSM / ref)
    end
end

function plotShellEigen(rows; modes=(1, 5), title="")
    fig = Figure(size=(1100, 450))
    for (ip, im) in enumerate(modes)
        ax = Axis(fig[1, ip]; xscale=log10, yscale=log10, xlabel="h (km)",
                  ylabel="|δΩ/Ω|", title="overtone n=$(im-1) " * title)
        for key in unique([(r.scheme, r.closure) for r in rows])
            rr = filter(r -> (r.scheme, r.closure) == key, rows)
            length(rr[1].relErr) >= im || continue
            h = getfield.(rr, :h_km)
            e = [abs(r.relErr[im]) for r in rr]
            ls = key[1] == "DSMconv" ? :solid : (occursin("traction", key[2]) ? :dash : :solid)
            scatterlines!(ax, h, max.(e, 1e-16); label="$(key[1]) $(key[2])", linestyle=ls)
        end
        ip == length(modes) && Legend(fig[1, ip+1], ax; framevisible=false, labelsize=10)
    end
    devSave(fig, "shell_eigen " * title)
    return fig
end

# max amplitude of the section = gain × station spacing (1.5–2 is the usual choice)
function plotRecordSection(traces, θdeg; keys=("DSMconv_ref", "DSMconv", "OPT3", "OPT5"), tmax=2000.0, title="", gain=1.5)
    fig = Figure(size=(1100, 900))
    Δθ = length(θdeg) > 1 ? minimum(diff(sort(collect(θdeg)))) : 10.0
    ax = Axis(fig[1, 1]; xlabel="time (s)", ylabel="epicentral distance (°)",
              title="SH (u_φ), f_y = 1 at 590 km, PREM " * title * "  [max amplitude = $(gain) × $(Δθ)°]")
    ref = traces[first(keys)]
    sel = ref.t .<= tmax
    scale = gain * Δθ / maximum(abs.(ref.traces[sel, :]))
    colors = (:black, :orange, :dodgerblue, :crimson)
    for (ik, k) in enumerate(keys)
        tr = traces[k]
        for (ir, θ) in enumerate(θdeg)
            lines!(ax, tr.t[sel], θ .+ scale .* tr.traces[sel, ir]; color=colors[ik], linewidth=ik == 1 ? 2 : 1,
                   label=ir == 1 ? k : nothing)
        end
    end
    axislegend(ax; position=:rb)
    devSave(fig, "record_section " * title); return fig
end

"""
    plotRecordSection3C(tracesByKey, θdeg; force, keys=..., tmax=2000.0, gain=1.5)

Three panels Z (u_r), R (u_θ), T (u_φ).  Per component one common scale for
all schemes: the largest |amplitude| of the reference (first key) is drawn as
`gain` × station spacing.  A component that is identically zero for this
source (e.g. T for f_z at the pole) is labelled and left empty.
"""
function plotRecordSection3C(tracesByKey, θdeg; force, keys=("DSMconv_ref", "DSMconv", "OPT3", "OPT5"),
                             tmax=2000.0, gain=2.0, normalize=:trace,
                             colors=(:black, :orange, :dodgerblue, :crimson),
                             sourceLabel="$(force) = 1 at 590 km", fileTag=String(force))
    # normalize = :trace  -> each station scaled by ITS OWN max |reference trace| (same
    #                         factor for every scheme at that station, so they stay comparable)
    # normalize = :global -> one factor per component (largest reference trace)
    # either way the scaled peak is gain × station spacing
    Δθ = length(θdeg) > 1 ? minimum(diff(sort(collect(θdeg)))) : 10.0
    fig = Figure(size=(1800, 1000))
    for (ic, (comp, cname)) in enumerate(((:Z, "vertical u_r"), (:R, "radial u_θ"), (:T, "transverse u_φ")))
        ref = getfield(tracesByKey[first(keys)], comp)
        sel = ref.t .<= tmax
        mx = maximum(abs.(ref.traces[sel, :]))
        ax = Axis(fig[1, ic]; xlabel="time (s)", ylabel=ic == 1 ? "epicentral distance (°)" : "",
                  title="$(cname), $(sourceLabel), PREM  [$(normalize) max = $(gain) × $(Δθ)°]")
        if mx == 0 || !isfinite(mx)
            text!(ax, 0.5, 0.5; text="identically zero for this source", space=:relative, align=(:center, :center))
            continue
        end
        stationMax = [maximum(abs.(@view ref.traces[sel, ir])) for ir in eachindex(θdeg)]
        scales = normalize === :trace ?
            [m > 1e-12 * mx ? gain * Δθ / m : 0.0 for m in stationMax] :   # nodal stations stay flat
            fill(gain * Δθ / mx, length(θdeg))
        for (ik, k) in enumerate(keys)
            tr = getfield(tracesByKey[k], comp)
            for (ir, θ) in enumerate(θdeg)
                lines!(ax, tr.t[sel], θ .+ scales[ir] .* tr.traces[sel, ir];
                       color=colors[mod1(ik, length(colors))], linewidth=ik == 1 ? 2 : 1,
                       label=(ir == 1 && ic == 1) ? k : nothing)
            end
        end
        ic == 1 && axislegend(ax; position=:rb)
    end
    devSave(fig, "record_section_3C_$(fileTag)_$(normalize)")
    return fig
end
