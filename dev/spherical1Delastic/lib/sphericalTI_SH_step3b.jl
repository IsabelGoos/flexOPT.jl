# helpers/sphericalTI_SH_step3b.jl
#
# Step 3b: SH (toroidal) synthetic seismograms in the time domain, DSM style,
# for a point force at the north pole and receivers on φ = 0.
#
#   u(r,θ,φ,ω) = Σ_l Σ_{|m|≤mmax} c_lm(ω) G_l(r,ω) C_lm(θ,φ)
#   C_lm = -k⁻¹ r̂×∇₁Y_lm = k⁻¹ [ θ̂ (im/sinθ) X_lm - φ̂ ∂θX_lm ] e^{imφ},
#          Y_lm = X_lm(θ) e^{imφ} fully normalised (Condon–Shortley), k = √(l(l+1))
#   c_lm = f · C*_lm(θs,φs)                       (point force f δ(x - xs))
#   radial: R[W_lm] = -c_lm δ(r - rs)/rs²  ⇒  Galerkin (ξ = r/R, Ω = ωR):
#           K_ξ(Ω,l) G = -X(ξs)/R,   K_ξ = Ω² T - H_L - (Λ-2) H_N
#   (|m| ≤ 2 is exact for a source on the polar axis: force |m| ≤ 1,
#    moment tensor |m| ≤ 2.)
#
# Global radial system: zones of the outer solid region (planet1D), one
# uniform h per zone, interfaces = one shared node where the two half-element
# rows of the adjacent zones are ADDED (continuity of W and of traction, the
# Galerkin/DSM way; flexOPT: native clipped p = 2 rows).  CMB and surface:
# traction-free (exact for SH).
#
# Linearity in the material fields is used to build, once per zone, the three
# basis matrices T, H_L, H_N; every (l, ω) system is then Ω²T - H_L - (Λ-2)H_N.
#
# Time domain: complex frequency ω - iω_I, ω_I = -ln(omegaiTlen)/tlen (DSM),
# flat source spectrum (fy = 1 at every ω), then inverse Fourier transform and
# e^{ω_I t} -- optionally with a cosine taper near fmax.
#
# Requires helpers/sphericalTI_SH*.jl (Step 1, 2, 3).

using LinearAlgebra, SparseArrays, Printf

# ---------------------------------------------------------------------------
# fully normalised associated Legendre functions, 0 ≤ m ≤ mmax
# ---------------------------------------------------------------------------

"""
    normalisedLegendre(lmax, mmax, θ)

X[l+1, m+1] = N_lm P_l^m(cos θ) (Condon–Shortley, ∫|Y_lm|²dΩ = 1) and its θ
derivative D, for 0 ≤ m ≤ mmax.  θ is clamped to [1e-5, π-1e-5] so that the
pole limits (X/sinθ, ∂θX) come out with ~1e-10 relative error.
Negative m: X_{l,-m} = (-1)^m X_{lm}.
"""
function normalisedLegendre(lmax::Integer, mmax::Integer, θ::Real)
    θc = clamp(θ, 1e-5, π - 1e-5)
    x, s = cos(θc), sin(θc)
    X = zeros(lmax + 1, mmax + 1)
    D = zeros(lmax + 1, mmax + 1)
    pmm = 1 / sqrt(4π)
    for m in 0:mmax
        m > 0 && (pmm *= -sqrt((2m + 1) / (2m)) * s)
        m <= lmax || break
        X[m+1, m+1] = pmm
        m + 1 <= lmax && (X[m+2, m+1] = sqrt(2m + 3) * x * pmm)
        for l in m+2:lmax
            a = sqrt((4l^2 - 1) / (l^2 - m^2))
            b = sqrt(((l - 1)^2 - m^2) / (4(l - 1)^2 - 1))
            X[l+1, m+1] = a * (x * X[l, m+1] - b * X[l-1, m+1])
        end
        for l in m:lmax
            prev = l - 1 >= m ? X[l, m+1] : 0.0
            c = l > m ? sqrt((2l + 1) / (2l - 1) * (l - m) * (l + m)) : 0.0
            D[l+1, m+1] = (l * x * X[l+1, m+1] - c * prev) / s
        end
    end
    return (X=X, D=D, sinθ=s)
end

"(C_θ, C_φ) of the normalised toroidal harmonic C_lm at (θ, φ), from a Legendre table."
function toroidalHarmonic(tab, l::Integer, m::Integer, φ::Real)
    l >= 1 || return (zero(ComplexF64), zero(ComplexF64))
    am = abs(m)
    sgn = (m < 0 && isodd(am)) ? -1.0 : 1.0
    X = sgn * tab.X[l+1, am+1]
    D = sgn * tab.D[l+1, am+1]
    k = sqrt(l * (l + 1))
    e = cis(m * φ)
    return ((im * m) * (X / tab.sinθ) * e / k, -D * e / k)
end

# ---------------------------------------------------------------------------
# global radial system (all zones of the outer solid region)
# ---------------------------------------------------------------------------

struct SHGlobalSystem
    ξ::Vector{Float64}
    T::SparseMatrixCSC{Float64,Int}
    HL::SparseMatrixCSC{Float64,Int}
    HN::SparseMatrixCSC{Float64,Int}
    label::String
end

"Global node list (shared interface nodes counted once) and per-zone offsets."
function _globalNodes(grid)
    ξ = Float64[]
    offsets = Int[]
    for (iz, seg) in enumerate(grid.segments)
        if iz == 1
            push!(offsets, 0); append!(ξ, seg.ξ)
        else
            abs(seg.ξ[1] - ξ[end]) <= 1e-12 || error("zones $(iz-1), $iz are not contiguous")
            push!(offsets, length(ξ) - 1); append!(ξ, seg.ξ[2:end])
        end
    end
    return ξ, offsets
end

function _scatterAdd!(G, Z, off)
    n = size(Z, 1)
    @inbounds for j in 1:n, i in 1:n
        v = Z[i, j]
        iszero(v) || (G[off+i, off+j] += v)
    end
    return G
end

# configs registry used by _flexOPTZoneBasis to find the OPT3 config
const _configsByName = Ref{Any}(Dict{String,Any}())

"Conventional DSM zone basis: T, H_L (Λ = 2), H_N (coefficient of Λ-2)."
function _dsmZoneBasis(model, seg)
    d1 = dsmSHConventional(model, seg, 1, 1.0)   # K = T - H(Λ=2)
    d2 = dsmSHConventional(model, seg, 2, 1.0)   # Λ = 6
    return (T=d1.T, HL=d1.H, HN=(d2.H .- d1.H) ./ 4)
end

"""
flexOPT zone basis.  The scheme's recipe for the interior rows, OPT3 rows
where the stencil cannot reach (nearRows = :fallback) or moving ν
(nearRows = :movingNu), native clipped (p = 2) half-element rows at both ends.
Rows are divided by their recipe scale (σh: σ = 2 for the strong-form FD3,
1 otherwise) so that the result is on the K_ξ scale of conventional DSM.
"""
function _flexOPTZoneBasis(config, recipes, model, seg; nearRows=:fallback, clipP=2)
    ξ = seg.ξ; h = seg.hξ; n = length(ξ)
    mats = [shZoneMaterials(model, seg.zone, x) for x in ξ]
    ρ = getfield.(mats, :ρ); L = getfield.(mats, :L); N = getfield.(mats, :N)
    z = zeros(n)
    sets = (T=Any[h^2 .* ρ .* ξ .^ 2, z, z, z, 1.0],
            HL=Any[z, L .* ξ .^ 2, h .* L .* ξ, h^2 .* L, 0.0],
            HN=Any[z, z, z, h^2 .* N, 0.0])
    rec = nearRows === :movingNu ?
        shRecipe(withRecipeExtras(config, Dict("nuGeometryMode" => :all))) :
        recipes[config.name]
    σ = config.orderBspace == -1 ? 2.0 : 1.0
    offs = recipeOffsets(recipes[config.name])
    bad = union(collect(1:(-minimum(offs))), collect((n - maximum(offs) + 1):n))
    inner = nearRows === :movingNu ? Int[] : setdiff(bad, [1, n])
    opt3 = _configsByName[]["OPT3"]
    clipRec = shRecipe(shClippedBoundaryConfig(opt3; p=clipP))
    rowScale = fill(σ * h, n)
    rowScale[inner] .= h; rowScale[[1, n]] .= h
    out = Dict{Symbol,Matrix{Float64}}()
    for (key, m) in pairs(sets)
        A = flexOPTSHMatrix(rec, m)
        if !isempty(inner)
            F = flexOPTSHMatrix(recipes["OPT3"], m)
            A[inner, :] .= F[inner, :]
        end
        top = flexOPTSHMatrix(clipRec, m)
        bot = flexOPTSHMatrix(clipRec, _mirroredModels(m))
        A[n, :] .= top[n, :]
        A[1, :] .= reverse(bot[n, :])
        out[key] = A ./ rowScale
    end
    # K = Ω²T - H  ⇔  A_T = s·T, A_HL = -s·H_L, A_HN = -s·H_N
    return (T=out[:T], HL=-out[:HL], HN=-out[:HN])
end

"""
    shGlobalSystem(model, grid; scheme="DSMconv", configs=shConfigs, recipes=nothing,
                   nearRows=:fallback)

Assemble T, H_L, H_N on the global node list (zones joined by adding the
half-element rows at the shared interface nodes).
"""
function shGlobalSystem(model, grid; scheme="DSMconv", configs=nothing, recipes=nothing,
                        nearRows=:fallback)
    ξ, offsets = _globalNodes(grid)
    nG = length(ξ)
    T = zeros(nG, nG); HL = zeros(nG, nG); HN = zeros(nG, nG)
    if scheme != "DSMconv"
        isnothing(configs) && error("flexOPT schemes need `configs`")
        _configsByName[] = Dict(c.name => c for c in configs)
        isnothing(recipes) && (recipes = Dict(c.name => shRecipe(c) for c in configs))
    end
    for (iz, seg) in enumerate(grid.segments)
        B = scheme == "DSMconv" ? _dsmZoneBasis(model, seg) :
            _flexOPTZoneBasis(_configsByName[][scheme], recipes, model, seg; nearRows=nearRows)
        _scatterAdd!(T, B.T, offsets[iz]); _scatterAdd!(HL, B.HL, offsets[iz])
        _scatterAdd!(HN, B.HN, offsets[iz])
    end
    label = scheme == "DSMconv" ? "DSMconv" : (nearRows === :movingNu ? scheme * "+mν" : scheme)
    return SHGlobalSystem(ξ, sparse(T), sparse(HL), sparse(HN), label)
end

"Hat-function load vector for a point source at ξs (must not sit on a node)."
function _hatLoad(ξ, ξs)
    j = findlast(<=(ξs), ξ)
    (isnothing(j) || j == length(ξ)) && error("source outside the grid")
    h = ξ[j+1] - ξ[j]
    b = zeros(length(ξ))
    b[j] = (ξ[j+1] - ξs) / h
    b[j+1] = (ξs - ξ[j]) / h
    return b
end

# ---------------------------------------------------------------------------
# frequency-domain synthetics and time series
# ---------------------------------------------------------------------------

"""
    shSynthetics(sys, R; θs=0.0, φs=0.0, rs_km, force=(0.0, 0.0, 1.0),
                 receiversθdeg=10:10:180, φr=0.0, tlen=3276.8, fmax=0.05,
                 omegaiTlen=1e-2, lmaxOf=ω -> ceil(Int, ω*R/3.0) + 30, lcap=1500, mmax=2)

Frequency-domain transverse/θ displacement at the surface for a point force
f = (f_r, f_θ, f_φ) at (rs, θs, φs) (default: f_φ = 1 at the pole, φs = 0, i.e.
fy = 1, flat in ω).  Returns frequencies, uθ[ω, receiver], uφ[ω, receiver].
"""
function shSynthetics(sys::SHGlobalSystem, R; θs=0.0, φs=0.0, rs_km,
                      force=(0.0, 0.0, 1.0), receiversθdeg=10:10:180, φr=0.0,
                      tlen=3276.8, fmax=0.05, omegaiTlen=1e-2,
                      lmaxOf=ω -> ceil(Int, ω * R / 3.0) + 30, lcap=1500, mmax=2,
                      verbose=true)
    ωI = -log(omegaiTlen) / tlen
    kmax = floor(Int, fmax * tlen)
    ωs = [2π * k / tlen for k in 0:kmax]
    θr = deg2rad.(collect(receiversθdeg))
    lmaxAll = min(lcap, maximum(lmaxOf(ω) for ω in ωs))
    recTabs = [normalisedLegendre(lmaxAll, mmax, θ) for θ in θr]
    srcTab = normalisedLegendre(lmaxAll, mmax, θs)
    # source coefficients c_lm = f·C*_lm(θs, φs) (toroidal: f_r does not couple)
    c = zeros(ComplexF64, lmaxAll + 1, 2mmax + 1)
    for l in 1:lmaxAll, m in -min(l, mmax):min(l, mmax)
        Cθ, Cφ = toroidalHarmonic(srcTab, l, m, φs)
        c[l+1, m+mmax+1] = force[2] * conj(Cθ) + force[3] * conj(Cφ)
    end
    # receiver harmonic sums Σ_m c_lm C_lm(θr, φr) per l
    Sθ = zeros(ComplexF64, lmaxAll + 1, length(θr))
    Sφ = zeros(ComplexF64, lmaxAll + 1, length(θr))
    for (ir, tab) in enumerate(recTabs), l in 1:lmaxAll, m in -min(l, mmax):min(l, mmax)
        cm = c[l+1, m+mmax+1]
        iszero(cm) && continue
        Cθ, Cφ = toroidalHarmonic(tab, l, m, φr)
        Sθ[l+1, ir] += cm * Cθ
        Sφ[l+1, ir] += cm * Cφ
    end
    activeL = [l for l in 1:lmaxAll if any(!iszero, @view Sφ[l+1, :]) || any(!iszero, @view Sθ[l+1, :])]
    b = -_hatLoad(sys.ξ, rs_km / R) ./ R
    nG = length(sys.ξ)
    uθ = zeros(ComplexF64, length(ωs), length(θr))
    uφ = zeros(ComplexF64, length(ωs), length(θr))
    t0 = time()
    for (iω, ω) in enumerate(ωs)
        Ωc = (ω - im * ωI) * R
        lmax = min(lmaxAll, lmaxOf(ω))
        for l in activeL
            l > lmax && break
            Λ = l * (l + 1)
            A = Ωc^2 .* sys.T .- sys.HL .- (Λ - 2) .* sys.HN
            G = lu(A) \ complex.(b)
            Gs = G[nG]
            @views uθ[iω, :] .+= Gs .* Sθ[l+1, :]
            @views uφ[iω, :] .+= Gs .* Sφ[l+1, :]
        end
        verbose && (iω % 20 == 0 || iω == length(ωs)) &&
            @printf("  %s: f = %.4f Hz (%d/%d), lmax = %d, %.1f s\n", sys.label,
                    ω / 2π, iω, length(ωs), lmax, time() - t0)
    end
    return (ω=ωs, ωI=ωI, tlen=tlen, θdeg=collect(receiversθdeg), uθ=uθ, uφ=uφ)
end

"""
    shTimeSeries(spec; np=4096, taper=0.2)

Inverse Fourier transform of the one-sided spectrum (frequencies 0..fmax,
zero above), with the DSM complex-frequency correction e^{ω_I t}.
`taper` = fraction of the band at the high end with a cosine taper
(0 = raw flat-spectrum impulse response).  Returns t and traces [t, receiver].
"""
function shTimeSeries(spec; np=4096, taper=0.2, component=:uφ)
    U = getfield(spec, component)
    nω = length(spec.ω)
    w = ones(nω)
    if taper > 0
        k0 = floor(Int, (1 - taper) * (nω - 1))
        for k in k0:nω-1
            w[k+1] = 0.5 * (1 + cos(π * (k - k0) / max(1, nω - 1 - k0)))
        end
    end
    dt = spec.tlen / np
    t = collect(0:np-1) .* dt
    Δω = 2π / spec.tlen
    traces = zeros(np, size(U, 2))
    for ir in 1:size(U, 2), (it, tt) in enumerate(t)
        acc = 0.5 * real(w[1] * U[1, ir])
        for k in 2:nω
            acc += real(w[k] * U[k, ir] * cis(spec.ω[k] * tt))
        end
        traces[it, ir] = acc * Δω / π * exp(spec.ωI * tt)
    end
    return (t=t, traces=traces)
end

# ---------------------------------------------------------------------------
# record section (needs a Makie backend loaded in Main, e.g. CairoMakie)
# ---------------------------------------------------------------------------

"""
    plotRecordSection(traces, θdeg; keys=..., tmax=2000.0, gain=1.5, title="", component="u_φ")

Record section of several schemes overlaid.  One common normalisation for all
traces and schemes: the largest |amplitude| of the reference (first key) over
the plotted window is drawn as `gain` × (station spacing in degrees), i.e.
1.5–2 inter-station distances as usual for record sections.
"""
function plotRecordSection(traces, θdeg; keys=("DSMconv_ref", "DSMconv", "OPT3", "OPT5"),
                           tmax=2000.0, gain=1.5, title="", component="u_φ",
                           colors=(:black, :orange, :dodgerblue, :crimson))
    Δθ = length(θdeg) > 1 ? minimum(diff(sort(collect(θdeg)))) : 10.0
    ref = traces[first(keys)]
    sel = ref.t .<= tmax
    scale = gain * Δθ / maximum(abs.(ref.traces[sel, :]))
    fig = Figure(size=(1100, 1000))
    ax = Axis(fig[1, 1]; xlabel="time (s)", ylabel="epicentral distance (°)",
              title="SH ($(component)), f_y = 1, PREM " * title *
                    "  [max amplitude = $(gain) × $(Δθ)°]")
    for (ik, k) in enumerate(keys)
        tr = traces[k]
        for (ir, θ) in enumerate(θdeg)
            lines!(ax, tr.t[sel], θ .+ scale .* tr.traces[sel, ir];
                   color=colors[mod1(ik, length(colors))], linewidth=ik == 1 ? 2 : 1,
                   label=ir == 1 ? k : nothing)
        end
    end
    axislegend(ax; position=:rb)
    display(fig)
    return fig
end
