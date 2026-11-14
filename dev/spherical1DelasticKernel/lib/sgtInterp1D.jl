# dev/spherical1DelasticKernel/lib/sgtInterp1D.jl
#
# On-the-fly fields from a stored disk: 2-D interpolation in (r, θ) inside
# each zone (block of consecutive disk radii of the same grid segment, so
# that no interpolant crosses a discontinuity), one frequency at a time.
#
#   method = :cubic | :linear   Interpolations.jl (cubic B-spline / linear) on
#                               the uniform (r, θ) grid of each block
#                               (blocks with < 4 radii fall back to linear);
#   method = :lagrange          built-in local Lagrange (order `rOrder` in r,
#                               `θOrder` in θ), no package needed.
# Default :lagrange (local): near the source the strain is singular and a
# global-in-the-zone cubic spline rings over the whole block (measured: 0.2
# vs 1e-2 at 216 km from a 600-km-deep source, both ~1e-4-1e-6 elsewhere).
# The accuracy is the user's choice through the disk sampling (`diskRadii`
# stride, θ step) and the method; `interpolationError` measures it against
# fields computed exactly at off-grid points.

using LinearAlgebra

const _HAS_INTERPOLATIONS = Ref(false)
try
    @eval using Interpolations
    _HAS_INTERPOLATIONS[] = true
catch
    @info "Interpolations.jl not available: sgtInterp1D uses method = :lagrange only"
end

"Blocks of consecutive disk radius indices sharing a grid segment (increasing r)."
function diskBlocks(disk)
    blocks = UnitRange{Int}[]
    i = 1; n = length(disk.r_km)
    while i <= n
        j = i
        while j < n && disk.seg[j+1] == disk.seg[i]
            j += 1
        end
        push!(blocks, i:j); i = j + 1
    end
    return blocks
end

_isUniform(v) = length(v) < 3 || maximum(abs.(diff(v) .- (v[end] - v[1]) / (length(v) - 1))) <=
                                  1e-9 * max(1.0, abs(v[end] - v[1]))
_asRange(v) = range(v[1], v[end]; length=length(v))

function _itp2(rv, θv, A, method)
    if method === :cubic && length(rv) >= 4 && _isUniform(rv) && _isUniform(θv)
        return Interpolations.cubic_spline_interpolation((_asRange(rv), _asRange(θv)), A;
                                                         extrapolation_bc=Interpolations.Line())
    elseif _isUniform(rv) && _isUniform(θv)
        return Interpolations.linear_interpolation((_asRange(rv), _asRange(θv)), A;
                                                   extrapolation_bc=Interpolations.Line())
    else
        return Interpolations.linear_interpolation((collect(rv), collect(θv)), A;
                                                   extrapolation_bc=Interpolations.Line())
    end
end

"""
    diskInterpolants(disk, iω; method=:lagrange, rOrder=3, θOrder=5)

Interpolants of every stored pair at frequency index iω (cheap; build one
per frequency and reuse it for all points Q).
"""
function diskInterpolants(disk, iω::Integer; method::Symbol=:lagrange, rOrder::Integer=3, θOrder::Integer=5)
    method in (:cubic, :linear, :lagrange) || error("method must be :cubic, :linear or :lagrange")
    (method !== :lagrange && !_HAS_INTERPOLATIONS[]) &&
        error("Interpolations.jl is not loaded: use method = :lagrange")
    blocks = diskBlocks(disk)
    θv = disk.θdeg
    out = []
    for b in blocks
        rv = disk.r_km[b]
        if method === :lagrange
            push!(out, (rows=b, r=rv, seg=disk.seg[first(b)]))
        else
            S = [_itp2(rv, θv, disk.strain[iω, ip, b, :], method) for ip in eachindex(disk.pairs.strain)]
            D = [_itp2(rv, θv, disk.disp[iω, ip, b, :], method) for ip in eachindex(disk.pairs.disp)]
            push!(out, (rows=b, r=rv, seg=disk.seg[first(b)], S=S, D=D))
        end
    end
    return (disk=disk, iω=iω, method=method, rOrder=rOrder, θOrder=θOrder, blocks=out)
end

"Block holding radius r (km); on a block boundary `side` chooses (:upper / :lower)."
function _blockOf(itp, r; side::Symbol=:upper)
    c = [k for (k, b) in enumerate(itp.blocks) if b.r[1] - 1e-9 <= r <= b.r[end] + 1e-9]
    isempty(c) && error("radius $(r) km is outside the disk radii")
    return itp.blocks[side === :upper ? last(c) : first(c)]
end

function _lagrangeRθ(itp, blk, A::Function, r, θdeg)
    rv = blk.r; n = length(rv)
    q = min(itp.rOrder, n - 1)
    j = clamp(searchsortedlast(rv, r), 1, n)
    lo = clamp(j - (q - 1) ÷ 2, 1, n - q)
    ridx = lo:lo+q
    wr = fornbergWeights(r, rv[ridx], 0)
    θidx, wθ = diskθWeights(itp.disk, θdeg; order=itp.θOrder)
    v = zero(ComplexF64)
    for (a, i) in enumerate(ridx), (b, j2) in enumerate(θidx)
        v += wr[a] * wθ[b] * A(blk.rows[i], j2)
    end
    return v
end

"""
    poleFieldAt(itp, a, r, θdeg, φ; which=:strain, side=:upper)

As `poleField` (lib/sgtDisk1D.jl) but at any radius r inside the disk
radii, from the interpolants `itp = diskInterpolants(disk, iω)`.
"""
function poleFieldAt(itp, a, r, θdeg, φ; which::Symbol=:strain, side::Symbol=:upper)
    disk = itp.disk
    b = rotateSourceAboutPole(disk.family, a, φ)
    blk = _blockOf(itp, r; side=side)
    pairs = which === :strain ? disk.pairs.strain : disk.pairs.disp
    out = zeros(ComplexF64, which === :strain ? 6 : 3)
    for (ip, (s, c)) in enumerate(pairs)
        iszero(b[s]) && continue
        v = if itp.method === :lagrange
            arr = which === :strain ? disk.strain : disk.disp
            _lagrangeRθ(itp, blk, (ir, iθ) -> arr[itp.iω, ip, ir, iθ], r, θdeg)
        else
            (which === :strain ? blk.S[ip] : blk.D[ip])(r, θdeg)
        end
        out[c] += b[s] * v
    end
    return which === :strain ? strainMatrix(out) : out
end

# ---------------------------------------------------------------------------
# all frequencies at once (geometry and interpolation weights do not depend on ω)
# ---------------------------------------------------------------------------

"""
    sgtAtAllω(disk, S, a, Q; side=:upper, rOrder=3, θOrder=5) -> (ε[nω, 3, 3], u[nω, 3])

Strain and displacement at Q (Q's local geographic basis) for every frequency
of the disk, source `a` at S (as `sgtAt`), local Lagrange interpolation in
(r, θ).  Much cheaper than calling `sgtAtInterp` per frequency.
"""
function sgtAtAllω(disk, S::PointGeo, a, Q::PointGeo; side::Symbol=:upper, rOrder::Integer=3, θOrder::Integer=5)
    abs(disk.rs_km - S.r) <= 1e-6 * S.r || error("disk computed for rs = $(disk.rs_km), point at $(S.r)")
    Rs, ap = sourceToPole(disk.family, a, S)
    θp, φp, T = frameOf(Rs, Q)
    b = rotateSourceAboutPole(disk.family, ap, φp)
    # radial block and weights
    c = [k for k in eachindex(disk.r_km) if true]
    blocks = diskBlocks(disk)
    cand = [bk for bk in blocks if disk.r_km[first(bk)] - 1e-9 <= Q.r <= disk.r_km[last(bk)] + 1e-9]
    isempty(cand) && error("radius $(Q.r) km is outside the disk radii")
    bk = side === :upper ? last(cand) : first(cand)
    rv = disk.r_km[bk]; n = length(rv)
    q = min(rOrder, n - 1)
    j = clamp(searchsortedlast(rv, Q.r), 1, n)
    lo = clamp(j - (q - 1) ÷ 2, 1, n - q)
    rows = bk[lo:lo+q]
    wr = fornbergWeights(Q.r, disk.r_km[rows], 0)
    θidx, wθ = diskθWeights(disk, θp; order=θOrder)
    nω = length(disk.ω)
    ε = zeros(ComplexF64, nω, 6); u = zeros(ComplexF64, nω, 3)
    for (ip, (s, cc)) in enumerate(disk.pairs.strain)
        iszero(b[s]) && continue
        for (ia, ir) in enumerate(rows), (ib, it) in enumerate(θidx)
            f = b[s] * wr[ia] * wθ[ib]
            @views ε[:, cc] .+= f .* disk.strain[:, ip, ir, it]
        end
    end
    for (ip, (s, cc)) in enumerate(disk.pairs.disp)
        iszero(b[s]) && continue
        for (ia, ir) in enumerate(rows), (ib, it) in enumerate(θidx)
            f = b[s] * wr[ia] * wθ[ib]
            @views u[:, cc] .+= f .* disk.disp[:, ip, ir, it]
        end
    end
    E = zeros(ComplexF64, nω, 3, 3); U = zeros(ComplexF64, nω, 3)
    for k in 1:nω
        Ek = T * strainMatrix(view(ε, k, :)) * T'
        E[k, :, :] .= Ek
        U[k, :] .= T * u[k, :]
    end
    return E, U
end
