# dev/raytheory/lib/rayGrid.jl
#
# Gridded 2-D media (Marmousi and any image/grid model) for the ray tracer,
# plus an independent first-arrival eikonal solver to check it.
#
# Ray theory needs a smooth model: Hamilton's equations use ∇v, and a model
# that is only piecewise constant or bilinear on a grid has ∇v discontinuous
# at every cell edge (rays then scatter numerically, not physically). So
#   1. the grid is smoothed with a Gaussian of chosen width σ (m or km),
#   2. v(x, z) and ∇v are evaluated with a cubic B-spline APPROXIMATION
#      (no prefilter: C², slightly smoother still, exact for linear v).
# σ is a physical choice — ray theory is valid when σ ≫ wavelength — and the
# ray field (caustics, multipathing) depends on it; it is kept explicit.
#
# Eikonal: first-order fast sweeping (Zhao 2005) for |∇T| = 1/v on the same
# smoothed model, with exact times near the source. It gives the FIRST
# arrival only, including head waves/diffractions into shadow zones, so it is
# compared to the earliest ray arrival only where rays exist.

"""
Regular 2-D grid: values[i, j] at x = x0 + (i-1) dx, z = z0 + (j-1) dz
(z = depth positive down, so j = 1 is the shallowest row).
"""
struct GridField2D
    values::Matrix{Float64}
    x0::Float64
    z0::Float64
    dx::Float64
    dz::Float64
end
gridX(g::GridField2D) = g.x0 .+ (0:size(g.values, 1)-1) .* g.dx
gridZ(g::GridField2D) = g.z0 .+ (0:size(g.values, 2)-1) .* g.dz
extentX(g::GridField2D) = (g.x0, g.x0 + (size(g.values, 1) - 1) * g.dx)
extentZ(g::GridField2D) = (g.z0, g.z0 + (size(g.values, 2) - 1) * g.dz)

"separable Gaussian smoothing, σx, σz in grid cells, edge values replicated"
function gaussianSmooth(A::AbstractMatrix, σx, σz)
    function smooth1(v, σ)
        σ <= 0 && return copy(v)
        r = ceil(Int, 3σ)
        w = [exp(-k^2 / (2σ^2)) for k in -r:r]; w ./= sum(w)
        n = length(v)
        [sum(w[k+r+1] * v[clamp(i + k, 1, n)] for k in -r:r) for i in 1:n]
    end
    B = similar(A, Float64)
    for j in axes(A, 2); B[:, j] = smooth1(A[:, j], σx); end
    C = similar(B)
    for i in axes(B, 1); C[i, :] = smooth1(B[i, :], σz); end
    return C
end

"""
    smoothGrid(values, x0, z0, dx, dz; σ) -> GridField2D

σ in the same length unit as dx, dz (0 = no smoothing).
"""
smoothGrid(values, x0, z0, dx, dz; σ=0.0) =
    GridField2D(gaussianSmooth(values, σ / dx, σ / dz), Float64(x0), Float64(z0), Float64(dx), Float64(dz))

# cubic B-spline basis and derivative, t in [0, 1)
@inline bsW(t) = ((1 - t)^3 / 6, (3t^3 - 6t^2 + 4) / 6, (-3t^3 + 3t^2 + 3t + 1) / 6, t^3 / 6)
@inline bsD(t) = (-(1 - t)^2 / 2, (9t^2 - 12t) / 6, (-9t^2 + 6t + 3) / 6, t^2 / 2)

"""
Value at node (i, j), linearly extrapolated outside the grid (ghost nodes),
so the B-spline stays exact for linear fields up to the edges; positions
further than one cell outside are clamped first.
"""
@inline function ghostValue(V, i, j, nx, nz)
    i = clamp(i, 0, nx + 1); j = clamp(j, 0, nz + 1)
    ic = clamp(i, 1, nx); jc = clamp(j, 1, nz)
    f = V[ic, jc]
    if i != ic
        i2 = i == 0 ? min(2, nx) : max(nx - 1, 1)
        f = 2f - V[i2, jc]
    end
    if j != jc
        j2 = j == 0 ? min(2, nz) : max(nz - 1, 1)
        fB = V[ic, j2]
        if i != ic
            i2 = i == 0 ? min(2, nx) : max(nx - 1, 1)
            fB = 2 * fB - V[i2, j2]
        end
        f = 2 * f - fB
    end
    return f
end

"(v, ∂v/∂x, ∂v/∂z) at (x, z) from the cubic B-spline approximation (linear ghost nodes at the edges)"
function gridEval(g::GridField2D, x, z)
    nx, nz = size(g.values)
    u = (x - g.x0) / g.dx; w = (z - g.z0) / g.dz
    i = floor(Int, u); j = floor(Int, w)
    tu = u - i; tw = w - j
    i += 1; j += 1                      # 1-based index of the node left of x
    Wu = bsW(tu); Du = bsD(tu); Ww = bsW(tw); Dw = bsD(tw)
    v = 0.0; vx = 0.0; vz = 0.0
    @inbounds for b in 1:4
        jj = j + b - 2
        for a in 1:4
            f = ghostValue(g.values, i + a - 2, jj, nx, nz)
            v += Wu[a] * Ww[b] * f
            vx += Du[a] * Ww[b] * f
            vz += Wu[a] * Dw[b] * f
        end
    end
    return v, vx / g.dx, vz / g.dz
end

"MRegion from gridded vp (and optionally vs; acoustic: vs = nothing → vs = vp/√3 placeholder)"
function gridRegion(gp::GridField2D, gs::Union{Nothing,GridField2D}=nothing)
    vpf(x) = gridEval(gp, x[1], x[2])[1]
    gvp(x) = (r = gridEval(gp, x[1], x[2]); [r[2], r[3]])
    if gs === nothing
        return MRegion(vpf, x -> vpf(x) / sqrt(3), gvp, x -> gvp(x) ./ sqrt(3))
    end
    return MRegion(vpf, x -> gridEval(gs, x[1], x[2])[1], gvp,
        x -> (r = gridEval(gs, x[1], x[2]); [r[2], r[3]]))
end

"""
    gridMedium2D(gp; gs=nothing) -> (medium, bbox)

One smooth region under a free surface at z = z0 (rays arrive there); rays
leaving the grid sideways or at the bottom die (pass `bbox` to the tracer).
"""
function gridMedium2D(gp::GridField2D; gs=nothing)
    med = PiecewiseMedium(Cartesian2D(), [PlaneInterface([0.0, gp.z0], [0.0, 1.0], :surface; freeSurface=true)],
        [gridRegion(gp, gs)])
    (xa, xb) = extentX(gp); (za, zb) = extentZ(gp)
    return med, [(xa, xb), (za - 1.0, zb)]
end

# ---------------------------------------------------------------- eikonal
"""
    eikonalFSM(slow, h, isrc, jsrc; v0=nothing, nsweep=8, rinit=5) -> T

First-order fast sweeping for |∇T| = slow on a square grid of spacing h,
source at node (isrc, jsrc) (fractional indices allowed). Nodes within
`rinit` cells of the source get the straight-ray time with the source
slowness (removes most of the source singularity error).
"""
function eikonalFSM(slow::AbstractMatrix, h, is, js; nsweep=8, rinit=5)
    nx, nz = size(slow)
    T = fill(Inf, nx, nz)
    fixed = falses(nx, nz)
    s0 = slow[clamp(round(Int, is), 1, nx), clamp(round(Int, js), 1, nz)]
    for j in 1:nz, i in 1:nx
        d = hypot(i - is, j - js)
        if d <= rinit
            T[i, j] = d * h * (slow[i, j] + s0) / 2; fixed[i, j] = true
        end
    end
    orders = ((1:nx, 1:nz), (nx:-1:1, 1:nz), (1:nx, nz:-1:1), (nx:-1:1, nz:-1:1))
    for _ in 1:nsweep
        for (I, J) in orders, j in J, i in I
            fixed[i, j] && continue
            a = min(i > 1 ? T[i-1, j] : Inf, i < nx ? T[i+1, j] : Inf)
            b = min(j > 1 ? T[i, j-1] : Inf, j < nz ? T[i, j+1] : Inf)
            (isinf(a) && isinf(b)) && continue
            f = slow[i, j] * h
            t = abs(a - b) >= f ? min(a, b) + f : (a + b + sqrt(2f^2 - (a - b)^2)) / 2
            t < T[i, j] && (T[i, j] = t)
        end
    end
    return T
end

"""
    eikonalOnGrid(gp, xs, zs; h=gp.dx/2) -> (T, x, z)

First arrivals on a uniform grid of spacing h resampled from the smooth
B-spline model gp (so ray tracer and eikonal see the same velocities).
"""
function eikonalOnGrid(gp::GridField2D, xs, zs; h=min(gp.dx, gp.dz), nsweep=8)
    (xa, xb) = extentX(gp); (za, zb) = extentZ(gp)
    xg = collect(xa:h:xb); zg = collect(za:h:zb)
    slow = [1 / gridEval(gp, x, z)[1] for x in xg, z in zg]
    T = eikonalFSM(slow, h, 1 + (xs - xa) / h, 1 + (zs - za) / h; nsweep=nsweep)
    return T, xg, zg
end
