# dev/raytheory/stages/plots.jl — figure helpers (CairoMakie), loaded by s1/s2/s3 only.
using CairoMakie
CairoMakie.activate!(type="png")

RT_COLORS = Makie.wong_colors()
rtColor(i) = RT_COLORS[mod1(i, length(RT_COLORS))]
RT_MODE_COLOR = Dict(:P => RT_COLORS[1], :S => RT_COLORS[6], :qP => RT_COLORS[1], :qS1 => RT_COLORS[6], :qS2 => RT_COLORS[2])

"draw a ray path; Cartesian (x, z) or (x, y, z) -> (x, z); colours by wave mode"
function plotRayPath!(ax, r::RayResult; linewidth=1.2, alpha=1.0)
    isempty(r.path) && return
    i0 = 1
    for i in 2:length(r.path)+1
        if i > length(r.path) || r.pathMode[i] != r.pathMode[i0]
            xs = [p[1] for p in r.path[i0:i-1]]; zs = [p[end] for p in r.path[i0:i-1]]
            if i <= length(r.path)                       # join to the next segment
                push!(xs, r.path[i][1]); push!(zs, r.path[i][end])
            end
            lines!(ax, xs, zs; color=(RT_MODE_COLOR[r.pathMode[i0]], alpha), linewidth=linewidth)
            i0 = i
        end
    end
end

"spherical ray (r, θ) in the disk (x = r sin θ, y = r cos θ)"
function plotSphericalRay!(ax, r::RayResult; linewidth=1.5)
    i0 = 1
    for i in 2:length(r.path)+1
        if i > length(r.path) || r.pathMode[i] != r.pathMode[i0]
            seg = r.path[i0:min(i, length(r.path))]
            lines!(ax, [p[1] * sin(p[2]) for p in seg], [p[1] * cos(p[2]) for p in seg];
                color=RT_MODE_COLOR[r.pathMode[i0]], linewidth=linewidth)
            i0 = i
        end
    end
end

"planet disk with its discontinuities"
function plotPlanetDisk!(ax, m::VerticalModel{Spherical1D})
    θ = range(0, 2π; length=721)
    for (k, b) in enumerate(m.boundaries)
        (k == 1) && continue
        (b.discontinuous || k == length(m.boundaries)) || continue
        r = b.zeta
        w = b.name in (:surface, :cmb, :icb) ? 1.5 : 0.5
        lines!(ax, r .* sin.(θ), r .* cos.(θ); color=(:gray30, b.name in (:surface, :cmb, :icb) ? 1.0 : 0.5), linewidth=w)
    end
    for k in 1:nzones(m)
        isFluid(m, k) || continue
        z = m.zones[k]
        poly!(ax, Makie.Point2f[(z.top * sin(t), z.top * cos(t)) for t in θ]; color=(:steelblue, 0.08))
    end
end
