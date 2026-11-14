# dev/spherical1DelasticKernel/stages/k3_born.jl
# K3: Born waveform kernels (lib/sgtKernels1D.jl) integrated over a shell vs brute
# force (central difference of pointSourceSynthetics of the perturbed models), in a
# layered model with a fluid core (planet1D layeredPSVModel):
#   δvs in a mantle shell (zone 3), δρ in the same shell, δvp in the whole fluid core;
# generic S and R (geometry/rotations included), components Z, R, T.
# Standalone run (DSMconv, h = 20 km, f ≤ 6 mHz): rel. diff ≤ 8e-4 (δvs), ≤ 1.4e-2 → 4e-3
# with h = 10 km (δρ, O(h²) of DSMconv), ≤ 8e-4 (fluid, 24 radial points); the same with
# disks stored on the grid and on-the-fly fields (k3Interp = :lagrange).

k3Schemes = ["DSMconv", "OPT3"]
k3h = 20.0
k3tlen = 1000.0; k3fmax = 0.006; k3Oi = 1e-2
k3Interp = nothing             # or :lagrange / :cubic (disks on the grid, fields on the fly)
k3R = 6371.0
k3Model(; dvs3=0.0, dvp1=0.0, drho3=0.0) = layeredPSVModel(; R=k3R,
    radii=[0.0, 3480.0, 4500.0, 5000.0, 6371.0],
    ρ=[[11.0, 0, -2.0, 0], 5.0, 5.0 * (1 + drho3), [4.0, 0.0, 0.0, 0.0]],
    Vp=[[10.0 * (1 + dvp1), 0, -1.5 * (1 + dvp1), 0], 13.0, 12.5, [9.0, 1.0, 0.0, 0.0]],
    Vs=[0.0, 7.0, 6.8 * (1 + dvs3), [5.0, 0.5, 0.0, 0.0]])
k3S = PointGeo(k3R - 600.0, deg2rad(50.0), deg2rad(20.0))
k3Rc = PointGeo(k3R, deg2rad(95.0), deg2rad(-30.0))
k3M = (1.0, 0.6, -0.8, -0.3, 0.4, -0.7)
xg, wg = gaussLegendre(6); k3rM = 4750.0 .+ 250.0 .* xg; k3wM = 250.0 .* wg
xf, wf = gaussLegendre(24); k3rF = 1740.0 .+ 1740.0 .* xf; k3wF = 1740.0 .* wf
k3m0 = k3Model()
k3Cases = [
    ("dvs_mantle", x -> k3Model(; dvs3=x), k3rM, k3wM,
        r -> (m = elasticMaterials(k3m0, 3, r / k3R); dL = 2e-3 * m.L; (L=dL, N=dL, F=-2dL))),
    ("drho_mantle", x -> k3Model(; drho3=x), k3rM, k3wM,
        r -> (m = elasticMaterials(k3m0, 3, r / k3R); e = 1e-3;
              (ρ=e * m.ρ, A=e * m.A, C=e * m.C, F=e * m.F, L=e * m.L, N=e * m.N))),
    ("dvp_core", x -> k3Model(; dvp1=x), k3rF, k3wF,
        r -> (m = elasticMaterials(k3m0, 1, r / k3R); dλ = 2e-3 * m.λ; (A=dλ, C=dλ, F=dλ))),
]
for scheme in k3Schemes, (name, bm, rq, wq, δ) in k3Cases
    params = Dict{String,Any}("stage" => "k3_born", "codeVersion" => CODE_VERSION_K["k3_born"],
        "scheme" => scheme, "h" => k3h, "case" => name, "tlen" => k3tlen, "fmax" => k3fmax,
        "omegaiTlen" => k3Oi, "interp" => string(k3Interp), "configs" => shConfigs)
    devRedo!("born", params, KER_CACHE_DIR, "k3_born")
    println("\n---- ", name, "  (", scheme, ", h = $(k3h) km)")
    rows = myProduceOrLoad(params, KER_CACHE_DIR, "k3_born") do cfg
        bornBruteForceCheck(; buildModel=bm,
            buildSystems=m -> (elasticGlobalSystem(:PSV, m, elasticRadialGrid(m; h=k3h, region=:all);
                                   scheme=scheme, configs=shConfigs, recipes=kerRecipes),
                               elasticGlobalSystem(:SH, m, elasticRadialGrid(m; h=k3h, region=:outerSolid);
                                   scheme=scheme, configs=shConfigs, recipes=kerRecipes)),
            S=k3S, M=k3M, Rcv=k3Rc, ε=1e-3, r_nodes=rq, r_weights=wq, δ=δ,
            tlen=k3tlen, fmax=k3fmax, omegaiTlen=k3Oi, nθ=128, nφ=256, iωs=[2, 4, 7],
            interp=k3Interp, p=sourceOrderOf(scheme)) |> rows -> Dict("rows" => rows)   # DrWatson saves Dicts only
    end["rows"]
    for r in rows
        @printf("  f=%.4f %s  rel.diff %.2e  (|δu|/|u| %.1e)\n", r.f, r.comp, r.rel, r.size)
    end
end
