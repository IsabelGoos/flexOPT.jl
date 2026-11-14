# dev/raytheory/lib/rayChristoffel.jl
#
# Layer 2 of the ray theory: the dispersion relation of a general
# (anisotropic) elastic medium,
#
#     det[ Γ(x, p) - I ] = 0 ,   Γ_ik = A_ijkl p_j p_l ,   A = C / ρ  (km²/s²)
#
# Its three eigen-branches G_m(x, p) (eigenvalues of Γ, sorted DESCENDING:
# m = 1 qP, 2 qS1 (fast shear), 3 qS2 (slow shear)) each give a Hamiltonian
#
#     H_m = ½ (G_m(x, p) - 1)     (homogeneous of degree 2 in p => σ = T)
#
# and, by Hellmann-Feynman with the unit polarisation g of branch m,
#
#     dx/dT = ∂H/∂p = A_ijkl g_i g_k p_l          (group velocity ∦ p !)
#     dp/dT = -∂H/∂x = -½ g_i g_k p_j p_l ∂A_ijkl/∂x
#
# Layer 3 (interface matching) for anisotropic media: with the tangential
# slowness t fixed, the normal component q_n solves det[Γ(t + q_n n) - I] = 0,
# a sextic, solved as a 6×6 generalized (Stroh-type) eigenproblem. Among the
# real roots on the requested branch, the outgoing one has its GROUP velocity
# leaving the interface on the right side (transmitted: v_g⋅n > 0,
# reflected: v_g⋅n < 0) — not its phase.
#
# Only Cartesian geometry here (Cartesian2D rays live in the x-z plane of a
# 3-D tensor; Cartesian3D is direct). Modes: :qP, :qS1, :qS2.
#
# Known limitation (classical ray theory): where qS1 and qS2 have equal
# velocities (shear singularities, e.g. kisses/crossings in VTI) the branch
# label by eigenvalue order is ill-defined and the ray may jump branches;
# quasi-isotropic / coupling ray theory would be needed there.

struct ChristoffelHamiltonian <: AbstractHamiltonian end

const VOIGT = [1 6 5; 6 2 4; 5 4 3]
branchIndex(mode::Symbol) = mode === :qP ? 1 : mode === :qS1 ? 2 : mode === :qS2 ? 3 :
    error("anisotropic modes are :qP, :qS1, :qS2 (got $mode)")

"""
Anisotropic region: `A(x)` returns the 6×6 density-normalised Voigt matrix
(km²/s²) and `dA(x)` the vector of its derivatives along each coordinate.
"""
struct ARegion{FA,FD}
    A::FA
    dA::FD
end
"homogeneous region"
function ARegion(A::AbstractMatrix; ndim=3)
    A6 = Matrix{Float64}(A)
    Z = [zeros(6, 6) for _ in 1:ndim]
    ARegion(x -> A6, x -> Z)
end
"A(x) = A0 + Σ_m (x_m - x0_m) G_m (G_m: 6×6 gradient matrices)"
function linearARegion(A0, x0, G::Vector)
    A0 = Matrix{Float64}(A0); x0 = collect(Float64, x0); G = [Matrix{Float64}(g) for g in G]
    ARegion(x -> A0 .+ sum((x[m] - x0[m]) .* G[m] for m in eachindex(G)), x -> G)
end

"isotropic density-normalised Voigt matrix"
function isotropicA(vp, vs)
    l2m, mu = vp^2, vs^2
    lam = l2m - 2mu
    A = zeros(6, 6)
    for i in 1:3, j in 1:3
        A[i, j] = i == j ? l2m : lam
    end
    A[4, 4] = A[5, 5] = A[6, 6] = mu
    return A
end

"""
    vtiA(vp0, vs0, ε, δ, γ) -> 6×6

VTI (symmetry axis z) from Thomsen (1986) parameters: vp0, vs0 along the
axis; C11 = C33(1+2ε), C66 = C44(1+2γ), (C13+C44)² = (C33-C44)² + 2δ C33 (C33-C44).
"""
function vtiA(vp0, vs0, ε, δ, γ)
    C33 = vp0^2; C44 = vs0^2
    C11 = C33 * (1 + 2ε); C66 = C44 * (1 + 2γ)
    C13 = sqrt((C33 - C44)^2 + 2δ * C33 * (C33 - C44)) - C44
    C12 = C11 - 2C66
    A = zeros(6, 6)
    A[1, 1] = A[2, 2] = C11; A[3, 3] = C33
    A[1, 2] = A[2, 1] = C12
    A[1, 3] = A[3, 1] = A[2, 3] = A[3, 2] = C13
    A[4, 4] = A[5, 5] = C44; A[6, 6] = C66
    return A
end

"full 3×3×3×3 tensor from Voigt"
function voigtToTensor(A)
    c = zeros(3, 3, 3, 3)
    for i in 1:3, j in 1:3, k in 1:3, l in 1:3
        c[i, j, k, l] = A[VOIGT[i, j], VOIGT[k, l]]
    end
    return c
end
function tensorToVoigt(c)
    A = zeros(6, 6)
    for i in 1:3, j in 1:3, k in 1:3, l in 1:3
        A[VOIGT[i, j], VOIGT[k, l]] = c[i, j, k, l]
    end
    return A
end
"rotate a Voigt matrix by the 3×3 rotation R (c'_ijkl = R_ia R_jb R_kc R_ld c_abcd), e.g. tilted TI"
function rotateA(A, R)
    c = voigtToTensor(A); d = zeros(3, 3, 3, 3)
    for i in 1:3, j in 1:3, k in 1:3, l in 1:3
        s = 0.0
        for a in 1:3, b in 1:3, cc in 1:3, e in 1:3
            s += R[i, a] * R[j, b] * R[k, cc] * R[l, e] * c[a, b, cc, e]
        end
        d[i, j, k, l] = s
    end
    return tensorToVoigt(d)
end
"rotation by angle θ (rad) about the y axis (tilts a z symmetry axis toward x)"
rotY(θ) = [cos(θ) 0 sin(θ); 0 1 0; -sin(θ) 0 cos(θ)]

embed3(::Cartesian2D, v) = [v[1], 0.0, v[2]]
embed3(::Cartesian3D, v) = collect(Float64, v)
project(::Cartesian2D, v) = [v[1], v[3]]
project(::Cartesian3D, v) = collect(Float64, v)

"Γ_ik = A_ijkl a_j b_l (symmetrised when a ≠ b is used through the caller)"
function christoffel(A, a, b=a)
    Γ = zeros(3, 3)
    for i in 1:3, k in 1:3
        s = 0.0
        for j in 1:3, l in 1:3
            s += A[VOIGT[i, j], VOIGT[k, l]] * a[j] * b[l]
        end
        Γ[i, k] = s
    end
    return Γ
end

"(eigenvalues descending, eigenvectors as columns in the same order)"
function branches(Γ)
    E = eigen(Symmetric(Γ))
    return reverse(E.values), E.vectors[:, 3:-1:1]
end

"group velocity of branch with polarisation g at slowness p (3-vectors)"
function groupVelocity(A, g, p)
    v = zeros(3)
    for j in 1:3
        s = 0.0
        for i in 1:3, k in 1:3, l in 1:3
            s += A[VOIGT[i, j], VOIGT[k, l]] * g[i] * g[k] * p[l]
        end
        v[j] = s
    end
    return v
end

stiffness(m::PiecewiseMedium, k::Int, x) = m.regions[k].A(x)
stiffnessGrad(m::PiecewiseMedium, k::Int, x) = m.regions[k].dA(x)

function rayRHS(::ChristoffelHamiltonian, geom::Union{Cartesian2D,Cartesian3D}, med, k, x, p, mode)
    A = stiffness(med, k, x)
    p3 = embed3(geom, p)
    λ, V = branches(christoffel(A, p3))
    g = V[:, branchIndex(mode)]
    dx = project(geom, groupVelocity(A, g, p3))
    dAs = stiffnessGrad(med, k, x)
    dp = zeros(length(x))
    for m in eachindex(dp)
        any(!iszero, dAs[m]) || continue
        dp[m] = -0.5 * dot(g, christoffel(dAs[m], p3) * g)
    end
    return dx, dp
end

function hamiltonianResidual(::ChristoffelHamiltonian, geom, med, k, x, p, mode)
    λ, _ = branches(christoffel(stiffness(med, k, x), embed3(geom, p)))
    return λ[branchIndex(mode)] - 1
end

"""
    christoffelRoots(A, t, n) -> Vector{(qn, slowness3, λ, V)}

All real normal slownesses q_n with det[Γ(t + q_n n) - I] = 0 (t ⊥ n,
3-vectors), from the companion linearisation of
(Q2 q² + Q1 q + Q0 - I) a = 0.
"""
function christoffelRoots(A, t, n)
    Q2 = christoffel(A, n)
    Q1 = christoffel(A, n, t) .+ christoffel(A, t, n)
    Q0 = christoffel(A, t) .- I(3)
    L = [zeros(3, 3) I(3); -Q0 -Q1]
    M = [I(3) zeros(3, 3); zeros(3, 3) Q2]
    E = eigen(L, M)
    out = []
    scale = 1 / sqrt(maximum(abs, Q2))
    for z in E.values
        isfinite(z) || continue
        abs(imag(z)) < 1e-7 * scale || continue
        qn = real(z)
        any(o -> abs(o[1] - qn) < 1e-10 * scale, out) && continue
        q = t .+ qn .* n
        λ, V = branches(christoffel(A, q))
        push!(out, (qn, q, λ, V))
    end
    return out
end

function interfaceSlowness(::ChristoffelHamiltonian, geom::Union{Cartesian2D,Cartesian3D}, med, kOut, x, q, n, mode, action)
    A = stiffness(med, kOut, x)
    q3 = embed3(geom, q); n3 = embed3(geom, n)
    t = q3 .- dot(q3, n3) .* n3
    want = branchIndex(mode)
    best = nothing
    for (qn, qq, λ, V) in christoffelRoots(A, t, n3)
        # which branch(es) does this root belong to (λ_m == 1)?
        for m in 1:3
            abs(λ[m] - 1) < 1e-7 || continue
            m == want || (abs(λ[want] - 1) < 1e-7) || continue   # degenerate shear: both labels ok
            vg = groupVelocity(A, V[:, m], qq)
            s = dot(vg, n3)
            ok = action === :transmit ? s > 1e-12 : s < -1e-12
            ok || continue
            # several (cusps): keep the smallest |q_n| (fastest) for now
            (best === nothing || abs(qn) < abs(best[1])) && (best = (qn, qq))
        end
    end
    best === nothing && return nothing
    return project(geom, best[2])
end

"""
    christoffelStart(med, xs, dir, mode) -> (x0, p0)

Start with phase normal `dir` (unit 2-D/3-D vector): p = dir / v_phase(dir).
"""
function christoffelStart(med, xs, dir, mode; region=nothing)
    geom = geometry(med)
    d = collect(Float64, dir) ./ norm(dir)
    k = region === nothing ? initialRegion(med, xs, d[end] > 0) : region
    λ, _ = branches(christoffel(stiffness(med, k, xs), embed3(geom, d)))
    return collect(Float64, xs), d ./ sqrt(λ[branchIndex(mode)])
end

"start from a phase-normal takeoff angle (from the downward vertical, toward +x)"
rayStart(::ChristoffelHamiltonian, med, xs, takeoff, mode; azimuth=0.0) =
    christoffelStart(med, xs, length(xs) == 2 ? [sin(takeoff), cos(takeoff)] :
        [sin(takeoff) * cos(azimuth), sin(takeoff) * sin(azimuth), cos(takeoff)], mode)
