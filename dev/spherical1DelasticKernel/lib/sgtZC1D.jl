# dev/spherical1DelasticKernel/lib/sgtZC1D.jl
#
# The Zhao & Chevrot (2011, GJI 185, 922, App. B) route, done right:
# the 10 RSGT functions h_i and 20 TSGT functions H_i from a pole disk, and
# the reconstruction of the strain at any Q with their (B10), (B11), (B14).
# Validated against the matrix-rotation route (sgtGeometry1D.jl) to ~1e-13:
# (B11) and (B14) are correct as printed; the old DSM-Kernel codes were not
# (see README.md, "audit").
#
# Angles (Dahlen & Tromp / ZC, as KernelMaker calculateSineCosine):
#   φ_XY = π - azimuth (clockwise from north) measured AT Y towards X,
# so φ_RQ (at Q towards R), φ_QR (at R towards Q), φ_QS (at S towards Q),
# φ_SQ (at Q towards S).
#
# Extraction geometry: the pole source/receiver and Q on the meridian φ = 0,
# where φ_RQ = π, φ_QR = 0 (receiver side) and φ_QS = 0, φ_SQ = π (source
# side) -- the old codes assumed the opposite for the terms odd in one angle.
#
# Contraction (B14 × M): H_ijpq M_pq over ALL p, q, i.e. off-diagonal moment
# components count twice (2M_rθ, 2M_rφ, 2M_θφ); `tsgt2h4freq` of KernelMaker
# used them once.

using LinearAlgebra

"Azimuth (radians, clockwise from north) at A towards B."
function azimuthAt(A::PointGeo, B::PointGeo)
    a = sphUnit(A.θ, A.φ); b = sphUnit(B.θ, B.φ); E = localBasis(A.θ, A.φ)
    t = b - dot(a, b) * a
    return atan(dot(t, E[:, 3]), dot(t, -E[:, 2]))
end

"ZC angle φ_XY = π - azimuth at Y towards X (in [0, 2π))."
zcAngle(X::PointGeo, Y::PointGeo) = mod(π - azimuthAt(Y, X), 2π)

"Epicentral distance (degrees)."
epicentralDeg(A::PointGeo, B::PointGeo) =
    rad2deg(acos(clamp(dot(sphUnit(A.θ, A.φ), sphUnit(B.θ, B.φ)), -1.0, 1.0)))

_unitM(i, j) = (A = zeros(3, 3); A[i, j] = 1.0; A[j, i] = 1.0; mtVector(A))

"""
    zcRSGT(disk, iω, ir, θdeg) -> h[1:10]

ZC's 10 receiver-side functions at (r = disk.r_km[ir], Δ_RQ = θdeg) from a
:force disk (receiver on the pole).  Z = r̂; L, T as (B10) with φ_QR = 0
at the extraction point: L̂ = -θ̂, T̂ = -φ̂ (pole-local basis).
"""
function zcRSGT(disk, iω, ir, θdeg)
    disk.family === :force || error("zcRSGT needs a :force disk")
    e(a) = poleField(disk, iω, a, ir, θdeg, 0.0; which=:strain)
    eZ = e((1.0, 0.0, 0.0)); eL = -e((0.0, 1.0, 0.0)); eT = -e((0.0, 0.0, 1.0))
    h = zeros(ComplexF64, 10)
    # (B11a) with φ_RQ = π
    h[1] = eZ[1, 1]; h[2] = eZ[1, 2]; h[3] = (eZ[2, 2] - eZ[3, 3]) / 2; h[4] = -(eZ[2, 2] + eZ[3, 3]) / 2
    # (B11b)
    h[5] = eL[1, 1]; h[6] = -(eL[2, 2] + eL[3, 3]) / 2; h[7] = (eL[3, 3] - eL[2, 2]) / 2; h[9] = -eL[1, 2]
    # (B11c)
    h[8] = -eT[2, 3]; h[10] = -eT[1, 3]
    return h
end

"""
    zcTSGT(disk, iω, ir, θdeg) -> H[1:20]

ZC's 20 source-side functions at (r = disk.r_km[ir], Δ_SQ = θdeg) from a
:moment disk (source on the pole), (B14) inverted at φ_QS = 0, φ_SQ = π.
"""
function zcTSGT(disk, iω, ir, θdeg)
    disk.family === :moment || error("zcTSGT needs a :moment disk")
    E = Dict(k => poleField(disk, iω, _unitM(ij...), ir, θdeg, 0.0; which=:strain)
             for (k, ij) in enumerate(((1, 1), (2, 2), (3, 3), (1, 2), (1, 3), (2, 3))))
    e(src, i, j) = E[src][i, j]               # src 1 rr 2 θθ 3 φφ 4 rθ 5 rφ 6 θφ (symmetric unit)
    H = zeros(ComplexF64, 20)
    H[1] = e(1, 1, 1)
    H[2] = -(e(2, 1, 1) + e(3, 1, 1)) / 2;  H[9] = (e(2, 1, 1) - e(3, 1, 1)) / 2
    H[5] = -e(4, 1, 1) / 2
    H[3] = -(e(1, 2, 2) + e(1, 3, 3)) / 2;  H[10] = (e(1, 2, 2) - e(1, 3, 3)) / 2
    H[7] = (e(4, 2, 2) + e(4, 3, 3)) / 4
    d = (e(4, 2, 2) - e(4, 3, 3)) / 4;  s = e(5, 2, 3) / 2          # H16 - H15, H16 + H15
    H[16] = (d + s) / 2; H[15] = (s - d) / 2
    H[6] = e(1, 1, 2)
    H[14] = (e(4, 1, 2) + e(5, 1, 3)) / 4; H[13] = (e(5, 1, 3) - e(4, 1, 2)) / 4
    H[8] = -(e(2, 1, 2) + e(3, 1, 2)) / 2
    H[17] = (e(2, 1, 2) - e(3, 1, 2)) / 4 - e(6, 1, 3) / 4
    H[18] = -(e(2, 1, 2) - e(3, 1, 2)) / 4 - e(6, 1, 3) / 4
    a, b, c, dd = e(2, 2, 2), e(3, 2, 2), e(2, 3, 3), e(3, 3, 3)
    H[4] = (a + b + c + dd) / 4; H[11] = (-a - b + c + dd) / 4; H[12] = (-a + b - c + dd) / 4
    sp = (a - b - c + dd) / 4; df = e(6, 2, 3) / 2                 # H20 + H19, H20 - H19
    H[20] = (sp + df) / 2; H[19] = (sp - df) / 2
    return H
end

"(B11): 3×3 strain at Q (Q's local basis) for receiver component comp ∈ (:Z, :L, :T) of (B10)."
function zcReceiverStrain(h, comp::Symbol, φRQ)
    c, s, c2, s2 = cos(φRQ), sin(φRQ), cos(2φRQ), sin(2φRQ)
    e = zeros(ComplexF64, 6)
    if comp === :Z
        e .= (h[1], h[3] * c2 - h[4], -h[3] * c2 - h[4], -h[2] * c, -h[2] * s, h[3] * s2)
    elseif comp === :L
        e .= (h[5], -h[7] * c2 - h[6], h[7] * c2 - h[6], h[9] * c, h[9] * s, -h[7] * s2)
    elseif comp === :T
        e .= (0, h[8] * s2, -h[8] * s2, -h[10] * s, h[10] * c, -h[8] * c2)
    else
        error("comp must be :Z, :L or :T")
    end
    return strainMatrix(e)
end

"""
    zcReceiverStrainSourceFrame(h, comp, S, Rcv, Q)

Receiver components defined w.r.t. the SOURCE (:Z, :R away from S, :T as
`receiverDirection`): rotation of (B10)'s (L, T) -- defined w.r.t. Q -- by
ψ = φ_SR - φ_QR at the receiver (missing in KernelMaker).
"""
function zcReceiverStrainSourceFrame(h, comp::Symbol, S::PointGeo, Rcv::PointGeo, Q::PointGeo)
    φRQ = zcAngle(Rcv, Q)
    comp === :Z && return zcReceiverStrain(h, :Z, φRQ)
    ψ = zcAngle(S, Rcv) - zcAngle(Q, Rcv)
    L = zcReceiverStrain(h, :L, φRQ); T = zcReceiverStrain(h, :T, φRQ)
    # (B10) as printed: L̂ = -θ̂ cos φ_QR - φ̂ sin φ_QR, T̂ = θ̂ sin φ_QR - φ̂ cos φ_QR (checked to 1e-14);
    # the source-referenced R̂ (away from S), T̂ (`receiverDirection`) are their rotation by ψ
    return comp === :R ? (cos(ψ) * L + sin(ψ) * T) : (-sin(ψ) * L + cos(ψ) * T)
end

"""
    zcSourceStrain(H, M, φQS, φSQ)

(B14) double-dotted with M = (Mrr, Mrθ, Mrφ, Mθθ, Mθφ, Mφφ) (S's local basis),
full contraction (off-diagonal components twice).  Returns the 3×3 strain at Q.
"""
function zcSourceStrain(H, M, φQS, φSQ)
    cqs, sqs, cqs2, sqs2 = cos(φQS), sin(φQS), cos(2φQS), sin(2φQS)
    csq, ssq, csq2, ssq2 = cos(φSQ), sin(φSQ), cos(2φSQ), sin(2φSQ)
    m = (M[1], M[4], M[6], 2M[2], 2M[3], 2M[5])           # rr θθ φφ, 2rθ 2rφ 2θφ
    T(i) = H[i]
    h01 = T(7)*cqs; h02 = T(7)*sqs; h03 = T(8)*csq; h04 = T(8)*ssq; h05 = T(9)*cqs2; h06 = T(10)*csq2
    h07 = T(11)*csq2; h08 = T(11)*ssq2; h09 = T(12)*cqs2; h10 = T(12)*sqs2
    h11 = T(13)*cos(φSQ + φQS); h12 = T(13)*sin(φSQ + φQS); h13 = T(14)*cos(φSQ - φQS); h14 = T(14)*sin(φSQ - φQS)
    h15 = T(15)*cos(2φSQ + φQS); h16 = T(15)*sin(2φSQ + φQS); h17 = T(17)*cos(φSQ + 2φQS); h18 = T(17)*sin(φSQ + 2φQS)
    h19 = T(19)*cos(2φSQ + 2φQS); h20 = T(19)*sin(2φSQ + 2φQS); h21 = T(20)*cos(2φSQ - 2φQS); h22 = T(20)*sin(2φSQ - 2φQS)
    h23 = T(16)*cos(2φSQ - φQS); h24 = T(16)*sin(2φSQ - φQS); h25 = T(18)*cos(φSQ - 2φQS); h26 = T(18)*sin(φSQ - 2φQS)
    e = zeros(ComplexF64, 6)
    e[1] = T(1)*m[1] + (h05-T(2))*m[2] + (-h05-T(2))*m[3] - T(5)*(cqs*m[4]+sqs*m[5]) + T(9)*sqs2*m[6]
    e[2] = (h06-T(3))*m[1] + (h21+h19-h07-h09+T(4))*m[2] + (-h21-h19-h07+h09+T(4))*m[3] + (h23-h15+h01)*m[4] + (-h24-h16+h02)*m[5] + (-h22+h20-h10)*m[6]
    e[3] = (-h06-T(3))*m[1] + (-h21-h19+h07-h09+T(4))*m[2] + (h21+h19+h07+h09+T(4))*m[3] + (-h23+h15+h01)*m[4] + (h24+h16+h02)*m[5] + (h22-h20-h10)*m[6]
    e[4] = -T(6)*csq*m[1] + (h25-h17+h03)*m[2] + (-h25+h17+h03)*m[3] + (-h13+h11)*m[4] + (h14+h12)*m[5] + (-h26-h18)*m[6]
    e[5] = -T(6)*ssq*m[1] + (h26-h18+h04)*m[2] + (-h26+h18+h04)*m[3] + (-h14+h12)*m[4] + (-h13-h11)*m[5] + (h25+h17)*m[6]
    e[6] = T(10)*ssq2*m[1] + (h22+h20-h08)*m[2] + (-h22-h20-h08)*m[3] + (h24-h16)*m[4] + (h23+h15)*m[5] + (h21-h19)*m[6]
    return strainMatrix(e)
end

"""
    checkZCRoute(dS, dR, iω; ntrials=40, rng)

Worst relative difference between the ZC route (zcRSGT/zcTSGT + B10/B11/B14)
and the matrix-rotation route (`sgtAt`) for random S (at dS.rs_km), R (at
dR.rs_km) and Q (at the first common disk radius).
"""
function checkZCRoute(dS, dR, iω; ntrials=40, rng=Random.MersenneTwister(1), verbose=true)
    rQ = first(filter(r -> r in dR.r_km, dS.r_km))
    irS = diskRadiusIndex(dS, rQ); irR = diskRadiusIndex(dR, rQ)
    wS = 0.0; wR = 0.0
    for _ in 1:ntrials
        S = PointGeo(dS.rs_km, acos(2rand(rng) - 1), 2π * rand(rng))
        Rc = PointGeo(dR.rs_km, acos(2rand(rng) - 1), 2π * rand(rng))
        Q = PointGeo(rQ, acos(2rand(rng) - 1), 2π * rand(rng))
        ΔS, ΔR = epicentralDeg(S, Q), epicentralDeg(Rc, Q)
        (min(ΔS, ΔR) < 2 || max(ΔS, ΔR) > 178) && continue
        M = Tuple(randn(rng, 6))
        εZC = zcSourceStrain(zcTSGT(dS, iω, irS, ΔS), M, zcAngle(Q, S), zcAngle(S, Q))
        wS = max(wS, norm(εZC - sgtAt(dS, iω, S, M, Q)) / norm(sgtAt(dS, iω, S, M, Q)))
        h = zcRSGT(dR, iω, irR, ΔR)
        for comp in (:Z, :R, :T)
            ε1 = zcReceiverStrainSourceFrame(h, comp, S, Rc, Q)
            ε2 = sgtAt(dR, iω, Rc, receiverDirection(comp, S, Rc), Q)
            wR = max(wR, norm(ε1 - ε2) / norm(ε2))
        end
    end
    verbose && @printf("ZC route (h1-h10, H1-H20 + B10/B11/B14) vs rotation route: receiver %.2e, source %.2e\n", wR, wS)
    return (receiver=wR, source=wS)
end
