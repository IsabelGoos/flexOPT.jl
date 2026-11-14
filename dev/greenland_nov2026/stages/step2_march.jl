# dev/greenland_nov2026/stages/step2_march.jl
# Runs ONLY on a cache miss of step2_opt.jl; its last expression is the Dict saved.
#
# Same predictor-corrector as src propagateLinearSystem(...; scheme=:predictor_corrector,
# predictor_mode=:lumped_future), written out here because in 3-D we cannot
# keep (a) the whole wavefield history nor (b) a source array over every grid
# point and time step. Instead:
#   b = -L u_known + R_src f_src     (only the 9 force columns of the source node)
#   u ← P⁻¹ b ;  repeat: u ← u + ω P⁻¹ (b - A u)     (P = lumped future blocks)
# and only the station nodes (every step) and surface maps (every snapshotEvery s)
# are kept.

Np = preparedOPT.NpointsSpace
NF = preparedOPT.NField
Afuture = preparedOPT.A_unknown
Lknown = preparedOPT.L_known
timedPredictor = @timed prepareLumpedFuturePredictor(Afuture, preparedOPT)
lumpedPredictor = timedPredictor.value
@show extrema(lumpedPredictor.condition_estimates)

wholeLI = LinearIndices(preparedOPT.spaceShape)
padLower = CartesianIndex(Tuple(optPadding[1, :]))
toWhole(I::CartesianIndex{3}) = wholeLI[I + padLower]
sourceLinear = toWhole(sourceIndex)
stationLinear = [toWhole(s.index) for s in stations]
surfaceLinear = [toWhole(CartesianIndex(i, j, topIndex[i, j])) for i in 1:Nx, j in 1:Ny]

# force columns of R for (component c, time slot s): node + (c-1) Np + (s-1) Np NF
NFP, NFF = preparedOPT.NforcePoints, preparedOPT.NForceField
sourceColumns = [sourceLinear + (c - 1) * NFP + (s - 1) * NFP * NFF for c in 1:3, s in 1:3]
Rsource = preparedOPT.R_force[:, vec(sourceColumns)]          # c fastest, then s
# body-force density F/Δx³ in the scaled PDE: f' = F Δt² / (ρ0 Δx³)  [metres]
sourceScale = SOURCE.peakForce * dtOPT^2 / (rho0 * dx^3)
forceValues = zeros(9)

knownField = zeros(Np, NF, 2)            # (point, field, slot): slot 2 = newest
u = zeros(Np * NF)                       # column = point + Np (field-1)
rhs = zeros(size(Afuture, 1))
residual = similar(rhs)
correction = similar(u)
nSta = length(stations)
traceRaw = zeros(Float32, nSta, 3, optSteps + 1)
snapshotStride = max(1, round(Int, snapshotEvery / dtOPT))
snapshots = Array{Float32,3}[]
snapshotTimes = Float64[]
iterationsPerStep = zeros(Int, optSteps)
residualPerStep = zeros(optSteps)
stoppedEarly = false
lastStep = optSteps
progressStride = max(1, optSteps ÷ 25)
marchStart = time()

for n in 1:optSteps
    # f at t = (n-2)Δt, (n-1)Δt, nΔt for the solve of u(nΔt) (Kirishima convention)
    for s in 1:3
        g = landslide_force_shape((n - 3 + s) * dtOPT; onset=SOURCE.onset, duration=SOURCE.duration)
        for c in 1:3
            forceValues[c + 3(s - 1)] = sourceScale * forceDirection[c] * g
        end
    end
    mul!(rhs, Lknown, vec(knownField), -1.0, 0.0)
    mul!(rhs, Rsource, forceValues, 1.0, 1.0)
    applyLumpedFuturePredictor!(u, lumpedPredictor, rhs)
    rhsNorm = max(norm(rhs), eps())
    used = 0
    for it in 1:optCorrectorIterations
        mul!(residual, Afuture, u)
        residual .= rhs .- residual
        norm(residual) / rhsNorm <= optCorrectorTolerance && break
        applyLumpedFuturePredictor!(correction, lumpedPredictor, residual)
        u .+= optCorrectorRelaxation .* correction
        used = it
    end
    mul!(residual, Afuture, u)
    residual .= rhs .- residual
    iterationsPerStep[n] = used
    residualPerStep[n] = norm(residual) / rhsNorm

    amplitude = maximum(abs, u)
    if !isfinite(amplitude) || amplitude > blowupLimit
        @warn "stopping: displacement blew up" n amplitude
        global stoppedEarly = true
        global lastStep = n - 1
        break
    end
    knownField[:, :, 1] .= @view knownField[:, :, 2]
    knownField[:, :, 2] .= reshape(u, Np, NF)
    for (is, p) in enumerate(stationLinear), c in 1:3
        traceRaw[is, c, n + 1] = u[p + (c - 1) * Np]
    end
    if n % snapshotStride == 0
        push!(snapshots, Float32[u[p + (c - 1) * Np] for c in 1:3, p in surfaceLinear])
        push!(snapshotTimes, n * dtOPT)
    end
    if n % progressStride == 0 || n == optSteps
        elapsed = time() - marchStart
        @printf("  step %6d / %d  t = %7.1f s  max|u| = %.2e m  P-C it %d res %.1e  (%.0f s, ETA %.0f s)\n",
            n, optSteps, n * dtOPT, amplitude, used, residualPerStep[n], elapsed,
            elapsed / n * (optSteps - n))
    end
end
marchTime = time() - marchStart

# box xyz -> Z, N, E at every station
rawTimes = (0:lastStep) .* dtOPT
traceBoxAll = traceRaw[:, :, 1:lastStep+1]
traceZNEAll = similar(traceBoxAll)
for (is, s) in enumerate(stations), k in axes(traceBoxAll, 3)
    v = SVector{3,Float64}(traceBoxAll[is, 1, k], traceBoxAll[is, 2, k], traceBoxAll[is, 3, k])
    traceZNEAll[is, 1, k] = dot(s.zLocal, v)
    traceZNEAll[is, 2, k] = dot(s.nLocal, v)
    traceZNEAll[is, 3, k] = dot(s.eLocal, v)
end
# resample to traceSampling (linear interpolation; dt << traceSampling not required)
outTimes = collect(0.0:traceSampling:last(rawTimes))
resampleTraces(T) = [linear_interpolation(rawTimes, Float64.(T[is, c, :]))(t)
                     for is in axes(T, 1), c in 1:3, t in outTimes]

Dict{String,Any}(
    "times" => outTimes,
    "tracesBox" => resampleTraces(traceBoxAll),
    "tracesZNE" => resampleTraces(traceZNEAll),
    "rawDt" => dtOPT,
    "rawTracesZNE" => traceZNEAll,
    "stationCodes" => [s.code for s in stations],
    "snapshots" => snapshots, "snapshotTimes" => snapshotTimes,
    "sourceScale" => sourceScale, "rho0" => rho0,
    "pcDiagnostics" => (iterations=iterationsPerStep[1:max(lastStep, 1)],
                        residuals=residualPerStep[1:max(lastStep, 1)]),
    "stoppedEarly" => stoppedEarly,
    "timings" => merge(operatorTimings, (predictor=timedPredictor.time, march=marchTime)),
    "wholeShape" => preparedOPT.spaceShape,
)
