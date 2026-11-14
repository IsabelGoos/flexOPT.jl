# dev/kirishima_oct2026/stages/step2_opt_propagate.jl
# Ported from notebooks/Elastic2DtopoBenchmark_oct2026.ipynb, cell step2-*-opt-propagation.
# Runs ONLY on a cache miss of step2_opt.jl; its last expression is the Dict saved.
optPadding = isnothing(bcOPT.cerjan) ? zeros(Int, 2, 2) : cerjan_padding(bcOPT.cerjan)
@assert optOperatorUnitConvention === optFormulation "restart the kernel and rerun the OPT operator cell"
xOPT = range(first(xOPTPhysical) - optPadding[1,1] * dxOPT;
    step=dxOPT, length=preparedOPT.spaceShape[1])
zOPT = range(first(zOPTPhysical) - optPadding[1,2] * dxOPT;
    step=dxOPT, length=preparedOPT.spaceShape[2])
sourceIXOPT = argmin(abs.(xOPTPhysical .- sourcePosition.x))
sourceIZOPT = argmin(abs.(zOPTPhysical .- sourcePosition.z))
sourceIndexOPT = CartesianIndex(
    sourceIXOPT + optPadding[1,1], sourceIZOPT + optPadding[1,2])
sourceLinearOPT = LinearIndices(preparedOPT.spaceShape)[sourceIndexOPT]
sourceSpecsOPT = spatial_source_specs(xOPTPhysical, zOPTPhysical)
sourceNodesOPT = [(linear=LinearIndices(preparedOPT.spaceShape)[
        CartesianIndex(argmin(abs.(xOPT .- s.x)),
            argmin(abs.(zOPT .- s.z)))], weight=s.weight)
    for s in sourceSpecsOPT]
@assert isapprox(sum(s.weight for s in sourceNodesOPT), 1.0; atol=1e-14)
optSteps = ceil(Int, duration / dtOPT)
optOutputStride = max(1, round(Int, outputSampling / dtOPT))
# Align the force stencil with (past..., future): for three time points,
# the first solve for u(+Δt) receives f(-Δt), f(0), f(+Δt).
nPastForceLevels = preparedOPT.timePointsUsedForOneStep - 1
sourceTimesOPT = ((1 - nPastForceLevels):optSteps) .* dtOPT
@assert sourceTimesOPT[preparedOPT.timePointsUsedForOneStep] == dtOPT
waveletOPT = rickerSource.(sourceTimesOPT)
# In SI, famousEquations uses f=ρa and the cell body-force density is
# sourceForce/(ΔxΔz). Under the scaled-coordinate PDE this becomes
# f'=f Δt²/ρ = sourceForce Δt²/(ρ ΔxΔz).
sourceScaleOPT = optFormulation === :scaled_coordinates ?
    sourceForce * dtOPT^2 / (rho0 * dxOPT^2) :
    sourceForce / dxOPT^2
@assert optFormulation !== :scaled_coordinates ||
    isapprox(sourceScaleOPT * rho0 * dxOPT^2 / dtOPT^2, sourceForce)
sourceOPT = zeros(Float64, preparedOPT.NforcePoints,
    preparedOPT.NForceField, length(sourceTimesOPT))
for sourceNode in sourceNodesOPT
    sourceOPT[sourceNode.linear, 2, :] .+=
        sourceNode.weight .* sourceScaleOPT .* waveletOPT
end
predictorSelection = optPredictorMode === :lumped_future ?
    (predictor_mode=:lumped_future,) :
    (predictor_mode=:prepared, predictor_prepared=preparedPredictorVolume)
optMarchingOptions = optTimeMarchingScheme === :direct ?
    (scheme=:direct,) : merge((
        scheme=:predictor_corrector,
        corrector_iterations=optCorrectorIterations,
        corrector_tolerance=optCorrectorTolerance,
        corrector_relaxation=optCorrectorRelaxation,
    ), predictorSelection)
propagationOPT = propagateLinearSystem(
    preparedOPT, optSteps, dtOPT;
    sourceFull=sourceOPT, output_stride=optOutputStride,
    blowup_limit=1e8,
    solver_name="OPT3",
    optMarchingOptions...,
)
@assert !propagationOPT.stopped_early
if !isnothing(propagationOPT.predictor_corrector)
    pc = propagationOPT.predictor_corrector
    @show extrema(pc.iterations_per_step) extrema(pc.relative_residuals)
    @show propagationOPT.timing
end

# Controlled Γ diagnostic. Preserve, for every temporal source slot, the
# signed column sum of R_force but place it only on the collocated equation
# row. Thus only the spatial redistribution is removed.
function collocated_source_operator(prepared, source_nodes; force_field=2)
    rows=Int[]; columns=Int[]; values=Float64[]
    nspace, nfield = prepared.NpointsSpace, prepared.NField
    for node in source_nodes, time_slot in 1:prepared.timePointsUsedForOneStep
        column = node.linear + (force_field-1)*nspace +
            (time_slot-1)*nspace*prepared.NForceField
        coefficient = sum(@view prepared.R_force[:, column])
        push!(rows, force_field + nfield*(node.linear-1))
        push!(columns, column); push!(values, coefficient)
    end
    sparse(rows, columns, values, size(prepared.R_force)...)
end
RForceCollocated = collocated_source_operator(preparedOPT, sourceNodesOPT)
preparedOPTWithoutGamma = runOPTWithoutGamma ?
    flexOPT._prepared_with_matrices(preparedOPT, preparedOPT.A_unknown,
        preparedOPT.L_known, RForceCollocated; source_mode=:collocated) : nothing
propagationOPTWithoutGamma = runOPTWithoutGamma ? propagateLinearSystem(
    preparedOPTWithoutGamma, optSteps, dtOPT; sourceFull=sourceOPT,
    output_stride=optOutputStride, blowup_limit=1e8,
    solver_name="OPT3 without Γ redistribution",
) : nothing
runOPTWithoutGamma && @assert !propagationOPTWithoutGamma.stopped_early
uxOPT = propagationOPT.history[:, :, 1, :]
uzOPT = propagationOPT.history[:, :, 2, :]
optTimes = propagationOPT.times
uxOPTWithoutGamma = runOPTWithoutGamma ?
    propagationOPTWithoutGamma.history[:, :, 1, :] : nothing
uzOPTWithoutGamma = runOPTWithoutGamma ?
    propagationOPTWithoutGamma.history[:, :, 2, :] : nothing

# Higher-order runs are gated by ElasticLHSOperatorAudit.ipynb.
if runOPT4
    sourceOPT4 = zeros(Float64, preparedOPT4.NforcePoints,
        preparedOPT4.NForceField, length(sourceTimesOPT))
    for sourceNode in sourceNodesOPT
        sourceOPT4[sourceNode.linear, 2, :] .+=
            sourceNode.weight .* sourceScaleOPT .* waveletOPT
    end
    propagationOPT4 = propagateLinearSystem(preparedOPT4, optSteps, dtOPT;
        sourceFull=sourceOPT4, output_stride=optOutputStride,
        blowup_limit=1e8, solver_name="OPT4")
    @assert !propagationOPT4.stopped_early
    uxOPT4 = propagationOPT4.history[:, :, 1, :]
    uzOPT4 = propagationOPT4.history[:, :, 2, :]
    opt4Times = propagationOPT4.times
else
    propagationOPT4 = uxOPT4 = uzOPT4 = opt4Times = nothing
end

# OPT5 remains available, but is skipped in the quick benchmark by default.
if runOPT5
    sourceOPT5 = zeros(Float64, preparedOPT5.NforcePoints,
        preparedOPT5.NForceField, length(sourceTimesOPT))
    for sourceNode in sourceNodesOPT
        sourceOPT5[sourceNode.linear, 2, :] .+=
            sourceNode.weight .* sourceScaleOPT .* waveletOPT
    end
    propagationOPT5 = propagateLinearSystem(
        preparedOPT5, optSteps, dtOPT;
        sourceFull=sourceOPT5, output_stride=optOutputStride,
        blowup_limit=1e8, solver_name="OPT5",
    )
    @assert !propagationOPT5.stopped_early
    uxOPT5 = propagationOPT5.history[:, :, 1, :]
    uzOPT5 = propagationOPT5.history[:, :, 2, :]
    opt5Times = propagationOPT5.times
else
    propagationOPT5 = uxOPT5 = uzOPT5 = opt5Times = nothing
end

# Run the less stable delta-test conventional recipe last, so an eventual
# late stop cannot discard the expensive OPT5 result from this session.
sourceConvFD = zeros(Float64, preparedConvFD.NforcePoints,
    preparedConvFD.NForceField, length(sourceTimesOPT))
for sourceNode in sourceNodesOPT
    sourceConvFD[sourceNode.linear, 2, :] .+=
        sourceNode.weight .* sourceScaleOPT .* waveletOPT
end
# orderB=-1 is an interior conventional stencil, not a stable closure for
# the abruptly truncated outer rows. Stop before the earliest source tail
# can excite those rows unless the user explicitly requests the diagnostic.
convFDBoundaryGuard = 0.5 / sourceFrequency
convFDBoundarySafeTime = max(0.0, earliestPEdgeArrival - convFDBoundaryGuard)
convFDDuration = runConvFDBeyondBoundarySafeTime ? duration :
    min(duration, convFDBoundarySafeTime)
convFDSteps = ceil(Int, convFDDuration / dtOPT)
propagationConvFD = propagateLinearSystem(
    preparedConvFD, convFDSteps, dtOPT;
    sourceFull=sourceConvFD, output_stride=optOutputStride,
    blowup_limit=1e8,
    solver_name="convFD3 (orderB=-1)",
)
uxConvFD = propagationConvFD.history[:, :, 1, :]
uzConvFD = propagationConvFD.history[:, :, 2, :]
convFDTimes = propagationConvFD.times
propagationConvFD.stopped_early && @warn(
    "convFD3 stopped early; its partial history remains available",
    final_time=last(convFDTimes), requested_duration=duration)
@show convFDBoundarySafeTime convFDDuration propagationConvFD.stopped_early
@show dtOPT size(uzOPT) maximum(abs, uxOPT) maximum(abs, uzOPT)
@show size(uzConvFD) maximum(abs, uxConvFD) maximum(abs, uzConvFD)
runOPT5 && @show size(uzOPT5) maximum(abs, uxOPT5) maximum(abs, uzOPT5)

# ---- returned to step2_opt.jl and saved by myProduceOrLoad ----
Dict{String,Any}(
    "xOPT" => xOPT,
    "zOPT" => zOPT,
    "dxOPT" => dxOPT,
    "dtOPT" => dtOPT,
    "optSteps" => optSteps,
    "uxOPT" => uxOPT,
    "uzOPT" => uzOPT,
    "optTimes" => optTimes,
    "uxOPTWithoutGamma" => uxOPTWithoutGamma,
    "uzOPTWithoutGamma" => uzOPTWithoutGamma,
    "uxOPT4" => uxOPT4,
    "uzOPT4" => uzOPT4,
    "opt4Times" => opt4Times,
    "uxOPT5" => uxOPT5,
    "uzOPT5" => uzOPT5,
    "opt5Times" => opt5Times,
    "uxConvFD" => uxConvFD,
    "uzConvFD" => uzConvFD,
    "convFDTimes" => convFDTimes,
    "convFDBoundarySafeTime" => convFDBoundarySafeTime,
    "convFDDuration" => convFDDuration,
    "operatorTimings" => operatorTimings,
    "pcDiagnostics" => isnothing(propagationOPT.predictor_corrector) ? nothing :
        (iterations_per_step=propagationOPT.predictor_corrector.iterations_per_step,
         relative_residuals=propagationOPT.predictor_corrector.relative_residuals),
    "optTiming" => propagationOPT.timing,
)
