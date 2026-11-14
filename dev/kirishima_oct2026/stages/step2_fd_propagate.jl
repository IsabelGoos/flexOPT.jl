# dev/kirishima_oct2026/stages/step2_fd_propagate.jl
# Ported from notebooks/Elastic2DtopoBenchmark_oct2026.ipynb, cell step2-*-fd-propagation.
# Runs ONLY on a cache miss of step2_fd.jl; its last expression is the Dict saved.
fdConfig = ElasticThreePointConfig2D(
    pointsInSpace=3, pointsInTime=3, supplementaryOrder=2, cfl=0.38,
)
fd = prepare_elastic_wave_2d(
    model, (dx, dx);
    material_mask=material,
    boundary_conditions=bcFD,
    config=fdConfig,
)
fdCoordinates = elastic_wave_coordinates(x, z, fd)
fdPadding = fd.padding
sourcePhysical = CartesianIndex(
    argmin(abs.(x .- sourcePosition.x)),
    argmin(abs.(z .- sourcePosition.z)),
)
sourceFD = sourcePhysical + CartesianIndex(Tuple(fdPadding[1, :]))
sourceSpecsFD = spatial_source_specs(fdCoordinates.x, fdCoordinates.z)
sourceNodesFD = [(index=CartesianIndex(
    argmin(abs.(fdCoordinates.x .- s.x)),
    argmin(abs.(fdCoordinates.z .- s.z))), weight=s.weight)
    for s in sourceSpecsFD]
fdSteps = ceil(Int, duration / fd.dt)
fdOutputStride = max(1, round(Int, outputSampling / fd.dt))
fdFramesX = Matrix{Float32}[copy(fd.ux)]
fdFramesZ = Matrix{Float32}[copy(fd.uz)]
fdTimes = Float64[0.0]
for step in 1:fdSteps
    step_elastic_wave_2d!(fd)
    for sourceNode in sourceNodesFD
        elasticWave2D.add_ricker_source!(
            fd, sourceNode.index;
            f0=sourceFrequency, t0=sourceDelay,
            amplitude=sourceForce * sourceNode.weight,
            component=:z, source_kind=:force,
        )
    end
    if step % fdOutputStride == 0 || step == fdSteps
        push!(fdFramesX, copy(fd.ux))
        push!(fdFramesZ, copy(fd.uz)); push!(fdTimes, fd.time)
    end
end
uxFD = cat(fdFramesX...; dims=3)
uzFD = cat(fdFramesZ...; dims=3)
@show fd.dt size(uzFD) maximum(abs, uxFD) maximum(abs, uzFD)

# ---- returned to step2_fd.jl and saved by myProduceOrLoad ----
Dict{String,Any}(
    "uxFD" => uxFD,
    "uzFD" => uzFD,
    "fdTimes" => fdTimes,
    "fdCoordinates" => fdCoordinates,
    "fdPadding" => fdPadding,
    "fdDt" => Float64(fd.dt),
)
