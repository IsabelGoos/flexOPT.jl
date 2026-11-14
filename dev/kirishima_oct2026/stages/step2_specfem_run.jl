# dev/kirishima_oct2026/stages/step2_specfem_run.jl
# Ported from notebooks/Elastic2DtopoBenchmark_oct2026.ipynb, cell step2-*-specfem (run part).
# Runs ONLY on a cache miss of step2_specfem.jl; its last expression is the Dict saved.
dtSPECFEM = dtOPT
# caseDirectoryTag keeps different steps' SPECFEM2D working directories
# apart (Step 2: "no_free_surface"; Step 3: "flat_free_surface", set in
# step3_config.jl) without needing a per-step copy of this file.
caseDirectoryTag = @isdefined(caseDirectoryTag) ? caseDirectoryTag : "no_free_surface"
caseDirectory = joinpath(flexopt_root, "data", "specfem2d_benchmarks",
    "$(caseLabel)_kirishima_$(caseDirectoryTag)")
solidZ = applyFreeSurface ? findall(z .<= 0.0) : collect(eachindex(z))
zSPECFEM = z[solidZ]
vp = vpField[:, solidZ]
vs = vsField[:, solidZ]
rho = rhoField[:, solidZ]
# Refining FD/OPT must not silently refine the SPECFEM spectral mesh.
specfemNxElements = round(Int, (last(x)-first(x)) / (4referenceDx))
specfemNzElements = round(Int, (last(zSPECFEM)-first(zSPECFEM)) / (4referenceDx))
specfemCase = prepare_specfem2d_case(
    caseDirectory, x, zSPECFEM, vp, vs, rho, fill(last(zSPECFEM), nx);
    source=sourcePosition, sources=sourceSpecsReference, receivers=receiverX,
    duration=duration, dt=dtSPECFEM, f0=sourceFrequency,
    source_factor=sourceForce, source_angle=0.0,
    source_time_function=rickerSource,
    nx_elements=specfemNxElements, nz_elements=specfemNzElements,
    receiver_z=receiverZ, free_surface=applyFreeSurface,
    snapshot_interval_steps=max(1, round(Int, outputSampling / dtSPECFEM)),
    snapshot_image_type=5, # vertical velocity, matching the seismograms
    output_wavefield_dumps=true,
    wavefield_dump_type=1, # displacement vector (uₓ, u_z)
    binary_wavefield_dumps=true,
)
if runSPECFEM2D
    specfemRun = run_specfem2d_case(specfemCase.case_directory)
else
    specfemRun = (output=specfemCase.output,)
end
verticalFiles = find_specfem2d_traces(specfemRun.output; component=:z)
@assert length(verticalFiles) == length(receiverX)
specfemTracesZ = [read_specfem2d_trace(file;
    time_shift=specfemCase.time_axis_shift) for file in verticalFiles]
horizontalFiles = find_specfem2d_traces(specfemRun.output; component=:x)
@assert length(horizontalFiles) == length(receiverX)
specfemTracesX = [read_specfem2d_trace(file;
    time_shift=specfemCase.time_axis_shift) for file in horizontalFiles]

# (the full wavefield dumps stay on disk in specfemRun.output; read them with
#  read_specfem2d_wavefield_dumps(specfemRun.output; ...) when needed)
specfemSnapshotVideo = make_specfem2d_snapshot_video(
    specfemRun.output;
    output_path=joinpath(specfemRun.output, "SPECFEM2D_vertical_velocity.mp4"),
    framerate=20,
)
solverLogText = read(joinpath(specfemCase.case_directory, "solver.log"), String)
specfemCFLMatch = match(r"Max CFL stability condition[^=]*=\s*([0-9.Ee+-]+)",
    solverLogText)
specfemCFL = isnothing(specfemCFLMatch) ? missing :
    parse(Float64, specfemCFLMatch.captures[1])

# ---- returned to step2_specfem.jl and saved by myProduceOrLoad ----
Dict{String,Any}(
    "specfemTracesX" => specfemTracesX,
    "specfemTracesZ" => specfemTracesZ,
    "specfemCase" => (case_directory=specfemCase.case_directory,
        time_axis_shift=specfemCase.time_axis_shift,
        nx_elements=specfemCase.nx_elements, nz_elements=specfemCase.nz_elements,
        output=specfemRun.output),
    "specfemCFL" => specfemCFL,
    "specfemSnapshotVideoPath" => specfemSnapshotVideo.path,
    "zSPECFEM" => zSPECFEM,
)
