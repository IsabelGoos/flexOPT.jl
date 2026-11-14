# application/Elastic_Spherical1D_Freq/run.jl -- 1-D spherical (TI, anelastic) elastic waves,
# frequency domain, DSM-style (SH toroidal, P-SV spheroidal, fluid core in Q = P/ω).
#
# Thin entry point: the code lives in dev/ while it is being benchmarked.
#   task "forward"  (default) -> dev/spherical1Delastic/run.jl        (seismograms, eigenfrequencies)
#   task "partials"           -> dev/spherical1DelasticKernel/run.jl  (strain Green's tensors, Fréchet kernels)
#
#   julia --project=. application/Elastic_Spherical1D_Freq/run.jl forward p5      # PREM 3-component seismograms
#   julia --project=. application/Elastic_Spherical1D_Freq/run.jl forward step3b  # PREM SH seismograms
#   julia --project=. application/Elastic_Spherical1D_Freq/run.jl partials k1     # kernel checks
#   julia --project=. application/Elastic_Spherical1D_Freq/run.jl partials k3     # Born vs brute force
# Every other word/option is passed unchanged to the dev run.jl (steps, --redo=..., --full).

let args = copy(@isdefined(APP_ARGS) ? APP_ARGS : ARGS)
    task = (!isempty(args) && args[1] in ("forward", "partials")) ? popfirst!(args) : "forward"
    global DEV_ARGS = args
    sub = task == "partials" ? "spherical1DelasticKernel" : "spherical1Delastic"
    println("Elastic_Spherical1D_Freq: task=$task -> dev/$sub ", args)
    Base.include(Main, joinpath(@__DIR__, "..", "..", "dev", sub, "run.jl"))
end
