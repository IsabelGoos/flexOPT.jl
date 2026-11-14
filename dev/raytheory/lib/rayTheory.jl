# dev/raytheory/lib/rayTheory.jl — loads the ray-theory core (stdlib only).
# Layers:  rayHamiltonian.jl (1: propagation), rayChristoffel.jl (2: dispersion
# relation / anisotropic branches), interfaceSlowness in both (3: interface
# matching); phases = event lists (rayPhases.jl); fast 1-D path rayTauP.jl.
include("rayModels.jl")
include("rayPhases.jl")
include("rayHamiltonian.jl")
include("rayTracer.jl")
include("rayChristoffel.jl")
include("rayTauP.jl")
include("rayShooting.jl")
include("rayGrid.jl")
