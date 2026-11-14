# dev/raytheory/reference/writeND.jl
# Write a flexOPT DSM .psv.model (isotropic Voigt reduction, as readDSMModel)
# as a TauP .nd file: every cubic zone sampled every 10 km (TauP interpolates
# linearly in between; effect on times < 1e-3 s), with the "mantle",
# "outer-core", "inner-core" labels TauP needs.
#   julia dev/raytheory/reference/writeND.jl dataInput/model/earth/PREM.psv.model premDSM.nd
include(joinpath(@__DIR__, "..", "lib", "rayTheory.jl"))
function writeND(m, path; perkm=0.1)
    open(path, "w") do io
        for k in nzones(m):-1:1
            z = m.zones[k]
            b = m.boundaries[k+1]
            b.name === :moho && println(io, "mantle")
            b.name === :cmb && println(io, "outer-core")
            b.name === :icb && println(io, "inner-core")
            n = max(2, ceil(Int, (z.top - z.bottom) * perkm) + 1)
            for ζ in range(z.top, z.bottom; length=n)
                @printf(io, "%10.4f %9.5f %9.5f %9.5f\n", m.surface - ζ, zoneVelocity(m, k, ζ, :P),
                    zoneVelocity(m, k, ζ, :S), polyval4(z.rho, ζ / m.L))
            end
        end
    end
end
if abspath(PROGRAM_FILE) == @__FILE__
    writeND(readDSMModel(ARGS[1]), ARGS[2])
    println("written ", ARGS[2])
end
