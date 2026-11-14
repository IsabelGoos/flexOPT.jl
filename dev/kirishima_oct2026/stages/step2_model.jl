# dev/kirishima_oct2026/stages/step2_model.jl
# Ported from notebooks/Elastic2DtopoBenchmark_oct2026.ipynb, cells step2-homog-model / step2-heter-model.
# One file for both sub-cases; run.jl sets step2Case = :homog or :heter.
caseLabel = step2Case === :homog ? "homogeneous" : "heterogeneous"
caseTitle = step2Case === :homog ? "Homogeneous" : "Heterogeneous"
if step2Case === :homog
    # Homogeneous sub-case: a single representative medium (the box median),
    # so FD3/OPT3/OPT4/OPT5 are compared against SPECFEM2D exactly as in
    # HomogeneousElastic2DBenchmark.ipynb, but using Kirishima-scale constants.
    vpField = fill(vp0, nx, nz)
    vsField = fill(vs0, nx, nz)
    rhoField = fill(rho0, nx, nz)
elseif step2Case === :heter
    # Heterogeneous sub-case: the real Kirishima NIED velocity/density field
    # over the Step 2 box, reused verbatim (km/s, g/cm^3 -> m/s, kg/m^3).
    vpField = seismicModel2DStep2.Vpv .* 1e3
    vsField = seismicModel2DStep2.Vsv .* 1e3
    rhoField = seismicModel2DStep2.ρ .* 1e3
else
    error("unknown step2Case=$step2Case")
end
model = (ρ=rhoField ./ 1e3, Vpv=vpField ./ 1e3, Vsv=vsField ./ 1e3)
vpMaximum, vsMinimum = maximum(vpField), minimum(vsField)
pointsPerSWavelength = (
    FD3=vsMinimum / (sourceFrequency * dx),
    OPT3=vsMinimum / (sourceFrequency * (dx / optSpatialRefinement)),
    OPT5=vsMinimum / (sourceFrequency * (dx / optSpatialRefinement)),
)
earliestPEdgeArrival = significantSourceStart +
    distanceToNearestPhysicalEdge / vpMaximum
@assert duration < earliestPEdgeArrival
@assert minimum(values(pointsPerSWavelength)) >=
    minimumAcceptedPointsPerSWavelength
@show vpMaximum vsMinimum pointsPerSWavelength earliestPEdgeArrival
# Content hash of the (vp, vs, rho) fields: part of every Step-2 cache key,
# so a changed model can never be confused with a cached one.
modelHash = string(hash((vpField, vsField, rhoField)); base=16)
@show caseLabel modelHash
