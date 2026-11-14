# dev/kirishima_oct2026/stages/step3_model.jl
# Analogous to step2_model.jl, for the Step 3 (flat free surface) box.
# vp/vs/rho and all derived safety checks are restricted to the solid
# (z<=0) part of the box: the void cap above the flat surface carries
# whatever getParamsAndTopo happened to return there, and it is masked
# out (never used) by both FD3's material_mask and OPT's solidOPT mask.
caseLabel = step3Case === :homog ? "homogeneous" : "heterogeneous"
caseTitle = step3Case === :homog ? "Homogeneous" : "Heterogeneous"
if step3Case === :homog
    # Homogeneous sub-case: a single representative medium (the box's
    # own solid-part median), compared against SPECFEM2D exactly as in
    # HomogeneousElastic2DBenchmark_FreeSurface.ipynb, but using
    # Kirishima-scale constants.
    vpField = fill(vp0, nx, nz)
    vsField = fill(vs0, nx, nz)
    rhoField = fill(rho0, nx, nz)
elseif step3Case === :heter
    # Heterogeneous sub-case: the real Kirishima NIED velocity/density
    # field over the Step 3 box, reused verbatim below the flat surface
    # (km/s, g/cm^3 -> m/s, kg/m^3); the void half (above z=0) is never
    # read back downstream.
    vpField = seismicModel2DStep3.Vpv .* 1e3
    vsField = seismicModel2DStep3.Vsv .* 1e3
    rhoField = seismicModel2DStep3.ρ .* 1e3
else
    error("unknown step3Case=$step3Case")
end
model = (ρ=rhoField ./ 1e3, Vpv=vpField ./ 1e3, Vsv=vsField ./ 1e3)
vpMaximum, vsMinimum = maximum(vpField[:, solidMaskStep3]), minimum(vsField[:, solidMaskStep3])
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
# Content hash of the (vp, vs, rho) fields: part of every Step-3 cache key,
# so a changed model can never be confused with a cached one.
modelHash = string(hash((vpField, vsField, rhoField)); base=16)
@show caseLabel modelHash
