# dev/spherical1Delastic/stages/configs.jl
# Ported from notebooks/SphericalTI_SH_oct2026.ipynb, cell sh-step2-03 (configs part).
# Used by every step from Step 2 on.
shConfigs = NamedTuple[]
push!(shConfigs, (name="convFD3", orderBtime=-1, orderBspace=-1, pointsInSpace=3, pointsInTime=3, supplementaryOrder=0,
    fieldItpl=(ptsSpace=1, ptsTime=1, offsetSpace=1.0, offsetTime=1, YorderBspace=-1, YorderBtime=-1),
    materItpl=(ptsSpace=1, ptsTime=1, offsetSpace=1.0, offsetTime=1, YorderBspace=-1, YorderBtime=-1)))
push!(shConfigs, (name="OPT3", orderBtime=1, orderBspace=1, pointsInSpace=3, pointsInTime=3, supplementaryOrder=2,
    fieldItpl=(ptsSpace=1, ptsTime=1, offsetSpace=1.0, offsetTime=1, YorderBspace=-1, YorderBtime=-1),
    materItpl=(ptsSpace=1, ptsTime=1, offsetSpace=1.0, offsetTime=1, YorderBspace=-1, YorderBtime=-1)))
push!(shConfigs, (name="OPT4", orderBtime=1, orderBspace=1, pointsInSpace=4, pointsInTime=1, supplementaryOrder=2,
    fieldItpl=(ptsSpace=1, ptsTime=1, offsetSpace=1.5, offsetTime=1, YorderBspace=-1, YorderBtime=1),
    materItpl=(ptsSpace=1, ptsTime=1, offsetSpace=1.5, offsetTime=1, YorderBspace=-1, YorderBtime=1)))
push!(shConfigs, (name="OPT5", orderBtime=1, orderBspace=1, pointsInSpace=5, pointsInTime=1, supplementaryOrder=2,
    fieldItpl=(ptsSpace=1, ptsTime=1, offsetSpace=2.0, offsetTime=1, YorderBspace=-1, YorderBtime=1),
    materItpl=(ptsSpace=1, ptsTime=1, offsetSpace=2.0, offsetTime=1, YorderBspace=-1, YorderBtime=1)))
