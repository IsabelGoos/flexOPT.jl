# dev/kirishima_oct2026/stages/shared_site.jl
# Ported from notebooks/Elastic2DtopoBenchmark_oct2026.ipynb, cells step1-box-model / step1-topography-plots.
# Geographic reference shared by Step 1 and Step 2 (so Step 2 does not
# need to load the big 3-D Step-1 model).
kirishimaSummit = GeoPoint(dmsToDecimal(31, 56, 03), dmsToDecimal(130, 51, 42))
p0 = kirishimaSummit

niedVelocitySource = DEFAULT_NIED_VELOCITY_SOURCE[]
velocitySourceKey = nied_velocity_source_key(niedVelocitySource)

# getParamsAndTopo uses ρAir=0.001 by default; 0.01 cleanly separates
# rock from air/void (same threshold as SimuKirishima.ipynb).
air_density_cutoff = 0.01
