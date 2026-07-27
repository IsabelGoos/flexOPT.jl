
function getTopoViaGMT(params::Dict)
    @unpack precision, region = params
    resolution=String(chopprefix(precision, "@earth_relief_"))
    if DEFAULT_PLANET[] === :Mars
        G = gmtread(remotegrid("mars", res = resolution))
    elseif DEFAULT_PLANET[] === :Earth
        # do nothing
    else
        @show "not yet coded for your planet"
    end

    topo = GMT.grdcut(precision, region=region)
    return @strdict(topo)
end



function makeAdHocSeismicModel(modelArray::Array, ρMin, ρMax, VpvMin, VpvMax, VsvMin, VsvMax; VphMin=VpvMin, VphMax=VpvMax, VshMin=VsvMin, VshMax=VsvMax, QμMin = 0.0, QμMax=0.0, QκMin = 0.0, QκMax = 0.0, QμPowerMin = 0.0, QμPowerMax = 0.0, QκPowerMin = 0.0, QκPowerMax = 0.0, ηMin = 1.0, ηMax = 1.0)
    #if 'modelArray' is described in 0.0-1.0 range
    seismicModel = initiateSeismicModel(modelArray)
    @. seismicModel.ρ = ρMin + modelArray*(ρMax-ρMin)
    @. seismicModel.Vpv = VpvMin + modelArray*(VpvMax-VpvMin)
    @. seismicModel.Vph = VphMin + modelArray*(VphMax-VphMin)
    @. seismicModel.Vsv = VsvMin + modelArray*(VsvMax-VsvMin)
    @. seismicModel.Vsh = VshMin + modelArray*(VshMax-VshMin)
    @. seismicModel.Qμ  = QμMin  + modelArray*(QμMax - QμMin)
    @. seismicModel.Qκ  = QκMin  + modelArray*(QκMax - QκMin)
    @. seismicModel.QμPower = QμPowerMin + modelArray*(QμPowerMax-QμPowerMin)
    @. seismicModel.QκPower = QκPowerMin + modelArray*(QκPowerMax-QκPowerMin)
    @. seismicModel.η = ηMin + modelArray*(ηMax-ηMin)
    return seismicModel
end

function initiateSeismicModel(newArray::Array)
    seismicModel=(ρ=zeros(Float64,size(newArray)...),Vpv=zeros(Float64,size(newArray)...),Vph=zeros(Float64,size(newArray)...),Vsv=zeros(Float64,size(newArray)...),Vsh=zeros(Float64,size(newArray)...),Qμ=zeros(Float64,size(newArray)...),Qκ=zeros(Float64,size(newArray)...),QμPower=zeros(Float64,size(newArray)...),QκPower=zeros(Float64,size(newArray)...),η=zeros(Float64,size(newArray)...))
    return seismicModel
end


function getParamsWithoutTopo(allGridsInGeoPoints,effectiveRadii;NradiusNodes=500,eps=10.0,VpWater=1.5,ρWater=1.0,VpAir=0.314,ρAir=0.001,hasAirModel=false)

    
    #@enum Couche Graine Noyau Manteau Océane Atmosphère Ionosphère Dehors

    #pointCharacter = Array{Couche}(undef,size(allGridsInGeoPoints)...)


    # some 'chelou' parameters that I need to control

    #NradiusNodes =500 # I don't know how to make this number reasonable for interpolation
    #eps = 10.0 # in metre "below" option should be enough but who knows
    #VpWater = 1500.0 # m/s
    #ρWater = 1000.0 # kg/m3
    #VpAir = 314.0
    #ρAir = 1.0



    #precision = GMTprecision(precisionInKm)  # this should be in Km


    #seismicModel=(ρ=zeros(Float64,size(allGridsInGeoPoints)...),Vpv=zeros(Float64,size(allGridsInGeoPoints)...),Vph=zeros(Float64,size(allGridsInGeoPoints)...),Vsv=zeros(Float64,size(allGridsInGeoPoints)...),Vsh=zeros(Float64,size(allGridsInGeoPoints)...),Qμ=zeros(Float64,size(allGridsInGeoPoints)...),Qκ=zeros(Float64,size(allGridsInGeoPoints)...),QμPower=zeros(Float64,size(allGridsInGeoPoints)...),QκPower=zeros(Float64,size(allGridsInGeoPoints)...),η=zeros(Float64,size(allGridsInGeoPoints)...))


    seismicModel = initiateSeismicModel(allGridsInGeoPoints)

    # get the extremeties in radius

    
    tmpNradiusNodes = NradiusNodes
    if NradiusNodes === 1
        tmpNradiusNodes = 1
    end

    ΔradiusIncrementInKm = (maximum(effectiveRadii)-minimum(effectiveRadii))/(tmpNradiusNodes-1) *1.e-3
    linearRadiiInKm =(collect(1:1:NradiusNodes) .- 1)*ΔradiusIncrementInKm .+ minimum(effectiveRadii)*1.e-3

    push!(linearRadiiInKm, planet1D.my1DDSMmodel.averagedPlanetRadiusInKilometer) 
    newRadii,params=planet1D.compute1DseismicParamtersFromPolynomialCoefficientsWithGivenRadiiArray(planet1D.my1DDSMmodel,linearRadiiInKm,"below")


    # get the extremeties in lat and lon


    for i in CartesianIndices(allGridsInGeoPoints)
        tmpPoint = allGridsInGeoPoints[i]
       

        if hasAirModel === false
            if 0.0 < tmpPoint.alt
                seismicModel.ρ[i]=ρAir
                seismicModel.Vpv[i]=VpAir
                seismicModel.Vph[i]=VpAir
                seismicModel.Vsv[i]=0.0
                seismicModel.Vsh[i]=0.0
            end
        end

    end

    seismicModel = interpolate_params(params, newRadii, effectiveRadii)



    return seismicModel

end



"""
    getParamsAndTopo(points, effectiveRadii, precisionInKm;
        velocity_model=:planet1D, ...)

Build the seismic parameters and apply GMT topography. `velocity_model` may be:

- `:planet1D`: retain the radial planet1D velocities (default);
- `:NIED` (or the `:NIED2023` alias): replace reliable in-domain `Vp` and `Vs` with the NIED ALJ
  three-dimensional model, while retaining planet1D density, attenuation and
  other parameters.

For NIED, `nied_source` may be either a downloadable ZIP URL or a local ZIP
path. Parsed data are cached and the official slowness interpolation convention
is used.
"""
function getParamsAndTopo(
    allGridsInGeoPoints,
    effectiveRadii,
    precisionInKm::Float64;
    NradiusNodes=500,
    eps=10.0,
    VpWater=1.5,
    ρWater=1.0,
    VpAir=0.314,
    ρAir=0.001,
    hasAirModel=false,
    velocity_model::Symbol=:planet1D,
    nied_source::AbstractString=DEFAULT_NIED_VELOCITY_SOURCE[],
    nied_cache_dir::AbstractString=DEFAULT_NIED_VELOCITY_CACHE[],
    nied_confidence_max::Union{Nothing,Real}=0.8,
    nied_outside::Symbol=:planet1D,
    nied_low_confidence::Symbol=:planet1D,
    nied_force_download::Bool=false,
    nied_force_parse::Bool=false,
)
    velocity_model in (:planet1D, :NIED, :NIED2023) ||
        throw(ArgumentError(
            "velocity_model must be :planet1D, :NIED, or :NIED2023",
        ))
    use_nied = velocity_model in (:NIED, :NIED2023)

    
    #@enum Couche Graine Noyau Manteau Océane Atmosphère Ionosphère Dehors

    #pointCharacter = Array{Couche}(undef,size(allGridsInGeoPoints)...)


    # some 'chelou' parameters that I need to control

    #NradiusNodes =500 # I don't know how to make this number reasonable for interpolation
    #eps = 10.0 # in metre "below" option should be enough but who knows
    #VpWater = 1500.0 # m/s
    #ρWater = 1000.0 # kg/m3
    #VpAir = 314.0
    #ρAir = 1.0



    precision = GMTprecision(precisionInKm)  # this should be in Km


    #seismicModel=(ρ=zeros(Float64,size(allGridsInGeoPoints)...),Vpv=zeros(Float64,size(allGridsInGeoPoints)...),Vph=zeros(Float64,size(allGridsInGeoPoints)...),Vsv=zeros(Float64,size(allGridsInGeoPoints)...),Vsh=zeros(Float64,size(allGridsInGeoPoints)...),Qμ=zeros(Float64,size(allGridsInGeoPoints)...),Qκ=zeros(Float64,size(allGridsInGeoPoints)...),QμPower=zeros(Float64,size(allGridsInGeoPoints)...),QκPower=zeros(Float64,size(allGridsInGeoPoints)...),η=zeros(Float64,size(allGridsInGeoPoints)...))


    # get the extremeties in radius

    
    tmpNradiusNodes = NradiusNodes
    if NradiusNodes === 1
        tmpNradiusNodes = 1
    end

    ΔradiusIncrementInKm = (maximum(effectiveRadii)-minimum(effectiveRadii))/(tmpNradiusNodes-1) *1.e-3
    linearRadiiInKm =(collect(1:1:NradiusNodes) .- 1)*ΔradiusIncrementInKm .+ minimum(effectiveRadii)*1.e-3

    push!(linearRadiiInKm, planet1D.my1DDSMmodel.averagedPlanetRadiusInKilometer) 
    newRadii,params=planet1D.compute1DseismicParamtersFromPolynomialCoefficientsWithGivenRadiiArray(planet1D.my1DDSMmodel,linearRadiiInKm,"below")


    # get the extremeties in lat and lon

    lats = [p.lat for p in allGridsInGeoPoints]
    lons = [p.lon for p in allGridsInGeoPoints]

    lat_min, lat_max = extrema(lats)
    lon_min, lon_max = extrema(lons)
    if lat_min === lat_max
        if lat_min > 1.0
            lat_min = lat_min - 0.5
        else
            lat_max = lat_max + 0.5
        end
    end

    if lon_min === lon_max
        if lon_min > 1.0
            lon_min = lon_min - 0.5
        else
            lon_max = lon_max + 0.5
        end
    end


    region = [lon_min,lon_max,lat_min,lat_max]

    paramsForGMT = @strdict precision region
    topo_out=myProduceOrLoad(getTopoViaGMT,paramsForGMT,"topoViaGMT")
    topo=topo_out["topo"]
    
    # below is a very strange behaviour from GMT and I need to clarify soon

    if size(topo.x)[1] !== size(topo.z)[2]
        x = (topo.x[1:end-1] .+ topo.x[2:end]) ./ 2
        y = (topo.y[1:end-1] .+ topo.y[2:end]) ./ 2
        z = topo.z
    else
        x=topo.x
        y=topo.y
        z=topo.z
    end

    #@show size(x), size(y), size(z)


    #topoInterpolater = interpolate((topo.y,topo.x),topo.z',Gridded(Linear())) #LinearInterpolation((x, y), z; extrapolation_bc=Flat())

    topoInterpolater=LinearInterpolation((x, y), z'; extrapolation_bc=Flat())
    #topoInterpolater = interpolate((topo.y,topo.x),topo,Gridded(Linear()))


    

    # Cache the topography because it is also needed after the radial
    # interpolation. Each threaded iteration reads and writes distinct cells.
    topography = similar(effectiveRadii, Float64)
    Threads.@threads for i in eachindex(allGridsInGeoPoints, effectiveRadii, topography)
        tmpPoint = allGridsInGeoPoints[i]
        topographyHere = topoInterpolater(tmpPoint.lon, tmpPoint.lat)
        topography[i] = topographyHere

        if 0.0 < tmpPoint.alt <= topographyHere
            # it might be very time-consuming if we do this for 3D Cartesian points ...
            effectiveRadii[i]=planet1D.my1DDSMmodel.averagedPlanetRadiusInKilometer*1.e3 - eps
     
        end
    end

    seismicModel = interpolate_params(params, newRadii, effectiveRadii)

    if use_nied
        nied_fields = apply_nied_velocity_model!(
            seismicModel,
            allGridsInGeoPoints;
            source=nied_source,
            cache_dir=nied_cache_dir,
            confidence_max=nied_confidence_max,
            outside=nied_outside,
            low_confidence=nied_low_confidence,
            force_download=nied_force_download,
            force_parse=nied_force_parse,
        )
        seismicModel = merge(
            seismicModel,
            nied_fields,
            (velocity_model=velocity_model,),
        )
    else
        seismicModel = merge(seismicModel, (velocity_model=:planet1D,))
    end

    Threads.@threads for i in eachindex(allGridsInGeoPoints, topography)
        tmpPoint = allGridsInGeoPoints[i]
        topographyHere = topography[i]

        # Reapply surface materials after interpolate_params, which creates a
        # fresh NamedTuple and would otherwise overwrite the earlier air mask.
        if !hasAirModel && topographyHere < tmpPoint.alt
            seismicModel.ρ[i] = ρAir
            seismicModel.Vpv[i] = VpAir
            seismicModel.Vph[i] = VpAir
            seismicModel.Vsv[i] = 0.0
            seismicModel.Vsh[i] = 0.0
            if use_nied
                seismicModel.nied_P_mask[i] = false
                seismicModel.nied_S_mask[i] = false
                seismicModel.nied_mask[i] = false
            end
        elseif topographyHere <= tmpPoint.alt <= 0.0
            seismicModel.ρ[i] = ρWater
            seismicModel.Vpv[i] = VpWater
            seismicModel.Vph[i] = VpWater
            seismicModel.Vsv[i] = 0.0
            seismicModel.Vsh[i] = 0.0
            if use_nied
                seismicModel.nied_P_mask[i] = false
                seismicModel.nied_S_mask[i] = false
                seismicModel.nied_mask[i] = false
            end
        end
    end


    return seismicModel

end

function interpolate_params(params, newRadii, effectiveRadii)
    (; (name => begin
        itp = LinearInterpolation(newRadii .* 1e3,
                                  getfield(params, name);
                                  extrapolation_bc = Flat())
        itp.(effectiveRadii)
    end for name in fieldnames(typeof(params)))...)
end

function GMTprecision(requiredResolutionInKmÀPeuPrès::Float64)
    if requiredResolutionInKmÀPeuPrès > 55.0
        return "@earth_relief_30m"
    elseif  18.0 < requiredResolutionInKmÀPeuPrès <= 55.0
        return "@earth_relief_10m"
    elseif  9.0 < requiredResolutionInKmÀPeuPrès <= 18.0
        return "@earth_relief_05m"
    elseif 3.6 < requiredResolutionInKmÀPeuPrès <= 9.0
        return "@earth_relief_02m"
    elseif 1.8 < requiredResolutionInKmÀPeuPrès <= 3.6
        return "@earth_relief_01m"
    elseif 0.45 < requiredResolutionInKmÀPeuPrès <= 1.8
        return "@earth_relief_15s"
    elseif 0.30 < requiredResolutionInKmÀPeuPrès <= 0.45
        return "@earth_relief_10s"
    elseif 0.15 < requiredResolutionInKmÀPeuPrès <= 0.30
        return "@earth_relief_05s"
    elseif 0.09 < requiredResolutionInKmÀPeuPrès <= 0.15
        return "@earth_relief_03s"
    elseif 0.06 < requiredResolutionInKmÀPeuPrès <= 0.09
        return "@earth_relief_02s"
    else
        return "@earth_relief_01s"
        # after this
        """
        # Load your own local DEM (GeoTIFF, etc.)
        topo = GMT.read("Copernicus_DSM_30m_Europe.tif")

        # Cut a region and maybe downsample
        region = [-5, 0, 43, 47]
        topo_sub = GMT.grdcut(topo, region=region)

        # Optional resampling (to ~10 m)
        topo_10m = GMT.grdsample(topo_sub, inc="10m")

        """
    end
end
