# some useful structures and functions to deal with planetary 3D coordinates

module GeoPoints
    # be careful, this module needs planet1D


    using ..planet1D
    using ..commonBatchs
    using Geodesy, StaticArrays,LinearAlgebra,GMT, UnPack
    using CSV, Downloads, JLD2, SHA, ZipFile
    using Interpolations
    using FileIO,CairoMakie
    using Colors, Images
    using Base.Threads

    import Base: +,-,/,*


    



    # specific functions: former GeoPoints
    include("GeoPoints/GeoPoints.jl")
    include("GeoPoints/niedVelocityModel.jl")
    include("GeoPoints/getSeismicParamTopo.jl")
    
    # specific functions: former imageReder
    include("imageReader/imageReader.jl")
    include("imageReader/modelMaker.jl")


    # functions GeoPoints
    export realToVirtualEarthCoords
    export DEFAULT_PLANET,DEFAULT_ELLIPSOID, planet_ellipsoid, set_default_planet!, localCoord2D, localCoord3D
    export p_ECEF_to_local, p_local_to_ECEF, p_ECEF_to_local2D
    export GeoPoint, +, -, /, *
    export GeoPointSet,GeoPoints_to_local
    export rotation_between_geopoints, transform_geopoints
    export theta_phi_shift, effective_cartesian_coordinates
    export PathProfile2D, creationPaths
    #export effectiveRadius, makeLocalCoordinateUnitVectors, makeLocalCoordinateUnitVectors, makeLocalCoordinateUnitVectors
    export constructLocalBox
    export getParamsAndTopo,getParamsWithoutTopo #, interpolate_params, GMTprecision
    export DEFAULT_NIED_VELOCITY_PAGE, DEFAULT_NIED_VELOCITY_SOURCE
    export DEFAULT_NIED_VELOCITY_CACHE, set_default_nied_velocity_source!
    export NIEDVelocityModel, load_nied_velocity_model
    export nied_velocity_source_key
    export makeAdHocSeismicModel,initiateSeismicModel # soon we need to make a structure of seismicModel (1D-3D)
    
    # functions imageReader
    export defineModel
    export read2DimageModel
    export adjustArray
end
