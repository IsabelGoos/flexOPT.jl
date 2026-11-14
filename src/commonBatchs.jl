module commonBatchs

    # batch files
    using Symbolics
    using StaticArrays
    using Interpolations
    using DIVAnd
  
    include("batchFiles/batchDrWatson.jl")
    include("batchFiles/batchImages.jl")
    include("batchFiles/batchUseful.jl")
    include("batchFiles/batchStagYY.jl")
    include("batchFiles/batchSymbolics.jl")
    include("batchFiles/batchInterpolations.jl")
    export dmsToDecimal
    export SVector, @SVector
    export interpolateField,getCartesianField
    export myProduceOrLoad,lazyProduceOrLoad, @strdict, @unpack
    export cacheKeyExtras
    export getColorPalette, regenerataionColorMap, color2Float
    export safeget, _paramget
    export var2string,reinterpolateArrayMembers,expandVectors,string_as_varname
    export car2vec,carDropDim,carAddDim,vec2car,car2svec,svec2car
    export flatten,deep_flatten,is_all_less_than_or_equal,distance2_point_to_box
    export myInv
    # batchStagYY
    export myListDir,readStagYYFiles,readStagYYFilesAverage
    export quarterDiskExtrapolationRawGrid!, quarterDiskExtrapolation
    export makeZOverALayers, layeredZOverA, getZOverAField
    export electronDensity, radialAverageAnomaly
    # symbolify
    export x,y,z,t
    export ∂x,∂y,∂z,∂t
    export ∂x²,∂y²,∂z²,∂t² 
    export Num2Float64,usefulPartials,myCoeff,mySolvefor,mySimplify,integrateTaylorPolynomials
    export mySimplify

end
