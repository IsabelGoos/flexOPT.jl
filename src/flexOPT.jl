module flexOPT

    using Symbolics,UnPack,Tullio
    using Dates
    using LinearAlgebra
    using KernelAbstractions
    using ..commonBatchs
    using Adapt

    # GPU backend (batchFiles/batchGPU.jl should be called in Main)
    backend = Main.backend
    makeGPUarray = Main.makeGPUarray


    # wrapping functions

    include("../src/motorsOPT/makeOPTrecette.jl")
    export TaylorOptions,makeOPTsemiSymbolic
    
    # deriving symbolic and semi-symbolic OPT operators
    include("../src/compactSymbolicFunctions/compactFunctionsArray.jl")
    include("../src/compactSymbolicFunctions/BsplineHelpers.jl")
    include("../src/CompactSymbolicFunctions/TaylorExpansionHelpers.jl")
    include("../src/CompactSymbolicFunctions/integralWYYKK.jl")


    # famous equations etc.

    include("../src/motorsOPT/others.jl")
    export timeDimensionString
    include("../src/motorsOPT/famousEquations.jl")
    include("../src/motorsOPT/famousSourceFunctions.jl")
    include("../src/motorsOPT/famousBoundaryConditions.jl")

    # semi-symbolics operators to fully numerical operators
    #include("fullyNumericalOPT/makeCostFunctions.jl")
    include("fullyNumericalOPT/makeCostFunctions.jl")
    include("fullyNumericalOPT/others.jl")
    export numericalOperatorConstruction,getModelPoints

    include("../src/numSolvers/linearSolution.jl")
    include("../src/numSolvers/timeMarchingSchemes.jl")
    include("../src/numSolvers/diffTools.jl")
    export prepareNumericalOperators,timeMarchingSchemePrepared,prepareLinearSystem,evaluateLinearSystem!,evaluateLinearSystem,timeMarchingSchemeLinear


    export quasiNumericalOperatorConstruction,constructingNumericalDiscretisedEquations

    export famousEquations, Ricker

end
