module flexOPT

    using Symbolics,UnPack,Tullio
    using Dates
    using LinearAlgebra
    using KernelAbstractions
    using ..commonBatchs
    using Adapt

    # Use the backend selected by batchFiles/batchGPU.jl when the caller has
    # loaded it.  A module import must nevertheless remain valid on its own:
    # recipe construction then falls back to the KernelAbstractions CPU.
    const backend = isdefined(Main, :backend) ?
        getfield(Main, :backend) : KernelAbstractions.CPU()
    const makeGPUarray = isdefined(Main, :makeGPUarray) ?
        getfield(Main, :makeGPUarray) :
        ((selected_backend, array) -> Adapt.adapt(selected_backend, array))

    """
    Global default formulation used by `makeOPTsemiSymbolic` when the recipe
    parameters do not contain `variationalForm`.

    Use `set_default_form!` rather than assigning a new value to this constant.
    """
    const DEFAULT_FORM = Ref{Symbol}(:weak)

    function set_default_form!(form::Symbol)
        normalised = Symbol(lowercase(String(form)))
        normalised in (:weak, :strong) ||
            throw(ArgumentError("form must be :weak or :strong"))
        DEFAULT_FORM[] = normalised
        return normalised
    end

    set_default_form!(form::AbstractString) = set_default_form!(Symbol(form))
    export DEFAULT_FORM, set_default_form!


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
    export timeDimensionString, WeakTerm, BoundaryFlux
    export naturalWeakForm, weakTermGroups
    include("../src/motorsOPT/nondimensionalisation.jl")
    export opt_nondimensionalization, isotropic_elasticity_tensor
    export nondimensionalize_elasticity_tensor, nondimensionalize_body_force
    export nondimensionalize_moment_tensor, prepare_nondimensional_recipe
    export prepare_nondimensional_elasticity, elasticity_component_models
    include("../src/motorsOPT/famousEquations.jl")
    include("../src/motorsOPT/famousSourceFunctions.jl")
    include("../src/motorsOPT/famousBoundaryConditions.jl")
    export CerjanBoundarySpec, BoundaryConditionSet, boundary_geometry
    export cerjan_padding, famousBoundaryConditions
    # moment theory of the test function W (Claude-authored, Nov 2026)
    include("../src/motorsOPT/testFunctionMoments.jl")
    using .TestFunctionMoments
    export bsplineTestFunction, hatTestFunction, centredBspline, mixTestFunctions,
           testFunctionMoments, predictedOrder, optimalTestFunctionWeights, flexOPTRefPoints

    # semi-symbolics operators to fully numerical operators
    #include("fullyNumericalOPT/makeCostFunctions.jl")
    include("fullyNumericalOPT/makeCostFunctions.jl")
    include("fullyNumericalOPT/others.jl")
    export numericalOperatorConstruction,getModelPoints

    include("../src/numSolvers/linearSolution.jl")
    include("../src/numSolvers/timeMarchingSchemes.jl")
    include("../src/numSolvers/diffTools.jl")
    export prepareNumericalOperators,timeMarchingSchemePrepared,prepareLinearSystem,evaluateLinearSystem!,evaluateLinearSystem,timeMarchingSchemeLinear,propagateLinearSystem,overlapBoundaryLinearSystem
    export prepareLumpedFuturePredictor,applyLumpedFuturePredictor!


    export quasiNumericalOperatorConstruction,constructingNumericalDiscretisedEquations

    export famousEquations, Ricker

end
