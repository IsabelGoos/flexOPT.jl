# Analytical computation of WYYKK integral (look at demo1DKernelIntegral.ipynb)

using .commonBatchs, UnPack, Symbolics, SHA, Dates


function short_array_string(v)
    join(string.(Float64.(v)), "_")
end

function report_slug(params)
    @unpack orderBspline1D, YorderBspline1Dμᶜ, YorderBspline1Dμ,
            μᶜs, μs, maxNode, ν, lᶜ_nᶜ_max, l_n_max = params

    readable = join([
        "W$(orderBspline1D)",
        "Yc$(YorderBspline1Dμᶜ)",
        "Y$(YorderBspline1Dμ)",
        "N$(maxNode)",
        "nν$(length(ν))",
        "nμc$(length(μᶜs))",
        "nμ$(length(μs))",
        "lc$(lᶜ_nᶜ_max)",
        "l$(l_n_max)",
    ], "__")

    raw_signature = join([
        "nu=" * short_array_string(ν),
        "muc=" * short_array_string(μᶜs),
        "mu=" * short_array_string(μs),
    ], "__")

    digest = bytes2hex(sha1(raw_signature))[1:10]

    return readable * "__" * digest
end



function WYYKKIntegralNumerical(params;ImakeReport=true)
    
    @unpack orderBspline1D, YorderBspline1Dμᶜ, YorderBspline1Dμ, μᶜs, μs, maxNode, ν, lᶜ_nᶜ_max, l_n_max, Δ = params
    νRefPoints = hasproperty(params, :νRefPoints) ? params.νRefPoints : collect(1:maxNode)
    WνOffset = hasproperty(params, :WνOffset) ? params.WνOffset : 0.0

    paramsForSymbolic = @strdict orderBspline1D YorderBspline1Dμᶜ YorderBspline1Dμ μᶜs μs maxNode ν νRefPoints WνOffset lᶜ_nᶜ_max l_n_max ImakeReport
    # v6 deliberately invalidates every pre-half-shift WYYKK result in tmp/.
    paramsForSymbolic["symbolic_endpoint_version"] = "segmentwise_definite_integrals_v6_shift_modes"

    output = myProduceOrLoad(WYYKKIntegralPureSymbolic,paramsForSymbolic,"WYYKKIntegralSymbolic")

    WYYKK_integral = output["WYYKK_integral"]

    coefWYYKK = Array{Float64, 5}(undef,l_n_max+1,lᶜ_nᶜ_max+1,length(μs),length(μᶜs),length(ν))

    nodes = WYYKK_integral.nodes
    numericalNodes = Δ .* nodes
    numNodes = WYYKK_integral.numberNodes
    x = WYYKK_integral.variables[1]
    Δx = WYYKK_integral.variables[2]

    for iν ∈ eachindex(ν), iμᶜ ∈ eachindex(μᶜs), iμ ∈ eachindex(μs), lᶜ_nᶜ ∈ 0:lᶜ_nᶜ_max, l_n ∈ 0:l_n_max
        l_n_slot=l_n+1
        lᶜ_nᶜ_slot = lᶜ_nᶜ+1
        tmpAntiDerivative=WYYKK_integral.data[:,1,l_n_slot,lᶜ_nᶜ_slot,iμ,iμᶜ,iν]
        tmpCoef = 0.0
        # since 
        for ι in 1:numNodes-1
            xLeft = numericalNodes[ι]
            xRight = numericalNodes[ι+1]
            expr = tmpAntiDerivative[ι]
            rightValue = Symbolics.value(Symbolics.substitute(expr, Dict(x => xRight, Δx => Δ)))
            leftValue = Symbolics.value(Symbolics.substitute(expr, Dict(x => xLeft, Δx => Δ)))
            tmpCoef += rightValue-leftValue
        end
        coefWYYKK[l_n_slot,lᶜ_nᶜ_slot,iμ,iμᶜ,iν] = tmpCoef
    end
    return coefWYYKK
end

function WYYKKIntegralPureSymbolic(params::Dict)


    # orders: -1 -> indicator function, 0 -> box car, >=1 -> B-spline

    # this computes the analytical value of the 1D integral between B-spline fns and weighted Taylor kernels
    # \int dx Bspline Y_μᶜ Y_μ  K_{lᶜ-nᶜ}(y-y_μᶜ) K_{l-n}(y-y_μ)

    # unlike the previous integralBsplineTaylorKernels1D, it computes for a specific ν
    # Cˡη;μ are computed for a specific geometry, so even though the boundary condition reduces
    # the number of available points, each Taylor expansion for K_{l-n}(y-y_μ) should be Ok

    @unpack orderBspline1D,YorderBspline1Dμᶜ,YorderBspline1Dμ,μᶜs,μs,maxNode,ν,lᶜ_nᶜ_max,l_n_max,ImakeReport = params
    νRefPoints = haskey(params, "νRefPoints") ? params["νRefPoints"] : collect(1:maxNode)
    WνOffset = haskey(params, "WνOffset") ? params["WνOffset"] : 0.0
    Wν = ν .+ WνOffset
    WνRefPoints = νRefPoints .+ WνOffset

    nodesFromOne = collect(1:1:maxNode) # ∈ Z like [1,2,3], an array of integers collect(1:1:N) (nothing else!!)

    allNodes = unique(sort(vcat(
        Float64.(nodesFromOne),
        Float64.(Wν),
        Float64.(WνRefPoints),
        Float64.(μs),
        Float64.(μᶜs),
    )))

    to_indices(xs, master) = searchsortedfirst.(Ref(master), Float64.(xs))

    idx_nodesFromOne = to_indices(nodesFromOne, allNodes)
    idx_ν            = to_indices(Wν, allNodes)
    idx_νRefPoints   = to_indices(WνRefPoints, allNodes)
    idx_μs           = to_indices(μs, allNodes)
    idx_μᶜs          = to_indices(μᶜs, allNodes)


    # for B-spline

    paramsBSν  = (maximumOrder=orderBspline1D, allNodes = allNodes, idx_nodesNum = idx_nodesFromOne, idx_refPoints = idx_νRefPoints, idx_selectedPoints = idx_ν)
    paramsBSμᶜ = (maximumOrder=YorderBspline1Dμᶜ, allNodes = allNodes, idx_nodesNum = idx_nodesFromOne, idx_refPoints = idx_μᶜs, idx_selectedPoints = idx_μᶜs)
    paramsBSμ  = (maximumOrder=YorderBspline1Dμ, allNodes = allNodes, idx_nodesNum = idx_nodesFromOne, idx_refPoints = idx_μs, idx_selectedPoints = idx_μs)
    # idx_nodesNum : an ordinary consecutive integer increment from 1 (the numerical nodes with Δy)
    # idx_refPoints_original : supporting nodes to construct the Bspline family
    # idx_selectedPoints  : the node addresses that user needs to take, should be a subset of idx_refPoints_original


    # for Taylor Expansions

    paramsTEμᶜ = (maxL_MINUS_N=lᶜ_nᶜ_max, allNodes = allNodes, idx_nodesNum = idx_nodesFromOne, idx_selectedPoints = idx_μᶜs)
    paramsTEμ  = (maxL_MINUS_N=l_n_max, allNodes = allNodes, idx_nodesNum = idx_nodesFromOne, idx_selectedPoints = idx_μs)

    # Computing Wν,Yμᶜ,Yμ 

    Wν = constructBsplineFamily(paramsBSν);
    Yμᶜ = constructBsplineFamily(paramsBSμᶜ);
    Yμ = constructBsplineFamily(paramsBSμ);

    # Computing Kμᶜ, Kμ

    Kμᶜ=constructTaylorExpansions(paramsTEμᶜ);
    Kμ=constructTaylorExpansions(paramsTEμ);

    # Computing WYYKK

    referenceCSF = Kμ.k
    WYYKK=CompactSymbolicFunctions(referenceCSF.nodes,1;auxDims=(l_n_max+1, lᶜ_nᶜ_max+1, length(μs),length(μᶜs),length(ν)),variables=referenceCSF.variables)

    orderSpecν = bspline_order_spec(orderBspline1D)
    orderSpecμᶜ = bspline_order_spec(YorderBspline1Dμᶜ)
    orderSpecμ = bspline_order_spec(YorderBspline1Dμ)

    # The public API still gives integer orders. From here on, the integral uses
    # semantic modes: :indicator has multiplicative factor 1; :bspline reads the
    # explicit B-spline slot. This keeps -1 out of array indexing.
    function bspline_factor(family, iFunction, derivSlot, spec)
        if spec.kind == :indicator
            is_indicator_family(family) || error("Expected :indicator family for order $(spec.order)")
            return 1
        end
        is_bspline_family(family) || error("Expected :bspline family for order $(spec.order)")
        return family.b.data[:, iFunction, derivSlot, bspline_order_slot(spec)]
    end

    derivSlot = 1 # no derivatives

    for iν ∈ eachindex(ν), iμᶜ ∈ eachindex(μᶜs), iμ ∈ eachindex(μs), lᶜ_nᶜ ∈ 0:lᶜ_nᶜ_max, l_n ∈ 0:l_n_max
        l_n_slot=l_n+1
        lᶜ_nᶜ_slot = lᶜ_nᶜ+1
        Wfactor = bspline_factor(Wν, iν, derivSlot, orderSpecν)
        YcFactor = bspline_factor(Yμᶜ, iμᶜ, derivSlot, orderSpecμᶜ)
        YFactor = bspline_factor(Yμ, iμ, derivSlot, orderSpecμ)
        WYYKK.data[:,1,l_n_slot, lᶜ_nᶜ_slot, iμ, iμᶜ,iν] = mySimplify(
            Wfactor .*
            YcFactor .*
            YFactor .*
            Kμᶜ.k.data[:,iμᶜ,lᶜ_nᶜ_slot] .*
            Kμ.k.data[:,iμ,l_n_slot]
        )
    end




    WYYKK_integral = integrate(WYYKK,x)
    # `WYYKKIntegralNumerical` evaluates F(right)-F(left) independently on
    # every polynomial segment, so additive integration constants cancel.
    # Enforcing continuity is unnecessary and fails in SymbolicUtils for some
    # higher-order B-spline expressions.
    
    reportDir = nothing

    if ImakeReport
        
        reportParentDir = joinpath(pwd(), "tmp/WYYKKIntegralPureSymbolic_report")
        reportDir = joinpath(reportParentDir, report_slug(params) * "__" * Dates.format(now(), "yyyymmdd_HHMMSS"))
        mkpath(reportDir)
        
        open(joinpath(reportDir, "params_summary.txt"), "w") do io
            println(io, "orderBspline1D = ", orderBspline1D)
            println(io, "YorderBspline1Dμᶜ = ", YorderBspline1Dμᶜ)
            println(io, "YorderBspline1Dμ = ", YorderBspline1Dμ)
            println(io, "maxNode = ", maxNode)
            println(io, "ν = ", ν)
            println(io, "μᶜs = ", μᶜs)
            println(io, "μs = ", μs)
            println(io, "lᶜ_nᶜ_max = ", lᶜ_nᶜ_max)
            println(io, "l_n_max = ", l_n_max)
        end


        # Page 1: B-spline families
        save_optional_bspline_report_plot(
            joinpath(reportDir, "01_Wnu.pdf"),
            Wν;
            derivOrder=0,
            order=orderBspline1D,
            N=80,
            title="Wν",
        )

        save_optional_bspline_report_plot(
            joinpath(reportDir, "02_Ymuc.pdf"),
            Yμᶜ;
            derivOrder=0,
            order=YorderBspline1Dμᶜ,
            N=80,
            title="Yμᶜ",
        )

        save_optional_bspline_report_plot(
            joinpath(reportDir, "03_Ymu.pdf"),
            Yμ;
            derivOrder=0,
            order=YorderBspline1Dμ,
            N=80,
            title="Yμ",
        )

        # Kμᶜ pages
        for iμᶜ in eachindex(μᶜs), lᶜ_nᶜ in 0:lᶜ_nᶜ_max
            lᶜ_nᶜ_slot = lᶜ_nᶜ + 1
            save_report_plot(
                joinpath(reportDir, "10_Kmuc_i$(iμᶜ)_l$(lᶜ_nᶜ).pdf"),
                Kμᶜ.k,
                (lᶜ_nᶜ_slot,);
                N=80,
                title="Kμᶜ, iμᶜ=$iμᶜ, lᶜ-nᶜ=$lᶜ_nᶜ",
                ylabel="Taylor kernel",
            )
        end

        # Kμ pages
        for iμ in eachindex(μs), l_n in 0:l_n_max
            l_n_slot = l_n + 1
            save_report_plot(
                joinpath(reportDir, "20_Kmu_i$(iμ)_l$(l_n).pdf"),
                Kμ.k,
                (l_n_slot,);
                N=80,
                title="Kμ, iμ=$iμ, l-n=$l_n",
                ylabel="Taylor kernel",
            )
        end

        # WYYKK pages
        for iν in eachindex(ν), iμᶜ in eachindex(μᶜs), iμ in eachindex(μs), lᶜ_nᶜ in 0:lᶜ_nᶜ_max, l_n in 0:l_n_max
            l_n_slot = l_n + 1
            lᶜ_nᶜ_slot = lᶜ_nᶜ + 1

            save_report_plot(
                joinpath(reportDir, "30_WYYKK_inu$(iν)_imuc$(iμᶜ)_imu$(iμ)_lc$(lᶜ_nᶜ)_l$(l_n).pdf"),
                WYYKK,
                (l_n_slot, lᶜ_nᶜ_slot, iμ, iμᶜ, iν);
                N=80,
                title="WYYKK, iν=$iν, iμᶜ=$iμᶜ, iμ=$iμ, lᶜ-nᶜ=$lᶜ_nᶜ, l-n=$l_n",
                ylabel="integrand",
            )
        end

        # WYYKK_integral pages
        for iν in eachindex(ν), iμᶜ in eachindex(μᶜs), iμ in eachindex(μs), lᶜ_nᶜ in 0:lᶜ_nᶜ_max, l_n in 0:l_n_max
            l_n_slot = l_n + 1
            lᶜ_nᶜ_slot = lᶜ_nᶜ + 1

            save_report_plot(
                joinpath(reportDir, "40_WYYKK_integral_inu$(iν)_imuc$(iμᶜ)_imu$(iμ)_lc$(lᶜ_nᶜ)_l$(l_n).pdf"),
                WYYKK_integral,
                (l_n_slot, lᶜ_nᶜ_slot, iμ, iμᶜ, iν);
                N=80,
                title="∫WYYKK, iν=$iν, iμᶜ=$iμᶜ, iμ=$iμ, lᶜ-nᶜ=$lᶜ_nᶜ, l-n=$l_n",
                ylabel="integral",
            )
        end
    end

    
    return @strdict(WYYKK_integral=WYYKK_integral, reportDir=reportDir)

end


function save_report_plot(
    filepath::AbstractString,
    csf::CompactSymbolicFunctions,
    slots::NTuple;
    N=80,
    Δxval=1.0,
    title="",
    ylabel="value",
)
    fig = plotCompactSymbolicFunctions(csf, slots; N=N, Δxval=Δxval, ylabel=ylabel)
    ax = content(fig[1, 1])
    ax.title = title
    save(filepath, fig)
    return filepath
end

function save_bspline_report_plot(
    filepath::AbstractString,
    csf::CompactSymbolicFunctions;
    derivOrder=0,
    order=0,
    N=80,
    Δxval=1.0,
    title="",
)
    slots = (derivOrder + 1, order + 1)
    ylabel = derivOrder == 0 ? "B-spline value" : "Derivative order $derivOrder"
    return save_report_plot(filepath, csf, slots; N=N, Δxval=Δxval, title=title, ylabel=ylabel)
end

function save_optional_bspline_report_plot(
    filepath::AbstractString,
    family;
    derivOrder=0,
    order=0,
    N=80,
    Δxval=1.0,
    title="",
)
    spec = bspline_order_spec(order)
    if spec.kind == :indicator
        is_indicator_family(family) || error("Expected :indicator family for report order $(spec.order)")
        txtpath = replace(filepath, r"\.pdf$" => ".txt")
        open(txtpath, "w") do io
            println(io, title)
            println(io, "family.kind = :indicator")
            println(io, "order = ", order)
            println(io, "This corresponds to the direct/indicator case.")
        end
        return txtpath
    end

    is_bspline_family(family) || error("Expected :bspline family for report order $(spec.order)")
    return save_bspline_report_plot(
        filepath,
        family.b;
        derivOrder=derivOrder,
        order=spec.order,
        N=N,
        Δxval=Δxval,
        title=title,
    )
end
