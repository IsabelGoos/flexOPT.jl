# application/Poisson_1D/run.jl -- 1-D heterogeneous Poisson, the template application
#
# From a terminal, in the repository root:
#   julia --project=. application/Poisson_1D/run.jl                       # forward, OPT3, case sameλ
#   julia --project=. application/Poisson_1D/run.jl forward scheme=convFD3 case=quadratic h=0.1
#   julia --project=. application/Poisson_1D/run.jl convergence           # convFD3 vs OPT3
#   julia --project=. application/Poisson_1D/run.jl convergence schemes=convFD3,convFD5,OPT3 case=twiceλ
#   julia --project=. application/Poisson_1D/run.jl forward scheme=OPT3_stag experimental=true
#   julia --project=. application/Poisson_1D/run.jl list                  # schemes and cases
#
# From a REPL (fastest: recipes and compiled code stay in memory):
#   julia> include("application/Poisson_1D/run.jl")          # runs the default once
#   julia> r = poisson1DForward(scheme="convFD5", case="λ_2", h=0.05); r.error
#   julia> poisson1D("convergence", "case=quadratic")         # same as the command line
#
# Outputs: application/Poisson_1D/outputs/<task>/ (PNG + CSV; git-ignored)

FLEXOPT_MODULES = [:flexOPT]
include(joinpath(@__DIR__, "..", "common", "appBootstrap.jl"))
using CairoMakie, DelimitedFiles
CairoMakie.activate!()
include(joinpath(@__DIR__, "app.jl"))

poisson1D(args...) = (global APP_ARGS = collect(String, string.(args)); include(joinpath(@__DIR__, "run.jl")))

let A = appArgs(@isdefined(APP_ARGS) ? APP_ARGS : ARGS;
                scheme="OPT3", case="sameλ", h=exp(-2.0), schemes=["convFD3", "OPT3"],
                logInvH=[0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0], experimental=false, precision="Float64")
    task = isempty(A.positional) ? "forward" : A.positional[1]
    precision = A.precision == "Float32" ? Float32 : Float64
    schemes = A.schemes isa AbstractVector ? String.(string.(A.schemes)) : [string(A.schemes)]
    println("Poisson_1D: task=$task  ", A)

    if task == "list"
        println("\nschemes:")
        for (k, s) in sort(collect(POISSON1D_SCHEMES); by=first)
            println("  ", rpad(k, 12), rpad(string(s.status), 14), s.version, "   ", s.from)
        end
        println("\ncases (κ(x), exact T = cos x on [0, 10π]): ", join(sort(collect(keys(POISSON1D_CASES))), ", "))

    elseif task == "forward"
        r = poisson1DForward(; scheme=string(A.scheme), case=string(A.case), h=Float64(A.h),
                             precision, experimental=A.experimental)
        println("Nx=$(r.Nx)  Δx=$(round(r.Δx, sigdigits=4))  error=$(r.error)  ($(round(r.seconds, digits=2)) s)")
        dir = appOutputDir(@__DIR__, "forward")
        fig = Figure(size=(800, 600))
        ax1 = Axis(fig[1, 1], ylabel="κ(x)", title="Poisson_1D  $(r.case)  $(r.scheme)  Nx=$(r.Nx)")
        lines!(ax1, r.X, r.κ)
        ax2 = Axis(fig[2, 1], ylabel="T(x)")
        lines!(ax2, r.X, r.exact, color=:black, label="exact cos x")
        scatter!(ax2, r.X, r.T, color=:red, markersize=4, label=r.scheme)
        axislegend(ax2)
        ax3 = Axis(fig[3, 1], xlabel="x", ylabel="T − exact")
        lines!(ax3, r.X, r.T .- r.exact, color=:red)
        appSave(fig, dir, "forward_$(r.case)_$(r.scheme)_Nx$(r.Nx)")
        writedlm(joinpath(dir, "forward_$(r.case)_$(r.scheme)_Nx$(r.Nx).csv"),
                 vcat(["x" "kappa" "f" "T" "exact"], hcat(r.X, r.κ, r.f, r.T, r.exact)), ',')
        global poisson1DLast = r

    elseif task == "convergence"
        rows = poisson1DConvergence(; schemes, case=string(A.case), logInvH=A.logInvH,
                                    precision, experimental=A.experimental)
        dir = appOutputDir(@__DIR__, "convergence")
        fig = Figure(size=(700, 500))
        ax = Axis(fig[1, 1], xscale=log10, yscale=log10, xlabel="Δx", ylabel="‖T − exact‖₂ / Nx",
                  title="Poisson_1D convergence, case $(A.case)")
        for s in schemes
            r = filter(x -> x.scheme == s, rows)
            scatterlines!(ax, [x.Δx for x in r], [x.error for x in r], label=s)
            println("observed orders ", rpad(s, 10), ": ", round.(poisson1DOrders(rows, s), digits=2))
        end
        axislegend(ax, position=:rb)
        appSave(fig, dir, "convergence_$(A.case)_" * join(schemes, "-"))
        writedlm(joinpath(dir, "convergence_$(A.case)_" * join(schemes, "-") * ".csv"),
                 vcat(["scheme" "case" "dx" "Nx" "error" "seconds"],
                      permutedims(hcat([[x.scheme, x.case, x.Δx, x.Nx, x.error, x.seconds] for x in rows]...))), ',')
        global poisson1DLast = rows
    else
        error("unknown task '$task'; use forward, convergence or list")
    end
end
