# 




using Symbolics, LinearAlgebra, Tullio


# famousEquation: we should also need to think how to treat equations with variable division


@doc raw"""
    famousEquations(name)

This is just a gallery of well-known PDEs and model construction methods to test.

you can call your coordinates (locally cartesian) as much as you want but time should be called t not T not τ!!!
# this is important because timeMarching scheme mode should be detected by that!

# Examples
```julia-repl
julia>  famousEquations("1Dacceleration")
```
"""

function famousEquations(name::AbstractString) 
    exprs,fields,vars,extexprs,extfields,extvars,coordinates,∂,∂²=famousEquation(Val(Symbol("eq_"*name)))
    equationCharacteristics=(exprs=exprs, fields=fields, vars=vars, coordinates=coordinates, ∂=∂, ∂²=∂²)
    equationCharacteristicsForce=(exprs=extexprs, fields=extfields, vars=extvars, coordinates=coordinates, ∂=∂, ∂²=∂²)

    return equationCharacteristics,equationCharacteristicsForce
end

function famousEquation(::Val{:eq_1Dacceleration})
     @variables u(t) ω
    exprs =∂t²(u)
    vars = ω
    fields = u
    extexprs=0
    extfields=0
    extvars=nothing
    coordinates =(t)
    ∂, ∂² = usefulPartials(coordinates)
    return exprs, fields, vars, extexprs, extfields, extvars, coordinates, ∂, ∂²
end


function famousEquation(::Val{:eq_1Dlaplacian})
    @variables u(x)
    expr = ∂x²(u)
    exprs = expr
    fields = u
    vars = 1
    extexprs=0
    extfields=0
    extvars=nothing
    coordinates =(x)
    ∂, ∂² = usefulPartials(coordinates)
    return exprs, fields, vars, extexprs, extfields, extvars, coordinates, ∂, ∂²
end

function famousEquation(::Val{:eq_1DsismoFreqHomo})
    @variables ρ μ ω u(x) f(x)
    expr = ρ*ω^2*u + μ*∂x²(u) 

    # declaration of physics
    exprs = expr
    #exts = f

    vars = ρ, μ, ω
    fields = u # for the moment we need to put like ux uy uz to distinguish different components

    extexprs=f
    extfields=f
    extvars=f
    coordinates =(x)
    ∂, ∂² = usefulPartials(coordinates)
    return exprs, fields, vars, extexprs, extfields, extvars, coordinates, ∂, ∂²
end

function famousEquation(::Val{:eq_1DsismoFreqHetero})
    # 1D SH in freq. problem
    @variables ρ(x) μ(x) ω u(x) f(x)
    expr = ρ*ω^2*u + ∂x(μ*∂x(u)) 

    # declaration of physics
    exprs = expr
    #exts = f

    vars = ρ, μ, ω
    fields = u # for the moment we need to put like ux uy uz to distinguish different components

    extexprs=f
    extfields=f
    extvars=1
    coordinates =(x)
    ∂, ∂² = usefulPartials(coordinates)
    return exprs, fields, vars, extexprs, extfields, extvars, coordinates, ∂, ∂²
end

function famousEquation(::Val{:eq_2DsismoSHFreq})
    # 1D SH in freq. problem
    @variables ρ(x,y) μ(x,y) ω u(x,y) f(x,y)
    expr = ρ*ω^2*u + μ*∂x²(u) + μ*∂y²(u)

    # declaration of physics
    exprs = expr
    #exts = f

    vars = ρ, μ, ω, f 
    fields = u # for the moment we need to put like ux uy uz to distinguish different components

    extexprs=f
    extfields=f
    extvars=1
    coordinates =(x,y)
    ∂, ∂² = usefulPartials(coordinates)
    return exprs, fields, vars, extexprs, extfields, extvars, coordinates, ∂, ∂²
end

function famousEquation(::Val{:eq_1DsismoTime})
    # 1D SH in time problem
        
    @variables ρ(x) μ(x) u(x,t) f(x,t) 
    expr = ρ*∂t²(u)- ∂x(μ*∂x(u)) # scalar PDE(s) to be solved

    # declaration of physics
    exprs = expr
    vars = ρ, μ
    fields = u # for the moment we need to put like ux uy uz to distinguish different components

    extexprs=f
    extfields=f
    extvars=f

    coordinates =(x,t)
    ∂, ∂² = usefulPartials(coordinates)
    return exprs, fields, vars, extexprs, extfields, extvars, coordinates, ∂, ∂²
end


function famousEquation(::Val{:eq_1DsismoTimeHomo})
    # 1D SH in time problem
        
    @variables ρ μ u(x,t) f(x,t) 
    expr = ρ*∂t²(u)- ∂x(μ*∂x(u)) # scalar PDE(s) to be solved

    # declaration of physics
    exprs = expr
    vars = ρ, μ
    fields = u # for the moment we need to put like ux uy uz to distinguish different components

    extexprs=f
    extfields=f
    extvars=f

    coordinates =(x,t)
    ∂, ∂² = usefulPartials(coordinates)
    return exprs, fields, vars, extexprs, extfields, extvars, coordinates, ∂, ∂²
end

function famousEquation(::Val{:eq_1DpoissonHetero})
     # 1D Poisson hetero
    @variables κ(x) T(x) f(x)
    expr = ∂x(κ*∂x(T))

    exprs = mySimplify(expr)
    vars = κ
    fields = T

    extexprs=f
    extfields=f
    extvars=f
    coordinates =(x)
    ∂, ∂² = usefulPartials(coordinates)
    return exprs, fields, vars, extexprs, extfields, extvars, coordinates, ∂, ∂²
end

function famousEquation(::Val{:eq_1DpoissonHomo})
     # 1D Poisson hetero
    @variables κ T(x) f(x)
    expr = ∂x(κ*∂x(T))

    exprs = mySimplify(expr)
    vars = κ
    fields = T

    extexprs=f
    extfields=f
    extvars=f
    coordinates =(x)
    ∂, ∂² = usefulPartials(coordinates)
    return exprs, fields, vars, extexprs, extfields, extvars, coordinates, ∂, ∂²
end

function famousEquation(::Val{:eq_2DpoissonHomo})
    # 2D Poisson homo

    @variables κ T(x,y) source(x,y)
    expr = ∂x²(κ*T)+∂y²(κ*T)

    exprs = expr
    vars = κ

    fields = T

    extexprs=source
    extfields=source
    extvars=nothing
    coordinates =(x,y)
    ∂, ∂² = usefulPartials(coordinates)
    return exprs, fields, vars, extexprs, extfields, extvars, coordinates, ∂, ∂²
end


function famousEquation(::Val{:eq_2DpoissonHetero})
    # 2D Poisson hetero

    @variables κ(x,y) T(x,y)  f(x,y) c
    #expr = ∂x²(κ*T)+∂y²(κ*T)
    expr = c*(∂x(κ*∂x(T))+∂y(κ*∂y(T)))
    exprs = mySimplify(expr)
    vars = κ, c
    exts=0
    fields = T

    extexprs=f
    extfields=f
    extvars=f

    coordinates =(x,y)
    ∂, ∂² = usefulPartials(coordinates)
    return exprs, fields, vars, extexprs, extfields, extvars, coordinates, ∂, ∂²
end

function famousEquation(::Val{:eq_2DsismoTimeIsoHeteroForce})
    @variables ρ(x,y) (C(x,y))[1:2,1:2,1:2,1:2] u(x,y,t)[1:2] F(x,y,t)[1:2]
    @variables λ(x,y) μ(x,y)

    δ = Matrix(I, 2, 2)
    @tullio C[i,j,k,l] := λ * δ[i,j]*δ[k,l] + μ*(δ[i,k]*δ[j,l] + δ[i,l]*δ[j,k])
    @tullio traction[i] := ∇₂[j](C[i,j,k,l] * ∇₂[l](u[k]))

    exprs = ρ * ∂t²(u[1]) - traction[1],
            ρ * ∂t²(u[2]) - traction[2]

    fields = u[1], u[2]
    vars = ρ, λ, μ

    extexprs = F[1], F[2]
    extfields = F[1], F[2]
    extvars = F[1], F[2]

    coordinates = (x,y,t)
    ∂, ∂² = usefulPartials(coordinates)
    return exprs, fields, vars, extexprs, extfields, extvars, coordinates, ∂, ∂²
end


function famousEquation(::Val{:eq_2DsismoTimeIsoHomo})

    # 2D wave equation with double couple source
    @variables ρ (C)[1:2,1:2,1:2,1:2] u(x,y,t)[1:2] (M(x,y,t))[1:2,1:2]
    @variables λ μ
    δ=Matrix(I, 2, 2) # Kronecker's delta
    @tullio C[i,j,k,l] := λ * δ[i,j]*δ[k,l]+μ*(δ[i,k]*δ[j,l]+δ[i,l]*δ[j,k])
    @tullio traction[i] := ∇₂[j](C[i,j,k,l]*∇₂[l](u[k]))
    @tullio derivMoment[i] := ∇₂[j](M[i,j])

    exprs = ρ* ∂t²(u[1]) - traction[1], ρ* ∂t²(u[2]) - traction[2]
    fields=u[1], u[2]
    vars = ρ, λ, μ

    extexprs = derivMoment[1], derivMoment[2]
    extfields = M[1,1], M[1,2], M[2,1], M[2,2]
    extvars = M[1,1], M[1,2], M[2,1], M[2,2]
    coordinates =(x,y,t)

    ∂, ∂² = usefulPartials(coordinates)
    return exprs, fields, vars, extexprs, extfields, extvars, coordinates, ∂, ∂²
end

function famousEquation(::Val{:eq_2DsismoTimeIsoHetero})

    # 2D wave equation with double couple source
    @variables ρ(x,y) (C(x,y))[1:2,1:2,1:2,1:2] u(x,y,t)[1:2] (M(x,y,t))[1:2,1:2]
    @variables λ(x,y) μ(x,y)
    δ=Matrix(I, 2, 2) # Kronecker's delta
    @tullio C[i,j,k,l] := λ * δ[i,j]*δ[k,l]+μ*(δ[i,k]*δ[j,l]+δ[i,l]*δ[j,k])
    @tullio traction[i] := ∇₂[j](C[i,j,k,l]*∇₂[l](u[k]))
    @tullio derivMoment[i] := ∇₂[j](M[i,j])

    exprs = ρ* ∂t²(u[1]) - traction[1], ρ* ∂t²(u[2]) - traction[2]

    fields=u[1], u[2]
    vars = ρ, λ, μ

    extexprs = derivMoment[1], derivMoment[2]
    extfields = M[1,1], M[1,2], M[2,1], M[2,2]
    extvars = M[1,1], M[1,2], M[2,1], M[2,2]
    coordinates =(x,y,t)
    ∂, ∂² = usefulPartials(coordinates)
    return exprs, fields, vars, extexprs, extfields, extvars, coordinates, ∂, ∂²
end


function famousEquation(::Val{:eq_2DsismoTimeIsoHeteroSingleForce})

    # 2D heterogeneous elastic wave equation with a Cartesian body force.
    @variables ρ(x,y) (C(x,y))[1:2,1:2,1:2,1:2] u(x,y,t)[1:2] (f(x,y,t))[1:2]
    @variables λ(x,y) μ(x,y)
    δ=Matrix(I, 2, 2) # Kronecker's delta
    @tullio C[i,j,k,l] := λ * δ[i,j]*δ[k,l]+μ*(δ[i,k]*δ[j,l]+δ[i,l]*δ[j,k])
    @tullio traction[i] := ∇₂[j](C[i,j,k,l]*∇₂[l](u[k]))
    exprs = ρ* ∂t²(u[1]) - traction[1], ρ* ∂t²(u[2]) - traction[2]

    fields=u[1], u[2]
    vars = ρ, λ, μ

    extexprs = f[1], f[2]
    extfields = f[1], f[2]
    extvars = f[1], f[2]
    coordinates =(x,y,t)
    ∂, ∂² = usefulPartials(coordinates)
    return exprs, fields, vars, extexprs, extfields, extvars, coordinates, ∂, ∂²
end

function famousEquation(::Val{:eq_2DacousticTime})
    # 2D wave equation with double couple source
    @variables v(x,y)  u(x,y,t) f(x,y,t)
    
    exprs =  ∂t²(u)- v^2 *(∂x²(u) + ∂y²(u))
    fields=u
    vars = v

    extexprs = f
    extfields = f
    extvars = f
    coordinates =(x,y,t)
    ∂, ∂² = usefulPartials(coordinates)
    return exprs, fields, vars, extexprs, extfields, extvars, coordinates, ∂, ∂²
end


function famousEquation(::Val{:eq_2DacousticHomoTime})
    # 2D wave equation with double couple source
    @variables v  u(x,y,t) f(x,y,t)
    
    exprs =  ∂t²(u)- v^2 *(∂x²(u) + ∂y²(u))
    fields=u
    vars = v

    extexprs = f
    extfields = f
    extvars = f
    coordinates =(x,y,t)
    ∂, ∂² = usefulPartials(coordinates)
    return exprs, fields, vars, extexprs, extfields, extvars, coordinates, ∂, ∂²
end

function famousEquation(::Val{:eq_3DsismoTimeIso})
    # 3D isotropic elastic wave equation with a time-dependent moment tensor.
    # M must depend on t: otherwise a three-point temporal recipe receives a
    # one-point force axis and cannot construct its temporal B-spline family.
    @variables ρ(x,y,z) (C(x,y,z))[1:3,1:3,1:3,1:3] u(x,y,z,t)[1:3] (M(x,y,z,t))[1:3,1:3]
    @variables λ(x,y,z) μ(x,y,z)
    #@variables ρ C[1:3,1:3,1:3,1:3] u(x,y,z,t)[1:3] f1,f2,f3 #(M(x,y,z))[1:3,1:3]
    #@variables λ μ
    δ=Matrix(I, 3, 3) # Kronecker's delta
    @tullio C[i,j,k,l] := λ * δ[i,j]*δ[k,l]+μ*(δ[i,k]*δ[j,l]+δ[i,l]*δ[j,k])
    
    @tullio traction[i] := ∇₃[j](C[i,j,k,l]*∇₃[l](u[k]))
    @tullio derivMoment[i] := ∇₃[j](M[i,j])

    exprs = ρ* ∂t²(u[1]) - traction[1], ρ* ∂t²(u[2]) - traction[2], ρ* ∂t²(u[3]) - traction[3]
    fields=u[1], u[2], u[3]
    vars = ρ, λ, μ

    extexprs = derivMoment[1], derivMoment[2], derivMoment[3]
    extfields = M[1,1], M[1,2], M[1,3], M[2,2], M[2,3],M[3,3]
    extvars = M[1,1], M[1,2], M[1,3], M[2,2], M[2,3],M[3,3]
    #extexprs = f1,f2,f3
    #extfields = f1,f2,f3
    #extvars = f1,f2,f3
    coordinates =(x,y,z,t)
    ∂, ∂² = usefulPartials(coordinates)
    return exprs, fields, vars, extexprs, extfields, extvars, coordinates, ∂, ∂²
end

function famousEquation(::Val{:eq_highSchoolProblem})
    
    @variables (X(t))[1:2] g m vh
    
    exprs= m*∂t²(X[2]) - m*g , ∂t(X[1])-vh
    fields = X[1], X[2]
    vars = g, m, vh

    extexprs = 0
    extfields = 0
    extvars = extfields

    coordinates =(t)

    ∂, ∂² = usefulPartials(coordinates)
    return exprs, fields, vars, extexprs, extfields, extvars, coordinates, ∂, ∂²
end
    
