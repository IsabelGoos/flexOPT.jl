using Symbolics

# this batchNewSymbolics (March 2025-) uses xyz-t coordinates since t should be the fourth coordinates for many reasons

#region some useful partial notations in Cartesian
@variables  x y z t # the physical equation should be a PDE w.r.t. these 4 coordinates in Cartesian (for the moment) 

∂x = Differential(x)
∂y = Differential(y)
∂z = Differential(z)
∂t = Differential(t)

∂x² = ∂x^2
∂y² = ∂y^2
∂z² = ∂z^2
∂t² = ∂t^2


# we can also use these vectors 
∂=[]
∂=push!(∂,∂t)
∂=push!(∂,∂x)
∂=push!(∂,∂y)
∂=push!(∂,∂z)

∇₃=[∂x ; ∂y ; ∂z]
∇₂=[∂x ; ∂y]

∂²=[]
∂²=push!(∂²,∂t²)
∂²=push!(∂²,∂x²)
∂²=push!(∂²,∂y²)
∂²=push!(∂²,∂z²)

#endregion

function Num2Float64(x)
    val = Symbolics.value(x)
    if val isa Num
        error("Cannot convert symbolic Num to Float64: $x")
    end
    return Float64(val)
end

function usefulPartials(coordinates)
    ∂=[]
    ∂²=[]
    for y in coordinates
        partial = Differential(y)
        ∂=push!(∂,partial)
        ∂²=push!(∂²,partial^2)
    end
    return ∂,∂²
end

function myCoeff(eq,x)
    result = Symbolics.coeff(eq,x)
    return result
end

function mySolvefor(eq,x)
    #@show typeof(eq)
    #print("mySolvefor called\n")
    result=Symbolics.symbolic_linear_solve(eq,x)
    #result=Symbolics.solve_for(eq,x)
    #@show typeof(result)
    result=mySimplify(result)
    #print("mySolvefor ended\n")
    #@show typeof(result)
    return result
end

function mySimplify(eq)
    # this function is a buffer for Symbolics.simplify that has some problems when there is a fraction
    # the initial procedure shares the same test with integrateTaylorPolynomials but I put this function independently
    # with the hope that one day I just declare mySimplify = Symbolics.simplify

    #neweq = simplify.(eq;expand=true,threaded=false,simplify_fractions=true)
    #print("mySimplify called with $(length(eq)) equations\n")
    #@time neweq=expand.(eq)
    #@time neweq=simplify_fractions.(neweq)
    #@time neweq=simplify.(eq;expand=true,threaded=false,simplify_fractions=true)
    neweq=expand.(eq)
    #neweq=simplify_fractions.(neweq)
    neweq = expand_derivatives.(neweq)
    #print("mySimplify ended\n")
    
    return neweq 
end


function integrateTaylorPolynomials(eq, x; simplify_expr=mySimplify)
    eq = simplify_expr(eq)

    num, den = if Symbolics.isdiv(eq)
        eq.num, eq.den
    else
        eq, 1
    end

    highestOrder = Symbolics.degree(num, x)
    highestOrder isa Integer || error("expression is not a polynomial in $x")
    highestOrder < 0 && return zero(num)

    # constant term
    # Recover the constant term without substitution. SymbolicUtils currently
    # fails to substitute x = 0 in some higher-order B-spline expressions, and
    # `x^0` may simplify to a plain number before reaching `coeff`.
    constant_term = num
    for i in 1:highestOrder
        constant_term -= Symbolics.coeff(num, x^i) * x^i
    end
    neweq = simplify_expr(constant_term) * x

    # positive powers
    for i in 1:highestOrder
        c = Symbolics.coeff(num, x^i)
        neweq += c * x^(i + 1) / (i + 1)
    end

    return neweq / den
end


function integrateTaylorPolynomials_slow(eq, x)
    # this function works only for (positive) polynomials 
    #
    # the function needs only Symbolics.jl
    #
    # Since the Symbolics.coeff does not support expression with division, lowestNegativeOrder::Int option should be implemented afterwards

    lowestNegativeOrder::Int=0
    eq = mySimplify(eq)

    # eqval=Symbolics.value(eq)
    eqval=eq
    #Symbolics.isdiv(eqval)

    # The above statement is to switch off the negative orders integral
   
    if Symbolics.isdiv(eqval)
        eq=eqval.num 
        denominator=eqval.den 
    else 
        denominator =1
    end
    
    highestOrder=Symbolics.degree(eq,x)


    tmpeq = eq
    neweq = 0
    if highestOrder>0
        for i in 1:highestOrder
            old_coef=Symbolics.coeff(eq,x^i)
            tmpeq-=old_coef*x^i
            new_coef=old_coef/(i+1)
            neweq+=new_coef*x^(i+1)        
        end
    end
    
    if lowestNegativeOrder<0
        old_coef=Symbolics.coeff(eq,x^(-1))
        tmpeq-=old_coef*x^(-1)
        new_coef=old_coef
        neweq+=new_coef*log(x)

        if lowestNegativeOrder<1
            for i in lowestNegativeOrder:-2
                old_coef=Symbolics.coeff(eq,x^i)
                tmpeq-=old_coef*x^i
                new_coef=old_coef//(i+1)
                neweq+=new_coef*x^(i+1)
            end
        end
    end

    # the rest of tmpeq should be the constant
    neweq+=tmpeq*x
    neweq/=denominator
    return mySimplify(neweq)
end
