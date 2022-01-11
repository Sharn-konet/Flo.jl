module Interface

function collectSymbols!(expression::Expr, symbols::Vector{Symbol})
    for arg in expression.args
        if hasproperty(arg, :args)
            collectSymbols!(arg, symbols)
        else
            push!(symbols, arg)
        end
    end
end

macro ODE(∂x::Expr, ∂y::Expr, ∂z::Expr)
            
    symbols = Vector{Symbol}()

    for expression in (∂x, ∂y, ∂z)
        collectSymbols!(expression, symbols)
    end

    symbols = filter!(!Meta.isoperator, symbols)

    diff_symbols = Set{Symbol}([expr.args[1] for expr in [∂x, ∂y, ∂z]])
    dependent_vars = Set{Symbol}([Symbol(string(symbol)[end]) for symbol in diff_symbols])
    constants = setdiff(symbols, diff_symbols, dependent_vars)

    # Check that ∂ is not used within the code
    @assert nor([occursin("∂", string(expr)) for expr in [∂x, ∂y, ∂z]]...) "Remove usage of ∂."

    quote
        function ODEFunc(t::Float64, u::Matrix{Float64}, $(constants...))

            x, y, z = u[1,:], u[2,:], u[3,:]
        
            $(@.(∂x)); $(@.(∂y)); $(@.(∂z))

            return $(diff_symbols...)
        end
    end 
end

macro ODE(derivatives)
    derivatives = filter(x -> !(x isa LineNumberNode), derivatives.args)
    # @assert length(derivatives) == 3 "Incorrected number of derivatives defined."
    @eval @ODE($(derivatives...))
end

end