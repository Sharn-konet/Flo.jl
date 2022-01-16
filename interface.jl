module Interface

using MacroTools: postwalk

function collectSymbols!(expression::Expr, symbols::Vector{Symbol})
    for arg in expression.args
        if hasproperty(arg, :args)
            collectSymbols!(arg, symbols)
        else
            push!(symbols, arg)
        end
    end
end

function isFunction(symbol::Symbol)
    if isdefined(Main, symbol)
        return isa(getfield(Main, symbol), Function)
    end
    return false
end

struct ODEFunction <: Function
    func::Function
    dimensions::Int8
end


# """ 
#     ODEFunc!(t, u, [constants...])


# Update u based on current position and system defined using @ODE and provided
# constants.

# # Example
# ``` julia-repl
# julia> lorenz = @ODE quote
#     dx = σ*(y-x) 
#     dy = x*(ρ-z) - y 
#     dz = x * y - β*z
# end

# julia> u = [
#     [1,2,3]
#     [2,3,4]
#     [4,5,6]
# ]

# julia> lorenz(0.1, u, σ = 10, ρ = 28, β = 8/3)
# ```
# """
macro ODE(∂x::Expr, ∂y::Expr, ∂z::Expr)
  
    symbols = Vector{Symbol}()

    for expression in (∂x, ∂y, ∂z)
        collectSymbols!(expression, symbols)
    end

    symbols = filter!(!Meta.isoperator, symbols)
    symbols = filter!(!isFunction, symbols)

    diff_symbols = Vector{Symbol}([expr.args[1] for expr in [∂x, ∂y, ∂z]])
    dependent_vars = Vector{Symbol}([Symbol(string(symbol)[end]) for symbol in diff_symbols])
    constants = setdiff(symbols, diff_symbols, dependent_vars)
    constants = [Expr(:(::), constant, :(Real)) for constant in constants]

    # @show dependent_vars = Expr(:tuple, dependent_vars[1:end-1]..., :($(dependent_vars[end]) = eachrow(u)))
    diff_var_initialisation = [:($symbol = Vector{Float64}(undef, size(u,2))) for symbol in diff_symbols]
    system = [var"@__dot__"(LineNumberNode(1), Main, ∂).args[1] for ∂ in [∂x, ∂y, ∂z]] # Apply broadcasting

    func = quote
        function ODEFunc(t::Real, u::Matrix{<:Real}; $(constants...))

            x, y, z = eachrow(u)

            $(diff_var_initialisation...)
            $(system...)

            return $(diff_symbols...)
        end

        ODEFunction(ODEFunc, 3)
    end

    @show func

    return func
end

# Method handles quote blocks being used to define ODEs
macro ODE(derivatives::Expr)
    derivatives = filter(x -> !(x isa LineNumberNode), derivatives.args)
    # @assert length(derivatives) == 3 "Incorrected number of derivatives defined."
    @eval @ODE($(derivatives...))
end

end