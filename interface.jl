module Interface

using MacroTools: postwalk, isline

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

function replaceVars(expr::Union{Expr, Symbol}, mapping::Dict)
    if expr in keys(mapping)
        return :(u[$(mapping[expr]), :])
    end
    
    return expr
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
macro ODE(dx⃗::Expr...)

    dim = length(dx⃗)
  
    symbols = Vector{Symbol}()

    for expression in dx⃗
        collectSymbols!(expression, symbols)
    end

    # Filter out operators and functions
    symbols = filter!(!Meta.isoperator, symbols)
    symbols = filter!(!isFunction, symbols)

    diff_symbols = Vector{Symbol}([expr.args[1] for expr in dx⃗])

    dependent_vars = Vector{Symbol}([Symbol(string(symbol)[end]) for symbol in diff_symbols])
    mapping_to_indices = Dict(symbol => index for (index, symbol) in enumerate(dependent_vars))

    constants = setdiff(symbols, diff_symbols, dependent_vars)
    constants = [Expr(:(::), constant, :(Real)) for constant in constants]

    diff_var_initialisation = [:($symbol = Vector{Float64}(undef, size(u, 2))) for symbol in diff_symbols]
    system = [var"@__dot__"(LineNumberNode(1), Main, expr).args[1] for expr in dx⃗] # Apply broadcasting
    system = [postwalk(x -> replaceVars(x, mapping_to_indices), expr) for expr in system] # Replace variables with indices

    func = quote
        function ODEFunc(t::Real, u::Matrix{<:Real}; $(constants...))::Matrix{Float64}
            $(diff_var_initialisation...)
            $(system...)

            return hcat($(diff_symbols...))
        end
    end
    @show func
    return func
end

# Method handles quote blocks being used to define ODEs
macro ODE(derivative_block::Expr)
    derivatives = filter(x -> !isline(x), derivative_block.args)
    @eval @ODE($(derivatives...))
end

end