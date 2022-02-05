module Interface

using MacroTools: postwalk, isline, inexpr

function collectSymbols!(expression::Expr, symbols::Vector{Symbol})
    for arg in expression.args
        if hasproperty(arg, :args)
            collectSymbols!(arg, symbols)
        elseif !(typeof(arg) <: Number)
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

function replaceVars(expr::Union{Expr,Symbol}, mapping::Dict)
    if expr in keys(mapping)
        return :(u[$(mapping[expr]), :])
    end

    return expr
end

# Handle usage of Integers/Floats within expressions (such as exponents)
replaceVars(expr::Number, mapping::Dict) = return expr

struct ODEFunction <: Function
    func::Function
    dimensions::Int8
end


""" 
    ODEFunc!(t, u, [constants...])


Update u based on current position and system defined using @ODE and provided
constants.

# Example
``` julia-repl
julia> lorenz = @ODE quote
    dx = σ*(y-x) 
    dy = x*(ρ-z) - y 
    dz = x * y - β*z
end

julia> u = [
    [1,2,3]
    [2,3,4]
    [4,5,6]
]

julia> lorenz(0.1, u, σ = 10, ρ = 28, β = 8/3)
```
"""
macro ODE(dx⃗::Expr...)

    dim = length(dx⃗)

    symbols = Vector{Symbol}()

    for expression in dx⃗
        collectSymbols!(expression, symbols)
    end

    # Filter out operators and functions
    filter!(!Meta.isoperator, symbols)
    filter!(!isFunction, symbols)
    filter!(!isequal(:t), symbols)

    diff_symbols = Vector{Symbol}([expr.args[1] for expr in dx⃗ if length(string(expr.args[1])) == 2])

    dependent_vars = Vector{Symbol}([Symbol(string(symbol)[end]) for symbol in diff_symbols])
    mapping_to_indices = Dict(symbol => index for (index, symbol) in enumerate(dependent_vars))

    constants = setdiff(symbols, diff_symbols, dependent_vars, [:t])
    constant_defaults = Dict{Symbol,Expr}([expr.args[1] => expr for expr in dx⃗ if (length(string(expr.args[1])) == 1) & (expr.args[1] in constants)])
    constants = [constant in keys(constant_defaults) ? Expr(:(kw), Expr(:(::), Symbol(constant), :(Real)), constant_defaults[constant].args[2]) : Expr(:(::), constant, :(Real)) for constant in constants]

    filter(expr -> expr.args[1] in keys(constant_defaults), dx⃗)

    diff_var_initialisation = [:($symbol = Vector{Float64}(undef, size(u, 2))) for symbol in diff_symbols]
    system = [var"@__dot__"(LineNumberNode(1), Main, expr).args[1] for expr in dx⃗ if !(expr in values(constant_defaults))] # Apply broadcasting
    system = [postwalk(x -> replaceVars(x, mapping_to_indices), expr) for expr in system] # Replace variables with indices

    if any(inexpr(expr, :t) for expr in dx⃗)
        func = quote
            function ODEFunc(t::Real, u::Matrix{<:Real}; $(constants...))::Matrix{Float64}
                $(diff_var_initialisation...)
                $(system...)
                return hcat($(diff_symbols...))
            end
        end
    else
        func = quote
            function ODEFunc(u::Matrix{<:Real}; $(constants...))::Matrix{Float64}
                $(diff_var_initialisation...)
                $(system...)
                return hcat($(diff_symbols...))
            end
        end
    end

    return func
end

# Method handles quote blocks being used to define ODEs
macro ODE(derivative_block::Expr)
    derivatives = filter(x -> !isline(x), derivative_block.args)
    @eval @ODE($(derivatives...))
end

end