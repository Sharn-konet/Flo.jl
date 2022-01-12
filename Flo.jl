using Revise
using GLMakie

struct FlowFunction <: Function
    fx::Expr
    fy::Expr
    fz::Expr
end

struct Swarm
    func::Function
    size::Int64
    dim::Int8
    tol::Float64
    positions::Matrix{Float64}
    step_size::Vector{Float64}
    error_history::Array{Float64}

    function Swarm(func::Function, dim::Int64, size::Int64, step_size::Float64)
        
        if dim != 3
            error("Number of dimensions not supported.")
        end
        
        new(func, size, 3, 1e-8,
            Matrix{Float64}(undef, 3, size),
            repeat([step_size], size),
            Array{Float64}(undef, 3, size, 0)
        )
    end
end