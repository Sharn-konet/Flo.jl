using Revise
using GLMakie

struct FlowFunction <: Function
    fx::Expr
    fy::Expr
    fz::Expr
end

Base.@kwdef struct Swarm
    func::Function
    dim::Int
    positions::Vector{Float64}
    step_size::Vector{Float64}
    tol::Float64
    
    error_history::Vector{Float64} = Vector{Float64}()
    Swarm(dim, positions, step_size, tol, func)
end