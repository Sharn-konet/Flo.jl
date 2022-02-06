module Swarms

include("Solvers.jl")

using TensorOperations

export Swarm, step!

"""
    Swarm(func; [size, dim, tol, step_size, positions])

Structure storing all information needed to describe a swarm of particles.
"""
struct Swarm
    size::Int64
    dim::Int8
    tol::Float64
    positions::Matrix{Float64}
    step_size::Vector{Float64}
    error_history::Array{Float64}

    function Swarm(; 
        size::Int64 = 1000, 
        dim::Int64 = 3, 
        tol::Float64 = 1e-8, 
        step_size::Vector{Float64} = repeat([1e-2], size), 
        positions::Matrix{Float64} = rand(Float64, dim, size) .*2 .-1)
        
        if dim != 3 error("Number of dimensions not supported.") end
        
        new(size, dim, tol, positions, step_size, Array{Float64}(undef, dim, size, 0))
    end
end

"""
    step!(swarm, solver)

Step the swarm one step using the provided solver.
"""
function step!(swarm::Swarm, func::Function, solver::RungeKuttaMethod)

    h = swarm.step_size

    dx⃗ = zeros(Float64, size(swarm.positions)..., length(solver.β))
    y_step = Matrix{Float64}(undef, size(swarm.positions)...)

    for i in 1:length(solver.β)
        weights = solver.γ[i,:]
        @tensor y_step[i,j] = dx⃗[i,j,k]*weights[k]
        dx⃗[:,:,i] = func(swarm.positions + h'.*y_step)'
    end

    @tensor y_step[i,j] = dx⃗[i,j,k]*solver.α[k]

    swarm.positions .+= h'.*y_step
end

end