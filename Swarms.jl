module Swarms

include("Solvers.jl")

using TensorOperations

export Swarm, step!

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
            rand(Float64, 3, size),
            repeat([step_size], size),
            Array{Float64}(undef, 3, size, 0)
        )
    end
end

function step!(swarm::Swarm, solver::RungeKuttaMethod)

    h = swarm.step_size
    sf = 0.9

    dx⃗ = zeros(Float64, size(swarm.positions)..., length(solver.β))
    y_step = Matrix{Float64}(undef, size(swarm.positions)...)

    for i in 1:length(solver.β)
        weights = solver.γ[i,:]
        @tensor y_step[i,j] = dx⃗[i,j,k]*weights[k]
        dx⃗[:,:,i] = swarm.func(swarm.positions + h'.*y_step)'
    end

    @tensor y_step[i,j] = dx⃗[i,j,k]*solver.α[k]

    swarm.positions .+= h'.*y_step

end

end