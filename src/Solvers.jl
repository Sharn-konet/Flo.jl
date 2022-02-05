abstract type ButcherTableau end

"""
    RungeKuttaMethod(α, β, γ, order)

Simple RungeKuttaMethod which can be used to solve a set of differential equations.
It's required that the number of rows and columns of γ match the dimensions of α and β respectively.
"""
struct RungeKuttaMethod <: ButcherTableau
    α::Vector{<:Real}
    β::Vector{<:Real}
    γ::Matrix{<:Real}
    order::Int
    RungeKuttaMethod(α, β, γ, order) = (length(α) == size(γ, 1)) & (length(β) == size(γ,2)) ? 
        new(α, β, γ, order) : error("Not a valid Butcher Tableau")
end

"""
    AdaptiveRungeKuttaMethod(α, β, γ, order)

Extended version of RungeKuttaMethod where methods of two different orders are used to approximate the 
error of each calculation. Allows for each step to be within a specified tolerance.
"""
struct AdaptiveRungeKuttaMethod <: ButcherTableau
    α::Vector{<:Real}
    β::Matrix{<:Real}
    γ::Matrix{<:Real}
    order::Int
    AdaptiveRungeKuttaMethod(α, β, γ, order) = length(α) == size(γ, 1) & size(β, 2) == size(γ,2) & size(β, 1) == 2 ? 
        new(α, β, γ, order) : error("Not a valid Butcher Tableau")
end

RK4 = RungeKuttaMethod(
    [1/6, 1/3, 1/3, 1/6], # α
    [0, 1/2, 1/2, 1], # β
    [[0, 1/2, 0, 0] [0, 0, 1/2, 0] [0, 0, 0, 1] [0, 0, 0, 0]], # γ
    4 # order
)