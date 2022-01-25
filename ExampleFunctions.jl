using Revise

include("Interface.jl")

module StrangeAttractors

using ..Interface: @ODE

export lorenz, TSUCS1, TSUCS2, yuwang, aizawa

function LORENZ(t::Float64, u::Matrix{Float64}; σ::Float64 = 10.0, ρ::Float64 = 28.0, β::Float64 = 8 / 3)

    x, y, z = u[1, :], u[2, :], u[3, :]

    ∂x = σ * (y - x)
    ∂y = x * (ρ - z) - y
    ∂z = x * y - β * z

    return [∂x ∂y ∂z]

end

lorenz = @ODE quote
    dx = σ * (y - x)
    dy = x * (ρ - z) - y
    dz = x * y - β * z

    σ = 10
    ρ = 28
    β = 8/3
end

TSUCS1 = @ODE quote
    dx = α * (y - x) + δ * x * z
    dy = ζ * y - x * z
    dz = β * z + x * y - ϵ * x * x
end

TSUCS2 = @ODE quote
    dx = α * (y - x) + δ * x * z
    dy = (ς * x) - (x * z) + ζ * y
    dz = (β * z) + (x * y) - (ϵ * x * x)

    α = 40
    δ = 0.16
    ς = 55
    ζ = 20
    β = 1.833
    ϵ = 0.65
end

yuwang = @ODE quote
    dx = α * (y - x)
    dy = β * x - σ * x * z
    dz = exp(x * y) - δ * z

    α = 10
    β = 40
    σ = 2
    δ = 2.5
end

aizawa = @ODE quote
    dx = (z - β) * x - δ * y
    dy = δ * x + (z - β) * y
    dz = σ + α * z - (z^3)/3 - (x^2 + y^2)*(1 + ϵ * z) + ζ * z * x^3

    α = 0.95
    β = 0.7
    σ = 0.6
    δ = 3.5
    ϵ = 0.25
    ζ = 0.1
end

lorenz_mod_2 = @ODE quote
    dx = -α*x + y^2 - z^2 + α*σ
    dy = x*(y-β*z) + δ
    dz = -z + x*(β * y + z)

    α = 0.9
    β = 5
    σ = 9.9
    δ = 1
end

end