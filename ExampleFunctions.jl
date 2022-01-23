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
end

# lorenz = u::Matrix{Float64} -> lorenz(u; σ = 10, ρ = 28, β = 8/3)

TSUCS1 = @ODE quote
    dx = α * (y - x) + δ * x * z
    dy = ζ * y - x * z
    dz = β * z + x * y - ϵ * x * x
end

TSUCS2 = @ODE quote
    dx = α * (y - x) + δ * x * z
    dy = (ς * x) - (x * z) + ζ * y
    dz = (β * z) + (x * y) - (ϵ * x * x)
end

yuwang = @ODE quote
    dx = α * (y - x)
    dy = β * x - σ * x * z
    dz = exp(x * y) - δ * z
end

aizawa = @ODE quote
    dx = (z - β) * x - δ * y
    dy = δ * x + (z - β) * y
    dz = σ + α * z - (z^3)/3 - (x^2 + y^2)*(1 + ϵ * z) + ζ * z * x^3
end

lorenz_mod_2 = @ODE quote
    dx = -α*x + y^2 - z^2 + α*σ
    dy = x*(y-β*z) + δ
    dz = -z + x*(β * y + z)
end

end