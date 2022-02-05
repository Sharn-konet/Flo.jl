using Revise

include("Interface.jl")

module Attractors

using ..Interface: @ODE

export lorenz, TSUCS1, TSUCS2, yuwang, aizawa, 
       rucklidge, genesio_tesi, finance, shimizu_morioka, 
       noose_hoover, liu_chen, arneodo, bouali,
       burke_shaw, chen_celikovsky, chen_lee, dequan_li,
       hadley, halvorsen, newton_leipnik

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

    α = 40
    β = 0.833
    δ = 0.5
    ϵ = 0.65
    ζ = 20
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

thomas = @ODE quote
    dx = -β * x + sin(y)
    dy = -β * y + sin(z)
    dz = -β * z + sin(x)

    β = 0.19
end

rucklidge = @ODE quote
    dx = -κ*x + α*y - y* z
    dy = x
    dz = -z + y^2

    κ = 2
    α = 6.7
end

genesio_tesi = @ODE quote
    dx = y
    dy = z
    dz = -σ * x - β * y - α * z + x^2

    α = 0.44
    β = 1.1
    σ = 1.0
end

finance = @ODE quote
    dx = ((1/β) - a)*x + z + x*y
    dy = -β*y - x^2
    dz = -x - σ*z

    α = 0.001
    β = 0.2
    σ = 1.1
end

shimizu_morioka = @ODE quote
    dx = y
    dy = (1-z)*x - α*y
    dz = x^2 - β*z

    α = 0.75
    β = 0.45
end

noose_hoover = @ODE quote
    dx = y
    dy = -x + y*z
    dz = α - y^2

    α = 1.5
end

liu_chen = @ODE quote
    dx = α*y + β*x + σ*y*z
    dy = δ*y - z + ϵ*x*z
    dz = ζ*z + η*x*y

    α = 2.4
    β = -3.78
    σ = 14
    δ = -11
    ϵ = 4
    ζ = 5.58
    η = 1
end

arneodo = @ODE quote
    dx = y
    dy = z
    dz = -α*x - β*y - z + σ*x^3

    α = -5.5
    β = 3.5
    σ = -1.0
end

bouali = @ODE quote
    dx = x*(4-y) + α*z
    dy = -y*(1-x^2)
    dz = -x*(1.5 - σ*z) - 0.05*z

    α = 0.3
    σ = 1.0
end

burke_shaw = @ODE quote
    dx = -σ*(x+y)
    dy = -y -σ*x*z
    dz = σ*x*y + ν

    σ = 10
    ν = 4.272
end

chen_celikovsky = @ODE quote
    dx = α*(y-x)
    dy = -x*z + σ*y
    dz = x*y - β*z

    α = 36
    β = 3
    σ = 20
end

chen_lee = @ODE quote
    dx = α*x - y*z
    dy = β*y + x*z
    dz = σ*z + x*(y/3)

    α = 5
    β = -10
    σ = -0.38
end

dequan_li = TSUCS2

hadley = @ODE quote
    dx = -y^2 -z^2 -α*x + α*ζ
    dy = x*y - β*x*z - y + η
    dz = β*x*y + x*z - z

    α = 0.2
    β = 4
    ζ = 8
    η = 1
end

halvorsen = @ODE quote
    dx = -α*x -4*y -4*z - y^2
    dy = -α*y - 4*z - 4*x - z^2
    dz = -α*z - 4*x - 4*y - x^2

    α = 1.4
end

newton_leipnik = @ODE quote
    dx = -α*x + y + 10*y*z
    dy = -x - 0.4*y + 5*x*z
    dz = β*z -5*x*y

    α = 0.4
    β = 0.175
end

end