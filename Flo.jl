using Revise

using GLMakie
using Makie.Colors

include("Swarms.jl")
include("ExampleFunctions.jl")

using .StrangeAttractors: lorenz, TSUCS2
using .Swarms: Swarm, step!, RK4

const FPS = 144

lorenz_attractor = Swarm(x -> lorenz(x, σ = 10, ρ = 28, β = 8/3), 3, 1000, 1e-2)
# TSUCS2_attractor = Swarm(x -> lorenz(x, α = 40, δ = 0.16, ς = 55, ζ = 20, β = 1.833, ϵ = 0.65), 3, 1000, 1e-4)

# set_theme!()

points = Observable(Point3f.(Vector{Float64}.([eachcol(lorenz_attractor.positions)...])))

fig, ax, l = scatter(points, axis = (; type = Axis3, protrusions = (0, 0, 0, 0),
viewmode = :fit, limits = (-30, 30, -30, 30, 0, 50)), markersize = 2000)

for i in 1:5000
    step!(lorenz_attractor, RK4)
    points[] = [eachcol(lorenz_attractor.positions)...]
    notify(points)
    sleep(1/FPS)
end


