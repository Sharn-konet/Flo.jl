using Revise

using GLMakie
using Makie.Colors

include("Swarms.jl")
include("ExampleFunctions.jl")

using .StrangeAttractors: lorenz
using .Swarms: Swarm, step!, RK4

const FPS = 144

lorenz_attractor = Swarm(x -> lorenz(x, σ = 10, ρ = 28, β = 8/3), 3, 5, 0.5)

# set_theme!()

points = Observable(Point3f[eachcol(lorenz_attractor.positions)...])

fig, ax, l = scatter(points)

for i in 1:100
    step!(lorenz_attractor, RK4)
    points[] = push!(points[], eachcol(lorenz_attractor.positions)...)
    notify(points)
    sleep(1/FPS)
end


