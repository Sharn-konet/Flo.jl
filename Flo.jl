using Revise

using GLMakie
using Makie.Colors

include("Swarms.jl")
include("ExampleFunctions.jl")
include("Solvers.jl")

using .StrangeAttractors: lorenz
using .Swarms: Swarm, step!
using .Solvers: RK4

const FPS = 144

lorenz_attractor = Swarm(x -> lorenz(x, σ = 10, ρ = 28, β = 8/3), 3, 5, 1e-6)

# set_theme!()

points = Observable(Point3f[])

fig, ax, l = scatter(points,
    colormap = :inferno, transparency = true,
    axis = (; type = Axis3, protrusions = (0, 0, 0, 0),
        viewmode = :fit, limits = (-30, 30, -30, 30, 0, 50)))

for i in 1:10000
    push!(points[], step!(lorenz_attractor, RK4))
end

notify.((points, colors))
