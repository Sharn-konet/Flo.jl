using Revise

using GLMakie
using Makie.Colors

include("Swarms.jl")
include("ExampleFunctions.jl")

using .StrangeAttractors: lorenz, TSUCS2, yuwang, aizawa
using .Swarms: Swarm, step!, RK4

const FPS = 144

lorenz_attractor = Swarm(x -> lorenz(x, σ = 10, ρ = 28, β = 8/3), 3, 1000, 1e-2)
TSUCS2_attractor = Swarm(x -> TSUCS2(x, α = 40, δ = 0.16, ς = 55, ζ = 20, β = 1.833, ϵ = 0.65), 3, 10000, 1e-2)
yuwang_attractor = Swarm(x -> yuwang(x, α = 10, β = 40, σ = 2, δ = 2.5), 3, 2500, 1e-2)
aizawa_atrractor = Swarm(x -> aizawa(x, α = 0.95, β = 0.7, σ = 0.6, δ = 3.5, ϵ = 0.25, ζ = 0.1), 3, 3000, 1e-2)

attractor = aizawa_atrractor

# set_theme!()

points = Observable(Point3f.(Vector{Float64}.([eachcol(attractor.positions)...])))

fig = Figure()
# ax = Axis3(fig[1,1], protrusions = (0, 0, 0, 0),
# viewmode = :fit, limits = (-3, 3, -5, 5, 0, 20))
ax = Axis3(fig[1,1], protrusions = (0, 0, 0, 0),
viewmode = :fit, limits = (-1.5, 1.5, -1.5, 1.5, -0.5, 2), perspectiveness = 0.5)
# ax = Axis3(fig[1,1], protrusions = (0, 0, 0, 0),
# viewmode = :fit, limits = (-100, 100, -100, 100, 0, 100))
fig[2,1] = buttongrid = GridLayout(tellwidth = false)
sim_state = Observable("Run")
run_button = buttongrid[1,1] = Button(fig, label = sim_state)

scatter!(points, markersize = 20)

run = Observable(false)

on(run_button.clicks) do clicks
    states = Dict(true => "Stop", false => "Run")
    run[] = !run[]
    sim_state[] = states[run[]]
    notify.((run, sim_state))
end

while true
    # Main sim loop
    while run.val
        step!(attractor, RK4)
        points[] = [eachcol(attractor.positions)...]
        notify(points)
        sleep(1/FPS)
    end
    sleep(0.01)
end

