using Revise

using GLMakie
using Makie.Colors

include("Swarms.jl")
include("ExampleFunctions.jl")

using .StrangeAttractors: lorenz, TSUCS2, yuwang, aizawa, lorenz_mod_2
using .Swarms: Swarm, step!, RK4

const FPS = 144

lorenz_attractor = Swarm(lorenz)
TSUCS2_attractor = Swarm(TSUCS2, (α = 40, δ = 0.16, ς = 55, ζ = 20, β = 1.833, ϵ = 0.65))
yuwang_attractor = Swarm(yuwang, (α = 10, β = 40, σ = 2, δ = 2.5))
aizawa_atrractor = Swarm(aizawa, (α = 0.95, β = 0.7, σ = 0.6, δ = 3.5, ϵ = 0.25, ζ = 0.1))
lorenz_mod_2_attractor = Swarm(lorenz_mod_2, (α = 0.9, β = 5, σ = 9.9, δ = 1))

attractor = lorenz_attractor

fig = Figure()
ax = Axis3(fig[1,1], protrusions = (0, 0, 0, 0),
viewmode = :fit, limits = (-3, 3, -5, 5, 0, 20))
# ax = Axis3(fig[1,1], protrusions = (0, 0, 0, 0),
# viewmode = :fit, limits = (-1.5, 1.5, -1.5, 1.5, -0.5, 2), perspectiveness = 0.5)

# ax = Axis3(fig[1,1], protrusions = (0, 0, 0, 0),
# viewmode = :fit, limits = (-100, 100, -100, 100, 0, 100))
fig[2,1] = buttongrid = GridLayout(tellwidth = false)
sim_state = Observable("Run")
run_button = buttongrid[1,1] = Button(fig, label = sim_state)

run = Observable(false)

on(run_button.clicks) do clicks
    run[] = !run[]
    sim_state[] = run.val ? "Stop" : "Run"
    notify.((run, sim_state))
end

lsgrid = labelslidergrid!(
    fig, 
    [string.(keys((σ = 10, ρ = 28, β = 8/3)))...],
    [values((σ = 10, ρ = 28, β = 8/3))...],
    width = 350,
    tellheight = false
)

fig[1,2] = lsgrid.layout

sliderobservables = [s.value for s in lsgrid.sliders]

for observable in sliderobservables
    on(observable) do val
        attractor = Swarm(lorenz, lorenz_defaults)
    end
end

bars = lift(sliderobservables...) do slvalues...
    [slvalues...]
end

solver = RK4

points = Observable(Point3f.(Vector{Float64}.([eachcol(attractor.positions)...])))

scatter!(points, markersize = 2000)

while true
    # Main sim loop
    while run.val
        step!(attractor, solver)
        points[] = [eachcol(attractor.positions)...]
        notify(points)
        sleep(0.001)
    end
    sleep(0.01)
end



