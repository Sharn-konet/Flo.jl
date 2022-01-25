using Revise

using GLMakie
using Makie.Colors

include("Swarms.jl")
include("ExampleFunctions.jl")

using .StrangeAttractors: lorenz, TSUCS2, yuwang, aizawa, lorenz_mod_2, TSUCS1
using .Swarms: Swarm, step!, RK4

const FPS = 144

attractor = Swarm(lorenz)

function createFigure(fig::Makie.Figure = Figure())
    ax = Axis3(fig[1,1], protrusions = (0, 0, 0, 0), viewmode = :fit,  aspect = :data)

    fig[2,1] = buttongrid = GridLayout(tellwidth = false)
    sim_state = Observable("Run")
    run_button = buttongrid[1,1] = Button(fig, label = sim_state)

    run = Observable(false)

    on(run_button.clicks) do clicks
        run[] = !run[]
        sim_state[] = run.val ? "Stop" : "Run"
        notify.((run, sim_state))
    end

    return fig, run
end

function createSliders!(fig, parameters::NamedTuple)
    lsgrid = labelslidergrid!(
        fig, 
        [string.(keys(parameters))...],
        [values(parameters)...],
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
end

figure, run_var = createFigure(); figure

# createSliders!(figure, (σ = 10, ρ = 28, β = 8/3))

solver = RK4

points = Observable(Point3f.(Vector{Float64}.([eachcol(attractor.positions)...])))

scatter!(points, markersize = 2000)

while true
    # Main sim loop
    while run_var[]
        step!(attractor, solver)
        points[] = [eachcol(attractor.positions)...]
        notify(points)
        sleep(1/FPS)
    end
    sleep(0.0001)
end



