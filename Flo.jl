using Revise

using GLMakie
using Makie.Colors
using LinearAlgebra: norm
using ElasticArrays

include("Swarms.jl")
include("ExampleFunctions.jl")

using .StrangeAttractors: lorenz, TSUCS2, yuwang, aizawa, lorenz_mod_2, TSUCS1
using .Swarms: Swarm, step!, RK4

const FPS = 144

ode_func = lorenz
attractor = Swarm(ode_func)

function createFigure(; fig::Makie.Figure = Figure(), limits)

    include("Theme.jl")

    ax = Axis3(fig[1,1], aspect = (1,1,1), autolimitaspect = true, limits = limits)

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

function findLimits(func::Function)::Tuple
    attractor = Swarm(func, size = 2, positions = (rand(Float64, 3, 2) .-1) .* 1e-6)
    prev_positions = attractor.positions
    history = ElasticArray{Float64}(undef, 3, 2, 0)
    append!(history, prev_positions)
    step!(attractor, RK4)
    while norm(attractor.positions[:,1] .- attractor.positions[:,2]) < 60
        prev_positions = attractor.positions
        append!(history, prev_positions)
        step!(attractor, RK4)
    end

    maxInHistory = history -> maximum(history, dims = 3)
    minInHistory = history -> minimum(history, dims = 3)

    maxInSwarm = swarm_positions -> maximum(swarm_positions, dims = 2)
    minInSwarm = swarm_positions -> minimum(swarm_positions, dims = 2)

    upper_limits = (maxInSwarm ∘ maxInHistory)(history)
    lower_limits = (minInSwarm ∘ minInHistory)(history)

    limits = Tuple([(zip(lower_limits, upper_limits)...)...]) .* 1.2
end

limits = findLimits(ode_func)

figure, run_var = createFigure(limits = limits); figure

# createSliders!(figure, (σ = 10, ρ = 28, β = 8/3))

solver = RK4

points = Observable(Point3f.(Vector{Float64}.([eachcol(attractor.positions)...])))

scatter!(points, markersize = 0.05, markerspace = SceneSpace)

while events(figure.scene).window_open.val
    # Main sim loop
    while run_var[] & events(figure.scene).window_open.val
        step!(attractor, solver)
        points[] = [eachcol(attractor.positions)...]
        notify(points)
        sleep(1/FPS - 0.000208)
    end
    sleep(1e-12)
end



