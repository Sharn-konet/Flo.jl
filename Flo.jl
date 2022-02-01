using Revise

using GLMakie
import GLMakie: lines!
using Makie.Colors
using LinearAlgebra: norm
using ElasticArrays
using Statistics: quantile, mean, median
using DataStructures: CircularBuffer
import Base: push!

include("Swarms.jl")
include("ExampleFunctions.jl")

using .Attractors
using .Swarms: Swarm, step!, RK4

const FPS = 144

function createFigure(; fig::Makie.Figure = Figure(), limits)

    include("Theme.jl")

    ax = Axis3(fig[1:2,2], aspect = (1,1,1), autolimitaspect = true, limits = limits[])

    fig[3,1] = buttongrid = GridLayout(tellwidth = false)
    sim_state = Observable("Run")
    run_button = buttongrid[1,1] = Button(fig, label = sim_state)

    run = Observable(false)

    on(run_button.clicks) do clicks
        run[] = !run[]
        sim_state[] = run.val ? "Stop" : "Run"
        notify.((run, sim_state))
    end

    return fig, ax, run
end

function createDropdown!(fig, available_functions::Vector{Symbol})
    menu = Menu(fig, options = sort(string.(available_functions)))

    fig[1, 1] = vgrid!(Label(fig, "Function", width = nothing), menu, tellheight = false, width = events(figure.scene).window_area[].widths[1]/3)

    return menu
end

function createSliders!(fig, parameters::NamedTuple)
    lsgrid = labelslidergrid!(
        fig, 
        [string.(keys(parameters))...],
        [values(parameters)...],
        width = 350,
        tellheight = false
    )

    fig[2,2] = lsgrid.layout

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

function findLimits(func::Function; q = 0.001)::Tuple
    num_particles = 100

    attractor = Swarm(func, size = num_particles, positions = (rand(Float64, 3, num_particles) .-0.5) .* 2)
    history = Array{Float64}(undef, 3, num_particles, 10000)

    dims, particles, steps = map(x -> 1:size(history, x), [1,2,3])

    # Run a time-limited simulation
    for i in steps
        history[:,:,i] = attractor.positions
        step!(attractor, RK4)
    end

    high_quantile = x -> quantile(x, 1 - q)
    low_quantile = x -> quantile(x, q)

    upper_limits = mean([high_quantile(history[i,j,:]) for i in dims, j in particles], dims = 2)
    lower_limits = mean([low_quantile(history[i,j,:]) for i in dims, j in particles], dims = 2)

    limits = (Tuple([(zip(lower_limits, upper_limits)...)...]) .* 1.5)
end

function push!(buffer::Observable{CircularBuffer{Vector{Float64}}}, args...)
    push!(buffer[], args...)
end

# function lines!(buffer::Observable{CircularBuffer{Vector{Float64}}}, args...)
#     lines!()
# end

ode_func = aizawa
attractor = Observable(Swarm(ode_func))

limits = Observable(findLimits(ode_func))

figure, ax, run_var = createFigure(limits = limits); figure

menu = createDropdown!(figure, names(Main.Attractors))

on(menu.selection) do func
    attractor[] = eval(:(Swarm($(Symbol(func)))))
    ax.limits[] = findLimits(eval(:($(Symbol(func)))))
end

# createSliders!(figure, (σ = 10, ρ = 28, β = 8/3))
solver = RK4

points = Observable(Point3f.(Vector{Float64}.([eachcol(attractor[].positions)...])))

history = [Observable(CircularBuffer{Vector{Float64}}(50)) for _ in 1:1000]

scatter!(points, markersize = 0.025, markerspace = SceneSpace)

for _ in 1:2
    push!.(history, [eachcol(attractor[].positions)...])
    step!(attractor[], solver)
end

for hist in history
    lines!(hcat(hist[]...), linewidth = 1)
end

while events(figure.scene).window_open.val
    # Main sim loop
    while run_var[] & events(figure.scene).window_open.val
        push!.(history, [eachcol(attractor[].positions)...])
        step!(attractor[], solver)
        points[] = [eachcol(attractor[].positions)...]
        notify(points)
        if length(history[1][][1]) > 2 notify.(history) end
        sleep(1/FPS - 0.000208)
    end
    sleep(1e-12)
end



