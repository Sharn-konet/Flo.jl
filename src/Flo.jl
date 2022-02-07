module Flo

using Revise

using GLMakie; GLMakie.activate!
using Makie.Colors
using LinearAlgebra: norm
using Statistics: quantile, mean, median, std
using Printf: @sprintf

include("Swarms.jl")
include("ExampleFunctions.jl")

using .Attractors
using .Swarms: Swarm, step!, RK4

include("FigureElements.jl")

const FPS = 144

"""
    findLimits(func, q)

Runs a short simulation to determine appropriate limits for displaying a particular function.
"""
function findLimits(func::Function; q = 0.001)::Tuple
    num_particles = 100

    attractor = Swarm(size = num_particles, positions = (rand(Float64, 3, num_particles) .-0.5) .* 2)
    history = Array{Float64}(undef, 3, num_particles, 5000)

    dims, particles, steps = map(x -> 1:size(history, x), [1,2,3])

    # Run a time-limited simulation
    for i in steps
        history[:,:,i] = attractor.positions
        step!(attractor, func, RK4)
    end

    high_quantile = x -> quantile(x, 1 - q)
    low_quantile = x -> quantile(x, q)

    upper_limits = mean([high_quantile(history[i,j,:]) for i in dims, j in particles], dims = 2)
    lower_limits = mean([low_quantile(history[i,j,:]) for i in dims, j in particles], dims = 2)

    limits = (Tuple([(zip(lower_limits, upper_limits)...)...]) .* 1.5)
end

"""
    findStats(func)

Simulates a function for a number of steps, gives statistics to normalise function on each dimension.
"""
function findStats(func::Function)::Tuple
    num_particles = 100

    attractor = Swarm(size = num_particles, positions = (rand(Float64, 3, num_particles) .-0.5) .* 2)
    history = Array{Float64}(undef, 3, num_particles, 5000)

    dims, particles, steps = map(x -> 1:size(history, x), [1,2,3])

    # Run a time-limited simulation
    for i in steps
        history[:,:,i] = attractor.positions
        step!(attractor, func, RK4)
    end

    # Reshape to create set of observations for each dim
    history = reshape(history, dims, length(particles)*length(steps))

    return (mean(history, dims = 2), std(history, dims = 2))
end

"""
    dropdownMapping(func_names)

Generates a Dictionary that cleans a set of function names and creates a mapping
from styled names => function names.
"""
function dropdownMapping(func_names::Vector{Symbol})
    dropdown_style = str -> replace(str, "_" => " ") |> titlecase |> str -> replace(str, " " => "-")
    styled_names = (dropdown_style ∘ string).(func_names)
    return (Dict ∘ zip)(styled_names, func_names)
end

"""
    getDefaultArgs(func)

Return all default arguments of a function and their associated variable names.
"""
function getDefaultArgs(func)
    arguments = code_lowered(func)[1].code[1].args
    arguments = filter(x -> isa(x, Real), arguments)

end

function refreshPoints(point::Real, limits::NTuple{2, Real})::Float64
    if (point < limits[1]) | (point > limits[2])
        return limits[1] + (limits[2] - limits[1])*rand(Float64)
    else
        return point
    end
end

function refreshPoints(positions::Vector{<:Real}, limits::NTuple{6, Real})
    for dim in 1:length(positions)
        positions[dim] = refreshPoints(positions[dim], limits[2*dim-1:2*dim])
    end

    return positions
end

function julia_main()::Cint
    μ = 0; σ = 1

    ode_func = Observable{Function}(Attractors.aizawa)
    attractor = Observable(Swarm(size = 2000, step_size = repeat([1e-2], 2000)))
    
    # Create main figure
    include("src/Theme.jl")
    limits = Observable(findLimits(ode_func[]))
    figure, ax, run_var = createFigure(limits = limits); display(figure)

    # Create slider to adjust speed
    speed_slider, _, _, layout = labelslider!(figure.scene, "Speed", exp10.(range(-4, -1, length = 400)), startvalue = 1e-2, format = x -> @sprintf("%.1e", x), snap =false)
    figure[3,2] = layout

    on(speed_slider.value) do speed
        attractor[].step_size .= speed
    end

    normalisation = Toggle(figure, active = false)
    label = Label(figure, lift(x -> x ? "Normalised" : "Unnormalised", normalisation.active))

    on(normalisation.active) do normalise
        if normalise
            μ, σ = findStats(eval(:($(ode_func[]))))
            ax.limits[] = (-3., 3., -3., 3., -3., 3.)
        else
            ax.limits[] = findLimits(eval(:($(ode_func[]))))
        end
    end

    figure[1:2, 3] = grid!(hcat(normalisation, label), tellheight = false)

    dropdown_dict = dropdownMapping(names(Flo.Attractors)[2:end])
    menu = createDropdown!(figure, [keys(dropdown_dict)...])

    on(menu.selection) do dropdown_item
        ode_func[] = eval(dropdown_dict[dropdown_item])
        if normalisation.active[]
            μ, σ = findStats(eval(:($(ode_func[]))))
        else
            ax.limits[] = findLimits(eval(:($(ode_func[]))))
        end
    end

    # createSliders!(figure, (σ = 10, ρ = 28, β = 8/3))
    solver = RK4

    points = Observable(Point3f.(Vector{Float64}.([eachcol(attractor[].positions)...])))

    scatter!(points, markersize = 0.025, markerspace = SceneSpace)

    while events(figure.scene).window_open.val
        # Main sim loop
        while run_var[] & events(figure.scene).window_open.val
            step!(attractor[], ode_func[], solver)
            if normalisation.active[]
                points[] = [(x -> refreshPoints(x, ax.limits[][2*row-1:2*row])).(((attractor[].positions .- μ) ./ σ)[row,:]) for row in 1:(size(attractor[].positions, 1))]
            else
                println(size([(x -> refreshPoints(x, ax.limits[][2*row-1:2*row])).(attractor[].positions[row, :]) for col in 1:(size(attractor[].positions, 2))]))
                println(size(points[]))
                points[] = [(x -> refreshPoints(x, ax.limits[][2*row-1:2*row])).(attractor[].positions[:, col]) for col in 1:(size(attractor[].positions, 2))]
            end
            sleep(1/FPS)
        end
        sleep(1e-12)
    end

    return 0
end

end