module Flo

using Revise

using GLMakie; GLMakie.activate!
using Makie.Colors
using LinearAlgebra: norm
using Statistics: quantile, mean, median
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
    history = Array{Float64}(undef, 3, num_particles, 10000)

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
    dropdownMapping(func_names)

Generates a Dictionary that cleans a set of function names and creates a mapping
from styled names => function names.
"""
function dropdownMapping(func_names::Vector{Symbol})
    dropdown_style = str -> replace(str, "_" => " ") |> titlecase |> str -> replace(str, " " => "-")
    styled_names = (dropdown_style ∘ string).(func_names)
    return (Dict ∘ zip)(styled_names, func_names)
end

function julia_main()::Cint

    ode_func = Observable{Function}(Attractors.aizawa)

    attractor = Observable(Swarm(size = 2000, step_size = repeat([1e-2], 2000)))

    limits = Observable(findLimits(ode_func[]))

    include("src/Theme.jl")

    figure, ax, run_var = createFigure(limits = limits); display(figure)

    speed_slider, _, _, layout = labelslider!(figure.scene, "Speed", exp10.(range(-4, -1, length = 50)), startvalue = 1e-2, format = x -> @sprintf("%.1e", x))

    figure[3,2] = layout

    on(speed_slider.value) do speed
        attractor[].step_size .= speed
    end

    dropdown_dict = dropdownMapping(names(Flo.Attractors)[2:end])
    menu = createDropdown!(figure, [keys(dropdown_dict)...])

    on(menu.selection) do dropdown_item
        ode_func[] = eval(dropdown_dict[dropdown_item])
        ax.limits[] = findLimits(eval(:($(ode_func[]))))
    end

    # createSliders!(figure, (σ = 10, ρ = 28, β = 8/3))
    solver = RK4

    points = Observable(Point3f.(Vector{Float64}.([eachcol(attractor[].positions)...])))

    scatter!(points, markersize = 0.025, markerspace = SceneSpace)

    while events(figure.scene).window_open.val
        # Main sim loop
        while run_var[] & events(figure.scene).window_open.val
            step!(attractor[], ode_func[], solver)
            points[] = [eachcol(attractor[].positions)...]
            notify(points)
            sleep(1/FPS)
        end
        sleep(1e-12)
    end

    return 0
end

end