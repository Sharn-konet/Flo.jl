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

    return ([mean(history, dims = 2)...], [std(history, dims = 2)...])
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

# function outsideLimits(point::Real, limits::NTuple{2, <:Real})::Float64
#     if (point < limits[1]) | (point > limits[2])
#         return limits[1] + (limits[2] - limits[1]) * rand(Float64)
#     else
#         return point
#     end
# end

function refreshPoints(position::Vector{<:Real}, limits::NTuple{6, <:Real})
    outsideLimits(point::Real, limits::NTuple{2, <:Real}) = (point < 2*limits[1]) | (point > 2*limits[2])

    dims = length(position)
    if any([outsideLimits(position[dim], limits[2*dim-1:2*dim]) for dim in 1:dims])
        new_position = rand(Float64, dims)
        for dim in 1:dims
            new_position[dim] = (limits[2*dim-1] + (limits[2*dim] - limits[2*dim-1]) * new_position[dim])*(3/4)
        end
        return new_position
    else
        return position
    end
end

function julia_main()::Cint

    ode_func = Observable{Function}(Attractors.aizawa)
    attractor = Observable(Swarm(size = 1000, step_size = repeat([1e-2], 1000)))
    
    # Create main figure
    include("src/Theme.jl")
    limits = Observable(findLimits(ode_func[]))
    figure, ax, run_var = createFigure(limits = limits); display(figure)
    # ax.limits[] = findLimits(eval(:($(ode_func[]))))

    # Create slider to adjust speed
    speed_slider, _, _, layout = labelslider!(figure.scene, "Speed", exp10.(range(-4, -1, length = 400)), startvalue = 1e-2, format = x -> @sprintf("%.1e", x), snap =false)
    figure[3,2] = layout

    on(speed_slider.value) do speed
        attractor[].step_size .= speed
    end

    normalisation = Toggle(figure, active = false)
    label = Label(figure, lift(x -> x ? "Normalised" : "Unnormalised", normalisation.active))

    μ = lift((func, norm) -> norm ? findStats(eval(:($(func))))[1] : [0., 0., 0.], ode_func, normalisation.active)
    σ = lift((func, norm) -> norm ? findStats(eval(:($(func))))[2] : [1., 1., 1.], ode_func, normalisation.active)
    
    limits = lift((func, norm) -> norm ? (-3.,3.,-3.,3.,-3.,3.) : findLimits(eval(:($(func)))), ode_func, normalisation.active)
    
    on(limits) do limit
        ax.limits[] = limit
    end

    figure[1:2, 3] = grid!(hcat(normalisation, label), tellheight = false)

    dropdown_dict = dropdownMapping(names(Flo.Attractors)[2:end])
    menu = createDropdown!(figure, [keys(dropdown_dict)...])

    on(menu.selection) do dropdown_item
        ode_func[] = eval(dropdown_dict[dropdown_item])
    end 

    points = @lift begin
        positions = $(normalisation.active) ? ($(attractor).positions .- $μ) ./ $σ : $(attractor).positions
        positions = ((x -> refreshPoints(x, $(ax.limits))) ∘ Vector{Float64}).([eachcol(positions)...])
        return Point3f.(positions)
    end

    scatter!(points, markersize = 0.025, markerspace = SceneSpace)

    while events(figure.scene).window_open.val
        # Main sim loop
        while run_var[] & events(figure.scene).window_open.val
            step!(attractor[], ode_func[], RK4)
            for j in 1:size(attractor[].positions, 2)
                attractor[].positions[:,j] = refreshPoints(attractor[].positions[:,j], ax.limits[])
            end
            notify.((points, attractor, limits))
            sleep(1/FPS)
        end
        sleep(1e-12)
    end

    return 0
end

end