using GLMakie: Observable, Figure

"""
    createFigure([; fig, limits])

Adds a 3D axis and a run button to the figure provided through `fig`.
If a figure is not supplied, a new figure is generated within the function.

#### Parameters:
* fig: Figure object to apply Axis3 to.
* limits: Used to set axis limits for Axis3 objects.
"""
function createFigure(; fig::Makie.Figure = Figure(), limits::Observable{NTuple{6, Float64}})

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


"""
    createDropdown!(fig, available_functions)

Creates a dropdown menu of all available functions and applies it to the provided figure.
Returns the menu as an object which can then be observed for selections.
"""
function createDropdown!(fig::Figure, available_functions::Vector{String})
    menu = Menu(fig, options = sort(available_functions))

    fig[1, 1] = vgrid!(Label(fig, "Function", width = nothing), menu, tellheight = false, width = events(fig.scene).window_area[].widths[1]/3)

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