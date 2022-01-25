using GLMakie: Theme, Axis3, RGBAf

flo_theme = Theme(
    Axis3 = (
        backgroundcolor = :black,
        protrusions = (0,0,0,0),
        titlealign = :left,
        viewmode = :fit,
        limits = (-100,100,-100,100,-100,100)
    )
)

set_theme!(flo_theme)