using Revise

include("ExampleFunctions.jl")
include("interface.jl")

using .Interface
using .StrangeAttractors: yuwang, TSUCS2, lorenz, TSUCS1

println(methods(TSUCS2))
lorenz.func(0., [[1.,2.,3.] [3.,4.,5.] [6.,7., 8.]], σ = 10, ρ = 28, β = 8/3)
TSUCS2.func(0., [[1,2,3] [3,4,5] [6,7, 8]], α = 1, δ = 2, ς = 4, ζ = 8, β = 5, ϵ = 6)
yuwang.func(0., [[1,2,3] [3,4,5] [6,7, 8]], α = 1, β = 2, σ = 4, δ = 8)

