using Revise
using Test: @testset, @test

include("Interface.jl")
include("ExampleFunctions.jl")

using .Interface: @ODE
using .StrangeAttractors: yuwang, TSUCS2, lorenz, TSUCS1

u = [[1., 2., 3.] [4., 5., 6.] [7., 8., 9.]]

lorenz(u)

@testset "Interface" begin
    @testset "Three-Dimensional Functions" begin

        @testset "Typing" begin
            u = [[1., 2., 3.] [4., 5., 6.] [7., 8., 9.]]

            @test lorenz(u) ≈ [[10., 10., 10.] [23., 83., 125.] [-6., 4., 32.]]
        end

        @testset "Time Functionality" begin
            u = [[1., 2., 3.] [4., 5., 6.] [7., 8., 9.]]

            time_ODE = Interface.@ODE quote
                dx = t
                dy = x + y
                dz = y + z
            end

            @test length(methods(time_ODE)) == 1

            @test time_ODE(0., u)[:,1] ≈ Vector{Float64}([0, 0, 0])
            @test time_ODE(2., u)[:,1] ≈ Vector{Float64}([2, 2, 2])
        end
    end

    @testset "Two-Dimensional Functions" begin
        
    end
end

lorenz([[1.,2.,3.] [3.,4.,5.] [6.,7., 8.]], σ = 10, ρ = 28, β = 8/3)
TSUCS2(0., [[1,2,3] [3,4,5] [6,7, 8]], α = 1, δ = 2, ς = 4, ζ = 8, β = 5, ϵ = 6)
yuwang(0., [[1,2,3] [3,4,5] [6,7, 8]], α = 1, β = 2, σ = 4, δ = 8)

