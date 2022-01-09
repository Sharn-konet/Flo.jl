using Revise

module Solvers

export RK4, DormandPrince

struct ButcherTableau
    α::Vector{Number}
    β::Vector{Number}
    γ::Matrix{Number}
    order::Int
    ButcherTableau(α, β, γ, order) = (length(α) == length(β)) & (size(γ, 1) == size(γ,2)) ? 
        new(α, β, γ, order) : error("Not a valid Butcher Tableau")
end

function Base.show(io::IO, bt::ButcherTableau)
    for i in 1:length(bt.β)
        print(string(bt.β[i], "|", bt.γ[i,:]...), "\n")
    end
    
    print("-+", ("-" for _ in size(bt.γ, 2))...)
    print(" |", bt.α...)
end

RK4 = ButcherTableau(
    [1/6, 1/3, 1/3, 1/6], # α
    [0, 1/2, 1/2, 1], # β
    [[0, 0,0, 0] [1/2, 0, 0, 0] [0, 1/2, 0, 0] [0, 0, 1, 0 ]], # γ
    4 # order
)

end