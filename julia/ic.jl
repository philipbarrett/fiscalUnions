#= ic.jl
Philip Barrett, pobarrett@gmail.com
20jul2016, Chicago

Computes and imposes the incentive compatibility constraints  =#

"""
    maxCont( surp::Float64, pdLoss::Vector{Float64}, b::Float64,
                        r::Float64, bgrid::LinSpace, W::Array{Polygon,2}, Q::Vector{Float64}, betta_hat::Vector{Float64},
                        dirs::Matrix{Float64}, outer::Bool=true )
Computes the maximum continuation discounted DWL following a particular joint action.  Returns a triple
"""
function maxCont( surp::Float64, pdLoss::Vector{Float64}, b::Float64,
                  r::Float64, bgrid::LinSpace, W::Array{Polygon,2}, Q::Vector{Float64}, betta_hat::Vector{Float64},
                  dirs::Matrix{Float64}, outer::Bool=true )

end
