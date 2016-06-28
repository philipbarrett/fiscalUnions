#= autarkyPlot.jl
Philip Barrett, pobarrett@gmail.com
25may2016, Chicago

Plots an AutarkySol object
=#

using Gadfly, Colors

"""
    plot_as( x::LinSpace, y::Matrix{Float64} )
Core plotting function for the autarky solutions
"""
function plot_as( x::LinSpace, y::Matrix{Float64} )
  color_vec = [ "magenta" "red" "blue" "black" "green" "cyan" "orange"  ]
  # Theme(default_color=)?
  Gadfly.plot( [ layer( x=x, y=y[:,i], Geom.line,
             Theme(default_color=color(parse(Colorant, color_vec[i%7+1]))) )
             for i in 1:(size(y)[2]) ]... )
end


"""
    plot_as( as::AutarkySol, part::AbstractString="V" )
Core plotting function for the autarky solutions
"""
function plot_as( as::AutarkySol, part::AbstractString="V" )
  if part == "V"
    y = as.V
  end
      # Value function
  if part == "g"
    y = as.g
  end
      # Government consumption
  if part == "b"
    y = as.bprime
  end
      # Continuation debt
  plot_as( as.am.bgrid, y )
      # Plot
end

# """
#     plot_sim( sim::Matrix, pds=1:100 )
# Plots (some periods of) a simulation
# """
# function plot_sim( sim::Matrix, pds=1:100 )
#
# end
