#= ploygonPlot.jl
Philip Barrett, pobarrett@gmail.com
27may2016, Chicago

Provides plotting functionality for the polygon class
=#

using Gadfly

function polyPlot( pts::Matrix )
  P = [ pts ; pts[1,:] ]
      # Wrap the matrix of points
  plot( P[:,1], P[:,2], lw=2 )
end
