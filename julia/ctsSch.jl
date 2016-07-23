#= ctsSch.jl
Philip Barrett, pobarrett@gmail.com
22jul2016, Chicago

Computes the direction search for the continuous action problem
=#

using Optim

function dirMax( bprime::Float64, b::Float64, chi::Vector, psi::Vector,
                  x1::Vector, x2::Vector, R1bds::Vector,
                  a::Vector, sumg::Float64,
                  rho::Float64, r::Float64, dir::Vector )
  obj(R1) = - ( dir[1] * rho *
                  w_eval( R1, chi[1], psi[1], x1[1], x2[1], a[1] ) +
               dir[2] * (1-rho) *
                  w_eval( (1+r)*b + sumg - bprime - R1, chi[2], psi[2], x1[2], x2[2], a[2] ) )
      # The objective function
  res = optimize( obj, R1bds[1], R1bds[2] )
      # The result
  R1 = res.minimum
  R2 = (1+r)*b + sumg - bprime - R1
      # Revenues
  dist = obj(R1)
      # Distance in search direction
  return R1, R2, dist
end
