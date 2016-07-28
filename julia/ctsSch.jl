#= ctsSch.jl
Philip Barrett, pobarrett@gmail.com
22jul2016, Chicago

Computes the direction search for the continuous action problem
=#

using Optim

function dirMax( bprime::Float64, b::Float64, chi::Vector, psi::Vector,
                  x1::Vector, x2::Vector, R1bds::Vector,
                  a::Vector, sumg::Float64,
                  rho::Float64, r::Float64, dir::Vector,
                  verbose::Bool=false )
  obj(R1) = - ( dir[1] * rho *
          w_eval( R1, chi[1], psi[1], x1[1], x2[1], a[1], rho, verbose ) +
               dir[2] * (1-rho) *
          w_eval( (1+r)*b + sumg - bprime - R1, chi[2], psi[2], x1[2],
          x2[2], a[2], 1-rho, verbose ) )::Float64
      # The objective function

  if verbose
    println("R1bds = ", R1bds )
    println("R2bds = ", (1+r)*b + sumg - bprime - R1bds )
    println("obj(R1[1]) = ", obj(R1bds[1]) )
    println("obj(R1[2]) = ", obj(R1bds[2]) )
  end

  res = optimize( obj, R1bds[1], R1bds[2] )
      # The result
  R1 = (res.minimum)::Float64
  R2 = ((1+r)*b + sumg - bprime - R1)::Float64
      # Revenues
  W1 = rho * w_eval( R1, chi[1], psi[1], x1[1], x2[1], a[1], rho )::Float64
  W2 =(1-rho) *  w_eval( R2, chi[2], psi[2], x1[2], x2[2], a[2], 1-rho )::Float64
      # Welfare
  dist = - obj(R1)::Float64
      # Distance in search direction
  return R1, R2, W1, W2, dist
end
