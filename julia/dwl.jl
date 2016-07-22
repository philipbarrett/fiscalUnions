#= dwl.jl
Philip Barrett, pobarrett@gmail.com
07jul2016, Chicago

Creates a deadweight loss object
=#

using Gadfly, Colors

type DWL
  x1::Vector{Float64}   # Lower bound on leisure
  x2::Vector{Float64}   # Upper bound on leisure
  x::Array{LinSpace{Float64},1}  # Leisure grids
  tau::Matrix{Float64}  # Tax rates
  R::Matrix{Float64}    # Revenue
  c::Matrix{Float64}    # Consumption
  W::Matrix{Float64}    # Deadweight loss
  blim::Float64         # Upper bound on debt
  nS::Int               # Number of exogenous states
  nR::Int               # Number of R values
end

function DWL( ; A::Vector=[0.0], g::Vector=[0.0], nR::Int=1,
              psi::Float64 = 1, chi::Float64 = 1, r::Float64 = .03,
              rho::Float64=1.0 )

  m::Int = size(A)[1]
      # Number of states
  if chi < maximum(A)
    error("H'(0) >= max_s A(s) not satisfied")
  end
    # Check that x=0 is never a solution

  x1 = [ 1 - ( A[i] / chi ) ^ psi for i in 1:m ]
      # Lower bound for leisure
  x2 = [ 1 - ( A[i] * psi / ( chi * ( 1 + psi ) ) ) ^ psi for i in 1:m ]
      # Upper bound for leisure
  x = [ linspace( x1[i], x2[i], nR ) for i in 1:m ]
      # The range of x values
  tau = [ max( 0, 1 - chi * ( 1 - x[i][j] ) ^ ( 1 / psi ) / A[i] )
            for i in 1:m, j in 1:nR ]
      # Taxes.  Should be same for all rows due to functional form.
      # Need max to prevent tiny tiny negative taxes.
  R = [ tau[i,j] * ( 1 - x[i][j] ) * rho * A[i] for i in 1:m, j in 1:nR ]
      # Revenue
  c = [ ( 1 - x[i][j] ) * A[i] - g[i] / rho for i in 1:m, j in 1:nR ]
      # Per capita consumption
  W = [ x[i][j] * rho * A[i] - rho * chi * ( 1 - x[i][j] ) ^
          ( 1 + 1 / psi ) / ( 1 + 1 / psi ) for i in 1:m, j in 1:nR ]
      # Total deadweight loss (i.e. not per capita)
  blim = 1 / r * minimum( R[:,nR] - g )
      # Nautral borrowing limit
  return DWL( x1, x2, x, tau, R, c, W, blim, m, nR )
      # The dwl object
end

function DWLplot( dw::DWL, chart="tR" )

  color_vec = [ "magenta" "red" "blue" "black" "green" "cyan" "orange"  ]

  if chart == "tR"
    x = dw.tau
    y = dw.R
    xlab = "tau"
    ylab = "Revenue"
  end

  if chart == "RW"
    x = dw.R
    y = dw.W
    xlab = "Revenue"
    ylab = "Deadweight loss"
  end

  if chart == "xc"
    y = dw.c
    x = [ dw.x[i][j] for i in 1:(size(y)[1]), j in 1:(size(y)[2]) ]
    xlab = "Leisure"
    ylab = "Consumption"
  end

  if chart == "xR"
    y = dw.R
    x = [ dw.x[i][j] for i in 1:(size(y)[1]), j in 1:(size(y)[2]) ]
    xlab = "Leisure"
    ylab = "Revenue"
  end

  Gadfly.plot( [ layer( x=x[i,:], y=y[i,:], Geom.line,
             Theme(default_color=color(parse(Colorant, color_vec[i%7+1]))) )
             for i in 1:(size(y)[1]) ]...,
               Guide.xlabel(xlab), Guide.ylabel(ylab) )

end
