#= dwl.jl
Philip Barrett, pobarrett@gmail.com
07jul2016, Chicago

Creates a deadweight loss object
=#

function dwl( A::Vector, g::Vector, n::Int, psi::Float64,
              chi::Float64, r::Float64, rho::Float64=1.0 )

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
  x = [ linspace( x1[i], x2[i], n ) for i in 1:m ]
      # The range of x values
  tau = [ max( 0, 1 - chi * ( 1 - x[i][j] ) ^ ( 1 / psi ) / A[i] )
            for i in 1:m, j in 1:n ]
      # Taxes.  Should be same for all rows due to functional form.
      # Need max to prevent tiny tiny negative taxes.
  R = [ tau[i,j] * ( 1 - x[i][j] ) * rho * A[i] for i in 1:m, j in 1:n ]
      # Revenue
  c = [ ( 1 - x[i][j] ) * A[i] - g[i] / rho for i in 1:m, j in 1:n ]
      # Per capita consumption
  W = [ x[i][j] * rho * A[i] - rho * chi * ( 1 - x[i][j] ) ^
          ( 1 + 1 / psi ) / ( 1 + 1 / psi ) for i in 1:m, j in 1:n ]
      # Total deadweight loss (i.e. not per capita)
  blim = 1 / r * minimum( R[:,n] - g )
      # Nautral borrowing limit

  return x1, x2, x, tau, R, c, W, blim
end

### TODO: ADD PLOTTING HERE: (tau,R), (x,c), (x,R), (R,W)
