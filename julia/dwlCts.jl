#= dwlCts.jl
Philip Barrett, pobarrett@gmail.com
22jul2016, Chicago

Creates a deadweight loss object for the continuous-action roblem
=#

using Gadfly, Colors, Roots

type DWLC

  nS::Int64                   # Number of exogenous states
  nb::Int64                   # Number of endogenous states

  chi::Float64                # Number of exogenous states
  psi::Float64                # Frisch elasticity of labor supply

  blim::Float64               # Upper bound on debt
  bgrid::LinSpace{Float64}    # Debt grid

  x1::Matrix{Float64}         # Lower bound on leisure
  x2::Vector{Float64}         # Upper bound on leisure
  R1::Matrix{Float64}         # Lower bound on revenue
  R2::Vector{Float64}         # Upper bound on revenue
  tau1::Matrix{Float64}       # Lower bound on taxes
  tau2::Vector{Float64}       # Upper bound on taxes

end

function dwlc( ; nS::Int=1, nb::Int=1, A::Vector=[0.0],
              g::Vector=[0.0],
              psi::Float64 = .75, chi::Float64 = 2.0,
              r::Float64 = .03, bmin::Float64=0.0 )

  # if chi < minimum(A)
  #   error("H'(1) < min_s A(s) not satisfied")
  # end
  #   # Check that x=1 is never a solution
  # NOT SURE WHAT THE RIGHT CHECK IS!!

  x2root = [ function(x)
              chi*x^(-1/psi) - A[i] + chi/psi*(1-x)*x^(-1-1/psi) end
                for i in 1:nS ]
      # Error on x2 equation
  x2 = [ fzero( x2root[i],0,1)::Float64 for i in 1:nS ]
      # Solve for the root
  R2 = ( A - chi * x2 .^ ( -1 / psi ) ) .* ( 1 - x2 )
      # Corresponding revenue
  blim = 1 / r * minimum( R2 - g )
      # The natural upper bound on debt
  bgrid = linspace( bmin, blim, nb )
      # The grid of debt levels
  R1 = [ (1+r)*bgrid[j] + g[i] - blim for i in 1:nS, j in 1:nb ]
      # Lowest revenue possible
  x1root = [ function(x)
              (A[i]-chi*x^(-1/psi))*(1-x) - R1[i,j] + 1e-14 end
              for i in 1:nS, j in 1:nb ]
      # Error on x2 equation.  Add a tiny amount to make
      # sure can find the root
  x1 = [ fzero( x1root[i,j],0,x2[i]) for i in 1:nS, j in 1:nb ]
      # Solve for the root
  tau2 = 1 - chi * x2.^( - 1 / psi ) ./ A
  tau1 = 1 - chi * x1.^( - 1 / psi ) ./ ( A * ones(1,nb) )
      # Range of taxes
  return DWLC( nS, nb, chi, psi, blim, bgrid, x1, x2, R1, R2,
                  tau1, tau2 )
end

function w_eval( R::Float64, chi::Float64, psi::Float64,
                    x1::Float64, x2::Float64, a::Float64 )
  x_err(x)=(a-chi*x^(-1/psi))*(1-x) - R + 1e-14
      # The error funciton on x
  x_rt = fzero( x_err, x1, x2 )
      # The root
  return x_rt * a - chi * x_rt ^ ( 1 - 1 / psi ) / ( 1 - 1 / psi )
end
