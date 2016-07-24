#= dwlCts.jl
Philip Barrett, pobarrett@gmail.com
22jul2016, Chicago

Creates a deadweight loss object for the continuous-action problem
=#

using Gadfly, Colors, Roots

type DWLC

  nS::Int64                   # Number of exogenous states
  nb::Int64                   # Number of endogenous states

  chi::Float64                # Number of exogenous states
  psi::Float64                # Frisch elasticity of labor supply

  blim::Float64               # Upper bound on debt
  bgrid::LinSpace{Float64}    # Debt grid

  x1::Array{Float64,3}        # Lower bound on leisure
  x2::Vector{Float64}         # Upper bound on leisure
  R1::Array{Float64,3}        # Lower bound on revenue
  R2::Vector{Float64}         # Upper bound on revenue
  tau1::Array{Float64,3}      # Lower bound on taxes
  tau2::Vector{Float64}       # Upper bound on taxes

end

type DWLC2

  nS::Int64                   # Number of exogenous states
  nb::Int64                   # Number of endogenous states

  blim::Float64               # Upper bound on debt
  bgrid::LinSpace{Float64}    # Debt grid
  bprimeposs::Array{BitArray{1},2}    # Poss. cont. debt
  nposs::Matrix{Int}          # Number of poss cont debt levels

  xlow::Array{Array{Float64,1},3}     # Lower bound on leisure
  xhigh::Matrix{Float64}              # Upper bound on leisure
  Rlow::Array{Array{Float64,1},3}     # Lower bound on revenue
  Rhigh::Matrix{Float64}              # Upper bound on revenue
  taulow::Array{Array{Float64,1},3}   # Lower bound on taxes
  tauhigh::Matrix{Float64}            # Upper bound on taxes

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

function dwlc2( ; nS::Int=1, nb::Int=1,
            A::Matrix=[0.0 0.0], g::Matrix=[0.0 0.0],
            psi::Vector = [.75, .75], chi::Vector = [2.0, 2.0],
            r::Float64 = .03, bmin::Float64=0.0, rho::Float64=.5 )

  vrho = [ rho, 1-rho ]
      # Relative sizes

  ## UPPER BOUNDS ##
  xhighroot = [ function(x)
              chi[j]*x^(-1/psi[j]) - A[i,j] +
                chi[j]/psi[j]*(1-x)*x^(-1-1/psi[j]) end
                for i in 1:nS, j in 1:2 ]
      # Error on xhigh equation
  xhigh = [ fzero( xhighroot[i,j],0,1)::Float64
                      for i in 1:nS, j in 1:2 ]
      # Solve for the root
  Rhigh = ( A - ( ones(nS,1) * chi' ) .* xhigh .^
            ( -1 ./ ( ones(nS,1) * psi' ) ) ) .*
            ( 1 - xhigh ) .* ( ones(nS) * vrho' )
      # Corresponding revenue
  tauhigh = [ 1 - chi[k] * xhigh[i,k]^( - 1 / psi[k] ) / A[i,k]
                for i in 1:nS, k in 1:2 ]
      # Upper bound on taxes

  ## DEBT BOUNDS AND POSSIBILTIES ##
  blim = 1 / r * minimum( sum(Rhigh - g, 2 ) ) - 1e-10
      # The natural upper bound on debt
  bgrid = linspace( bmin, blim, nb )
      # The grid of debt levels
  bprimeposs = [ ((1+r)*bgrid[j] + sum(g[i,:]) - bgrid + 1e-14 .<=
                      sum(Rhigh[i,:]))::BitArray{1}
                  for i in 1:nS, j in 1:nb ]
      # Possible bprime from (s,b)
  nposs = [ sum(bprimeposs[i,j])::Int for i in 1:nS, j in 1:nb ]
      # Number of possible actions

  ## LOWER BOUNDS ##
  Rlow = [ ((1+r)*bgrid[j] + sum(g[i,:]) - bgrid[bprimeposs[i,j]] -
              Rhigh[i,3-k])::Vector{Float64}
              for i in 1:nS, j in 1:nb, k in 1:2 ]
      # Lowest revenue possible *given* a continuation debt
      # level.  Only returns the lower limit when the proposed
      # continuation debt level is possible.  Remember to
      # subtract the *other* country's revenue from the
      # outstanding debt.
  xlowroot = [
      [ function(x)
        (A[i,k]-chi[k]*x^(-1/psi[k]))*(1-x)*vrho[k] -
        Rlow[i,j,k][l] + 1e-14 end for l in 1:nposs[i,j] ]
              for i in 1:nS, j in 1:nb, k in 1:2 ]
      # Error on x1 equation.  Add a tiny amount to make
      # sure can find the root
  xlow = [ [ fzero( xlowroot[i,j,k][l],0,xhigh[i,k])::Float64
                for l in 1:nposs[i,j] ]
                for i in 1:nS, j in 1:nb, k in 1:2 ]
      # Solve for the root
  taulow = [ 1 - chi[k]*xlow[i,j,k].^(-1/psi[k]) / A[i,k]
                for i in 1:nS, j in 1:nb, k in 1:2 ]
      # Range of taxes

  return DWLC2( nS, nb, blim, bgrid, bprimeposs, nposs,
            xlow, xhigh, Rlow, Rhigh, taulow, tauhigh )
end


function w_eval( R::Float64, chi::Float64, psi::Float64,
                    xlow::Float64, xhigh::Float64, a::Float64 )
  x_err(x) = (a-chi*x^(-1/psi))*(1-x) - R + 1e-14
      # The error function on x
# println( "xlow=", xlow )
# println( "xhigh=", xhigh )
# println( "x_err(xlow)=", x_err(xlow) )
# println( "x_err(xhigh)=", x_err(xhigh) )
  x_rt = fzero( x_err, xlow, xhigh )
      # The root
  return x_rt * a - chi * x_rt ^ ( 1 - 1 / psi ) / ( 1 - 1 / psi )
end
