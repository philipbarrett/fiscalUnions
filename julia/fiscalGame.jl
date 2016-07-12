#= fiscalGame.jl
Philip Barrett, pobarrett@gmail.com
11jul2016, Chicago

Sets up the game objects  =#

using Polygons

"""
Type defining the fiscal game.  Contains:
  * r: Interest rate
  * betta: Discount factor
  * delta: (vector of) discount govt factors
  * psi: Vector of labor supply elasticities
  * P: Transition matrix for states
  * T: Vector of tax levels
  * bgrid: The linspace grid of values for the debt
"""
type fiscalGame

  # Deep parameters
  r::Float64      # Interest rate
  betta::Float64  # Discount factor
  delta::Vector{Float64}  # Extra government discount factor
  psi::Vector{Float64}    # Frisch elasticity of labor supply
  chi::Vector{Float64}    # Leisure utility scale parameter

  # Stochastic process
  A::Matrix{Float64}      # Labor productivity (nSx2)
  g::Matrix{Float64}      # Government spending (nSx2)
  P::Matrix{Float64}      # Transition probabilities (nSxnS)
  nS::Int                 # Number of States

  # B grid
  bgrid::LinSpace         # Debt grid
  nb::Int                 # Length of b grid

  # DWL object
  dw::Vector{DWL}         # Deadweight loss

  # Game formulation
  nR::Int                    # Number of actions
  actions::Array{Float64,3}  # Matrix of action pairs (nR^2x2xnS)
  pdLoss::Array{Float64,3}   # Matrix of payoff pairs (nR^2x2xnS)
  surp::Matrix{Float64}       # Vector of primary defecits (nR^2xnS)

  # Algorithm details
  dirs::Matrix{Float64}        # Search objects (nDx2)
end

function fiscalGame( ; r=0.04, delta=[.95, .95], psi=[75, .75], chi=[7.0, 7.0],
                      rho=.5, A=-1, g=-1, P=-1, nR::Int=20, nb=150,
                      bmin=0, ndirs=32 )

  if ( A[1] < 0 || g[1] < 0 || P[1] < 0 )
    Atemp, gtemp, P = defaultStates()
    A = hcat( Atemp, Atemp )
    g = hcat( gtemp, gtemp )
  end
  nS = size(A)[1]
  vrho = [ rho 1-rho ]
      # Vector of weights

  dw =[ DWL( A=A[:,i], g=g[:,i], nR=nR, psi=psi[i],
                chi=chi[i], r=r, rho=vrho[i] ) for i in 1:2 ]
      # Deadwight loss object

  blim = 1 / r * minimum( dw[1].R[:,nR] - g[:,1] +
                            dw[2].R[:,nR] - g[:,2] )
      # Maximum joint borrowing
  bgrid = linspace( bmin, blim * .975, nb )

  actions = zeros( nR ^ 2, 2, nS )
  pdLoss =  zeros( nR ^ 2, 2, nS )
      # Initialize
  for i in 1:nS
    for j in 1:nR
      for k in 1:nR
        actions[ (j-1)*nR+k, :, i ] = [ dw[1].R[i,k] dw[2].R[i,j] ]
        pdLoss[ (j-1)*nR+k, :, i ] = [ dw[1].W[i,k] dw[2].W[i,j] ]
      end
    end
  end

  surp = [ sum( actions[i,:,j] ) - sum( g[j,:] ) for i in 1:(nR^2),
                                                          j in 1:nS ]
      # Primary defecits
  dirs = hcat( [ cos((i-1)/(ndirs-1)) for i in 1:ndirs ],
               [ sin((i-1)/(ndirs-1)) for i in 1:ndirs ] )

  fiscalGame( r, 1/(1+r), delta, psi, chi, A, g, P, nS, bgrid, nb,
              dw, nR, actions, pdLoss, surp, dirs )
end
