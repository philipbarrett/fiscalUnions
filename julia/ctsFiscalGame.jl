#= ctsFiscalGame.jl
Philip Barrett, pobarrett@gmail.com
26jul2016, Chicago

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
type CtsFiscalGame

  # Deep parameters
  r::Float64      # Interest rate
  betta::Float64  # Discount factor
  delta::Vector{Float64}  # Extra government discount factor
  psi::Vector{Float64}    # Frisch elasticity of labor supply
  chi::Vector{Float64}    # Leisure utility scale parameter
  rho::Float64            # Relative size of country 1

  # Stochastic process
  A::Matrix{Float64}      # Labor productivity (nSx2)
  g::Matrix{Float64}      # Government spending (nSx2)
  P::Matrix{Float64}      # Transition probabilities (nSxnS)
  nS::Int                 # Number of States

  # B grid
  bgrid::LinSpace         # Debt grid
  nb::Int                 # Length of b grid

  # DWL object
  dw::DWLC2               # Deadweight loss

  # Useful summations
  gSum::Vector{Float64}       # Joint expenditure (nS)

  # Algorithm details
  ndirs::Int                  # No. search directions
  par::Bool                   # Parallel execution flag

end

function ctsfiscalgame( ; r=0.04, delta=[.95, .95], psi=[75, .75],
                      chi=[7.0, 7.0], rho=.5, A=-1, g=-1, P=-1,
                      nb=150, bmin=0, ndirs=24, par=false )

  # if ( A[1] < 0 || g[1] < 0 || P[1] < 0 )
  #   Atemp, gtemp, P = defaultStates()
  #   A = hcat( Atemp, Atemp )
  #   g = hcat( gtemp, gtemp )
  # end
  nS = size(A)[1]
  vrho = [ rho 1-rho ]
      # Vector of weights

  dw = dwlc2( nS=nS, nb=nb, A=A, g=g, psi=psi, chi=chi,
                r=r, bmin=bmin, rho=rho )
      # Deadwight loss object

  blim = dw.blim
      # Maximum joint borrowing
  bgrid = dw.bgrid
      # Borrowing grid
  gSum = vec( sum( g, 2 ) )
      # Sum of expenditure

  CtsFiscalGame( r, 1/(1+r), delta, psi, chi, rho, A, g, P, nS, bgrid,
                  nb, dw, gSum, ndirs, par )
end

function initGame( cfg::CtsFiscalGame )

  vrho = [ rho, 1-rho ]
      # Vector of rhos
  Rmin = [ minimum( cfg.dw.Rlow[i,j,k] )
              for i in 1:cfg.nS, j in 1:cfg.nb, k in 1:2 ]
      # Minimum revenue
  xMin = [ minimum( cfg.dw.xlow[i,j,k] )
              for i in 1:cfg.nS, j in 1:cfg.nb, k in 1:2 ]
      # Minimum leisure
  wMin = vrho .*
          [ minimum(
            [ w_eval( Rmin[iS, ib, i], cfg.chi[i], cfg.psi[i],
              xMin[iS,ib,i], cfg.dw.xhigh[iS,i], cfg.A[iS,i],
              vrho[i] )
              for iS in 1:cfg.nS, ib in 1:cfg.nb ] ) for i in 1:2 ]
  wMax = vrho .*
    maximum( [ w_eval( cfg.dw.Rhigh[i,j], cfg.chi[j],
          cfg.psi[j], minimum(xMin[i,:]), cfg.dw.xhigh[i,j],
          cfg.A[i,j], vrho[j] )::Float64
          for i in 1:cfg.nS, j in 1:2 ], 1)
      # Maximum welfare
  U = Polygon( pts=[ wMin[1] wMax[2]
                      wMin[1] wMin[2]
                      wMax[1] wMin[2]
                      wMax[1] wMax[2] ] )
      # The set f period values U
  W = [ U::Polygon for i in 1:cfg.nS, j in 1:cfg.nb ]
      # The period payoff polygons.
  # ndirs = [ cfg.ndirsl::Int for i in 1:cfg.nS, j in 1:cfg.nb ]
      # Matrix of number of directions for each state
  # return W, ndirs
  return W
end

function pdPayoffs( cfg::CtsFiscalGame, dirs::Matrix, outer::Bool=true )

  ndirs = size(dirs)[1]
      # Number of search directions
  ibprimeidx = [ find( cfg.dw.bprimeposs[iS,ib] )
                  for iS in 1:cfg.nS, ib in 1:cfg.nb ]
      # Array of indices for bprime

  if outer

    # for iS in 1:cfg.nS, ib in 1:cfg.nb
    #   for ibprime in 1:cfg.dw.nposs[iS,ib], idir in 1:ndirs
    #     println( "(iS,ib,ibprime,idir)=", (iS,ib,ibprime,idir) )
    #     dm = dirMax( cfg.bgrid[ibprimeidx[iS,ib][ibprime]],
    #       cfg.bgrid[ib],
    #       [ cfg.dw.Rlow[iS,ib,i][ibprime] for i in 1:2 ],
    #       vec( cfg.dw.Rhigh[iS, :] ), cfg.gSum[iS], cfg.rho, cfg.r,
    #       vec(dirs[idir,:]),
    #       [ cfg.dw.apx_coeffs[iS,i]::Array{Float64,1} for i in 1:2 ],
    #       [ cfg.dw.apx_N[iS,i]::Array{Float64,1} for i in 1:2 ],
    #       false )[5]
    #   end
    # end

    dists = [ [ [
              dirMax( cfg.bgrid[ibprimeidx[iS,ib][ibprime]],
                cfg.bgrid[ib],
                [ cfg.dw.Rlow[iS,ib,i][ibprime] for i in 1:2 ],
                vec( cfg.dw.Rhigh[iS, :] ), cfg.gSum[iS], cfg.rho, cfg.r,
                vec(dirs[idir,:]),
                [ cfg.dw.apx_coeffs[iS,i]::Array{Float64,1} for i in 1:2 ],
                [ cfg.dw.apx_N[iS,i]::Array{Float64,1} for i in 1:2 ]
                )[5] for idir in 1:ndirs ]
                for ibprime in 1:cfg.dw.nposs[iS,ib] ]
                for iS in 1:cfg.nS, ib in 1:cfg.nb ]
                # vec( cfg.dw.xhigh[iS,:]),
                # [ cfg.dw.Rlow[iS,ib,1][ibprime], cfg.dw.Rhigh[iS,1] ],
                # vec(cfg.A[iS,:]), cfg.gSum[iS], cfg.rho, cfg.r,
                # vec(dirs[idir,:]) )[5] for idir in 1:ndirs ]
    ret = [ [ Polygon( dirs=dirs, dists=dists[iS,ib][ibprime] )::Polygon
                for ibprime in 1:cfg.dw.nposs[iS,ib] ]::Array{Polygon,1}
                for iS in 1:cfg.nS, ib in 1:cfg.nb ]
  else
    ptsArray = [ [ [
              dirMax( cfg.bgrid[ibprimeidx[iS,ib][ibprime]],
                cfg.bgrid[ib],
                [ cfg.dw.Rlow[iS,ib,i][ibprime] for i in 1:2 ],
                vec( cfg.dw.Rhigh[iS, :] ), cfg.gSum[iS], cfg.rho, cfg.r,
                vec(dirs[idir,:]),
                [ cfg.dw.apx_coeffs[iS,i]::Array{Float64,1} for i in 1:2 ],
                [ cfg.dw.apx_N[iS,i]::Array{Float64,1} for i in 1:2 ]
                )[3:4] for idir in 1:ndirs ]
                for ibprime in 1:cfg.dw.nposs[iS,ib] ]
                for iS in 1:cfg.nS, ib in 1:cfg.nb ]
    pts =  [ [ [ ptsArray[iS,ib][ibprime][idir][i]::Float64
                    for idir in 1:ndirs, i in 1:2 ]
                    for ibprime in 1:cfg.dw.nposs[iS,ib] ]
                    for iS in 1:cfg.nS, ib in 1:cfg.nb ]
    ret = [ [ Polygon( pts=pts[iS,ib][ibprime] )::Polygon
                for ibprime in 1:cfg.dw.nposs[iS,ib] ]::Array{Polygon,1}
                for iS in 1:cfg.nS, ib in 1:cfg.nb ]
  end

  return ret

end
