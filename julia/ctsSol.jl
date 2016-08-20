#= ctsSol.jl
Philip Barrett, pobarrett@gmail.com
28jul2016, Chicago

Solves the continuous action problem
=#

"""
    expectedVals(  conts::Array{Polygon, 2}, P::Matrix{Float64},
                    dirs::Matrix{Float64}, outer::Bool=true  )
Computes the set of expected values for each continuation debt level
"""
function expectedVals( conts::Array{Polygon, 2}, P::Matrix{Float64},
                        dirs::Matrix{Float64}, outer::Bool=true )
  nS, nb = size(conts)
  return [ weightedSum( conts[:,ibprime], vec(P[iS,:]), dirs, outer )
              for iS in 1:nS, ibprime in 1:nb ]
end

"""
    expectedValsExtrema(  conts::Array{Polygon, 2}, P::Matrix{Float64},
                            dirs::Matrix{Float64}, outer::Bool=true )
Computes the extremal distances/points in each of the search directions of
the expected continuation values.
"""
function expectedValsExtrema( conts::Array{Polygon, 2}, P::Matrix{Float64},
                              dirs::Matrix{Float64}, outer::Bool=true )
  nS, nb = size(conts)
  ndirs = size(dirs)[1]
  td = transpose(dirs)
      # Do this just once
  if outer
    dists = [ maximum( conts[i,j].pts * td, 1 ) for i in 1:nS, j in 1:nb ]
    return [ [ ( P[iS,:] * [ dists[kS,j][l] for kS in 1:nS ] )[1]::Float64
                for l in 1:ndirs ] for i in 1:nS, j in 1:nb ]
  else
    pts =  [ conts[i,j].pts[[indmax( conts[i,j].pts * td[:,k] )
              for k in 1:ndirs ],:]::Array{Float64,2} for i in 1:nS, j in 1:nb ]
    return [ [ ( P[iS,:] * [ pts[kS,j][l,r]::Float64 for kS in 1:nS ] )[1]::Float64
                  for l in 1:ndirs, r in 1:2 ]::Matrix{Float64}
                    for i in 1:nS, j in 1:nb ]::Array{Matrix{Float64},2}
  end
end


function valsUpdate( conts::Array{Polygon, 2},
                      pdPayoffs::Array{Array{Polygon,1},2},
                      P::Matrix{Float64}, dirs::Matrix{Float64},
                      dw::DWLC2, gSum::Vector, rho::Float64, r::Float64,
                      betta::Float64, outer::Bool=true, evbar=false )
  nS, nb = size(conts)
  ibprimeidx = [ find( dw.bprimeposs[iS,ib] )
                  for iS in 1:nS, ib in 1:nb ]
      # The indices of the continuation debt levels
  econt = expectedVals( conts, P, dirs, outer )
      # The continuation values
  if evbar == false
    actionsets = [ [ setSum( (1-betta) * pdPayoffs[iS,ib][ibprime],
                            betta * econt[iS,ibprimeidx[iS,ib][ibprime]],
                            dirs, outer )::Polygon
                          # Add period and continuation values
                    for ibprime in 1:dw.nposs[iS,ib] ]::Array{Polygon,1}
                    for iS in 1:nS, ib in 1:nb ]
        # Add the period payoffs to the appropriate continuation sets
    out = [ union(actionsets[iS,ib])::Polygon for iS in 1:nS, ib in 1:nb ]
  else
    outputtext = outer ? "dists" : "pts"
        # Text selector for output
    Tp = outer? Vector{Float64} : Matrix{Float64}
    actionsets = [ [ search_ic( dw.bgrid[ibprimeidx[iS,ib][ibprime]],
                          dw.bgrid[ib], [ dw.Rlow[iS,ib,k][ibprime] for k in 1:2 ],
                          vec( dw.Rhigh[iS, : ] ), gSum[iS], rho, r, betta,
                          vec(dw.apx_coeffs[iS,:]), vec(dw.apx_N[iS,:]),
                          dirs, econt[iS,ib], vec(evbar[iS,:]), 0,
                          outputtext )::Tp
                          # Add period and continuation values
                    for ibprime in 1:dw.nposs[iS,ib] ]
                    for iS in 1:nS, ib in 1:nb ]
        # Incentive compatibility
    if outer
      ndir = size(dirs)[1]
      dists = [ [ maximum( [ actionsets[iS,ib][ibprime][i]
                                  for ibprime in 1:dw.nposs[iS,ib] ] )::Float64
                                    for i in 1:ndir ]
                    for iS in 1:nS, ib in 1:nb ]
      out = [ Polygon(dirs=dirs, dists=dists[iS,ib])::Polygon
                for iS in 1:nS, ib in 1:nb ]
          # Formulate union as maximum in each search direction
    else
      out = [ Polygon( pts = vcat(actionsets[iS,ib]...))::Polygon
                        for iS in 1:nS, ib in 1:nb ]
          # Take the union over the action-specific sets
    end
  end

  return out
end

function eqm( cfg::CtsFiscalGame, maxiter::Int=200, outer::Bool=true,
                      IC=false, pd=false, w_in=false, tol::Float64=1e-04 )

  dirs = hcat( [ cos(i*2*pi/cfg.ndirs )::Float64 for i in 1:cfg.ndirs ] ,
               [ sin(i*2*pi/cfg.ndirs )::Float64 for i in 1:cfg.ndirs ] )
      # Search directions

  if pd == false
    println("** Computing period payoffs **")
    pd = pdPayoffs( cfg, dirs, outer )
  end
  println("** Computing initial continuation values **")
  if w_in == false
    w_old = initGame(cfg)
  else
    println("      Using user-supplied continuations")
    w_old = w_in
  end
  if IC
    println("** Computing deviating values **")
    vdev = devCont( cfg.dw, cfg.P, cfg.betta, cfg.rho, cfg.A )
  end

  println("** Main iteration **")
  hdiff = tol * 2
  it = 1
  while (it <= maxiter) & (hdiff > tol)
    println("    Iteration ", it )
    if IC
      ##### THIS PART NOT WORKNG YET
      w_new = valsUpdate( w_old, pd, cfg.P, dirs, cfg.dw, cfg.betta, outer, vdev )
      ############################33
    else
      w_new = valsUpdate( w_old, pd, cfg.P, dirs, cfg.dw, cfg.gSum,
                                cfg.rho, cfg.r, cfg.betta )
    end
    hdiff = maximum( hausdorff( w_old, w_new ) )
    println("      Normalized Hausdorff distance = ", hdiff )
    it += 1
    w_old=w_new
  end

  return w_old

end


function devCont( dw::DWLC2, P::Matrix, betta::Float64, rho::Float64,
                  A::Matrix, vbar_flag::Bool=false )

  vrho = [rho, 1-rho ]
  nS, nb = dw.nS, dw.nb

  # First, compute the state-dependent payoff at the debt limit
  wlim = [ vrho[j] * w_eval( dw.Rhigh[iS,j], chi[j], psi[j],
            dw.xlow[iS,nb,j][dw.nposs[iS,nb]], dw.xhigh[iS,j],
            A[iS,j], vrho[j] ) for iS in 1:nS, j in 1:2 ]
      # Period payoff
  vbar = inv( eye(nS) - betta * P ) * ( ( 1 - betta ) * wlim )
      # Solves v* = (1-beta).ulim + beta.P.v*
  if vbar_flag
    return vbar
  end
  wdev = [ vrho[j] * w_eval( minimum(dw.Rlow[iS,ib,j]), chi[j], psi[j],
            minimum(dw.xlow[iS,ib,j]), dw.xhigh[iS,j], A[iS,j], vrho[j] )
            for iS in 1:nS, ib in 1:nb, j in 1:2 ]
      # Deviating period payoff
  vdev = (1-betta) * wdev + betta *
            [ vbar[iS,j] for iS in 1:nS, ib in 1:nb, j in 1:2 ]
        # The minmax payoff
  return vdev
end
