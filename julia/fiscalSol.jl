#= fiscalSol.jl
Philip Barrett, pobarrett@gmail.com
12jul2016, Chicago

Solves the fiscalGame model  =#

"""
    function bpi( bgrid::linspace, bprime )
Returns a tuple (bl, bu, pl) where bl and bu are adjacent nodes on bgrid with bl < bprime < bu and bprime = pl * bl + (1-pl) * bu.
"""
function bpi( bgrid::LinSpace, bprime::Float64 )

  bmin = minimum(bgrid)
  bmax = maximum(bgrid)
  n = length(bgrid)

  if ( bprime < bmin )
    println( "bprime = ", bprime )
    println( "bmin = ", bmin )
    error("bprime must be greater than bmin")
  end

  if ( bprime > bmax )
    println( "bprime = ", bprime )
    println( "bmax = ", bmax )
    error("bprime must be less than bmax")
  end

  if bprime == bmin
    return 1, 2, bgrid[1], bgrid[2], 1
  elseif bprime== bmax
    return n-1, n, bgrid[n-1], bgrid[n], 0
  end

  nl = maximum( find( bgrid .< bprime ) )
  nu = nl + 1
  bl = bgrid[nl]
  bu = bgrid[nu]
      # Bounding indices and values
  pl = ( bu - bprime ) / ( bu - bl )

  return nl, nu, bl, bu, pl
end

"""
    uncSetUpdate( surp::Float64, pdloss::Float64, b::Float64, r::Float64,
                    bgrid::LinSpace,  W::Array{Polygon,2}, Q::Vector{Float64} )
Updates a particular polygon without any incentive-compatibility constraints for
a given surplus and period loss function.  Here Q is a vector of transition
probabilities for the *given* exogenous state
"""
function uncSetUpdate( surp::Float64, pdLoss::Vector{Float64}, b::Float64,
                        r::Float64, bgrid::LinSpace, W::Array{Polygon,2}, Q::Vector{Float64}, betta_hat::Vector{Float64},
                        dirs::Matrix{Float64}, outer::Bool=true )
## TODO: Add handling of empty sets

  bprime = (1+r) * b - surp
      # Continuation debt

  nl, nu, bl, bu, pl = bpi( bgrid, bprime )
      # The bounding debt choices and randomiation probability

# println("  nl=", nl)
# println("  nu=", nu)

  idx = (Q .> 0)
      # Because weighted sum does not play well with zero polygons
  wContl = weightedSum( W[idx,nl], Q[idx], dirs, outer )
  wContu = weightedSum( W[idx,nu], Q[idx], dirs, outer )
      # The bounding debt continuation value sets
  if ( pl == 0 )
    wCont = wContu
  elseif ( pl == 1 )
    wCont = wContl
  else
    wCont = weightedSum( [ wContl, wContu ], [pl, 1-pl], dirs, outer )
  end
      # The continuation value set: Need exceptions when pl is 0 or 1

  return ( ( 1 - betta_hat ) .* pdLoss ) + betta_hat .* wCont
end

"""
    function uncSetUpdate( fg::FiscalGame, W::Array{Polygon,2}, outer::Bool=true )
Updates the set using the unconstrained mapping for a full fiscal game
"""
function uncSetUpdate( fg::FiscalGame, W::Array{Polygon,2}, outer::Bool=true )
  nS, nb, nR, surp, pdLoss, bgrid, r, P, dirs, potFeas =
        fg.nS, fg.nb, fg.nR, fg.surp, fg.pdLoss, fg.bgrid, fg.r, fg.P, fg.dirs, fg.potFeas
      # Unpack fg
  betta_hat = fg.betta * fg.delta
      # Discuont rate
  out = copy(W)
      # initialize output
  for i in 1:nS, j in 1:nb
    println( "(i,j)=(", i, ",", j, ")" )
    out[i, j] = union( [ uncSetUpdate( surp[k,i], vec(pdLoss[k,:,i]), bgrid[j],
                                      r, bgrid, W, vec(P[i,:]), betta_hat,
                                      dirs, outer )::Polygon
                          for k in find(potFeas[i,j]) ] )
  end
  return out
end
