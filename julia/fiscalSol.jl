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
    function dirscreate( N::Int )
Creates the N search directions
"""
function dirscreate( N::Int )
  hcat( [ cos(2*pi*(i-1)/N) for i in 1:N ],
               [ sin(2*pi*(i-1)/N) for i in 1:N ] )
end

"""
    function uncSetUpdate( fg::FiscalGame, W::Array{Polygon,2}, outer::Bool=true )
Updates the set using the unconstrained mapping for a full fiscal game
"""
function uncSetUpdate( fg::FiscalGame, W::Array{Polygon,2},
                        ndirs::Matrix{Int}, outer::Bool=true )
  nS, nb, nR, surp, pdLoss, bgrid, r, P, potFeas, par =
        fg.nS, fg.nb, fg.nR, fg.surp, fg.pdLoss, fg.bgrid, fg.r, fg.P, fg.potFeas, fg.par
      # Unpack fg
  betta_hat = fg.betta * fg.delta
      # Discount rate
  out = copy(W)
      # initialize output

  if par
      # Parallel execution
    # @spawn r, bgrid, W, betta_hat, outer
    #     # Export variables - I wish that this worked! :(
    surp_ar = vec( [surp[:,i]::Vector{Float64}
                        for i in 1:nS, j in 1:nb] )
    pdLoss_ar = vec( [ pdLoss[:,:,i]::Matrix{Float64}
                        for i in 1:nS, j in 1:nb] )
    b_ar = vec( [ bgrid[j]::Float64
                        for i in 1:nS, j in 1:nb] )
    r_ar = vec( [ r::Float64 for i in 1:(nS*nb) ] )
    bgrid_ar = vec( [ bgrid::LinSpace{Float64} for i in 1:(nS*nb) ] )
    W_ar = vec( [ W::Array{Polygon,2} for i in 1:(nS*nb) ] )
    Q_ar = vec( [ vec(P[i,:])::Vector{Float64}
                        for i in 1:nS, j in 1:nb ] )
    bh_ar = vec( [ betta_hat::Vector{Float64} for i in 1:(nS*nb) ] )
    dirs_ar = vec( [ dirscreate(ndirs[i,j])::Matrix{Float64}
                        for i in 1:nS, j in 1:nb] )
    idx_ar = vec( [ find(potFeas[i,j])::Vector{Int}
                        for i in 1:nS, j in 1:nb] )
    outer_ar = vec( [ outer::Bool for i in 1:(nS*nb) ] )
        # Setting up the array inputs
    temp = pmap( (s,p,b,r,bg,WW,Q,bh,d,idx,out) ->
              union( [ uncSetUpdate( s[k], vec(p[k,:]), b, r, bg,
                                      WW, Q, bh, d, out )::Polygon
                         for k in idx ] )::Polygon,
          surp_ar, pdLoss_ar, b_ar, r_ar, bgrid_ar, W_ar, Q_ar,
          bh_ar,dirs_ar, idx_ar, outer_ar )
#    return temp
#    out = [ temp[i + j * (nS-1)]::Polygon for i in 1:nS, j in 1:nb ]
    out = reshape( temp, nS, nb )
    out = [ out[i,j]::Polygon for i in 1:nS, j in 1:nb ]
  else
      # Serial execution
    for i in 1:nS, j in 1:nb
      println( "(i,j)=(", i, ",", j, ")" )
      dirs = dirscreate( ndirs[i,j] )
          # The search directions
      out[i, j] = union( [ uncSetUpdate( surp[k,i], vec(pdLoss[k,:,i]), bgrid[j],
                                        r, bgrid, W, vec(P[i,:]), betta_hat,
                                        dirs, outer )::Polygon
                            for k in find(potFeas[i,j]) ] )
    end # for i in 1:nS, j in 1:nb
  end # if par

  return out
end

function ndirsUpdate( ndirs, hd, hddirs::Float64 = 1e-6, ndirsu=2^8 )
  N, M = size(ndirs)
      # Dimensions
  return [ hd[i,j] < hddirs ? min( 2 * ndirs[i,j], ndirsu ) : ndirs[i,j] for i in 1:N, j in 1:M ]
end

function uncSol( fg::FiscalGame, W::Array{Polygon,2},
                        ndirs::Matrix{Int}, outer::Bool=true,
                        maxiter::Int=200, tol=1e-05,
                        saveloc = "" )
  hd = [tol*2]
  it = 0
    # Initialize loop variables
  while ( it < maxiter && maximum(hd) > tol )
    it += 1
    println("*** it = ", it , " ***")
    W_new = uncSetUpdate( fg, W, ndirs )
    hd = hausdorff( W, W_new )
    println("   ** mean & max hausdorff dist = ", ( mean(hd), maximum(hd) ), "**" )
    W = copy(W_new)
    ndirs = ndirsUpdate( ndirs, hd, fg.hddirs, fg.ndirsu )
    println("   ** mean # search dirs = ", mean(ndirs), "**" )
    println("   ** extremal # search dirs = ", extrema(ndirs), "**\n" )
  end

  if saveloc != ""
    jldopen(saveloc, "w") do file
      addrequire(file, Polygons)
      write(file, "W", W)
      write(file, "hd", hd)
      write(file, "ndirs", ndirs )
    end
  end

  return W, ndirs, hd
end
