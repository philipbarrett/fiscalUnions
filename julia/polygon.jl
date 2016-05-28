#= ploygon.jl
Philip Barrett, pobarrett@gmail.com
26may2016, Chicago

Defines polygon type and various associated methods
=#

import Base: *

"""
    type polygon
Defines the polygon type with three entries
"""
type polygon
  pts::Matrix{Float64}     # The points representation
  dirs::Matrix{Float64}    # The directions representation
  dists::Vector{Float64}   # The distances associated with the vectors
end

"""
    polygon( ; pts=[ NaN NaN ], dirs=[ NaN NaN ], dists=[NaN]  )
Constructor for poylgon.  Maitained assumptions are that dirs is ordered
clockwise and that pts is already a convex hull.
"""
function polygon( ; pts=[ NaN NaN ], dirs=[ NaN NaN ], dists=[NaN]  )

  if( !isnan( pts[1] ) && isnan( dirs[1] ) )
    dirs, dists = ptsToDirs( pts )
  elseif( !isnan( dirs[1] ) && isnan( pts[1] ) )
    pts = dirsToPts( dirs, dists )
  end
  return polygon( pts, dirs, dists )
end

"""
    dirsToPts( dirs::Matrix, dists::Vector )
Given a normal/distance description of a set, returns a set of points mZ which
lie at the vertices.  This assumes that the vector of normals is ordered
clockwise already.
"""
function dirsToPts( dirs::Matrix{Float64}, dists::Vector{Float64} )
  nPts = length(dists)
      # Number of points to be computed
  pts = zeros(Float64, 2, nPts)
      # Initialize the (transpose) output vector
  counter = 1
      # Counter: needed in cases where multiple boundary lines intersect at only
      # one point

  # Create the extended matrices (put the first case at the end)
  dirsExt = [ dirs ; dirs[1,:] ]
  distsExt = [ dists; dists[1,:] ]

  # Main loop: Over all points
  for i in 1:nPts
    A = dirsExt[i:(i+1),:]
        # The directional vectors to be intersected
    b = distsExt[i:(i+1)]
        # And their distances
    sol = A \ b
        # The candidate solution

    # Now compute the error: Detects multiple intersections
    err = norm( A * sol - b )
        # Inversion error
    relerr = ( norm(b) == 0 ) ? err : err / norm(b)
        # Relative error (where possible)

    # If there is no error, then assign to the output
    if( relerr < 1e-14 )
      pts[:,counter] = sol
      counter += 1
    end
  end
  return( transpose(pts) )
end

"""
    ptsToDirs( pts::Matrix )
Given a set of points, returns the normal-distance representation of the edges
Assumes that the points are ordered clockwise
clockwise already.
"""
function ptsToDirs( pts::Matrix{Float64} )

  nPts = size(pts)[1]
      # Number of points
  flip = [ 0 1; -1 0 ]
      # Matrix to convert grdient to normal
  neighbors = [ pts[ 2:nPts, : ] ; pts[1,:] ]
      # The neighboring points
  unscaled = ( neighbors - pts ) * flip
      # The unscaled vaules
  dirs::Matrix{Float64} = unscaled ./ sqrt( sum( neighbors .* neighbors, 2 ) )
      # The directions
  dists::Vector{Float64} = vec(sum( pts .* dirs, 2))
      # The distances in each direction
  return dirs, dists
end

"""
    add( poly1::polygon, poly2::polygon, dirs, outer=true )
Addition functionality.  Provides either an inner or out approximation
"""
function add( poly1::polygon, poly2::polygon, dirs, outer=true )
  dists1 = maximum( poly1.pts * dirs', 1 )
  dists2 = maximum( poly2.pts * dirs', 1 )
      # The distances in each direction
  dists = vec( dists1 + dists2 )

  if( outer )
    return polygon( dirs=dirs, dists=dists )
  end
  return polygon( pts=dirsToPts( dirs, dists ) )
end

"""
    add( poly1::Array{polygon,1}, poly2::Array{polygon,1}, dirs, outer=true )
Adds an array of polygons
"""
function add( polys::Array{polygon,1}, dirs, outer=true )

  N = length(polys)
      # Number of polygons to add
  dists = zeros( size(dirs)[1] )
      # Initate the distances
  for( i in 1:N )
    dists += maximum( polys[i].pts * dirs', 1 )'
        # Sum the distances
  end

  if( outer )
    return polygon( dirs=dirs, dists=vec( dists ) )
  end
  return polygon( pts=dirsToPts( dirs, vec( dists ) ) )
end

"""
    times( poly1::Array{polygon,1}, poly2::Array{polygon,1}, dirs, outer=true )
Adds an array of polygons
"""
function (*)( k::Number, poly::polygon )
  return polygon( k * poly.pts, poly.dirs, k * poly.dists )
end

"""
    wtdSum( polys::Array{polygon,1}, wts::Vector dirs, outer=true )
Computes a weighted sum of polygons
"""
function wtdSum( polys::Array{polygon,1}, wts::Vector, dirs, outer=true )
  N = length( polys )
      # Number of polygons
  polys2 = polys
      # Initiate scaled polygons with input array
  for( i in 1:N )
    polys2[i] = wts[i] * polys[i]
  end
  return( add( polys2, dirs, outer ) )
end

# """
#     whichmin( x::Vector )
# Uncovers the location of the minimum of a vector
# """
# function whichmin( x::Vector )
#   i = 1
#   min_x=minimum(x)
#   while( x[i] > min_x )
#     i+=1
#   end
#   return i
# end

"""
    acw( p1, p2, p3 )
Positive if p1 -> p2 -> p3 is anti-clockwise
"""
function acw( p1, p2, p3 )
  return (p2[1] - p1[1]) * (p3[2] - p1[2]) - (p2[2] - p1[2]) * (p3[1] - p1[1])
end

"""
    gScan( pts::Matrix )
Computes a convex hull using the Graham Scan algortithm
"""

function gScan( pts::Matrix )

  N = size(pts)[1]
      # Number of points
  lex = sortrows( pts[ :, [2, 1 ] ] )[ :, [ 2, 1 ] ]
      # Sorts the points lexicographically by y then x
  P = lex[ 1, : ]
      # The initial point
  otherpts = lex[ 2:end, : ]
      # The other points
  orderpts = [ P; otherpts ; P ]
      # Initialize the anti-clockwise ordered points
  cah = zeros( N - 1 )
      # The cosine of the angle between each point and P.  Will be in
      # [-1,1] because of choice of P

  ## Order the points anti-clockwise ##
  for( i in 1:(N-1) )
    diff = otherpts[i,:] - P
    cah[i] = diff[1] / norm(diff)
        # Cosine is adjacent over hypotenuse
  end
  orderpts[2:N,:] = otherpts[ sortperm(cah, rev=true), : ]
      # Order the points by the angle measure

  ## Create the output ##
  out = similar( pts )
  out[ 1:2, : ] = orderpts[1:2,:]
      # Initialize
  M = 2
      # Counts number of rows in convex hull
  for( i in 3:N )
    while( acw( out[ M-1, : ], out[ M, : ], orderpts[ i, : ] ) <= 0 )
      if( M > 2 )
        M -= 1
      elseif( i == N )
        break
      else
        i += 1
      end
    end
    M += 1
        # Increment counter
    out[ M, : ] = orderpts[ i, : ]
        # Add to the output if we have an acw angle
  end
  return out[ 1:M, : ]
end
