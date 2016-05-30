#= polygon.jl
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
    add( poly::polygon, pt::Vector )
Adds a vector to all points in a polygon
"""
function add( poly::polygon, u::Vector )

  pts = similar( poly.pts )
  pts[:,1] = poly.pts[:,1] + u[1]
  pts[:,2] = poly.pts[:,2] + u[2]
      # Shift the points
  return polygon( pts = pts )
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
Scalar multiplication
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

"""
    crop( poly::polygon, dim::Int, dist, upper=true )
Returns a polygon cropped in either the x or y dimension, with upper = true
retaining the part of the polygon above the chop, and if false the part below.
The polygon should already be oreinted counterclockwise
"""
## TODO: DEBUG ME!
## Strategy here is all wrong.  Better: Identify all the points that are to the
## left (right) of the slice.  Then throw away them (and their associated
## normals).  Finally, reintroduce the slice boundary into the dirs/dist
## framework
function crop( poly::polygon, dim::Int, dist, upper=true )
  N = size(poly.pts)[1]

  # Compute the direction of the chop
  if( dim == 1 && upper )
    dir = [ -1, 0 ]
    dist = - dist
  elseif( dim == 1 && !upper )
    dir = [ 1, 0 ]
  elseif( dim==2 && upper )
    dir = [ 0, -1 ]
    dist = - dist
  else
    dir = [ 0, 1 ]
  end
  c = dir[1]
  s = dir[2]
      # Sign and cosine of the dir vector

  # The chop doesn't bind
  l = minimum( poly.pts[:,dim] )
  u = maximum( poly.pts[:,dim] )
  if( sum( dist * dir ) < l || sum( dist * dir ) > u )
    return(poly)
  end
  if( sum( dist * dir ) < u || sum( dist * dir ) > l )
    return(polygon())
  end

  # Find the points which are to be chopped
  # init = 1, term = 1
  #     # Initial and terminal chopping points
  # for( i in 1:N )
  #   if(  )
  # end

# 
#
#
#
#
#
#
#
#
#
#
#
#   # Now loop over the direction vectors
#   i = 1
#       # Counter
#   c1 = poly.dirs[1,1] / norm( poly.dirs[1,:] )
#   s1 = poly.dirs[1,2] / norm( poly.dirs[1,:] )
#       # Initiate the next point's sign and cos
#   while( i<(N-1) )
#     c0 = c1
#     s0 = c1
#         # Old is now new. How sad, how true.
#     c1 = poly.dirs[i+1,1] / norm( poly.dirs[i+1,:] )
#     s1 = poly.dirs[i+1,2] / norm( poly.dirs[i+1,:] )
#         # Cosine and sine of the next direction vector
#     if( min( s1, s0 ) > 0  )
#       if( c0 >= c > c1 )
#         break
#       end
#     elseif( max( s1, s0 ) < 0 )
#       if( c0 <= c < c1 )
#         break
#       end
#     elseif( s0 > 0 > s1 )
#       if( c <= min( c0, c1 ) )
#         break
#       end
#     else
#       if( c >= max( c0, c1 ))
#         break
#       end
#     end
#     i += 1
#   end
#       # So dir fits in between the ith and (i+1)th direction vectors
#   newdirs = [ poly.dirs[1:i, :] ; dir' ; poly.dirs[ (i+1):end, : ] ]
#   newdists = [ poly.dists[1:i] ; dist ; poly.dists[ (i+1):end ] ]
#       # Create the new directions and distances
#
#   ## Special cases for exact dir in poly.dirs ##
#   if( poly.dirs[i,:] == dir)
#     newdirs = poly.dirs
#     newdists = [ poly.dists[1:(i-1)] ; dist ; poly.dists[ (i+1):end ] ]
#   end
#   if( poly.dirs[i+1,:] == dir)
#     newdirs = poly.dirs
#     newdists = [ poly.dists[1:i] ; dist ; poly.dists[ (i+2):end ] ]
#   end
#
# println( "i=", i )
# println( "newdirs\n", newdirs )
# println( "newdists\n", newdists )
#
#   return( polygon( dirs = newdirs, dists = newdists ) )

end


























# Comment
