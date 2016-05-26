#= ploygon.jl
Philip Barrett, pobarrett@gmail.com
26may2016, Chicago

Defines polygon type and various associated methods
=#
"""
    type polygon
Defines the polygon type with three entries
"""
type polygon
  pts::Matrix     # The points representation
  dirs::Matrix    # The directions representation
  dists::Vector   # The distances associated with the vectors
end

"""
    polygon( ; pts::Matrix=Matrix[], dirs::Matrix=Matrix[], dists=[] )
Constructor for poylgon.  Maitained assumptions are that dirs is ordered
clockwise and that pts is already a convex hull.
"""
function polygon( ; pts::Matrix=Matrix[], dirs::Matrix=Matrix[], dists=[] )

  if( isempty( pts ) )
    pts = dirsToPts( dirs, dists )
  end

  if( isempty(dirs) )
    dirs, dists = ptsToDirs( pts )
  end

end

"""
    dirsToPts( dirs::Matrix, dists::Vector )
Given a normal/distance description of a set, returns a set of points mZ which
lie at the vertices.  This assumes that the vector of normals is ordered
clockwise already.
"""
function dirsToPts( dirs::Matrix, dists::Vector )
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
ptsToDirs( pts::Matrix )

  nPts = size(pts)[1]
      # Number of points
end
