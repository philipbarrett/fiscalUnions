#=== Things to keep ===#

# Include function definitions
include("../julia/autarky.jl")
include("../julia/autarkySim.jl")


### 1. Solving the model ###
am = AutarkyModel()
as = solve_am(am)
asim = sim_am(as)

prs_m = prsModel( )
prs_s = solve_am( prs_m )
prs_sim = sim_am( prs_s )

plot( asim[1:200,2] )
plot( asim[1:200,4] )
scatter( asim[1:200,2 ], asim[1:200,4 ] )
scatter( asim[1:200,2 ], asim[1:200,4 ] )

cor( asim[1:(end-1),2], asim[2:(end),2])
cor( asim[1:(end-1),4], asim[2:(end),4])

mean( asim, 1 )
mean( prs_sim, 1 ) / 2

### 2. Polygon stuff ###
dirs  = Matrix{Float64}( [ -1 0 ; 0 1 ; 1 0 ; 0 -1 ] )
dists = Vector{Float64}( [ 1, 1, 1, 1 ] )
z = dirsToPts( dirs, dists )

dirs  = Matrix{Float64}( [ -1 0 ; 0 1 ; 1 1 ; 1 0 ; 0 -1 ] )
dists = Vector{Float64}( [ 1, 1, 1, 1, 1 ] )
z = dirsToPts( dirs, dists )

dirs_ret, dist_ret = ptsToDirs( z )
z_ret = dirsToPts( dirs_ret, dist_ret )

s = polygon( z, dirs, dists )
q = polygonP( pts = z )
r = polygonD( dirs = dirs, dists = dists )

q.pts == s.pts == r.pts

dirs2  = Matrix{Float64}( [ -1 0 ; 0 1 ; 1 0 ; 0 -1 ] )
dists 2= Vector{Float64}( [ 1, 1, 1, 1 ] )
z = dirsToPts( dirs2, dists2 )

t = polygonP( pts = z )
e = add( s, s, dirs )
e = add( s, s, dirs )

add( t, s )
