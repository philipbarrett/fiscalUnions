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
dirs  = [ -1 0 ; 0 1 ; 1 0 ; 0 -1 ]
dists = [ 1 ; 1 ; 1 ; 1 ]
z = dirsToPts( dirs, dists )

dirs  = [ -1 0 ; 0 1 ; 1 1 ; 1 0 ; 0 -1 ]
dists = [ 1 ; 1 ; 1 ; 1 ; 1 ]
z = dirsToPts( dirs, dists )
