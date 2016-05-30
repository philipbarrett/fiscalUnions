#=== Things to keep ===#

# Include function definitions
include("../julia/autarky.jl")
include("../julia/autarkySim.jl")
include("../julia/polygon.jl")
include("../julia/polygonUtils.jl")
include("../julia/polygonPlot.jl")

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
dirs  = Matrix{Float64}( [ -1 0 ; 0 -1 ; 1 0 ; 0 1 ] )
dists = Vector{Float64}( [ 1, 1, 1, 1 ] )
z = dirsToPts( dirs, dists )

dirs  = Matrix{Float64}( [ -1 0 ; 0 -1 ; 1 0 ; 1 1 ; 0 1 ] )
dists = Vector{Float64}( [ 1, 1, 1, 1, 1 ] )
z = dirsToPts( dirs, dists )

dirs_ret, dist_ret = ptsToDirs( z )
z_ret = dirsToPts( dirs_ret, dist_ret )

s = polygon( z, dirs, dists )
q = polygon( pts = z )
r = polygon( dirs = dirs, dists = dists )

q.pts == s.pts == r.pts

dirs2  = Matrix{Float64}( [ -1 0 ; 0 1 ; 1 0 ; 0 -1 ] )
dists2= Vector{Float64}( [ 1, 1, 1, 1 ] )
z = dirsToPts( dirs2, dists2 )

t = polygon( pts = z )
ee = add( s, s, dirs )
f = add( s, s, dirs, false )
g = add( s, s, dirs2 )
h = add( s, s, dirs2, false )
# e = add( s, s, dirs )
plot( [ s, add( s, [1,1] ) ] )

add( t, s, dirs2 )
apoly = [ t, ee, r, q ]
jj = add( apoly, dirs2 )

hh = 4 * s

ll = wtdSum( apoly, [1:4], dirs2 )

ff = [ 1 0 ; 0 1 ; 0 0 ; .25 .25 ]
vv = gScan(ff)

gg = [ 0 1; .5 .5 ; .9 .2 ; 1 1 ; 1 0; .3 .6 ; 0 0 ]
nn = gScan(gg)

hh = [ gg ; .5 1 ]
mm = gScan(hh)

jj = [ .6 2 ; hh ; .5 1 ]
oo = gScan(jj)

kk = [ 0 1; .5 .5 ; .3 .6 ; 0 0 ; .6 2 ;
      .3 -1 ; .6 0 ; .8 1 ; .7 1.2 ]
pp = gScan(kk)
plot( layer( x=kk[:,1], y=kk[:,2], Geom.point),
      layer( x=pp[:,1], y=pp[:,2], Geom.path,
          Theme(default_color=colorant"red") ) )


ll = randn( 100, 2 )
@time qq = gScan( ll )
plot( layer( x=ll[:,1], y=ll[:,2], Geom.point),
      layer( x=qq[:,1], y=qq[:,2], Geom.path,
          Theme(default_color=colorant"red") ) )


# ll = randn( 100, 2 )
# @time qq = gScan( ll )
