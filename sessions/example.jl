using Polygons
using JLD

include("../julia/dwlCts.jl")
include("../julia/ctsFiscalGame.jl")
include("../julia/wApx.jl")
include("../julia/ctsSch.jl")
include("../julia/ctsSol.jl")

nS = 3
nb = 50
A = [ 3.01 3.0
      3.0  3.0
      3.0  3.01 ]
g =  [ .1  .1
       .0  .0
       .1  .1 ]
    # Need a low-g state so that blim can be escaped somehow
psi = [ .75, .75 ]
chi = [ 8.0, 8.0 ]
rho = .5
r = .04
bmin = 0.0
P = [ .6 .2 .2
      .2 .6 .2
      .2 .2 .6 ]
betta = 1 / ( 1 + r )
ndir = 48

cfg = ctsfiscalgame( r=0.04, delta=[.95, .95], psi=psi,
                          chi=chi, rho=.5, A=A, g=g, P=P,
                          nb=nb, bmin=0.0, ndirs=ndir,
                          par=false )

eq = eqm( cfg, 150 )



# jldopen("/home/philip/Dropbox/data/2016/fiscalUnions/unc.jld", "w") do file
#     addrequire(file, Polygons)
#     write(file, "W", W)
# end
#
# # Can then acecss W with:
# d = load("/home/philip/Dropbox/data/2016/fiscalUnions/unc.jld")
# W = d["W"]
