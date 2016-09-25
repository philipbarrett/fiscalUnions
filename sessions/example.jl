using Polygons
using JLD

include("../julia/dwlCts.jl")
include("../julia/ctsFiscalGame.jl")
include("../julia/wApx.jl")
include("../julia/ctsSch.jl")
include("../julia/ctsSol.jl")

nS = 3
nb = 30
A = [ 3.09 2.91
      3.0  3.0
      2.91  3.09 ]
g =  [ .1  .2
       .05  .05
       .2  .1 ]
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
ndir = 24

cfg = ctsfiscalgame( r=r, delta=[.95, .95], psi=psi,
                          chi=chi, rho=.5, A=A, g=g, P=P,
                          nb=nb, bmin=0.0, ndirs=ndir,
                          par=false )


# Try one iteration
dirs = hcat( [ cos(i*2*pi/ndir )::Float64 for i in 1:ndir ] ,
             [ sin(i*2*pi/ndir )::Float64 for i in 1:ndir ] )
W = initGame( cfg )
    # Initiate payoffs
pdout = pdPayoffs( cfg, dirs, true )
pdin = pdPayoffs( cfg, dirs, false )
    # Period payoffs.  Outer approx now good! :)
updateout = valsUpdate( W, pdout, cfg.P, dirs, cfg.dw, cfg.gSum,
                          cfg.rho, cfg.r, cfg.betta )
updatein = valsUpdate( W, pdin, cfg.P, dirs, cfg.dw, cfg.gSum,
                          cfg.rho, cfg.r, cfg.betta, false )

# polyPlot()  # Some stuff here

eqout = eqm( cfg, 200, true, false, pdout, W )
eqin = eqm( cfg, 20, false, false, pdin, eqout )
    # Woooooot solutions!

### Now: Apply the incentive compatibility condition
# eqoutIC = eqm( cfg, 5, true, true, pdout, eqout )
# eqinIC = eqm( cfg, 5, false, true, pdin, eqin )
eqoutIC1 = eqm( cfg, 1, true, true, pdout, eqout )
eqhybIC1 = eqm( cfg, 1, true, true, pdin, eqin )
  # Hybrid case: Outer apx using inner ctn
eqhybIC1X = eqm( cfg, 1, false, true, pdin, eqout )
  # Hybrid case 2: Inner apx using outer ctn
eqinIC1 = eqm( cfg, 1, false, true, pdin, eqin )

save("/home/philip/Dropbox/data/2016/fiscalUnions/eg.jld", "cfg", cfg,
      "eqout", eqout, "eqin", eqin, "eqoutIC", eqoutIC,
      "eqinIC", eqinIC, "pdin", pdin, "pdout", pdout )

# jldopen("/home/philip/Dropbox/data/2016/fiscalUnions/unc.jld", "w") do file
#     addrequire(file, Polygons)
#     write(file, "W", W)
# end
#
# # Can then acecss W with:
# d = load("/home/philip/Dropbox/data/2016/fiscalUnions/unc.jld")
# W = d["W"]
