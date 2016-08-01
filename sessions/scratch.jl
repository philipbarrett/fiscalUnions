using BenchmarkTools

nS = 3
nb = 50
A = [ 3.0 3.0
      3.0 3.0
      3.0 3.0 ]
g =  [ .05  .0
       .0   .0
       .0   .05 ]
    # Need a low-g state so that blim can be escaped somehow
psi = [ .75, .75 ]
chi = [ 1.75, 1.75 ]
rho = .5
r = .03
bmin = 0.0
P = [ .6 .2 .2
      .2 .6 .2
      .2 .2 .6 ]

dw2 = dwlc2( nS=nS, nb=nb, A=A, g=g, psi=psi, chi=chi,
              r=r, bmin=bmin, rho=rho )

ndir = 24
dirs = hcat( [ cos(i*2*pi/ndir )::Float64 for i in 1:ndir ] ,
             [ sin(i*2*pi/ndir )::Float64 for i in 1:ndir ] )
idir = 1
dir = vec(dirs[ idir, : ])
    # Search direction
iS, ib, ibprimeidx = 1,1,3
    # NB: ibprimeidx is on the grid 1:(nposs[iS,ib])
ibprime = find( dw2.bprimeposs[iS,ib] )[ibprimeidx]
    # Index of bprime in grid, so bprime = bgrid[ibprime]
nR = 100
vho = [ rho, 1-rho ]

# RR = linspace( dw2.Rlow[iS,ib,1][ibprimeidx], dw2.Rhigh[iS,1], nR )
# WW = [ w_eval( RR[i], chi[1], psi[1], dw2.xlow[iS,ib,1][ibprime],
#                   dw2.xhigh[iS,1], A[iS,1], rho ) for i in 1:nR ]
# UU = [ w_eval( (1+r)*dw2.bgrid[ib] + sum(g[iS,:]) -
#                   dw2.bgrid[ibprime] - RR[i],
#                   chi[2], psi[2], dw2.xlow[iS,ib,2][ibprime],
#                   dw2.xhigh[iS,2], A[iS,2], rho ) for i in 1:nR ]
# plot( layer( x=RR, y=WW, Geom.line ),
#       layer( x=RR, y=UU, Geom.line ) )

cfg = ctsfiscalgame( r=0.04, delta=[.95, .95], psi=psi,
                          chi=chi, rho=.5, A=A, g=g, P=P,
                          nb=nb, bmin=0.0, ndirs=ndir,
                          par=false )

ibprimeidx = [ find( cfg.dw.bprimeposs[iS,ib] )
                for iS in 1:cfg.nS, ib in 1:cfg.nb ]

R1, R2, W1, W2, dist =
      dirMax( cfg.bgrid[ibprimeidx[iS,ib][ibprime]],
        cfg.bgrid[ib], chi, psi,
        [ cfg.dw.xlow[iS,ib,i][ibprime] for i in 1:2 ],
        vec( cfg.dw.xhigh[iS,:]),
        [ cfg.dw.Rlow[iS,ib,1][ibprime], cfg.dw.Rhigh[iS,1] ],
        vec(cfg.A[iS,:]), cfg.gSum[iS], cfg.rho, cfg.r,
        vec(dirs[idir,:]), false )

bm = @benchmark R1, R2, W1, W2, dist =
  dirMax( cfg.bgrid[ibprimeidx[iS,ib][ibprime]],
    cfg.bgrid[ib], chi, psi,
    [ cfg.dw.xlow[iS,ib,i][ibprime] for i in 1:2 ],
    vec( cfg.dw.xhigh[iS,:]),
    [ cfg.dw.Rlow[iS,ib,1][ibprime], cfg.dw.Rhigh[iS,1] ],
    vec(cfg.A[iS,:]), cfg.gSum[iS], cfg.rho, cfg.r,
    vec(dirs[idir,:]), false )
println(bm)

W = initGame( cfg )

pdout = pdPayoffs( cfg, dirs, true )
pdin = pdPayoffs( cfg, dirs, false )

updateout = valsUpdate( W, pdout, cfg.P, dirs, cfg.dw, cfg.betta )
updatein = valsUpdate( W, pdin, cfg.P, dirs, cfg.dw, cfg.betta, false )

ib = 19
polyPlot([W[1,ib], pdout[1,ib]..., pdout[2,ib]..., updateout[1,ib],
          updateout[2,ib], updatein[1,ib], updatein[2,ib]])

hh = hausdorff( W, updateout ) / ( 1 - cfg.betta )
fb = eqm( cfg, 100, true, false, pdout, W )
fbin = valsUpdate(fb, pdin, cfg.P, dirs, cfg.dw, cfg.betta, false )

dev = devCont( cfg.dw, cfg.P, cfg.betta, cfg.rho, cfg.A )

fbIC1 = valsUpdate( fb, pdout, P, dirs, cfg.dw, cfg.betta, false, dev )
                      # pdPayoffs::Array{Array{Polygon,1},2},
                      # P::Matrix{Float64}, dirs::Matrix{Float64},
                      # dw::DWLC2,
                      # betta::Float64, outer::Bool=true, vdev=NaN )
fbIC = eqm( cfg, 20, true, true, pdout, fb )

# using AudioIO
# play([sin(x) for x=0:0.03*pi:441])
