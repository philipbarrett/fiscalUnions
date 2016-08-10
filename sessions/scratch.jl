using BenchmarkTools

nS = 3
nb = 30
A = [ 3.0 3.0
      3.0 3.0
      3.0 3.0 ]
g =  [ .2  .1
       .0  .0
       .1  .2 ]
    # Need a low-g state so that blim can be escaped somehow
psi = [ .75, .75 ]
chi = [ 8.0, 8.0 ]
rho = .5
r = .03
bmin = 0.0
P = [ .6 .2 .2
      .2 .6 .2
      .2 .2 .6 ]
betta = 1 / ( 1 + r )

dw2 = dwlc2( nS=nS, nb=nb, A=A, g=g, psi=psi, chi=chi,
              r=r, bmin=bmin, rho=rho )

ndir = 24
dirs = hcat( [ cos(i*2*pi/ndir )::Float64 for i in 1:ndir ] ,
             [ sin(i*2*pi/ndir )::Float64 for i in 1:ndir ] )
idir = 1
dir = vec(dirs[ idir, : ])
    # Search direction
iS, ib, ibprimeidx1 = 1,1,3
    # NB: ibprimeidx is on the grid 1:(nposs[iS,ib])
ibprime = find( dw2.bprimeposs[iS,ib] )[ibprimeidx1]
    # Index of bprime in grid, so bprime = bgrid[ibprime]
nR = 100
vho = [ rho, 1-rho ]

RR = linspace( dw2.Rlow[iS,ib,1][ibprimeidx1], dw2.Rhigh[iS,1], nR )
WW = [ w_eval( RR[i], chi[1], psi[1], dw2.xlow[iS,ib,1][ibprime],
                  dw2.xhigh[iS,1], A[iS,1], rho ) for i in 1:nR ]
NN, coeff, err = w_eval_fit( dw2.Rhigh[iS,1], chi[1], psi[1],
                  dw2.xlow[iS,ib,1][ibprime],
                  dw2.xhigh[iS,1], A[iS,1], rho )
UU = [ w_eval( (1+r)*dw2.bgrid[ib] + sum(g[iS,:]) -
                  dw2.bgrid[ibprime] - RR[i],
                  chi[2], psi[2], dw2.xlow[iS,ib,2][ibprime],
                  dw2.xhigh[iS,2], A[iS,2], rho ) for i in 1:nR ]
VV = [ w_eval_apx( RR[i], coeff, dw2.Rhigh[iS,1], NN ) for i in 1:nR ]
plot( layer( x=RR, y=WW, Geom.line, Theme(default_color=colorant"purple") ),
      layer( x=RR, y=VV, Geom.line, ) )

############################################################

using JuMP
using Polygons
using Ipopt

iS, ib, ibprimeidx1 = 1,1,3
ibprime = find( dw2.bprimeposs[iS,ib] )[ibprimeidx1]
idir = 16
vevbar = [ 2.3, 2.3 ]

ibprimeidx = [ find( dw2.bprimeposs[iS,ib] ) for iS in 1:nS, ib in 1:nb ]

dists = [ (1 + i / (10*ndir))::Float64 for i in 1:ndir ]
V = Polygon( dirs=dirs, dists=dists)

pd_1( R ) = w_eval_apx( R, dw2.apx_coeffs[iS,1], dw2.Rhigh[iS,1],
                            dw2.apx_N[iS,1] )::Float64
pd_2( R ) = w_eval_apx( R, dw2.apx_coeffs[iS,2], dw2.Rhigh[iS,2],
                            dw2.apx_N[iS,2] )
pd_1_d( R ) = w_eval_apx_d( R, dw2.apx_coeffs[iS,1], dw2.Rhigh[iS,1],
                            dw2.apx_N[iS,1] )
pd_2_d( R ) = w_eval_apx_d( R, dw2.apx_coeffs[iS,2], dw2.Rhigh[iS,2],
                            dw2.apx_N[iS,2] )
pd_1_d2( R ) = w_eval_apx_d2( R, dw2.apx_coeffs[iS,1], dw2.Rhigh[iS,1],
                            dw2.apx_N[iS,1] )
pd_2_d2( R ) = w_eval_apx_d2( R, dw2.apx_coeffs[iS,2], dw2.Rhigh[iS,2],
                            dw2.apx_N[iS,2] )
    # Period objective functions
JuMP.register(:pd_1, 1, pd_1, pd_1_d, pd_1_d2)
JuMP.register(:pd_2, 1, pd_2, pd_2_d, pd_2_d2)
JuMP.register(:convex_basis, 1, convex_basis, convex_basis_d,
                                convex_basis_d2)


t = Model(solver=Ipopt.IpoptSolver(print_level=5))
    # Test model
@NLparameter( t, a1 == dirs[idir,1] )
@NLparameter( t, a2 == dirs[idir,2] )
@variable( t, dw2.Rlow[iS,ib,1][ibprime] <= R1 <= dw2.Rhigh[iS,1] - 1e-9 )
@variable( t, dw2.Rlow[iS,ib,2][ibprime] <= R2 <= dw2.Rhigh[iS,2] - 1e-9 )
@variable( t, v1 )
@variable( t, v2 )
# @NLobjective(t, Max,
#   (1-betta) * ( a1 * pd_1( R1 ) + a2 * pd_2( R2 ) ) +
#    betta * ( a1 * v1 + a2 * v2 ) )
@NLobjective(t, Max,
  (1-betta) * ( a1 * pd_1( R1 )
              + a2 * pd_2( R2 ) ) +
   betta * ( a1 * v1 + a2 * v2 ) )
@constraint(t, R1 + R2 + dw2.bgrid[ibprimeidx[iS,ib][ibprime]] ==
                  dw2.bgrid[ib]*(1+r) + sum(g[iS,:]) )
@constraint(t, dirs * [ v1, v2 ] .<= V.dists )
@NLconstraint( t, (1-betta)*(pd_1(R1)-pd_1(dw2.Rlow[iS,ib,1][ibprime])) <=
                    betta*(vbar[1]-v1))
@NLconstraint( t, (1-betta)*(pd_2(R2)-pd_2(dw2.Rlow[iS,ib,2][ibprime])) <=
                    betta*(vbar[2]-v2))

print(t)

status = solve(t)

println("Objective value: ", getobjectivevalue(t))
println("R1 = ", getvalue(R1))
println("R2 = ", getvalue(R2))
println("v1 = ", getvalue(v1))
println("v2 = ", getvalue(v2))

# setvalue( a1, dirs[1,1] )
# setvalue( a2, dirs[1,2] )
# status = solve(t)
#
# println("Objective value: ", getobjectivevalue(t))
# println("R1 = ", getvalue(R1))
# println("R2 = ", getvalue(R2))
# println("v1 = ", getvalue(v1))
# println("v2 = ", getvalue(v2))

m = Model(solver=Ipopt.IpoptSolver(print_level=5))
    # Test model
@variable( m, dw2.Rlow[iS,ib,1][ibprime] <= R1 <= dw2.Rhigh[iS,1] - 1e-9 )
@variable( m, dw2.Rlow[iS,ib,2][ibprime] <= R2 <= dw2.Rhigh[iS,2] - 1e-9 )
@variable( m, v1 )
@variable( m, v2 )
@objective(m, Max, R1 )
    # The case where we want to maximize just one of the revenues at the
    # expense of the eother
@constraint(m, R1 + R2 + dw2.bgrid[ibprimeidx[iS,ib][ibprime]] ==
                  dw2.bgrid[ib]*(1+r) + sum(g[iS,:]) )
@constraint(m, dirs * [ v1, v2 ] .<= V.dists )
@NLconstraint( m, (1-betta)*(pd_1(R1)-pd_1(dw2.Rlow[iS,ib,1][ibprime])) <=
                    betta*(vbar[1]-v1))
@NLconstraint( m, (1-betta)*(pd_2(R2)-pd_2(dw2.Rlow[iS,ib,2][ibprime])) <=
                    betta*(vbar[2]-v2))

status = solve(m)

println("R1 = ", getvalue(R1))
println("R2 = ", getvalue(R2))
println("v1 = ", getvalue(v1))
println("v2 = ", getvalue(v2))


dists, W1, W2, R1, R2, v1, v2 = search_ic( dw2.bgrid[ibprimeidx[iS,ib][ibprime]],
                  dw2.bgrid[ib], [ dw2.Rlow[iS,ib,k][ibprime] for k in 1:2 ],
                  vec( dw2.Rhigh[iS, : ] ), sum(g[iS,:]), rho, r, betta,
                  vec(dw2.apx_coeffs[iS,:]), vec(dw2.apx_N[iS,:]),
                  dirs, V, vevbar, 0 )


############################################################

cfg = ctsfiscalgame( r=0.04, delta=[.95, .95], psi=psi,
                          chi=chi, rho=.5, A=A, g=g, P=P,
                          nb=nb, bmin=0.0, ndirs=ndir,
                          par=false )


R1, R2, W1, W2, dist =
dirMax( cfg.bgrid[ibprimeidx[iS,ib][ibprime]],
  cfg.bgrid[ib], [ cfg.dw.Rlow[iS,ib,i][ibprimeidx[iS,ib][ibprime]] for i in 1:2 ],
  vec( cfg.dw.Rhigh[iS, :] ), cfg.gSum[iS], cfg.rho, cfg.r, vec(dirs[idir,:]),
  [ cfg.dw.apx_coeffs[iS,i]::Array{Float64,1} for i in 1:2 ],
  [ cfg.dw.apx_N[iS,i]::Array{Float64,1} for i in 1:2 ], true )

bm = @benchmark R1, R2, W1, W2, dist =
  dirMax( cfg.bgrid[ibprimeidx[iS,ib][ibprime]],
    cfg.bgrid[ib], [ cfg.dw.Rlow[iS,ib,i][ibprimeidx[iS,ib][ibprime]] for i in 1:2 ],
    vec( cfg.dw.Rhigh[iS, :] ), cfg.gSum[iS], cfg.rho, cfg.r, vec(dirs[idir,:]),
    [ cfg.dw.apx_coeffs[iS,i]::Array{Float64,1} for i in 1:2 ],
    [ cfg.dw.apx_N[iS,i]::Array{Float64,1} for i in 1:2 ] )
println(bm)

W = initGame( cfg )
vbar = devCont( cfg.dw, cfg.P, cfg.betta, cfg.rho, cfg.A, true)
evbar = P * vbar

pdout = pdPayoffs( cfg, dirs, true )

pdin = pdPayoffs( cfg, dirs, false )

updateout = valsUpdate( W, pdout, cfg.P, dirs, cfg.dw, cfg.gSum,
                          cfg.rho, cfg.r, cfg.betta )


econt = expectedVals( W, P, dirs )
VV = search_ic( dw2.bgrid[ibprimeidx[iS,ib][ibprime]],
                      dw2.bgrid[ib], [ dw2.Rlow[iS,ib,k][ibprime] for k in 1:2 ],
                      vec( dw2.Rhigh[iS, : ] ), cfg.gSum[iS], rho, r, betta,
                      vec(dw2.apx_coeffs[iS,:]), vec(dw2.apx_N[iS,:]),
                      dirs, econt[iS,ib], vec(evbar[iS,:]), 5,
                      "outer" )
updateout_ic = valsUpdate( W, pdout, cfg.P, dirs, cfg.dw, cfg.gSum,
                            cfg.rho, cfg.r, cfg.betta, true, evbar )
updatein_ic = valsUpdate( W, pdout, cfg.P, dirs, cfg.dw, cfg.gSum,
                            cfg.rho, cfg.r, cfg.betta, false, evbar )
updatein = valsUpdate( W, pdin, cfg.P, dirs, cfg.dw, cfg.gSum, cfg.betta, false )

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
