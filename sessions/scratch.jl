A = [ 2.95, 2.975, 3, 3.025, 3.05 ]
g = [ .205, .2125, .21, .2075, .215 ]
A_jt = [ 2.95 2.95
         3.05 3.05
         3.05 2.95
         2.95 3.05 ]
g_jt = [ .21   .21
         .21   .21
         .215   .205
         .205   .215 ]
P = [ .4 .2 .2 .2
      .2 .4 .2 .2
      .2 .2 .4 .2
      .2 .2 .2 .4 ]
nR = 20
psi = .75
chi = 7.0
nb=10
r = .03
delta= 1.0


# dw = DWL( A=A, g=g, nR=nR, psi=psi, chi=chi, r=r )
# am = AutarkyModel( r=r, delta=delta, psi=psi, chi=chi, nR=nR )
#
# DWLplot( am.dw, "tR" )
# DWLplot( am.dw, "RW" )
# DWLplot( am.dw, "xR" )
#
# V, bprime, R, x = vbg_init( am )
# vOut = copy(V)
# bellman_operator!(am, V,vOut, bprime, R, x )
#
# as = solve_am(am, tol=1e-5, maxiter=200 )
# asim = sim_am(as)
#
# pm = prsModel( r=r, delta=delta, psi=[ psi psi ], chi=[ chi chi ],
#                 A=A_jt, g=g_jt, P=P, nR=[nR nR], rho=.5, lam=.5 )
# V, bprime, R1, R2, x1, x2 = vbg_init( pm )
# vOut = copy(V)
# bellman_operator!(pm, V,vOut, bprime, R1, R2, x1, x2 )

fg = FiscalGame( r=r, delta=[delta, delta], psi=[ psi, psi ], chi=[ chi, chi ],
                A=A_jt, g=g_jt, P=P, nR=nR, rho=.5, nb=nb )
init = initGame(fg)
polyPlot(init[:,1])

betta_hat = fg.betta * fg.delta
ff = uncSetUpdate( fg.surp[3,2], vec(fg.pdLoss[3,:,2]), fg.bgrid[4],
                        r, fg.bgrid, init, vec(P[2,:]), betta_hat,
                        fg.dirs )
gg = uncSetUpdate( fg.surp[3,2], vec(fg.pdLoss[3,:,2]), fg.bgrid[4],
                        r, fg.bgrid, init, vec(P[2,:]), betta_hat,
                        fg.dirs, false )
polyPlot([ff, gg])

hh = uncSetUpdate( fg.surp[1,1], vec(fg.pdLoss[1,:,1]), fg.bgrid[1],
                  r, fg.bgrid, init, vec(P[1,:]), betta_hat, fg.dirs )
# is = 1
# ib = 10
# for k in find(fg.potFeas[is,ib])
#   println( "k=", k )
#   jj = uncSetUpdate( fg.surp[k,is], vec(fg.pdLoss[k,:,is]), fg.bgrid[ib],
#                                     fg.r, fg.bgrid, init, vec(fg.P[is,:]),
#                                     betta_hat, fg.dirs )
# end
# kk = [ uncSetUpdate( fg.surp[k,is], vec(fg.pdLoss[k,:,is]), fg.bgrid[ib],
#                                   fg.r, fg.bgrid, init, vec(fg.P[is,:]), betta_hat,
#                                   fg.dirs )::Polygon
#                       for k in find(fg.potFeas[is,ib]) ]


W = uncSetUpdate( fg, init )

for it in 1:10
  println("*** it = ", it , " ***")
  W_new = uncSetUpdate( fg, W )
  W = copy(W_new)
end

using JLD

jldopen("/home/philip/Dropbox/data/2016/fiscalUnions/unc.jld", "w") do file
    addrequire(file, Polygons)
    write(file, "W", W)
end

# Can then acecss W with:
# d = load("/home/philip/Dropbox/data/2016/fiscalUnions/unc.jld")
# W = d["w"]
