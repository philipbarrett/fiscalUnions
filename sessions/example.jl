using Polygons
using CHull2D
using JLD

include("../julia/dwl.jl")
include("../julia/fiscalGame.jl")
include("../julia/fiscalSol.jl")

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
nb = 20
r = .03
delta = 1.0
ndirsl = 8

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
              A=A_jt, g=g_jt, P=P, nR=nR, rho=.5, nb=nb, ndirsl=ndirsl )
init, ndirs = initGame(fg)
polyPlot(init[:,1])

betta_hat = fg.betta * fg.delta
thisndirs = 20
dirs = hcat( [ cos(2*pi*(i-1)/thisndirs)::Float64 for i in 1:thisndirs ],
             [ sin(2*pi*(i-1)/thisndirs)::Float64 for i in 1:thisndirs ] )

# ff = uncSetUpdate( fg.surp[3,2], vec(fg.pdLoss[3,:,2]), fg.bgrid[4],
#                          r, fg.bgrid, init, vec(P[2,:]), betta_hat,
#                          dirs )
# gg = uncSetUpdate( fg.surp[3,2], vec(fg.pdLoss[3,:,2]), fg.bgrid[4],
#                         r, fg.bgrid, init, vec(P[2,:]), betta_hat,
#                         dirs, false )
# polyPlot([ff, gg])
#
# hh = uncSetUpdate( fg.surp[1,1], vec(fg.pdLoss[1,:,1]), fg.bgrid[1],
#                   r, fg.bgrid, init, vec(P[1,:]), betta_hat, dirs )
# iS = 1
# ib = 10
# for k in find(fg.potFeas[is,ib])
#   println( "k=", k )
#   jj = uncSetUpdate( fg.surp[k,is], vec(fg.pdLoss[k,:,is]), fg.bgrid[ib],
#                                     fg.r, fg.bgrid, init, vec(fg.P[is,:]),
#                                     betta_hat, fg.dirs )
# end
#
# kk = [ uncSetUpdate( fg.surp[k,iS], vec(fg.pdLoss[k,:,iS]), fg.bgrid[ib],
#                                   fg.r, fg.bgrid, init, vec(fg.P[iS,:]), betta_hat,
#                                   dirs )::Polygon
#                       for k in find(fg.potFeas[iS,ib]) ]
#     # Compilation
# using BenchmarkTools
# bmk = @benchmark kk = [ uncSetUpdate( fg.surp[k,iS], vec(fg.pdLoss[k,:,iS]), fg.bgrid[ib],
#                                   fg.r, fg.bgrid, init, vec(fg.P[iS,:]), betta_hat,
#                                   dirs )::Polygon
#                                   for k in find(fg.potFeas[iS,ib]) ]
#     # Averages about 3.5ms (0.0035s)
#     # Following aggresive typing in Polygons, this went *up*
# println(bmk)

# W = uncSetUpdate( fg, init, ndirs )
W = init

it = 0
hd = 10.0

while it < 80 && maximum(hd) > 1e-05
  it += 1
  println("*** it = ", it , " ***")
  W_new = uncSetUpdate( fg, W, ndirs )
  hd = hausdorff( W, W_new )
  println("   ** mean & max hausdorff dist = ", ( mean(hd), maximum(hd) ), "**" )
  W = copy(W_new)
  ndirs = ndirsUpdate( ndirs, hd, fg.hddirs, fg.ndirsu )
  println("   ** mean # search dirs = ", mean(ndirs), "**" )
  println("   ** extremal # search dirs = ", extrema(ndirs), "**\n" )
end




# jldopen("/home/philip/Dropbox/data/2016/fiscalUnions/unc.jld", "w") do file
#     addrequire(file, Polygons)
#     write(file, "W", W)
# end

# Can then acecss W with:
# d = load("/home/philip/Dropbox/data/2016/fiscalUnions/unc.jld")
# W = d["W"]
