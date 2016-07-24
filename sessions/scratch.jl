nS = 3
nb = 50
A = [ 2.95 3.05
      3.00 3.00
      3.05 2.95 ]
g =  .05 * ones( nS, 2 )
psi = [ .75, .75 ]
chi = [ 1.75, 1.75 ]
rho = .5
r = .03
bmin = 0.0

dw2 = dwlc2( nS=nS, nb=nb, A=A, g=g, psi=psi, chi=chi,
              r=r, bmin=bmin, rho=rho )

iS, ib, ibprimeidx = 1,6,5
    # NB: ibprimeidx is on the grid 1:(nposs[iS,ib])
ibprime = find( dw2.bprimeposs[iS,ib] )[ibprimeidx]
    # Index of bprime in grid, so brpime = bgrid[ibprime]
# ibprimepossidx = sum( dw2.bprimeposs[iS,ib][1:ib] )
#     # Index of bprime in possible
R1, R2, dist =
  dirMax( dw2.bgrid[ibprime], dw2.bgrid[ib], chi, psi,
    [ dw2.xlow[iS,ib,i][ibprimeidx] for i in 1:2 ],
    vec( dw2.xhigh[iS,:]),
    [ dw2.Rlow[iS,ib,1][ibprimeidx], dw2.Rhigh[iS,1] ],
    vec(A[iS,:]), sum(g[iS,:]), rho, r, [ -1, -1 ] )


    dirMax( bprime::Float64, b::Float64, chi::Vector, psi::Vector,
                      x1::Vector, x2::Vector, R1bds::Vector,
                      a::Vector, sumg::Float64,
