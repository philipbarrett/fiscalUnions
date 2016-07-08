A = [ 2.95, 2.975, 3, 3.025, 3.05 ]
g = [ .2, .205, .21, .215, .22 ]
P = [ .4 .3 .2 .1 .0
      .2 .4 .2 .1 .1
      .1 .2 .4 .2 .1
      .1 .1 .2 .4 .2
      .0 .1 .2 .3 .4 ]
nR = 15
psi = .75
chi = 7.0
r = .03
delta= .97


dw = DWL( A=A, g=g, nR=nR, psi=psi, chi=chi, r=r )
am = AutarkyModel( r=r, delta=delta, psi=psi, chi=chi, nR=nR )

DWLplot( am.dw, "tR" )
DWLplot( am.dw, "RW" )
DWLplot( am.dw, "xR" )

V, bprime, R = vbg_init( am )
vOut = copy(V)
bellman_operator!(am, V, vOut, bprime, R )

as = solve_am(am, maxiter=200 )
