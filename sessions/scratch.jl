A = 3 + .025 * [ i for i in -6:6]
g = .25 +  .01 * [ i for i in -6:6]
n = 20
psi = .75
chi = 7.0
r = .03

x1, x2, x, tau, R, c, W, blim = dwl( A, g, n, psi, chi, r)
