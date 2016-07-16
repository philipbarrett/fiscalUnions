#= parEg.jl
Philip Barrett, pobarrett@gmail.com
14jul2016, Chicago

Example with parallel calculation of sets=#

# Idea for a command line call:
#   julia -p <n> -L file1.jl -L file2.jl driver.jl
# Loads file1.jl, file2.jl etc

nprocs = 70 
using JLD, ClusterManagers
addprocs_pbs(nprocs)

@everywhere wd = "/home/pobarrett/code/2016/fiscalUnions"

@everywhere include("$wd/julia/dwl.jl")
@everywhere include("$wd/julia/fiscalGame.jl")
@everywhere include("$wd/julia/fiscalSol.jl")

@everywhere using Polygons, CHull2D

saveloc = "/home/pobarrett/data/2016/fiscalUnions/unc_32.jld"

A = [ 2.95, 2.975, 3, 3.025, 3.05 ]
g = [ .205, .2125, .21, .2075, .215 ]
A_jt = [ 2.95 2.95
         3.05 3.05
         3.05 2.95
         2.95 3.05 ]
g_jt = [ .21   .21
         .21   .21
         .205   .215
         .215   .205 ]
P = [ .4 .2 .2 .2
      .2 .4 .2 .2
      .2 .2 .4 .2
      .2 .2 .2 .4 ]
nR = 40
psi = .75
chi = 7.0
nb = 30
r = .03
delta = 1.0
ndirsl = 32
ndirsu = 32
par = true 

fg_par = FiscalGame( r=r, delta=[delta, delta], psi=[ psi, psi ], chi=[ chi, chi ],
              A=A_jt, g=g_jt, P=P, nR=nR, rho=.5, nb=nb, ndirsl=ndirsl, ndirsu=ndirsu, par=par )
# fg_ser = FiscalGame( r=r, delta=[delta, delta], psi=[ psi, psi ], chi=[ chi, chi ],
#              A=A_jt, g=g_jt, P=P, nR=nR, rho=.5, nb=nb, ndirsl=ndirsl, par=!par )
    # Serial and parallel fiscal games

init, ndirs = initGame(fg_par)
    # Initialize


# @time W_ser = uncSetUpdate( fg_ser, init, ndirs )
W_par = uncSetUpdate( fg_par, init, ndirs )
    # Compilation - also gets the warnings out of the way! :)
@time W_par = uncSetUpdate( fg_par, init, ndirs )

# hd = hausdorff( W_par, W_ser )
    # Just check that this is all ok

# println( "max(hd) = ", maximum(hd) )

W_sol, ndirs, hd = uncSol( fg_par, init, ndirs, true, 500, 1e-05, saveloc )

