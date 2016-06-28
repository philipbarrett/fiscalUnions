#= calSol.jl
Philip Barrett, pobarrett@gmail.com
24june2016, Chicago

Calculates the calibrated solution to the autarky and first best solutions  =#

### 0. Set up ###

## 0.1 Include function definitions ##
include("../julia/autarky.jl")
include("../julia/autarkySim.jl")
include("../julia/autarkyPlot.jl")

## 0.2 Libraries ##
using DataFrames, Polygons

## 0.3 model settings ##
param_file = "/home/philip/Dropbox/data/2016/fiscalUnions/taxCoeffs.rdata"
    # The location of the parameters


### 1. Read in the calibrated parameters ###
params_R = read_rda(param_file)
    # The raw parameters
nT = convert( Int64, params_R["n.ar1"].data[1] )
nT_jt = convert( Int64, params_R["n.Z"].data[1] )
    # The number of tax levels for hte infdividual and joint processes
gam = params_R["gam"].data[1]
nn = params_R["nn"].data[1]
rr = params_R["rr"].data[1]
    # Direct parameters: Trend revenue growth, population growth, and the real interest rate
debtTaxRatio = params_R["mu.debt.t"].data
sdGovExp = params_R["sd.gc.t"].data
    # Calibration targets: Average debt and standard deviation of governemnt expenditure
taxes = [ exp( params_R["l.indiv"].data[i].data[1].data )
            for i in 1:2 ]
taxes_jt = zeros( nT_jt , 2 )
[ taxes_jt[i] = params_R["T.vals.temp"].data[i] for i in 1:(2*nT_jt) ]
    # Individual and joint taxes
trans = [ zeros( nT , nT ) for i in 1:2 ]
for( i in 1:2)
  [ trans[i][j] = params_R["l.indiv"].data[i].data[2].data[j]
                      for j in 1:(nT^2) ]
end
trans_jt = zeros( nT_jt , nT_jt )
[ trans_jt[j] = params_R["T.p"].data[j] for j in 1:(nT_jt^2) ]
    # Individual and joint

### 1. Solving the model ###
sig = [ 8 8 ]
betta = [ 1.02 1.02 ]
gbar = [ .85 .85 ]
nb = 40
tol = 1e-5
maxiter=600
    # Start with a coarse solution
betta_hat = betta .* ( ( 1 + gam ) / ( 1 + nn ) ) .^ ( 1 - sig )
betta_lim = ( 1 + gam ) / ( 1 + rr )
    # For borrowing need betta_hat < betta_lim
indiv_am = [ AutarkyModel( r=rr, betta=betta[i], gam=gam, sig=sig[i],
                            gbar=gbar[i], nn=nn, T=taxes[i],
                            P=trans[i], nb=nb ) for i in 1:2 ]
indiv_as = [ solve_am(indiv_am[i], tol=tol, maxiter=maxiter ) for i in 1:2 ]
indiv_sim = [ sim_am(indiv_as[i]) for i in 1:2 ]

sim_mu = [ mean( indiv_sim[i], 1 ) for i in 1:2 ]
sim_sd = [ std( indiv_sim[i], 1 ) for i in 1:2 ]

sim_debt_max = [ maximum( indiv_sim[i][:,2] ) for i in 1:2 ]
grid_debt_max = [ maximum(indiv_as[i].am.bgrid ) for i in 1:2 ]
grid_ok = all( sim_debt_max .< grid_debt_max )

# prs_m = prsModel( )
# prs_s = solve_am( prs_m )
# prs_sim = sim_am( prs_s )
#
# plot( asim[1:200,2] )
# plot( asim[1:200,4] )
# scatter( asim[1:200,2 ], asim[1:200,4 ] )
# scatter( asim[1:200,2 ], asim[1:200,4 ] )
#
# cor( asim[1:(end-1),2], asim[2:(end),2])
# cor( asim[1:(end-1),4], asim[2:(end),4])
#
# mean( asim, 1 )
# mean( prs_sim, 1 ) / 2
