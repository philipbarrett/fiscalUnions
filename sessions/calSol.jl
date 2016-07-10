#= calSol.jl
Philip Barrett, pobarrett@gmail.com
24june2016, Chicago

Calculates the calibrated solution to the autarky and first best solutions  =#

### 0. Set up ###

## 0.1 Include function definitions ##
include("../julia/dwl.jl")
include("../julia/autarky.jl")
include("../julia/autarkySim.jl")
include("../julia/autarkyPlot.jl")
# include("../julia/prs.jl")
# include("../julia/prsSim.jl")

## 0.2 Libraries ##
using DataFrames, Polygons

## 0.3 model settings ##
param_file = "/home/philip/Dropbox/data/2016/fiscalUnions/cal.rdata"
    # The location of the parameters


### 1. Read in the calibrated parameters ###
params_R = read_rda(param_file)
    # The raw parameters
nS = convert( Int64, params_R["n.Z.joint"].data[1] )
    # The number of tax levels for hte infdividual and joint processes
rr = params_R["rr"].data[1]
    # Direct parameters: Trend revenue growth, population growth, and the real interest rate
debtTaxRatio = params_R["mu.debt"].data
    # Calibration target: Average debt
rho = params_R["rho"].data[1]
    # The scale paameters
##########
states = [ zeros( nS , 2 ) for i in 1:2 ]
states[1][:,1] = params_R["T.vals"].data[1:nS]
states[2][:,1] = params_R["T.vals"].data[nS + (1:nS)]
states[1][:,2] = params_R["T.vals"].data[2*nS + (1:nS)]
states[2][:,2] = params_R["T.vals"].data[3*nS + (1:nS)]
    # States for the individual problem
states_jt = zeros( nS , 4 )
[ states_jt[i] = params_R["T.vals.joint"].data[i] for i in 1:(4*nS) ]
    # Individual and joint taxes
trans = zeros( nS , nS )
[ trans[j] = params_R["T.p.joint"].data[j] for j in 1:(nS^2) ]
    # Transition probabilities
trans = trans ./ ( sum(trans, 2) * ones( 1, nS ) )
    # Just make sure that the rows sum to one for sure

### 1. Solving the model: One country ###
psi = 0.75
nR = 50
chi = [ 13.0 13.0 ]
delta= [ 1.0005 1.001 ]
nb = 120
tol = 1e-5
maxiter=140
    # Start with a coarse solution
indiv_am = [ AutarkyModel( r=rr, delta=delta[i], psi=psi, chi=chi[i],
                A=states[i][:,1], g=states[i][:,2], P=trans, nR=nR,
                nb=nb ) for i in 1:2 ]
# indiv_as = [ solve_am(indiv_am[i], tol=tol, maxiter=maxiter )
#                 for i in 1:2 ]
indiv_as = [ solve_am(indiv_am[i], indiv_as[i], tol=tol, maxiter=maxiter )
                for i in 1:2 ]
indiv_sim = [ sim_am(indiv_as[i]) for i in 1:2 ]

sim_mu = hcat( [ mean( indiv_sim[i], 1 )' for i in 1:2 ]... )
sim_sd = hcat( [ std( indiv_sim[i], 1 )' for i in 1:2 ]... )
sim_sd_log_xyc = hcat( [ std( log( indiv_sim[i][:, [ 6, 8, 9 ] ] ), 1 )'
                    for i in 1:2 ]... )

# play([sin(x) for x=0:0.03*pi:441])

# gc_sd = [ sim_sd[i][4] for i in 1:2 ]
# debt_mu = [ sim_mu[i][2] for i in 1:2 ]
#
# sim_debt_max = [ maximum( indiv_sim[i][:,2] ) for i in 1:2 ]
# grid_debt_max = [ maximum(indiv_as[i].am.bgrid ) for i in 1:2 ]
# grid_ok = all( sim_debt_max .< grid_debt_max )
#
# tax_sum::Vector = taxes_jt[:,1] + chi * taxes_jt[:,2]
# prs_m = prsModel( r=rr, betta=betta[1], gam=gam, sig=sig[1], gbar=gbar,
#                   nn=nn, T=tax_sum, P=trans_jt, lam=.5, chi=chi,
#                   rho=rho, nb=nb )
# prs_s = solve_pm( prs_m )
# prs_sim = sim_pm( prs_s )
#
# mu_prs = mean( prs_sim, 1 ) ./ [ 1+chi 1+chi 1+chi 1 1 ]
# sd_prs = std( prs_sim, 1 ) ./ [ 1+chi 1+chi 1+chi 1 1 ]




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
