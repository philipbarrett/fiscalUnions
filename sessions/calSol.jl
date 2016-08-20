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
include("../julia/prs.jl")
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
delta= [ 1 1 ] # [.999875 .999875] # [ 1.000005 1.000015 ]
nb = 120
tol = 1e-5
maxiter=150
    # Start with a coarse solution
indiv_am = [ AutarkyModel( r=rr, delta=delta[i], psi=psi, chi=chi[i],
                A=states[i][:,1], g=states[i][:,2], P=trans, nR=nR,
                nb=nb ) for i in 1:2 ]
indiv_as = [ solve_am(indiv_am[i], tol=tol, maxiter=maxiter )
                for i in 1:2 ]
# indiv_as = [ solve_am(indiv_am[i], indiv_as[i], tol=tol, maxiter=maxiter )
#                 for i in 1:2 ]
indiv_sim = [ sim_am(indiv_as[i]) for i in 1:2 ]

sim_mu = hcat( [ mean( indiv_sim[i], 1 )' for i in 1:2 ]... )
sim_sd = hcat( [ std( indiv_sim[i], 1 )' for i in 1:2 ]... )
sim_sd_log_xyc = hcat( [ std( log( indiv_sim[i][:, [ 6, 8, 9 ] ] ), 1 )'
                    for i in 1:2 ]... )
sim_cor = [ cor(indiv_sim[i]) for i in 1:2 ]
sim_check = [ indiv_sim[i][:,2] - indiv_sim[i][:,5] -
      indiv_sim[i][:,4] + ( 1 + indiv_am[i].r ) *
      indiv_sim[i][:,3] for i in 1:2 ]

### 2. Two countries, rational expectations equilibrium ###
ndir=24
cfg = ctsfiscalgame( r=rr, delta=[1/rr 1/rr], psi=[psi, psi],
                          chi=vec(chi), rho=rho, A=states_jt[:,1:2], g=states_jt[:,3:4],
                          P=trans, nb=nb, bmin=0.0, ndirs=ndir,
                          par=false )

eq = eqm( cfg, 150 )


### 3. Two countries, SPE ###
