#= calSol.jl
Philip Barrett, pobarrett@gmail.com
24june2016, Chicago

Calculates the calibrated solution to the autarky and first best solutions  =#

### 0. Set up ###

## 0.1 Include function definitions ##
include("../julia/autarky.jl")
include("../julia/autarkySim.jl")

## 0.2 Libraries ##
using DataFrames, Polygons

## 0.3 model settings ##
param_file = "/home/philip/Dropbox/data/2016/fiscalUnions/taxCoeffs.rdata"
    # The location of the parameters


### 1. Read in the calibrated parameters ###
params_R = read_rda(param_file)
    # The raw parameters
nT = params_R["n.ar1"].data[1]
nT_jt = params_R["n.Z"].data[1]
    # The number of tax levels for hte infdividual and joint processes





### 1. Solving the model ###
am = AutarkyModel()
as = solve_am(am)
asim = sim_am(as)

prs_m = prsModel( )
prs_s = solve_am( prs_m )
prs_sim = sim_am( prs_s )

plot( asim[1:200,2] )
plot( asim[1:200,4] )
scatter( asim[1:200,2 ], asim[1:200,4 ] )
scatter( asim[1:200,2 ], asim[1:200,4 ] )

cor( asim[1:(end-1),2], asim[2:(end),2])
cor( asim[1:(end-1),4], asim[2:(end),4])

mean( asim, 1 )
mean( prs_sim, 1 ) / 2
