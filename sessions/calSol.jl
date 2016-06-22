#= calSol.jl
Philip Barrett, pobarrett@gmail.com
22jun2016, Chicago

Calculates the calibrated solution to the model =#

# 0. Define objects
taxfile = "/home/philip/Dropbox/data/2016/fiscalUnions/taxCoeffs.rdata"

# 1. Libraries
using DataFrames

# 2. Include other files
include("../julia/autarky.jl")
include("../julia/autarkySim.jl")

# 3. Read in the model parameters
tax = read_rda(taxfile)
    # Read up the R data file
n_joint = convert( Int32, tax["n.Z"].data[1] )
n_indiv = convert( Int32, tax["n.ar1"].data[1] )
    # The metaparameters governing the dimensions of the
    # discretized tax matrices
p_joint = zeros( Float64, ( n_joint, n_joint ) )
T_joint = zeros( Float64, ( n_joint, 2 ) )
    # Initialize the joint tax processes
[ p_joint[i] = tax["T.p"].data[i] for i in 1:(n_joint)^2 ]
[ T_joint[i] = tax["T.vals.temp"].data[i] for i in 1:(2*n_joint) ]
    # The joint tax values and transition probabilities
p_indiv = [ zeros( Float64, ( n_indiv, n_indiv ) ) for i in 1:2 ]
T_indiv = [ zeros( Float64, n_indiv ) for i in 1:2 ]
    # Initialize the individual tax processes
[ p_indiv[j][i] = tax["l.indiv"].data[j].data[2].data[i]
                              for i in 1:(n_indiv)^2, j in 1:2 ]
[ T_indiv[j][i] = tax["l.indiv"].data[j].data[1].data[i]
                              for i in 1:n_indiv, j in 1:2 ]
    # The individual tax values and transition probabilities
