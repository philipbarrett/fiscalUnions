#=== Things to keep ===#

# Include function definitions
include("../julia/autarky.jl")
include("../julia/autarkySim.jl")

# Solving the model
am = AutarkyModel()
as = solve_am(am)
asim = sim_am(as)

plot( asim[1:200,2] )
plot( asim[1:200,4] )
scatter( asim[1:200,2 ], asim[1:200,4 ] )

cor( asim[1:(end-1),2], asim[2:(end),2])
cor( asim[1:(end-1),4], asim[2:(end),4])
