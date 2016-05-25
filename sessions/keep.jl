#=== Things to keep ===#

# Include function definitions
include("../julia/autarky.jl")
include("../julia/autarkySim.jl")

# Solving the model
am = AutarkyModel()
as = solve_am(am)
asim = sim_am(as)
