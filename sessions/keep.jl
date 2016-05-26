#=== Things to keep ===#

# Include function definitions
include("../julia/autarky.jl")
include("../julia/autarkySim.jl")

# Solving the model
am = AutarkyModel()
as = solve_am(am)
asim = sim_am(as)

prs_m = prsModel( )
prs_s = solve_am( prs_m )
prs_sim = sim_am( prs_s )

# plot( asim[1:200,2] )
# plot( asim[1:200,4] )
# scatter( asim[1:200,2 ], asim[1:200,4 ] )

cor( asim[1:(end-1),2], asim[2:(end),2])
cor( asim[1:(end-1),4], asim[2:(end),4])

println( mean( asim, 1 ) )
println( mean( prs_sim, 1 ) )

println( std( asim, 1 ) )
println( std( prs_sim, 1 ) )
