#= autarkySim.jl
Philip Barrett, pobarrett@gmail.com
25may2016, Chicago

Simulates an AutarkySol object
=#

using QuantEcon: MarkovChain, simulate

"""
    sim_am( as::AutarkySol, n::Int = 1e7, burn::Int = 1e5 )
Simulates an AutarkySol object for n periods, discarding burn periods
"""
function sim_am( as::AutarkySol, n::Int = 1_000_000, burn::Int = 10_000,
                  bInit::Float64=0.0 )

  # Set up
  bInit = 0
  taxMC = QuantEcon.MarkovChain( as.am.P )
  taxIdx = QuantEcon.simulate(taxMC, n+burn)
  taxSim = as.am.T[taxIdx]
  out = ones( n + burn, 4 )
  out[:,1] = taxSim
  out[1,2] = bInit
      # Initialize output matrix: (Taxes, b, bprime, g)
  bprimefn = interpolate( as.am, as.bprime )
  gfn = interpolate( as.am, as.g )

  for i = 1:(n+burn)
    if( i > 1 )
      out[i,2] = out[i-1,3]
        # Carry over b from previous period
    end
    out[i,3] = bprimefn[ out[i,2], taxIdx[i] ]
    out[i,4] = gfn[ out[i,2], taxIdx[i] ]
        # Spending and next period debt
  end
  return out[ (burn+1):end, : ]
      # Drop the burn
end
