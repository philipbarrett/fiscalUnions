#= prsSim.jl
Philip Barrett, pobarrett@gmail.com
30jun2016, Chicago

Simulates a prsSol object
=#

using QuantEcon: MarkovChain, simulate

"""
    sim_pm( ps::prsSol, n::Int = 1e7, burn::Int = 1e5 )
Simulates a prsSol object for n periods, discarding burn periods
"""
function sim_pm( ps::prsSol, n::Int = 1_000_000, burn::Int = 10_000,
                  bInit::Float64=0.0 )

  # Set up
  bInit = 0
  taxMC = QuantEcon.MarkovChain( ps.pm.P )
  taxIdx = QuantEcon.simulate(taxMC, n+burn)
  taxSim = ps.pm.T[taxIdx]
  out = ones( n + burn, 5 )
  out[:,1] = taxSim
  out[1,2] = bInit
      # Initialize output matrix: (Taxes, b, bprime, g)
  bprimefn = interpolate( ps.pm, ps.bprime )
  g1fn = interpolate( ps.pm, ps.g1 )
  g2fn = interpolate( ps.pm, ps.g2 )

  for i = 1:(n+burn)
    if( i > 1 )
      out[i,2] = out[i-1,3]
        # Carry over b from previous period
    end
    out[i,3] = bprimefn[ out[i,2], taxIdx[i] ]
    out[i,4] = g1fn[ out[i,2], taxIdx[i] ]
    out[i,5] = g2fn[ out[i,2], taxIdx[i] ]
        # Spending and next period debt
  end
  return out[ (burn+1):end, : ]
      # Drop the burn
end
