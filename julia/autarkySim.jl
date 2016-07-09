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
  sMC = QuantEcon.MarkovChain( as.am.P )
  sIdx = QuantEcon.simulate(sMC, n+burn)
  aSim = as.am.A[sIdx]
  gSim = as.am.g[sIdx]
  out = ones( n + burn, 9 )
  out[:,1] = aSim
  out[:,2] = gSim
  out[1,3] = bInit
      # Initialize output matrix: (Taxes, b, bprime, g)
  bprimefn = interpolate( as.am, as.bprime )
  Rfn = interpolate( as.am, as.R )
  xfn = interpolate( as.am, as.x )

  for i = 1:(n+burn)
    if( i > 1 )
      out[i,3] = out[i-1,4]
        # Carry over b from previous period
    end
    out[i,4] = bprimefn[ out[i,3], sIdx[i] ]    # Next pd debt
    out[i,5] = Rfn[ out[i,3], sIdx[i] ]         # Revenue
    out[i,6] = xfn[ out[i,3], sIdx[i] ]         # Leisure
    out[i,7] = out[i,5] - out[i,2]              # Surplus
    out[i,8] = out[i,1] * ( 1 - out[i,7] )      # Output
    out[i,9] = out[i,8] - out[i,2]              # Consumption
  end
  return out[ (burn+1):end, : ]
      # Drop the burn
end
