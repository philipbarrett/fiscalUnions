#= prs.jl
Philip Barrett, pobarrett@gmail.com
29june2016, Chicago

Calculates the perfect risk sharing model =#

using Interpolations, Optim

"""
Type defining PRS solution.  Contains:
  * r: Interest rate
  * betta: Discount factor
  * gam: (Constant) growth rate of the economy
  * sig: CRRA curvature coefficient
  * gbar: Minimum required expenditure
  * nn: Rate of population growth
  * P: Transition matrix for taxes
  * T: Vector of tax levels
  * lam: Pareto weight
  * chi: Relative income
  * rho: Relative populations
  * bgrid: The linspace grid of values for the debt
"""
type prsModel
  # Deep parameters
  r::Float64      # Interest rate
  betta::Float64  # Discount factor
  gam::Float64    # Growth rate
  sig::Float64    # CRRA parameter
  gbar::Vector   # Non-discretionary spending
  nn::Float64     # Rate of population growth

  # Tax process
  P::Matrix       # Transition matrix
  T::Vector       # Taxes: Need to sum with chi beforehand
  nT::Int         # Number of Taxes

  # Size parameters
  lam::Float64    # Pareto weight on country 1
  chi::Float64    # Relative income of country 2
  rho::Float64    # Relative population of country 2

  # B grid
  bgrid::LinSpace # Debt grid
  nb::Int         # Length of b grid
end


"""
    function prsModel(; r=0.04, betta=0.9, gam=0.02, T=-1, P=-1, nb=150, bmin=0)
Computes the perfect risk sharing model and formats it as an autarkic model
"""
function prsModel(; r=0.04, betta=0.9, gam=0.02, sig=1,
                      gbar= .7 * ones(2), nn=0, T=-1, P=-1,
                      lam=.5, chi=1, rho=1, nb=150, bmin=0)
  if( T[1] < 0 || P[1] < 0 )
    T1, P1 = defaultTaxes()
    nT1 = length(T1)
    P = kron(P1,P1)
    T =  [ T1[i] + chi * T1[j] for i in 1:nT1, j in 1:nT1 ][:]
  end

  bgrid = linspace( bmin, ( minimum(T) - gbar[1] - chi * gbar[2] )
                        / ( r - gam ) * .975, nb )

  return prsModel( r, betta, gam, sig, gbar, nn, P, T, length(T),
                    lam, chi, rho, bgrid, nb )
end

"""
Type defining the PRS solution.  Includes matrices for V,
bprime and the two gs, as well as the number of iterations required and
the convergence distance.  Also contains the model definition
"""
type prsSol

  # The model
  pm::prsModel

  # The solution objects
  V::Matrix
  bprime::Matrix
  g1::Matrix
  g2::Matrix

  # Algorithm reports
  iter::Int
  dist::Float64
end

"""
    Interpolations.interpolate(pm::prsModel, x::AbstractMatrix)
Given a matrix of size `(pm.nb, pm.nT)`, construct an interpolation
object that does linear interpolation in the asset dimension and has
a lookup table in the z dimension
"""
function Interpolations.interpolate(pm::prsModel, x::AbstractMatrix)
  bt = (pm.nb, pm.nT)
  if size(x) != bt
    msg = "x must have dimensions $(bt)"
    throw(DimensionMismatch(msg))
  end

  itp = interpolate(x, (BSpline(Linear()), NoInterp()), OnGrid())
  scale(itp, pm.bgrid, 1:bt[2])
end

"""
Initialize the matrices for V, bprime, and g
"""
function vbg_init( pm::prsModel )
  rho, chi, gbar, lam, sig = pm.rho, pm.chi, pm.gbar, pm.lam, pm.sig[1]
      # Extract parameters.  Assume common sigma across the two countries
  E = mean( pm.T )
      # Expenditure
  zet = ( lam / ( 1 - lam ) ) ^ ( - 1 / sig )
      # Relative weight thing
  g1 = ones( pm.nb, pm.nT ) * ( ( E - chi * gbar[2] + rho * zet * gbar[1] ) /
                                ( 1 + rho * zet ) )
  g2 = ( E - g1 ) / chi
      # Government consumptions
  bprime = ( ( 1 + pm.r ) * pm.bgrid * ones( 1, pm.nT ) + E -
                ones( pm.nb, 1 ) * pm.T' )/ (1+gam)
  pd1 = ( pm.sig == 1 ) ? log(g1-pm.gbar[1]) : (g1-pm.gbar[1]) .^ (1-pm.sig) ./ (1-pm.sig)
  V = lam * ( 1 + (rho/chi) ^ ( 1 - sig ) * zet ) * pd1
      # Period payoff from aggregate utillity function
  return V, bprime, g1, g2
end

"""
    bellman_operator!(pm::prsModel, V::Matrix, vOut::Matrix,
                      bOut::Matrix, g1Out::Matrix, g2Out::Matrix )
Apply the Bellman operator for a given model and initial value.
##### Arguments
- `pm::prsModel` : Instance of `prsModel`
- `v::Matrix`: Current guess for the value function
- `vOut::Matrix` : Storage for output value function
- `bOut::Matrix` : Storage for output policy function
- `g1Out::Matrix` : Storage for output policy function
- `g2Out::Matrix` : Storage for output policy function
##### Returns
None: `vOut`, `bOut`, `g1Out` and `g2Out` are updated in place.
"""
function bellman_operator!(pm::prsModel, V::Matrix, vOut::Matrix,
                        bOut::Matrix, g1Out::Matrix, g2Out::Matrix )
    # simplify names, set up arrays
  r, betta, gam, sig, gbar, P, T, nT, lam =
      pm.r, pm.betta, pm.gam, pm.sig, pm.gbar, pm.P, pm.T, pm.nT, pm.lam
  bgrid, nb = pm.bgrid, pm.nb
  bmin = minimum(bgrid)
  bmax = maximum(bgrid)

  T_idx = 1:length(T)
  vf = interpolate(pm, V)

  zet = ( lam / ( 1 - lam ) ) ^ ( - 1 / sig )
      # Scaled Pareto weight

  # solve for RHS of Bellman equation
  for (iT, thisT) in enumerate(T), (ib, thisb) in enumerate(bgrid)

    betta_hat = betta * ( (1+gam) / (1+nn) ) ^ (1-sig)
        # Effective discount rate
    function obj(bprime)
      cont = 0.0
      for j in T_idx
          cont += vf[bprime, j] * P[iT, j]
      end
      E = thisT + (1+gam) * bprime - (1+r)*thisb
          # Gov expenditure
      g1 = ( E - chi * gbar[2] + rho * zet * gbar[1] ) /
                                    ( 1 + rho * zet )
          # Gov cons in cty 1
      u1 = ( sig == 1 ) ? log(g1-gbar[1]) : (g1-gbar[1]) ^ (1-sig) / (1-sig)
          # Period utility
      util = lam * ( 1 + (rho/chi) ^ ( 1 - sig ) * zet ) * u1
          # Aggregate utility
      return -( (1-betta_hat)*util + betta_hat * cont )
    end

    opt_lb = ( (1+r) * thisb - thisT + gbar[1] + chi * gbar[2] ) / ( 1 + gam )
        # Debt cannot be so low that the government expenditure
        # is below gbar

    # println("iT", "=", iT )
    # println("thisT", "=", thisT )
    # println("opt_lb", "=", opt_lb )
    # println("bmax", "=", bmax, "\n" )

    res = optimize(obj, opt_lb, bmax + 1e-7 )
    bprime_star = res.minimum

    vOut[ib, iT] = - obj(bprime_star)
    bOut[ib, iT] = bprime_star
    E_star = thisT + (1+gam) * bprime_star - (1+r) * thisb
    g1Out[ib, iT] = ( E_star - chi * gbar[2] + rho * zet * gbar[1] ) /
                                  ( 1 + rho * zet )
    g2Out[ib, iT] = ( E_star - g1Out[ib, iT] ) / chi
  end
  return nothing
end

"""
    solve_pm(pm::prsModel)
Solves the autarky model using value function iteration
"""
function solve_pm(pm::prsModel; tol=1e-6, maxiter=500 )

    V, bprime, g1, g2 = vbg_init(pm)
        # Initialize the output matrices
    Vprime = similar(V)

    dist = 2*tol
    it = 0

    while (tol < dist) && (it < maxiter)
      it += 1
      bellman_operator!( pm, V, Vprime, bprime, g1, g2 )
      dist = maxabs(V - Vprime) / mean( abs(V) )
      copy!(V, Vprime)
      mod(it, 5) == 0 ? println(it, "\t", dist) : nothing
    end
    return prsSol( pm, Vprime, bprime, g1, g2, it, dist )
end
