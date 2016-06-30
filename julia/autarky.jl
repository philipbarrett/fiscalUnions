#= autarky.jl
Philip Barrett, pobarrett@gmail.com
24may2016, Chicago

Computes the government's optimal expenditure problem in response to uncertain
tax income.  Based on Spencer Lyon's ifp.jl file on quantecon =#
using Interpolations, Optim

# # utility and marginal utility functions
# u(x) = log(x)

"""
Type defining Autarky solution.  Contains:
  * r: Interest rate
  * betta: Discount factor
  * gam: (Constant) growth rate of the economy
  * sig: CRRA curvature coefficient
  * gbar: Minimum required expenditure
  * nn: Rate of population growth
  * P: Transition matrix for taxes
  * T: Vector of tax levels
  * bgrid: The linspace grid of values for the debt
"""
type AutarkyModel
  # Deep parameters
  r::Float64      # Interest rate
  betta::Float64  # Discount factor
  gam::Float64    # Growth rate
  sig::Float64    # CRRA parameter
  gbar::Float64   # Non-discretionary spending
  nn::Float64     # Rate of population growth

  # Tax process
  P::Matrix       # Transition matrix
  T::Vector       # Taxes
  nT::Int         # Number of Taxes

  # B grid
  bgrid::LinSpace # Debt grid
  nb::Int         # Length of b grid
end

"""
    AutarkyModel( r, betta, gam, T, P, nb, bmin )
Constructor for the autarky model object
"""
function AutarkyModel( ; r=0.04, betta=0.9, gam=0.02, sig=1, gbar=.7,
                        nn=0, T=-1, P=-1, nb=150, bmin=0)

  if ( T[1] < 0 || P[1] < 0 )
    T, P = defaultTaxes()
  end

  bgrid = linspace( bmin, ( minimum(T) - gbar ) / ( r - gam ) *
                              .975, nb )

  return AutarkyModel( r, betta, gam, sig, gbar, nn, P, T, length(T), bgrid, nb )
end

"""
A default choice for taxes.  Less ugly than putting in the arguments
of the constructor
"""
function defaultTaxes()

  T = [ .02, .04, .06, .08 ]
  P = [ .6 .2 .1 .1
        .2 .4 .3 .1
        .1 .3 .4 .2
        .1 .1 .2 .6 ]

  return T, P
end

"""
Type defining the autarky solution.  Includes matrices for V,
bprime and g, as well as the number of iterations required and
the convergence distance.  Also contains the model definition
"""
type AutarkySol

  # The model
  am::AutarkyModel

  # The solution objects
  V::Matrix
  bprime::Matrix
  g::Matrix

  # Algorithm reports
  iter::Int
  dist::Float64
end

"""
    AutarkySol( am, V, bprime, g, iter, dist )
Constructor for the autarky solution object
"""
function AutarkySol( am, V, bprime, g, iter, dist )
  return AutarkySol( am, V, bprime, g, iter, dist )
end

"""
    Interpolations.interpolate(am::AutarkyModel, x::AbstractMatrix)
Given a matrix of size `(am.nb, am.nT)`, construct an interpolation
object that does linear interpolation in the asset dimension and has
a lookup table in the z dimension
"""
function Interpolations.interpolate(am::AutarkyModel, x::AbstractMatrix)
  bt = (am.nb, am.nT)
  if size(x) != bt
    msg = "x must have dimensions $(bt)"
    throw(DimensionMismatch(msg))
  end

  itp = interpolate(x, (BSpline(Linear()), NoInterp()), OnGrid())
  scale(itp, am.bgrid, 1:bt[2])
end

"""
Initialize the matrices for V, bprime, and g
"""
function vbg_init( am::AutarkyModel )
  g = ones( am.nb, am.nT ) * mean( am.T )
      # Constant govt consumption (at avg)
  bprime = ( ( 1 + am.r ) * am.bgrid * ones( 1, am.nT ) + g -
                ones( am.nb, 1 ) * am.T' )/ (1+gam)
  pd = ( am.sig == 1 ) ? log(g-am.gbar) : (g-am.gbar) .^ (1-am.sig) ./ (1-am.sig)
  V = pd
  return V, bprime, g
end

"""
    bellman_operator!(am::AutarkyModel, V::Matrix,
                      vOut::Matrix, bOut::Matrix, gOut::Matrix )
Apply the Bellman operator for a given model and initial value.
##### Arguments
- `am::AutarkyModel` : Instance of `AutarkyModel`
- `v::Matrix`: Current guess for the value function
- `vOut::Matrix` : Storage for output value function
- `bOut::Matrix` : Storage for output policy function
- `gOut::Matrix` : Storage for output policy function
##### Returns
None: `vOut`, `bOut` and `gOut` are updated in place.
"""
function bellman_operator!(am::AutarkyModel, V::Matrix,
                  vOut::Matrix, bOut::Matrix, gOut::Matrix )
    # simplify names, set up arrays
  r, betta, gam, sig, gbar, P, T, nT =
        am.r, am.betta, am.gam, am.sig, am.gbar, am.P, am.T, am.nT
  bgrid, nb = am.bgrid, am.nb
  bmin = minimum(bgrid)
  bmax = maximum(bgrid)

  T_idx = 1:length(T)
  vf = interpolate(am, V)

  # solve for RHS of Bellman equation
  for (iT, thisT) in enumerate(T), (ib, thisb) in enumerate(bgrid)

    betta_hat = betta * ( (1+gam) / (1+nn) ) ^ (1-sig)
        # Effective discount rate
    function obj(bprime)
      cont = 0.0
      for j in T_idx
          cont += vf[bprime, j] * P[iT, j]
      end
      g = thisT + (1+gam) * bprime - (1+r)*thisb
          # Gov expenditure
      util = ( sig == 1 ) ? log(g-gbar) : (g-gbar) ^ (1-sig) / (1-sig)
          # Period utility

      return -( (1-betta_hat)*util + betta_hat * cont )
    end

    opt_lb = ( (1+r) * thisb - thisT + gbar ) / ( 1 + gam )
        # Debt cannot be so low that the government expenditure
        # is below gbar

    res = optimize(obj, opt_lb, bmax + 1e-10 )
    bprime_star = res.minimum

    vOut[ib, iT] = - obj(bprime_star)
    bOut[ib, iT] = bprime_star
    gOut[ib, iT] = thisT + (1+gam) * bprime_star - (1+r) * thisb
  end
  return nothing
end

"""
    solve_am(am::AutarkyModel)
Solves the autarky model using value function iteration
"""
function solve_am(am::AutarkyModel; tol=1e-6, maxiter=500 )

    V, bprime, g = vbg_init(am)
        # Initialize the output matrices
    Vprime = similar(V)

    dist = 2*tol
    it = 0

    while (tol < dist) && (it < maxiter)
      it += 1
      bellman_operator!( am, V, Vprime, bprime, g )
      dist = maxabs(V - Vprime) / mean( abs(V) )
      copy!(V, Vprime)
      mod(it, 10) == 0 ? println(it, "\t", dist) : nothing
    end
    return AutarkySol( am, Vprime, bprime, g, it, dist )
end

function solve_am(am::AutarkyModel, am_init::AutarkyModel; tol=1e-6, maxiter=500 )

    V, bprime, g = am_init.V, am_init.bprime, am_init.g
        # Initialize the output matrices
    Vprime = similar(V)

    dist = 2*tol
    it = 0

    while (tol < dist) && (it < maxiter)
      it += 1
      bellman_operator!( am, V, Vprime, bprime, g )
      dist = maxabs(V - Vprime) / mean( abs(V) )
      copy!(V, Vprime)
      mod(it, 10) == 0 ? println(it, "\t", dist) : nothing
    end
    return AutarkySol( am, Vprime, bprime, g, it, dist )
end

#
# hh = AutarkyModel()
# kk = solve_am(hh)
# plot( hh.bgrid, kk.V )
# plot( hh.bgrid, [ kk.bprime hh.bgrid ], xlims=(1.8,2), ylims=(1.8,2) )
# plot( hh.bgrid, kk.g, xlims=(1.5,2), ylims=(0,.15) )
