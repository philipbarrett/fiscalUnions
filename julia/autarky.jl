#= autarky.jl
Philip Barrett, pobarrett@gmail.com
24may2016, Chicago

Computes the government's optimal expenditure problem in response to uncertain
tax income.  Based on Spencer Lyon's ifp.jl file on quantecon =#
using Interpolations, Optim

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
  delta::Float64  # Extra government discount factor
  psi::Float64    # Firsch elasticity of labor supply
  chi::Float64    # Leisure utility scale parameter

  # Stochastic process
  A::Vector       # Labor productivity
  g::Vector       # Government spending
  P::Matrix       # Transition probabilities
  nS::Int         # Number of States

  # B grid
  bgrid::LinSpace # Debt grid
  nb::Int         # Length of b grid

  # DWL object
  dw::DWL         # Deadweight loss
end

"""
    AutarkyModel( r, betta, gam, T, P, nb, bmin )
Constructor for the autarky model object
"""
function AutarkyModel( ; r=0.04, delta=0.95, psi=0.75, chi=7,
                        A=-1, g=-1, P=-1, nR=20, nb=150, bmin=0)

  if ( A[1] < 0 || g[1] < 0 || P[1] < 0 )
    A, g, P = defaultStates()
  end
  nS = length(A)

  dw =DWL( A=A, g=g, nR=nR, psi=psi, chi=chi, r=r )
      # The deadweight object

  bgrid = linspace( bmin, dw.blim * .975, nb )

  return AutarkyModel( r, 1 / ( 1 +r ), delta, psi, chi, A, g, P, nS,
                        bgrid, nb, dw )
end

"""
A default choice for the states.  Less ugly than putting in the arguments
of the constructor
"""
function defaultStates()

  A = [ 2.95, 2.975, 3, 3.025, 3.05 ]
  g = [ .2, .205, .21, .215, .22 ]
  P = [ .4 .3 .2 .1 .0
        .2 .4 .2 .1 .1
        .1 .2 .4 .2 .1
        .1 .1 .2 .4 .2
        .0 .1 .2 .3 .4 ]
  return A, g, P
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
  R::Matrix

  # Algorithm reports
  iter::Int
  dist::Float64
end

"""
    AutarkySol( am, V, bprime, g, iter, dist )
Constructor for the autarky solution object
"""
function AutarkySol( am, V, bprime, R, iter, dist )
  return AutarkySol( am, V, bprime, R, iter, dist )
end

"""
    Interpolations.interpolate(am::AutarkyModel, x::AbstractMatrix)
Given a matrix of size `(am.nb, am.nS)`, construct an interpolation
object that does linear interpolation in the asset dimension and has
a lookup table in the z dimension
"""
function Interpolations.interpolate(am::AutarkyModel, x::AbstractMatrix)
  bt = (am.nb, am.nS)
  if size(x) != bt
    msg = "x must have dimensions $(bt)"
    throw(DimensionMismatch(msg))
  end

  itp = interpolate(x, (BSpline(Linear()), NoInterp()), OnGrid())
  scale(itp, am.bgrid, 1:bt[2])
end

"""
    function idxLoc( x::Float64, v::Vector )
Finds the non-integer locataton of a value in a vector
"""
function idxLoc( x::Float64, v::Vector )

  v1 = minimum(v)
  vn = maximum(v)
  n = length(v)

  if ( x < v1 || x > vn )
    error("x must be within the range of v")
  end

  if x == v1
    return 1.0
  elseif x== vn
    return Convert( n, Float64 )
  end

  n1 = maximum( find( v .< x ) )
  n2 = n1 + 1
      # Bounding indices
  return n1 + ( x - v[n1] ) / ( v[n2] - v[n1] ) * ( n2 - n1 )
end

"""
Initialize the matrices for V, bprime, and g
"""
function vbg_init( am::AutarkyModel )
  R::Matrix{Float64} = ones( am.nb, am.nS ) * mean( am.g )
      # Constant govt consumption (at avg)
  bprime::Matrix{Float64} = ( 1 + am.r ) * am.bgrid * ones( 1, am.nS ) - R +
                ones( am.nb, 1 ) * am.g'
  bprime[bprime.<am.bgrid[1]] = am.bgrid[1]
  wf = interpolate(am.dw.W, (NoInterp(), BSpline(Linear())), OnGrid())
  Ridx = [ idxLoc( mean(am.g), am.dw.R'[:,i] ) for i in 1:am.nS ]
      # Index of R in interpolation:
  V::Matrix{Float64} = ones( am.nb ) * [ wf[ i, Ridx[i] ] for i in 1:am.nS ]'
      # Copy average vaule across the rows
  return V, bprime, R
end

"""
    bellman_operator!(am::AutarkyModel, V::Matrix,
                      vOut::Matrix, bOut::Matrix, ROut::Matrix )
Apply the Bellman operator for a given model and initial value.
##### Arguments
- `am::AutarkyModel` : Instance of `AutarkyModel`
- `V::Matrix`: Current guess for the value function
- `vOut::Matrix` : Storage for output value function
- `bOut::Matrix` : Storage for output policy function
- `ROut::Matrix` : Storage for output policy function
##### Returns
None: `vOut`, `bOut` and `gOut` are updated in place.
"""
function bellman_operator!(am::AutarkyModel, V::Matrix,
                  vOut::Matrix, bOut::Matrix, ROut::Matrix )
    # simplify names, set up arrays
  r, betta, delta, W, Rgrid, P, A, g, nS, nR =
                    am.r, am.betta, am.delta, am.dw.W, am.dw.R, am.P,
                    am.A, am.g, am.nS, am.dw.nR
  bgrid, nb = am.bgrid, am.nb
  bmin = minimum(bgrid)
  bmax = maximum(bgrid)

  s_idx = 1:length(A)
  vf = interpolate(am, V)
  wf = interpolate(W, (NoInterp(), BSpline(Linear())), OnGrid())
  RR = Rgrid'

  # solve for RHS of Bellman equation
  for iS in s_idx, (ib, thisb) in enumerate(bgrid)

    betta_hat = betta * delta
        # Effective discount rate
    function obj(bprime)
      cont = 0.0
      for j in s_idx
          cont += vf[bprime, j] * P[iS, j]
      end
      R = g[iS] + ( 1 + r ) * thisb - bprime
          # Required revenue
      Ridx = idxLoc( R, RR[:,iS] )
          # Fractional index of revenue
      loss = wf[ iS, Ridx ]
          # Period deadweight loss
      return ( 1 - betta_hat ) * loss + betta_hat * cont
    end

    opt_lb = max( bmin, (1+r) * thisb + g[iS] - Rgrid[iS, nR] )
        # Debt cannot be so low that the government expenditure
        # is below gbar
    opt_ub = min( bmax, (1+r) * thisb + g[iS] - Rgrid[iS, 1] )

    res = optimize(obj, opt_lb, opt_ub )
    bprime_star = res.minimum

    vOut[ib, iS] = obj(bprime_star)
    bOut[ib, iS] = bprime_star
    ROut[ib, iS] = g[iS] + ( 1 + r ) * thisb - bprime_star
  end
  return nothing
end

"""
    solve_am(am::AutarkyModel)
Solves the autarky model using value function iteration
"""
function solve_am(am::AutarkyModel; tol=1e-6, maxiter=500 )

    V, bprime, R = vbg_init(am)
        # Initialize the output matrices
    Vprime = similar(V)

    dist = 2*tol
    it = 0

    while (tol < dist) && (it < maxiter)
      it += 1
      bellman_operator!( am, V, Vprime, bprime, R )
      dist = maxabs(V - Vprime) / mean( abs(V) )
      copy!(V, Vprime)
      mod(it, 10) == 0 ? println(it, "\t", dist) : nothing
    end
    return AutarkySol( am, Vprime, bprime, R, it, dist )
end

# function solve_am(am::AutarkyModel, am_init::AutarkyModel; tol=1e-6, maxiter=500 )
#
#     V, bprime, g = am_init.V, am_init.bprime, am_init.g
#         # Initialize the output matrices
#     Vprime = similar(V)
#
#     dist = 2*tol
#     it = 0
#
#     while (tol < dist) && (it < maxiter)
#       it += 1
#       bellman_operator!( am, V, Vprime, bprime, g )
#       dist = maxabs(V - Vprime) / mean( abs(V) )
#       copy!(V, Vprime)
#       mod(it, 10) == 0 ? println(it, "\t", dist) : nothing
#     end
#     return AutarkySol( am, Vprime, bprime, g, it, dist )
# end
#
# #
# # hh = AutarkyModel()
# # kk = solve_am(hh)
# # plot( hh.bgrid, kk.V )
# # plot( hh.bgrid, [ kk.bprime hh.bgrid ], xlims=(1.8,2), ylims=(1.8,2) )
# # plot( hh.bgrid, kk.g, xlims=(1.5,2), ylims=(0,.15) )
