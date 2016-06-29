#= prs.jl
Philip Barrett, pobarrett@gmail.com
29june2016, Chicago

Calculates the perfect risk sharing model =#

include("../julia/autarky.jl")

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
  gam::Vector    # Growth rate
  sig::Vector    # CRRA parameter
  gbar::Vector   # Non-discretionary spending
  nn::Float64     # Rate of population growth

  # Tax process
  P::Matrix       # Transition matrix
  T::Vector       # Taxes
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
function prsModel(; r=0.04, betta=0.9, gam=0.02, sig=1, gbar=.7,
                        nn=0, T=-1, P=-1, chi=1, rho=1, nb=150, bmin=0)
  if( T[1] < 0 || P[1] < 0 )
    T1, P1 = defaultTaxes()
    nT1 = length(T1)
    P = kron(P1,P1)
    T =  [ T1[i] + chi * T1[j] for i in 1:nT1, j in 1:nT1 ][:]
  end

  bgrid = linspace( bmin, ( minimum(T) - gbar[1] - chi * gbar[2] )
                        / ( r - gam ) * .975, nb )

  return prsModel( r, betta, gam, sig, gbar, nn, P, T, length(T), chi,
                    rho, bgrid, nb )
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
