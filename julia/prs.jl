#= prs.jl
Philip Barrett, pobarrett@gmail.com
29june2016, Chicago

Calculates the perfect risk sharing model =#

using Interpolations, Optim

include("dwl.jl")

"""
Type defining PRS solution.  Contains:
  * r: Interest rate
  * betta: Discount factor
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
  delta::Float64  # Extra government discount factor
  psi::Vector{Float64}    # Frisch elasticity of labor supply
  chi::Vector{Float64}    # Leisure utility scale parameter

  # Stochastic process
  A::Matrix{Float64}      # Labor productivity
  g::Matrix{Float64}      # Government spending
  P::Matrix{Float64}      # Transition probabilities
  nS::Int                 # Number of States

  # B grid
  bgrid::LinSpace # Debt grid
  nb::Int         # Length of b grid

  # DWL object
  dw::Array{DWL,1}        # Deadweight loss

  # Size parameters
  lam::Float64    # Pareto weight on country 1
  rho::Float64    # Population of country 1 (cty 2 = 1-rho)

end


"""
    function prsModel(; r=0.04, betta=0.9, gam=0.02, T=-1, P=-1, nb=150, bmin=0)
Computes the perfect risk sharing model and formats it as an autarkic model
"""
function prsModel(; r=0.04, delta=0.95, psi=[0.75 .75], chi=[7 7], A=-1,
                    g=-1, P=-1, lam=.5, nR=[20 20], rho=.5, nb=150, bmin=0)

  if ( A[1] < 0 || g[1] < 0 || P[1] < 0 )
    Atemp, gtemp, P = defaultStates()
    A = hcat( Atemp, Atemp )
    g = hcat( gtemp, gtemp )
  end
  nS = size(A)[1]
  vrho = [ rho 1-rho ]
      # Vector of weights

  dw =[ DWL( A=A[:,i], g=g[:,i], nR=nR[i], psi=psi[i],
                chi=chi[i], r=r, rho=vrho[i] ) for i in 1:2 ]
      # Deadwight loss object

  blim = 1 / r * minimum( dw[1].R[:,nR[1]] - g[:,1] +
                            dw[2].R[:,nR[2]] - g[:,2] )
      # Maximum joint borrowing

  bgrid = linspace( bmin, blim * .975, nb )

  return prsModel( r, 1 / ( 1 +r ), delta, squeeze(psi, 1),
                      squeeze(chi, 1), A, g, P, nS, bgrid,
                      nb, dw, lam, rho )
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
  R1::Matrix
  R2::Matrix
  x1::Matrix
  x2::Matrix

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
  bt = (pm.nb, pm.nS)
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
  rho, lam = pm.rho, pm.lam
      # Extract parameters.  Assume common sigma across the two countries
  vrho = [ rho, 1-rho ]
  vlam = [ lam, 1-lam ]
  R1::Matrix{Float64} = rho * ones( pm.nb, pm.nS ) * mean( pm.g )
  R2::Matrix{Float64} = ( 1 - rho ) * ones( pm.nb, pm.nS ) * mean( pm.g )
      # Revenues
  bprime::Matrix{Float64} = ( 1 + pm.r ) * pm.bgrid * ones( 1, pm.nS ) -
            R1 - R2 + ones( pm.nb, 1 ) * sum( pm.g, 2 )'
      # Continuation debt
  bprime[bprime.<pm.bgrid[1]] = pm.bgrid[1]
      # Remove anything that has snuck below the lower bound
  wf = [ interpolate(pm.dw[i].W, (NoInterp(), BSpline(Linear())), OnGrid())
            for i in 1:2 ]
      # Loss functions
  xf = [ interpolate(pm.dw[j].x[i], BSpline(Linear()), OnGrid())
              for i in 1:pm.nS, j in 1:2 ]
      # Interpolates the grid for x for each state
  Ridx = [ idxLoc( vrho[j] * mean(pm.g), pm.dw[j].R'[:,i] )
                        for i in 1:pm.nS, j in 1:2 ]
          # Index of R in interpolation:
  V::Matrix{Float64} = ones( pm.nb ) *
                  ( [ wf[j][ i, Ridx[i] ] for i in 1:pm.nS, j in 1:2 ] *
                    vlam )'
      # Copy average value across the rows
  x1::Matrix{Float64} = ones( pm.nb ) *
              [ xf[i,1][ Ridx[i] ]::Float64 for i in 1:pm.nS ]'
  x2::Matrix{Float64} = ones( pm.nb ) *
              [ xf[i,2][ Ridx[i] ]::Float64 for i in 1:pm.nS ]'
      # The associated x value
  return V, bprime, R1, R2, x1, x2
end

"""
    bellman_operator!(pm::prsModel, V::Matrix, vOut::Matrix,
                      bOut::Matrix, R1Out::Matrix, R2Out::Matrix,
                      x1Out::Matrix, x2Out::Matrix )
Apply the Bellman operator for a given model and initial value.
##### Arguments
- `am::AutarkyModel` : Instance of `AutarkyModel`
- `V::Matrix`: Current guess for the value function
- `vOut::Matrix` : Storage for output value function
- `bOut::Matrix` : Storage for output policy function
- `R1Out::Matrix` : Storage for output policy function
- `R2Out::Matrix` : Storage for output policy function
- `x1Out::Matrix` : Storage for output leisure effort
- `x2Out::Matrix` : Storage for output leisure effort
##### Returns
None: `vOut`, `bOut`, `ROut` and `xOut` are updated in place.
"""
function bellman_operator!(pm::prsModel, V::Matrix,
      vOut::Matrix, bOut::Matrix, R1Out::Matrix, R2Out::Matrix,
      x1Out::Matrix, x2Out::Matrix )
    # simplify names, set up arrays
  r, betta, delta, P, A, g, nS, lam =
            pm.r, pm.betta, pm.delta, pm.P, pm.A, pm.g, pm.nS, pm.lam
  W = [ pm.dw[i].W for i in 1:2 ]
  Rgrid = [ pm.dw[i].R for i in 1:2 ]
  RR = [ pm.dw[i].R' for i in 1:2 ]
  nR = [ pm.dw[i].nR for i in 1:2 ]
  bgrid, nb = pm.bgrid, pm.nb
  bmin = minimum(bgrid)
  bmax = maximum(bgrid)

  s_idx = 1:size(A)[1]
  vf = interpolate(pm, V)
  wf = [ interpolate(W[i], (NoInterp(), BSpline(Linear())), OnGrid()) for i in 1:2 ]
  xf = [ interpolate(pm.dw[j].x[i], BSpline(Linear()), OnGrid())
                    for i in 1:nS, j in 1:2 ]

  # solve for RHS of Bellman equation
  for iS in s_idx, (ib, thisb) in enumerate(bgrid)

    betta_hat = betta * delta
        # Effective discount rate

    opt_lb = max( bmin, (1+r) * thisb + sum(g[iS,:]) -
                Rgrid[1][iS, nR[1]] - Rgrid[2][iS, nR[2]] )
        # Debt cannot be so low that the government revenue is above
        # its max
    opt_ub = min( bmax, (1+r) * thisb + sum(g[iS,:]) -
                Rgrid[1][iS, 1] - Rgrid[2][iS, 1] )
        # Debt cannot be so high that the government revenue is below
        # its min

    function obj(x::Vector)

      bprime, R1 = x[1], x[2]
          # Extract from the vector
      R2 = sum(g[iS,:]) + ( 1 + r ) * thisb - bprime - R1
          # Required revenue

      if ( bprime < opt_lb )
        return - 1e5 * ( bprime - opt_lb )
      end
      if ( bprime > opt_ub )
        return 1e5 * ( bprime - opt_lb )
      end
      if ( R1 < Rgrid[1][iS, 1] )
        return - 1e5 * ( R1 - Rgrid[1][iS, 1] )
      end
      if ( R1 > Rgrid[1][iS, nR[1]] )
        return 1e5 * ( R1 - Rgrid[1][iS, nR[1]] )
      end
      if ( R2 < Rgrid[2][iS, 1] )
        return - 1e5 * ( R2 - Rgrid[2][iS, 1] )
      end
      if ( R2 > Rgrid[2][iS, nR[2]] )
        return 1e5 * ( R2 - Rgrid[2][iS, nR[2]] )
      end
          # Enforce the bounds

      cont = 0.0
      for j in s_idx
        # println( "vf[bprime, j] = ", vf[bprime, j] )
        cont += vf[bprime, j] * P[iS, j]
      end

      R1idx = idxLoc( R1, RR[1][:,iS] )
      R2idx = idxLoc( R2, RR[2][:,iS] )
          # Fractional index of revenue
      loss = lam * wf[1][ iS, R1idx ] + ( 1 - lam ) * wf[2][ iS, R2idx ]
          # Period deadweight loss
      return ( 1 - betta_hat ) * loss + betta_hat * cont
    end

    println("iS, ib", "=", iS, ", ", ib )
    # println("thisT", "=", thisT )
    # println("opt_lb", "=", opt_lb )
    # println("bmax", "=", bmax, "\n" )
# println("Derp?")
    # f0 = obj([bOut[ib, iS], R1Out[ib, iS]])
    # println( "f0 = " , f0 )

    res = optimize(obj, [bOut[ib, iS], R1Out[ib, iS]], #opt_lb, opt_ub,
                    method = GradientDescent() ) # )
    x_star = res.minimum

    vOut[ib, iS] = obj(x_star)
    bOut[ib, iS] = x_star[1]
    R1Out[ib, iS] = x_star[2]
    R2Out[ib, iS] = sum(g[iS,:]) + ( 1 + r ) * thisb - sum(x_star)
    x1Out[ib, iS] = xf[iS,1][R1Out[ib, iS]]
    x1Out[ib, iS] = xf[iS,2][R2Out[ib, iS]]
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
