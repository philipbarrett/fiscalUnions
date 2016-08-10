#= ctsSch.jl
Philip Barrett, pobarrett@gmail.com
22jul2016, Chicago

Computes the direction search for the continuous action problem
=#

using Optim
using JuMP
using Polygons
using Ipopt

# function dirMax( bprime::Float64, b::Float64, chi::Vector, psi::Vector,
#                   x1::Vector, x2::Vector, R1bds::Vector,a::Vector,
#                   sumg::Float64, rho::Float64, r::Float64, dir::Vector,
#                   verbose::Bool=false )
  # obj(R1) = - ( dir[1] * rho *
  #         w_eval( R1, chi[1], psi[1], x1[1], x2[1], a[1], rho, verbose ) +
  #               dir[2] * (1-rho) *
  #         w_eval( (1+r)*b + sumg - bprime - R1, chi[2], psi[2], x1[2],
  #         x2[2], a[2], 1-rho, verbose ) )::Float64
      # The objective function

function dirMax( bprime::Float64, b::Float64,
                  Rlow::Vector, Rhigh::Vector, sumg::Float64,
                  rho::Float64, r::Float64, dir::Vector,
                  apx_coeffs::Array{Array{Float64,1},1},
                  apx_N::Array{Array{Float64,1},1}, verbose::Bool=false )
  obj(R1) = - ( dir[1] * rho *
          w_eval_apx( R1, apx_coeffs[1], Rhigh[1], apx_N[1] ) +
                dir[2] * (1-rho) *
          w_eval_apx( (1+r)*b + sumg - bprime - R1,
                          apx_coeffs[2], Rhigh[2], apx_N[2] ) )
      ## Also try providing derivatives?

  if verbose
    println("Rlow = ", Rlow )
    println("Rhigh = ", Rhigh )
    println("obj(Rlow[1]) = ", obj(Rlow[1]) )
    println("obj(Rhigh[2]) = ", obj(Rhigh[2]) )
  end

  ### FIND OPTIMAL REVENUES FOR CTY 1 ###
  if (dir[1] < 0) & (dir[2] < 0)
    # Minimize joint deadweight losses
    # Interior convex optimization problem
    res = optimize( obj, Rlow[1], Rhigh[2] )
        # The result
    R1 = (res.minimum)::Float64
  elseif (dir[1] > 0) & (dir[2] > 0)
    # Maximize joint deadweight losses
    # Pick whichever endpoint gives a greater loss
    objl = obj(Rlow[1])
    obju = obj(Rhigh[1])
        # Objective at upper and lower bounds
    ibd = objl < obju ? 1 : 2
        # The index of the greater bound.
        # NB: Looks wrong way round, but remember that there is a negative
        # *in* the objective function
    R1 = (ibd==1) ? Rlow[1]::Float64 : Rhigh[2]::Float64
  elseif (dir[1] > 0) & (dir[2] < 0)
    # Minimize cty 2 losses, maximize cty 1
    R1 = Rhigh[1]::Float64
        # Upper bound on revenue
  else
    # Minimize cty 1 losses, maximize cty 2
    R1 = Rlow[1]::Float64
        # Lower bound on revenue
  end

  R2 = ((1+r)*b + sumg - bprime - R1)::Float64
      # Country 2 revenue
  W1 = rho *   w_eval_apx( R1, apx_coeffs[1], Rhigh[1], apx_N[1] )::Float64
  W2 =(1-rho) *   w_eval_apx( R2, apx_coeffs[2], Rhigh[2], apx_N[2] ) ::Float64
      # Welfare
  dist = - obj(R1)::Float64
      # Distance in search direction

  ### CHECK THAT BOUNDS ARE NOT SUPERIOR CHOICES

  return R1, R2, W1, W2, dist
end


function search_ic( bprime::Float64, b::Float64,
                    Rlow::Vector, Rhigh::Vector, sumg::Float64,
                    rho::Float64, r::Float64, betta::Float64,
                    apx_coeffs::Array{Array{Float64,1},1},
                    apx_N::Array{Array{Float64,1},1}, dirs::Matrix,
                    V::Polygon, evbar::Vector, print_level::Int=0,
                    output::ASCIIString="all" )

  ## 1. Create and register the period payoff functions
  pd_1( R ) = w_eval_apx( R, apx_coeffs[1], Rhigh[1], apx_N[1] )::Float64
  pd_2( R ) = w_eval_apx( R, apx_coeffs[2], Rhigh[2], apx_N[2] )::Float64
  pd_1_d( R ) = w_eval_apx( R, apx_coeffs[1], Rhigh[1], apx_N[1] )::Float64
  pd_2_d( R ) = w_eval_apx( R, apx_coeffs[2], Rhigh[2], apx_N[2] )::Float64
  pd_1_d2( R ) = w_eval_apx( R, apx_coeffs[1], Rhigh[1], apx_N[1] )::Float64
  pd_2_d2( R ) = w_eval_apx( R, apx_coeffs[2], Rhigh[2], apx_N[2] )::Float64
      # The period payoffs and derivatives
  try
    JuMP.register(:pd_1, 1, pd_1, pd_1_d, pd_1_d2)
  catch e
    if e.msg == "Operator pd_1 has already been defined"
      ind = pop!( ReverseDiffSparse.univariate_operator_to_id, :pd_1)
      deleteat!( ReverseDiffSparse.univariate_operators, ind)
      pop!( ReverseDiffSparse.user_univariate_operator_f, ind)
      pop!( ReverseDiffSparse.user_univariate_operator_fprime, ind)
      pop!( ReverseDiffSparse.user_univariate_operator_fprimeprime, ind)
      JuMP.register(:pd_1, 1, pd_1, pd_1_d, pd_1_d2);
    end
  end
  try
    JuMP.register(:pd_2, 1, pd_2, pd_2_d, pd_2_d2)
  catch e
    if e.msg == "Operator pd_2 has already been defined"
      ind = pop!( ReverseDiffSparse.univariate_operator_to_id, :pd_2)
      deleteat!( ReverseDiffSparse.univariate_operators, ind)
      pop!( ReverseDiffSparse.user_univariate_operator_f, ind)
      pop!( ReverseDiffSparse.user_univariate_operator_fprime, ind)
      pop!( ReverseDiffSparse.user_univariate_operator_fprimeprime, ind)
      JuMP.register(:pd_2, 1, pd_2, pd_2_d, pd_2_d2)
    end
  end
      # Register the functions

  ## 2. Create and fill the uni-directional solutions
  ##    (Used when one or both directions is positive)
  uni_sols = zeros(2,4)
      # The matrix of R1, R2, v1, v2, for the two solutions
  uni_mod = [ Model(solver=Ipopt.IpoptSolver(print_level=print_level)) for i in 1:2 ]
      # Initialize the model
  for i in 1:2
    @variable( uni_mod[i], Rlow[1] <= R1 <= Rhigh[1],
                            start = (i==1) ? Rhigh[1] : Rlow[1] )
    @variable( uni_mod[i], Rlow[2] <= R2 <= Rhigh[2],
                            start = (i==2) ? Rhigh[2] : Rlow[2] )
    @variable( uni_mod[i], v1 )
    @variable( uni_mod[i], v2 )
        # The variables
    if i == 1
      @NLobjective(uni_mod[i], Min,
          ( 1 - rho ) * ( betta * pd_2( R2 ) + ( 1 - betta ) * v2 ) )
    else
      @NLobjective(uni_mod[i], Min,
          rho * ( betta * pd_1( R1 ) + ( 1 - betta ) * v1 ) )
    end
        # Objective is to minimize the loss for R2 or R1
    @constraint(uni_mod[i], R1 + R2 + bprime == ( 1 + r ) * b + sumg )
        # Budget constraint
    @constraint(uni_mod[i], V.dirs * [ v1, v2 ] .<= V.dists )
        # Continuation set
    @NLconstraint( uni_mod[i],
          ( 1 - betta ) * rho * ( pd_1(R1) - pd_1(Rlow[1]) ) <=
                                      betta * ( evbar[1] - v1 ) )
    @NLconstraint( uni_mod[i],
          ( 1 - betta ) * ( 1 - rho ) * ( pd_2(R2) - pd_2(Rlow[2]) ) <=
                                      betta * ( evbar[2] - v2 ) )
        # Incentive compatibility constraints
    status = solve(uni_mod[i], suppress_warnings=false)
        # Solve the model
    if status != :Optimal
      R1_check = getvalue(R1)
      R2_check = getvalue(R2)
      v1_check = getvalue(v1)
      v2_check = getvalue(v2)
          # Get the values
      check = [ true, true, true, true ]
      check[1] = ( abs( R1_check + R2_check + bprime - ( 1 + r ) * b + sumg ) < 1e-10 )
      check[2] = all( V.dirs * [ v1_check, v2_check ] .<= V.dists )
      check[3] = ( ( 1 - betta ) * rho * ( pd_1(R1_check) - pd_1(Rlow[1]) ) <=
                                  betta * ( evbar[1] - v1_check ) )
      check[4] = ( ( 1 - betta ) * rho * ( pd_2(R2_check) - pd_2(Rlow[2]) ) <=
                                  betta * ( evbar[2] - v2_check ) )
      if all(check)
        warn("Possible sub-optimality in unidirectional problem #", i,
             "\nConstraints all hold ok at solution though.")
      else
        minv = minimum(V.pts,1)
        warn("Unidirectional problem #", i, " failing")
        println( "check(s) ", find(check), " failing")
        println("a = ", dirs[i,:])
        println("b = ", b)
        println("bprime = ", bprime)
        println("sumg = ", sumg)
        println("evbar = ", evbar)
        println("sum(Rhigh) = ", sum(Rhigh) )
        println("sum(Rlow) = ", sum(Rlow) )
        println("b + sumg - bprime = ", b + sumg - bprime )
        println("betta * ( evbar[1] - min(v1) ) = ",
                      betta * ( evbar[1] - minv[1] ) )
        println("( 1 - betta ) * rho * ( pd_1(Rhigh[1]) - pd_1(Rlow[1]) ) = ",
                      ( 1 - betta ) * rho * ( pd_1(Rhigh[1]) - pd_1(Rlow[1]) ) )
        println("betta * ( evbar[2] - min(v2) ) = ",
                      betta * ( evbar[2] - minv[2] ) )
        println("( 1 - betta ) * rho * ( pd_1(Rhigh[2]) - pd_1(Rlow[2]) ) = ",
                      ( 1 - betta ) * rho * ( pd_1(Rhigh[2]) - pd_1(Rlow[2]) ) )
        print(uni_mod[i])
      end
    end
        # Check that the model solves
    uni_sols[i,1] = getvalue(R1)
    uni_sols[i,2] = getvalue(R2)
    uni_sols[i,3] = getvalue(v1)
    uni_sols[i,4] = getvalue(v2)
        # Fill in the solutions
  end

  ## 3. Create the convex solution model
  convex_mod = Model(solver=Ipopt.IpoptSolver(print_level=print_level))
      # The model
  @NLparameter( convex_mod, a1 == dirs[1,1] )
  @NLparameter( convex_mod, a2 == dirs[1,2] )
      # Initialize the serach directions
  Rstart = ( Rlow + Rhigh ) / 2
  @variable( convex_mod, Rlow[1] <= R1 <= Rhigh[1],
                              start = Rstart[1] )
  @variable( convex_mod, Rlow[2] <= R2 <= Rhigh[2],
                              start = Rstart[2] )
  @variable( convex_mod, v1 )
  @variable( convex_mod, v2 )
      # The variables
  @NLobjective(convex_mod, Max,
    (1-betta) * ( a1 * rho * pd_1( R1 ) + a2 * ( 1 - rho ) * pd_2( R2 ) ) +
                    betta * ( a1 * v1 + a2 * v2 ) )
      # The objective function
  @constraint(convex_mod, R1 + R2 + bprime == ( 1 + r ) * b + sumg )
      # Budget constraint
  @constraint(convex_mod, V.dirs * [ v1, v2 ] .<= V.dists )
      # Continuation set
  @NLconstraint( convex_mod,
        ( 1 - betta ) * rho * ( pd_1(R1) - pd_1(Rlow[1]) ) <=
                                                    betta * ( evbar[1] - v1 ) )
  @NLconstraint( convex_mod,
        ( 1 - betta ) * ( 1 - rho ) * ( pd_2(R2) - pd_2(Rlow[2]) ) <=
                                                    betta * ( evbar[2] - v2 ) )
      # Incentive compatibility constraints

  ## 4. Initiate and fill the output
  ndir = size(dirs)[1]
      # Number of search directions
  # objout = zeros( ndir )
  R1out = zeros( ndir )
  R2out = zeros( ndir )
  v1out = zeros( ndir )
  v2out = zeros( ndir )
      # Initiate the output
  for i in 1:ndir

    if (dirs[i,1] < 0) & (dirs[i,2] < 0)
        # Interior maximization
      setvalue( a1, dirs[i,1] )
      setvalue( a2, dirs[i,2] )
          # Change the search direction in the JuMP problem
      status = solve(convex_mod, suppress_warnings=false )
          # Solve the model
      if status != :Optimal
        R1_check = getvalue(R1)
        R2_check = getvalue(R2)
        v1_check = getvalue(v1)
        v2_check = getvalue(v2)
            # Get the values
        check = [ false, false, false, false ]
        check[1] = ( abs( R1_check + R2_check + bprime - ( 1 + r ) * b - sumg ) < 1e-10 )
        check[2] = all( V.dirs * [ v1_check, v2_check ] - V.dists .<= 1e-08 )
        check[3] = ( ( 1 - betta ) * rho * ( pd_1(R1_check) - pd_1(Rlow[1]) ) <=
                                    betta * ( evbar[1] - v1_check ) )
        check[4] = ( ( 1 - betta ) * (1-rho) * ( pd_2(R2_check) - pd_2(Rlow[2]) ) <=
                                    betta * ( evbar[2] - v2_check ) )
        if all(check)
          # warn("Possible sub-optimality in convex problem #", i,
          #      "\nConstraints all hold at solution though.")
        else
          minv = minimum(V.pts,1)
          warn("Convex problem #", i, " failing")
          println( "\ncheck(s) ", find(!check), " failing")
          println("*********************************************")
          println("Variables:")
          println("*********************************************")
          println("b = ", b)
          println("bprime = ", bprime)
          println("sumg = ", sumg)
          println("evbar = ", evbar)
          println("[R1,R2,v1,v2] = ", (R1_check, R2_check, v1_check, v2_check) )
          println("*********************************************")
          println("Extreme values")
          println("*********************************************")
          println("Rhigh = ", Rhigh )
          println("Rlow = ", Rlow )
          println("sum(Rhigh) = ", sum(Rhigh) )
          println("sum(Rlow) = ", sum(Rlow) )
          println("( 1 - betta ) * rho * ( pd_1(Rhigh[1]) - pd_1(Rlow[1]) ) = ",
                        ( 1 - betta ) * rho * ( pd_1(Rhigh[1]) - pd_1(Rlow[1]) ) )
          println("betta * ( evbar[1] - min(v1) ) = ",
                        betta * ( evbar[1] - minv[1] ) )
          println("( 1 - betta ) * (1-rho) * ( pd_1(Rhigh[2]) - pd_1(Rlow[2]) ) = ",
                        ( 1 - betta ) * (1-rho) * ( pd_1(Rhigh[2]) - pd_1(Rlow[2]) ) )
          println("betta * ( evbar[2] - min(v2) ) = ",
                        betta * ( evbar[2] - minv[2] ) )
          println("*********************************************")
          println("Check 1")
          println("*********************************************")
          println("R1 + R2 = ", R1_check + R2_check )
          println("(1+r)*b + sumg - bprime = ", (1+r)*b + sumg - bprime )
          println("      err = ", R1_check + R2_check - ((1+r)*b + sumg - bprime))
          println("*********************************************")
          println("Check 2")
          println("*********************************************")
          v_err = V.dirs * [ v1_check, v2_check ] - V.dists
              # Should be negative
          println("Max err on v = ", maximum(v_err) )
          i_max = indmax(v_err)
          println("Err on boundary with normal = ", vec(V.dirs[i_max,:]) )
          println("*********************************************")
          println("Check 3")
          println("*********************************************")
          println("( 1 - betta ) * rho * ( pd_1(R1) - pd_1(Rlow[1]) ) = ",
                        ( 1 - betta ) * rho * ( pd_1(R1_check) - pd_1(Rlow[1]) ) )
          println("betta * ( evbar[1] - v1 ) = ",
                        betta * ( evbar[1] - v1_check ) )
          println("*********************************************")
          println("Check 4")
          println("*********************************************")
          println("( 1 - betta ) * rho * ( pd_2(R2) - pd_1(Rlow[2]) ) = ",
                        ( 1 - betta ) * rho * ( pd_2(R2_check) - pd_1(Rlow[2]) ) )
          println("betta * ( evbar[2] - v2 ) = ",
                        betta * ( evbar[2] - v2_check ) )
          # print(convex_mod)
        end
      end
          # Check that the model solves
      # objout[i] = getobjectivevalue(convex_mod)
      R1out[i] = getvalue(R1)
      R2out[i] = getvalue(R2)
      v1out[i] = getvalue(v1)
      v2out[i] = getvalue(v2)

          # Fill in the solutions
    elseif (dirs[i,1] >= 0) & (dirs[i,2] < 0)
        # "Punishment" for country 1
      R1out[i] = uni_sols[1,1]
      R2out[i] = uni_sols[1,2]
      v1out[i] = uni_sols[1,3]
      v2out[i] = uni_sols[1,4]
      # objout[i] = (1-betta) * ( dirs[i,1] * rho * pd_1( R1out[i] ) +
      #                         dirs[i,2] * ( 1 - rho ) * pd_2( R2out[i] ) ) +
      #                 betta * ( dirs[i,1] * v1out[i] + dirs[i,2] * v2out[i] )
    elseif (dirs[i,1] < 0) & (dirs[i,2] >= 0)
        # "Punishment" for country 2
      R1out[i] = uni_sols[2,1]
      R2out[i] = uni_sols[2,2]
      v1out[i] = uni_sols[2,3]
      v2out[i] = uni_sols[2,4]
      # objout[i] =
    else
        # "Punishing" both
      Wtest = [ (1-betta) * ( dirs[i,1] * rho * pd_1( uni_sols[j,1] ) +
                              dirs[i,2] * ( 1 - rho ) *
                                    pd_2( uni_sols[j,2] ) ) + betta *
                        ( dirs[i,1] * uni_sols[j,3] +
                          dirs[i,2] * uni_sols[j,4] ) for j in 1:2 ]
          # Try out the two extremes for W
      idx = indmax(Wtest)
          # Fid the best of the two in the current search direction
      R1out[i] = uni_sols[idx,1]
      R2out[i] = uni_sols[idx,2]
      v1out[i] = uni_sols[idx,3]
      v2out[i] = uni_sols[idx,4]
      # objout[i] = Wtest[idx]
    end
  end

  objout = [ ( (1-betta) * ( dirs[i,1] * rho * pd_1( R1out[i] ) +
                      dirs[i,2] * ( 1 - rho ) * pd_2( R2out[i] ) ) +
                    betta * ( dirs[i,1] * v1out[i] +
                              dirs[i,2] * v2out[i] ) )::Float64
                                  for i in 1:ndir ]
  W1out = [ ((1-betta) * rho * pd_1( R1out[i] ) + betta * v1out[i])::Float64
                    for i in 1:ndir ]
  W2out = [ ((1-betta) * rho * pd_2( R2out[i] ) + betta * v2out[i])::Float64
                    for i in 1:ndir ]
  if output == "outer"
    return Polygon( dirs=dirs, dists=objout )
  elseif output == "inner"
    return Polygon( pts = hcat( W1out, W2out ) )
  elseif output == "dists"
    return objout
  elseif output == "pts"
    return hcat( W1out, W2out )
  else
    return objout, W1out, W2out, R1out, R2out, v1out, v2out
  end
end
