using JuMP
using Polygons
using Ipopt

m = Model()
@variable(m, 0 <= x <= 2 )
@variable(m, 0 <= y <= 30 )

@objective(m, Max, 5x + 3*y )
@constraint(m, 1x + 5y <= 3.0 )

print(m)

status = solve(m)
@time status = solve(m)

println("Objective value: ", getobjectivevalue(m))
println("x = ", getvalue(x))
println("y = ", getvalue(y))



rbds = [ 0.0, 1.0 ]
B = [ .1, .5, 1 ]
s = 1
k = 1.0
betta = .75
vbar = [ .8, .8 ]
    # Model parameters
ndir = 24
dirs = hcat( [ cos(i*2*pi/ndir )::Float64 for i in 1:ndir ] ,
             [ sin(i*2*pi/ndir )::Float64 for i in 1:ndir ] )
dists = [ (1 + i / (10*ndir))::Float64 for i in 1:ndir ]
V = Polygon( dirs=dirs, dists=dists)
polyPlot(V)
    # The continuation set
a = [ -.1, -.9 ] / sqrt(.81+.01)
# a = [ -1.0, -1.0 ] / sqrt(2)

pd_i( Ri, Rj ) = Ri ^ 2 + B[s] * exp( Rj - Ri )
function pd_i_d( g, Ri, Rj )
  g[1] = 2 * Ri - B[s] * exp( Rj - Ri )
  g[2] = B[s] * exp( Rj - Ri )
end
    # Period objective function
JuMP.register(:pd_i, 2, pd_i, pd_i_d)

t = Model(solver=Ipopt.IpoptSolver(print_level=0))
    # Test model
@NLparameter( t, a1 == a[1] )
@NLparameter( t, a2 == a[2] )
@variable( t, 0 <= R1 <= 1 )
@variable( t, 0 <= R2 <= 1 )
@variable( t, v1 )
@variable( t, v2 )
@NLobjective(t, Max,
  (1-betta) * a1 * pd_i( R1, R2 ) + a2 * pd_i( R2, R1 ) +
   betta * ( a1 * v1 + a2 * v2 ))

@constraint(t, R1 + R2 == k )
@constraint(t, dirs * [ v1, v2 ] .<= V.dists )
@NLconstraint( t, (1-betta)*(pd_i(R1,R2)-pd_i(rbds[1],R2)) <=
                    betta*(vbar[1]-v1))
@NLconstraint( t, (1-betta)*(pd_i(R2,R1)-pd_i(rbds[1],R1)) <=
                    betta*(vbar[2]-v2))

print(t)

status = solve(t)

println("Objective value: ", getobjectivevalue(t))
println("R1 = ", getvalue(R1))
println("R2 = ", getvalue(R2))
println("v1 = ", getvalue(v1))
println("v2 = ", getvalue(v2))

a_new = [ -.9, -.1 ] / sqrt(.81+.01)
setvalue( a1, a_new[1] )
setvalue( a2, a_new[2] )
status = solve(t)
# @time status = solve(t)

println("Objective value: ", getobjectivevalue(t))
println("R1 = ", getvalue(R1))
println("R2 = ", getvalue(R2))
println("v1 = ", getvalue(v1))
println("v2 = ", getvalue(v2))


#############################################################

using JuMP
using Ipopt

# EnableNLPResolve()

f = (x) -> exp( A * x ) - x
f1 = (x) -> A * exp( A * x ) - 1.0
f2 = (x) -> A * A * exp( A * x )
    # Period objective function
JuMP.register(:f, 1, f, f1, f2)

A = 1.0
mod = Model(solver=Ipopt.IpoptSolver(print_level=0))
@variable(mod, - Inf <= x <= Inf )
@NLobjective(mod, Min, f(x) )
status=solve(mod)
println("x = ", getvalue(x))

A = 2.0
mod = Model(solver=Ipopt.IpoptSolver(print_level=0))
@variable(mod, - Inf <= x <= Inf )
@NLobjective(mod, Min, f(x) )
status=solve(mod)
println("x = ", getvalue(x))

f = (x) -> exp( 2 *A * x ) - x
f1 = (x) -> 2 * A * exp( 2 * A * x ) - 1.0
f2 = (x) -> 4 * A * A * exp( 2 * A * x )

A = 1.0
mod = Model(solver=Ipopt.IpoptSolver(print_level=0))
@variable(mod, - Inf <= x <= Inf )
@NLobjective(mod, Min, f(x) )
status=solve(mod)
println("x = ", getvalue(x))

##### New Session #####
using JuMP
using Ipopt

function set_A_sol( A )
  f = (x) -> exp( A * x ) - x
  f1 = (x) -> A * exp( A * x ) - 1.0
  f2 = (x) -> A * A * exp( A * x )
  # Local redefinition of f
  try
    JuMP.register(:f, 1, f, f1, f2)
  catch e
    if e.msg == "Operator f has already been defined"
      ind = pop!( ReverseDiffSparse.univariate_operator_to_id, :f);
      deleteat!( ReverseDiffSparse.univariate_operators, ind);
      pop!( ReverseDiffSparse.user_univariate_operator_f, ind);
      pop!( ReverseDiffSparse.user_univariate_operator_fprime, ind);
      pop!( ReverseDiffSparse.user_univariate_operator_fprimeprime, ind);
      JuMP.register(:f, 1, f, f1, f2);
    end
  end
  mod = Model(solver=Ipopt.IpoptSolver(print_level=0))
  @variable(mod, - Inf <= x <= Inf )
  @NLobjective(mod, Min, f(x) )
  status=solve(mod)
  return getvalue(x)
end

ans1 = set_A_sol(0.5)
ans2 = set_A_sol(1.0)
ans3 = set_A_sol(2.0)
