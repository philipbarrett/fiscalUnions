using JuMP
using Polygons

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
a = [ -1.0, -1.0 ] / sqrt(2)
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

pd_i( Ri, Rj ) = Ri ^ 2 + B[s] * exp( Rj - Ri )
function pd_i_d( g, Ri, Rj )
  g[1] = 2 * Ri - B[s] * exp( Rj - Ri )
  g[2] = B[s] * exp( Rj - Ri )
end
    # Period objective function
JuMP.register(:pd_i, 2, pd_i, pd_i_d)

t = Model(solver=Ipopt.IpoptSolver(print_level=3))
    # Test model
@variable( t, 0 <= R1 <= 1 )
@variable( t, 0 <= R2 <= 1 )
@variable( t, v1 )
@variable( t, v2 )
@NLobjective(t, Max,
  (1-betta) * a[1] * pd_i( R1, R2 ) + a[2] * pd_i( R2, R1 ) +
   betta * ( a[1] * v1 + a[2] * v2 ))

## NEED TO REGISTER MY OWN FUNCTION

@constraint(t, R1 + R2 == k )
@constraint(t, dirs * [ v1, v2 ] .<= V.dists )
@NLconstraint( t, (1-betta)*(pd_i(R1,R2)-pd_i(rbds[1],R2)) <=
                    betta*(vbar[1]-v1))
@NLconstraint( t, (1-betta)*(pd_i(R2,R1)-pd_i(rbds[1],R1)) <=
                    betta*(vbar[2]-v2))

print(t)

status = solve(t)
# @time status = solve(t)

println("Objective value: ", getobjectivevalue(t))
println("R1 = ", getvalue(R1))
println("R2 = ", getvalue(R2))
println("v1 = ", getvalue(v1))
println("v2 = ", getvalue(v2))


function jumpTest( rbd::Vector{Float64}, B::Vector{Float64},
                    s::Int, a::Vector{Float64}, k::Float64,
                    betta::Float64, vbar::Vector{Float64},
                    V::Polygon )
  t = Model(solver=Ipopt.IpoptSolver(print_level=0) )
      # Test model
  pd_i( Ri, Rj ) = Ri ^ 2 + B[s] * ( Rj - Ri ) ^ 3
      # Period objective function

  @variable( t, rbds[1] <= R1 <= rbds[2] )
  @variable( t, rbds[1] <= R2 <= rbds[2] )
  @variable( t, v1 <= 100 )
  @variable( t, v2 <= 100)
  @NLobjective(t, Min,
    ( ( (1-betta) * [ pd_i( R1, R2 ) pd_i( R2, R1 ) ] +
        betta * [ v1 v2 ] )
         * a )[1] )
  # @constraint(t, R1 + R2 == k )
  # @constraint(t, V.dirs * [ v1, v2 ] .<= V.dists )
  @NLconstraint( t, (1-betta)*(pd_i(R1,s)-pd_i(rbds[1],s)) <=
                      betta*(vbar[1]-v1))
  @NLconstraint( t, (1-betta)*(pd_i(R2,s)-pd_i(rbds[1],s)) <=
                      betta*(vbar[2]-v2))

  status = solve(t)

  obj = getobjectivevalue(t)
  R1 = getvalue(R1)
  R2 = getvalue(R2)
  W1 = (1-betta) * pd_i( R1, R2 ) + betta * v1
  W2 = (1-betta) * pd_i( R2, R1 ) + betta * v2
  v1 = getvalue(v1)
  v2 = getvalue(v2)

  return obj, W1, W2, v1, v2, R1, R2
end

@time x1 = jumpTest( [.1, 1], [ .5, 1 ], 2, [ -1, -1 ] / sqrt(2),
                      1.0, .5, [ .6, .6 ], V )
