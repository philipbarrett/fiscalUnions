xx = []

typeof(xx)

mm = Matrix[]

isempty(mm)

ff = [ 1 2 ; 3 4 ]

isempty(ff)
isempty(xx)

qq = zeros(Float64, 3, 2)

yy = [1; 5 ]
zz = ff \ yy

function testFun( )
  return [0,0,0]
end

dd=testFun()
typeof(dd)

scatter(ll[:,1], ll[:,2])
plot!( qq[:,1], qq[:,2], lw=2, col='r')
polyPlot(qq)

scatter(kk[:,1], kk[:,2])
plot!( pp[:,1], pp[:,2], lw=2, col='r')
oo = gScan(pp)
plot!( oo[:,1], oo[:,2], lw=2, col='g')
