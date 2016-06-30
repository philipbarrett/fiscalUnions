# xx = []

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

plot( layer( x=ll[:,1], y=ll[:,2], Geom.point),
      layer( x=qq[:,1], y=qq[:,2], Geom.line ) )
plot!( qq[:,1], qq[:,2], lw=2, col='r')
polyPlot(qq)

scatter(kk[:,1], kk[:,2])
plot!( pp[:,1], pp[:,2], lw=2, col='r')
oo = gScan(pp)
plot!( oo[:,1], oo[:,2], lw=2, col='g')

plot( layer( x=taxes_jt[:,1], y=taxes_jt[:,2], Geom.point,
                  Theme(default_color=color( "red") ) ),
      layer( x=taxes_jt[:,1], y=tax_sum/(1+chi), Geom.point,
                        Theme(default_color=color( "blue") ) ) )
