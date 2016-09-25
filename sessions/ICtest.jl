d = load("/home/philip/Dropbox/data/2016/fiscalUnions/eg.jld")
eqin = d["eqin"]
eqout = d["eqout"]
cfg = d["cfg"]
pdin = d["pdin"]


iS, ib, ibprime = 1, 1, 3

outer=false
    # Taking the inner continuations
ibprimeidx = [ find( dw.bprimeposs[iS,ib] )
                for iS in 1:cfg.nS, ib in 1:cfg.nb ]

econt = expectedVals( eqin, cfg.P, dirs, outer )
    # Continuations
vdev = devCont( cfg.dw, cfg.P, cfg.betta, cfg.rho, cfg.A, cfg.chi,
                cfg.psi, true )
    # Return the deviating *continuation* only
evdev = cfg.P * vdev

testin_pts_v = [ search_ic( cfg.dw.bgrid[ibprimeidx[iS,ib][ibprime]],
                          dw.bgrid[ib], [ dw.Rlow[iS,ib,k][ibprime] for k in 1:2 ],
                          vec( dw.Rhigh[iS, : ] ), cfg.gSum[iS], cfg.rho, cfg.r,
                          cfg.betta, vec(dw.apx_coeffs[iS,:]), vec(dw.apx_N[iS,:]),
                          dirs, econt[iS,ib], vec(evdev[iS,:]), 0, "pts" )::Matrix
                          for ibprime in 1:dw.nposs[iS,ib] ]
    # Inner approx
testout_dists_v = [ search_ic( cfg.dw.bgrid[ibprimeidx[iS,ib][ibprime]],
                          dw.bgrid[ib], [ dw.Rlow[iS,ib,k][ibprime] for k in 1:2 ],
                          vec( dw.Rhigh[iS, : ] ), cfg.gSum[iS], cfg.rho, cfg.r,
                          cfg.betta, vec(dw.apx_coeffs[iS,:]), vec(dw.apx_N[iS,:]),
                          dirs, econt[iS,ib], vec(evdev[iS,:]), 0, "dists" )::Vector
                          for ibprime in 1:dw.nposs[iS,ib] ]
    # Outer approx

testin_v_poly = [ Polygon( pts = z )::Polygon for z in testin_pts_v ]
testout_v_poly = [ Polygon( dirs=dirs, dists=d )::Polygon for d in testout_dists_v ]

testin_u_poly = Polygon( pts=vcat(testin_pts_v...) )
testout_u_poly = Polygon( dirs=dirs, dists=vec(maximum(hcat(testout_dists_v...),2)) )

# vv = valsUpdate( conts::Array{Polygon, 2},
#                       pdPayoffs::Array{Array{Polygon,1},2},
#                       P::Matrix{Float64}, dirs::Matrix{Float64},
#                       dw::DWLC2, gSum::Vector, rho::Float64, r::Float64,
#                       betta::Float64, outer::Bool=true, evbar=false )
#   nS, nb = size(conts)
## HERE: Do the valsUpdate calculation.  Then compare again to the examples
