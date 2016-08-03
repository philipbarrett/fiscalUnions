#= wApx.jl
Philip Barrett, pobarrett@gmail.com
02aug2016, Chicago

Code to approximate the deadweight loss function
=#

function w_eval_apx( R::Float64, coeff::Vector{Float64}, Rlim::Float64, NN=NaN )
  NN = (isnan(NN[1])) ? (2:length(coeff)) : NN
  if R < 0.0
    return coeff[1]
  else
    out = coeff[1]
        # Intercept
    counter = 1
    for n in NN
      counter += 1
      out += coeff[counter] * ( Rlim ^ ( 1 / n ) * ( 1 - R / ( n * Rlim ) ) -
                              ( max( Rlim - R, 0 ) ) ^ ( 1 / n ) )
    end
    return out
  end
end

function w_eval_apx_d( R::Float64, coeff::Vector{Float64}, Rlim::Float64, NN=NaN )
  NN = (isnan(NN[1])) ? (2:length(coeff)) : NN
  if R < 0.0
    return 0.0
  else
    out = 0.0
        # Intercept
    counter = 1
    for n in NN
      counter += 1
      out += coeff[counter] * ( - Rlim ^ ( 1 / n - 1 ) +
                                ( max( Rlim - R, 0 ) ) ^ ( 1 / n - 1 ) ) / n
    end
    return out
  end
end

function w_eval_apx_d2( R::Float64, coeff::Vector{Float64}, Rlim::Float64, NN=NaN )
  NN = (isnan(NN[1])) ? (2:length(coeff)) : NN
  if R < 0.0
    return 0.0
  else
    out = 0.0
        # Intercept
    counter = 1
    for n in NN
      counter += 1
      out += coeff[counter] * ( 1 - 1 / n ) / n *
                            ( max( Rlim - R, 0 ) ) ^ ( 1 / n - 2 )
    end
    return out
  end
end

function w_eval_fit( Rhigh::Float64, chi::Float64, psi::Float64,
                  xlow::Float64, xhigh::Float64, a::Float64, rho::Float64,
                  nR::Int = 40, tol::Float64=1e-05, max_n::Int=5 )
  RR = linspace( 0, Rhigh, nR )
      # The grid on which to assess the fit
  NN = [ 2.0, 3.0, 1.5 ]
      # The initial vector of orders for approximation
  X = [ (Rhigh ^ ( 1 / NN[j] ) * ( 1 - RR[i] / ( NN[j] * Rhigh ) ) -
          ( max( Rhigh - RR[i], 0 ) ) ^ ( 1 / NN[j] ))::Float64 for i in 1:nR, j in 1:3 ]
      # The data for the regression
  Y = [ w_eval( RR[i], chi, psi, xlow, xhigh, a, rho )::Float64 for i in 1:nR ]
      # Create the data for the function
  coeff = [ Y[1] ; (inv( X' * X ) * ( X' * (Y-Y[1]) ))]::Vector{Float64}
      # Regress the function on the data
  res = ([ w_eval_apx( RR[i], coeff, Rhigh, NN )::Float64 for i in 1:nR ]
          - Y)::Vector{Float64}
      # Residuals
  err = maximum(abs(res))
      # Approximation error
  counter = 2
  while ( err > tol ) & (counter < max_n)
    NN = [ NN; [ 2 ^ counter, 1 + 1 / ( 1 + counter ) ] ]
    X = [ (Rhigh ^ ( 1 / NN[j] ) * ( 1 - RR[i] / ( NN[j] * Rhigh ) ) -
            ( max( Rhigh - RR[i], 0 ) ) ^ ( 1 / NN[j] ))::Float64
            for i in 1:nR, j in (2*counter):(2*counter+1) ]
    coeff = [ coeff ; (inv( X' * X ) * ( X' * res ) )::Vector{Float64} ]
    res = [ w_eval_apx( RR[i], coeff, Rhigh, NN )::Float64 for i in 1:nR ] - Y
    err = maximum(abs(res))
        # Repeat fitting until error small enough
    counter += 1
  end
  return NN, coeff, err
end
