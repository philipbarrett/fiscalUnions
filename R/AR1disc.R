#############################################################################
# AR1disc.R
#
# Functions to discretize a univariate AR(1) process
# Philip Barrett, Chicago
# Created: 10jun2016
#############################################################################


ar1.disc <- function( n, rho, eps, ext=FALSE ){
  # Creates the discretized version of the
  # auto-regressive process
  
  ### Approximate the AR(1) using Adda & Cooper's (2003) method
  
  sig.z <- eps / sqrt( 1 - rho ^ 2 )
  # The variance of the stationary distribution
  if( ext ){
    pciles <- c( 0, sapply( 1:(n-1), 
                            function(i) .5 - .5 * cos( (2*i-1) / (2*(n-1)) * pi ) ), 1 )
    pciles[2] <- pciles[2] / 4
    pciles[length(pciles)-1] <- 1 - pciles[2]
    # Yet more stretching of the tail points
  }else{
    pciles <- 0:n / n
  } # If need an "extended" sample, overwieght the ends using the Chebychev nodes
  
  m <- sig.z * qnorm( pciles, 0 )
  # The equi-probable ranges of z
  #   v.z <- - sig.z * n *( dnorm( qnorm( pciles[-1] ) ) - dnorm( qnorm( pciles[-(n+1)] ) ) )
  v.z <- qnorm( .5 * ( pnorm( m[-1], 0, sig.z ) + 
                         pnorm( m[-(n+1)], 0, sig.z ) ) , 0, sig.z )
  m.p <- matrix( NA, nrow=n, ncol=n )
  # The transition probabilities
  
  integrand <- function(j){
    # Returns the function in the integrand for element i,j of the transition matrix
    return( function(z) ( pnorm( m[j+1], rho * z, eps ) - pnorm( m[j], rho * z, eps ) ) 
            * dnorm( z, 0, sig.z ) )
  }
  
  p.i.j <- function( i, j ){
    # Computes the i,j-th element of the transition matrix
    num <- integrate( integrand(j), m[i], m[i+1] )
    dom <- pnorm( m[i+1], 0, sig.z) - pnorm( m[i], 0, sig.z) 
    out <- num$value / dom
    # Create the integral
    if( num$message == 'OK' )
      return( out )
    #       return( n * num$value )
    stop('Numerical integration failed')
    # Check that the integration is ok
    
    #     out <- integrate( integrand(j), m[i], m[i+1] )
    #     # Create the integral
    #     if( out$message == 'OK' )
    #       return( n * out$value )
    #     stop('Numerical integration failed')
    #     # Cheack that the integration is ok
  }
  
  m.p <- sapply( 1:n, function(j) sapply( 1:n, function(i) p.i.j(i,j) ) )
  # The matrix of transition probabilities
  
  rowsums <- apply( m.p, 1, sum )
  disc.err <- max( abs( 1 - rowsums ) )
  n.disc.err <- sum( abs( 1 - rowsums ) > 1e-12 )
  if( n.disc.err != 0 ){
    message( 'Warning: AR(1) discretization has maximum error of ', 
             round( disc.err, 7 ), '  on ', n.disc.err, ' rows.' )
  }
  
  m.p <- m.p / apply( m.p, 1, sum )
  # Only need this for epsilon super-small
  
  return( list( z=v.z, p=m.p, m=m ) )
}

ar.test <- function( rho, eps, n.disc=41, n.sim=10000, ext=FALSE, plot=FALSE ){
  # Tests the discretization of the AR(1) process for theta
  
  ar1 <- ar1.disc( n.disc, rho, eps, ext )
  # The discretized MC
  set.seed( 42 )
  # Set the seed
  mc.theta <- new( "markovchain", states = as.character( ar1$z ),
                   transitionMatrix = ar1$p, name = "theta" )
  # Create the Markov chain objects
  sim.theta <- as.numeric( rmarkovchain( n=n.sim, object=mc.theta ) )
  # The simulations
  ar1.lm <- lm( sim.theta[-1] ~ sim.theta[-n.sim] )
  # The autoregression
  if( plot==TRUE ){
    plot( c(1,200), range(ar1$z), type='n', xlab='t', ylab='y(t)' )
    lines( 1:200, sim.theta[ (n.sim-199):n.sim ], lwd=2 )
    abline( h=ar1$z, lwd=.5, lty=2 )
  }
  
  return( rbind( sd=c( exact=eps/sqrt(1-rho^2), sim=sd(sim.theta) ), 
                 rho=c( exact=rho, sim=ar1.lm$coef[2] ) ) )
}