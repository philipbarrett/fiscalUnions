#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
arma::mat lyapunovEq( arma::mat phi, arma::mat sigma, 
                      int maxit = 1000, double tol = 1e-08, bool printit = false ) {
// Solves the Lyapunov matrix equation lambda = phi * lambda * phi.t + sigma
// via the iteration lambda = sum_{k=0}^inf phi^k * lambda * phi.t^k
  mat lam = sigma ;
      // Initalize the output
  mat term = sigma ;
      // The current term in the summation
  mat phiT = phi.t() ;
      // Store this transpose in memory
  int i = 0 ;
  double diff = 2 * tol ;
      // Loop counters
  while( i < maxit & diff > tol ){
    term = phi * term * phiT ;
        // Increment the summand
    lam = lam + term ;
        // Add to output
    diff = norm( term, "inf" ) ;
        // The difference measure - the maximum norm
    i++ ;
        // Increment counter
  }
  
  if(printit){
    Rcout << "iter = " << i << std::endl ;
    Rcout << "diff = " << diff << std::endl ;
  }
  return lam ;
}
