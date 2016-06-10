#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
arma::mat lyapunovEq( arma::mat& phi, arma::mat& sigma, 
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

// [[Rcpp::export]]
int find_nearest( arma::rowvec& X, arma::mat& Z, int nZ ){
// Find the nearest point to X in Z according to the Eucildian norm
  mat diffs = ones( nZ ) * X - Z ;
      // The matrix of differences
  vec dists = sum( diffs % diffs, 1 ) ;
      // The vector of distances
  uvec idx = find( dists == min(dists) ) ;
      // The vector of minimum indices
  return idx[0] ;
      // Return the location of the minimum distance
}

// [[Rcpp::export]]
arma::rowvec p_hat( arma::rowvec& eX, arma::mat& eps, arma::mat& Z ){
// Computes the vector of probabilities that X' is nearest to the corresponding
// element of Z

  int nZ = Z.n_rows ;
      // Number of points in Z
  int n_eps = eps.n_rows ;
  double inc_eps = 1 / double(eps.n_rows) ;
      // Number of shocks
  rowvec p = zeros<rowvec>(nZ) ;
      // Initialize the output
  int idx = 0 ;
      // Location of nearest point
  rowvec Xprime ;
      // value of X prime
  for( int i = 0 ; i < n_eps ; i++ ){
    Xprime = eX + eps.row(i) ;
    idx = find_nearest( Xprime, Z, nZ ) ;
        // Location of closest point to X'
    p[idx] += inc_eps ;
        // Increment the probability vector
//        Rcout << "i = " << i << std::endl ;
//        Rcout << "idx = " << idx << std::endl ;
//        Rcout << "p: \n" << p << std::endl << std::endl ;
  }
  return p ;
}

// [[Rcpp::export]]
arma::mat p_trans( arma::mat& eX, arma::mat& eps, arma::mat& Z ){
// Computes by Monte Carlo integration the 
  
}