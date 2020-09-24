#include <Rcpp.h>
using namespace Rcpp;

//' get_nb_cpp
//' Get negative binomial probability density for double 'count', given
//' size and mu parameters.
//' @param mu NumericVector
//' @param count double
//' @param size NumericVector
//' @export
// [[Rcpp::export]]
NumericVector get_nb_cpp(NumericVector mu, double count, NumericVector size) {
  int muSize = mu.length();
  NumericVector nb_ll(muSize);
  for (int i = 0; i < muSize; i++) {
    nb_ll[i] = R::dnbinom_mu(count, size[i], mu[i], 1);
  }
  // Rcpp::dnbinom_mu( x, size, mu, log = false )
  return nb_ll;
}

//' get_pois_countVec
//' Get poisson probability density for multiple counts.
//' @param count_vector NumericVector
//' @param mu double
//' @export
// [[Rcpp::export]]
NumericVector get_pois_countVec(NumericVector count_vector, double mu) {
  // if any are negative, stop.
  // if any counts aren't ints, stop.
  NumericVector poiss_ll = Rcpp::dpois(count_vector, mu, true);
  // Rcpp::as<type>()
  // Rcout << "val of count input:" << count_vector[1] << "\n";
  // Rcout << "val of mu:" << mu << "\n";
  // Rcout << "val of countVec poissll is:" << poiss_ll[1] << "\n";
  return poiss_ll;
}

//' get_pois_muVec
//' @param count doubl
//' @param mu_vector NumericVector
//' @export
// [[Rcpp::export]]
NumericVector get_pois_muVec(double count, NumericVector mu_vector) {
  int muSize = mu_vector.length();
  NumericVector poiss_ll(muSize);
  for (int i = 0; i < muSize; i++) {
    poiss_ll[i] = R::dpois(count, mu_vector[i], 1);
  }
  return poiss_ll;
}

//' get_norm_countVec
//' @param count_vector NumericVector
//' @param mu double
//' @param var double
//' @export
// [[Rcpp::export]]
NumericVector get_norm_countVec(NumericVector count_vector, double mu,
                                double var) {
  NumericVector norm_ll = Rcpp::dnorm(count_vector, mu, var);
  return norm_ll;
}
  