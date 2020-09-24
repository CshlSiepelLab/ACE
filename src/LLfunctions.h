#ifndef LLFUNCTIONS_H
#define LLFUNCTIONS_H
#include <Rcpp.h>
using namespace Rcpp;
NumericVector get_nb_cpp(NumericVector mu, double count, NumericVector size);
NumericVector get_pois_countVec(NumericVector count_vector, double mu);
NumericVector get_pois_muVec(double count, NumericVector mu_vector);
NumericVector get_norm_countVec(NumericVector count_vector, double mu,
                                double var);
#endif