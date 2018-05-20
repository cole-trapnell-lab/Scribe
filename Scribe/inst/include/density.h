#include <Rcpp.h>
using namespace Rcpp;

NumericMatrix kde_cpp(NumericMatrix data, int k = 1, int b = 2, int pdf = 1, int density_sample_type = 1);