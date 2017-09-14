// Copyright 2015 <Jeremy Yee> <jeremyyee@outlook.com.au>
// Radius search
////////////////////////////////////////////////////////////////////////////////

// #include <vector>
// #include <RcppArmadillo.h>
// #include <string>
//#include "flann/flann.hpp"

Rcpp::NumericVector Neighbour(Rcpp::NumericMatrix ref_query_,
                     // Rcpp::NumericMatrix ref_,
                     int k,
                     std::string build,
                     int cores,
                     int checks);