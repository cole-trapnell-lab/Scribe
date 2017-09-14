// Copyright 2015 <Jeremy Yee> <jeremyyee@outlook.com.au>
// Radius search
////////////////////////////////////////////////////////////////////////////////
// 
// #include <vector>
// #include <RcppArmadillo.h>
// #include <string>

Rcpp::NumericVector RadiusSearch(Rcpp::NumericMatrix ref_query_,
                        // Rcpp::NumericMatrix ref_,
                        Rcpp::NumericVector radius,
                        // int max_neighbour,
                        std::string build,
                        int cores,
                        int checks); 