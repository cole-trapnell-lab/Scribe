// Copyright 2015 <Jeremy Yee> <jeremyyee@outlook.com.au>
// Radius search
////////////////////////////////////////////////////////////////////////////////

#include <vector>
#include <RcppArmadillo.h>
#include <string>
#include "flann/flann.hpp"
#include "../inst/include/radius_search.h"

// //[[Rcpp::export]]
Rcpp::NumericVector RadiusSearch(Rcpp::NumericMatrix ref_query_,
                     // Rcpp::NumericMatrix ref_,
                     Rcpp::NumericVector radius,
                     // int k,
                     std::string build,
                     int cores,
                     int checks) {
  float order = std::numeric_limits<float>::infinity(); //infinity order for the minkowski norm 

  const std::size_t n_dim = ref_query_.ncol();
  const std::size_t n_ref_query = ref_query_.nrow();
  // const std::size_t n_ref = ref_.nrow();
  // Column major to row major
  arma::mat ref_query(n_dim, n_ref_query);
  {
    arma::mat temp_q(ref_query_.begin(), n_ref_query, n_dim, false);
    ref_query = arma::trans(temp_q);
  }
  flann::Matrix<double> rq_flann(ref_query.memptr(), n_ref_query, n_dim);
  // arma::mat ref(n_dim, n_ref);
  // {
  //   arma::mat temp_r(ref_.begin(), n_ref, n_dim, false);
  //   ref = arma::trans(temp_r);
  // }
  // flann::Matrix<double> ref_flann(ref.memptr(), n_ref, n_dim);
  // Setting the flann index params
  flann::IndexParams params;
  if (build == "kdtree") {
    params = flann::KDTreeSingleIndexParams(1);
  } else if (build == "kmeans") {
    params = flann::KMeansIndexParams(2, 10, flann::FLANN_CENTERS_RANDOM, 0.2);
  } else if (build == "linear") {
    params = flann::LinearIndexParams();
  }
  // Perform the radius search
  flann::Index<flann::MaxDistance<double> > index(rq_flann, params); // change to infinity norm (MinkowskiDistance)
  index.buildIndex();
  std::vector< std::vector<int> > indices_flann; //(max_neighbour)
  std::vector< std::vector<double> > dists_flann; // (max_neighbour) (n_ref_query, std::vector<double>)
  flann::SearchParams search_params;
  search_params.cores = cores;
  search_params.checks = checks;
  // search_params.max_neighbors = max_neighbour;

  // return only the k-th distance: 
  Rcpp::NumericVector points_in_radius(n_ref_query); 
  for(int i = 0; i < n_ref_query; i ++)
  {
    Rcpp::NumericMatrix ref_query_subset_ = ref_query_(Rcpp::Range(i, i), Rcpp::_);
    arma::mat ref_query_subset(n_dim, 1);
    {
    arma::mat temp_q(ref_query_subset_.begin(), 1, n_dim, false);
    ref_query_subset = arma::trans(temp_q);
    }
    
    flann::Matrix<double> rq_flann_subset(ref_query_subset.memptr(), 1, n_dim); 
    index.radiusSearch(rq_flann_subset, indices_flann, dists_flann, radius[i], // subset of rq_flann matrix 
                     search_params);

    Rcpp::Rcout << "Current i is" << i << "Current length of radius is " << indices_flann.size() << std::endl; 

    points_in_radius[i] = indices_flann[i].size(); 
  }
  return points_in_radius;   
  // return Rcpp::List::create(Rcpp::Named("indices") = indices_flann,
  //                           Rcpp::Named("distances") = dists_flann);
}

// // Export to R
// RcppExport SEXP rflann_RadiusSearch(SEXP query_SEXP,
//                                     SEXP ref_SEXP,
//                                     SEXP radiusSEXP,
//                                     SEXP max_neighbourSEXP,
//                                     SEXP buildSEXP,
//                                     SEXP coresSEXP,
//                                     SEXP checksSEXP) {
// BEGIN_RCPP
//     Rcpp::RObject __result;
//     Rcpp::RNGScope __rngScope;
//     Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type
//         query_(query_SEXP);
//     Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type
//         ref_(ref_SEXP);
//     Rcpp::traits::input_parameter< double >::type
//         radius(radiusSEXP);
//     Rcpp::traits::input_parameter< int >::type
//         max_neighbour(max_neighbourSEXP);
//     Rcpp::traits::input_parameter< std::string >::type
//         build(buildSEXP);
//     Rcpp::traits::input_parameter< int >::type
//         cores(coresSEXP);
//     Rcpp::traits::input_parameter< int >::type
//         checks(checksSEXP);
//     __result = Rcpp::wrap(RadiusSearch(query_, ref_, radius,
//                                        max_neighbour, build, cores, checks));
//     return __result;
// END_RCPP
// }
