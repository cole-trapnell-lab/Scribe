// Copyright 2015 <Jeremy Yee> <jeremyyee@outlook.com.au>
// Nearest Neighbours
////////////////////////////////////////////////////////////////////////////////

#include <RcppArmadillo.h>
#include <string>
#include "flann/flann.hpp"
#include "../inst/include/neighbour.h"

// //[[Rcpp::export]]
Rcpp::NumericVector Neighbour(Rcpp::NumericMatrix ref_query_,
                     // Rcpp::NumericMatrix ref_,
                     int k,
                     std::string build,
                     int cores,
                     int checks) {

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

  Rcpp::Rcout << "pass setting params.." << std::endl;

  // Finding the nearest neighbours
  flann::Index<flann::L2<double> > index(rq_flann, params); // change to infinity norm (MinkowskiDistance) 
  Rcpp::Rcout << "before knnSearch .." << std::endl;
  index.buildIndex();
  flann::Matrix<int> indices_flann(new int[n_ref_query * k], n_ref_query, k);
  flann::Matrix<double> dists_flann(new double[n_ref_query * k], n_ref_query, k);
  flann::SearchParams search_params;
  search_params.cores = cores;
  search_params.checks = 100000; //checks

  Rcpp::Rcout << "before knnSearch 1.1.." << std::endl;
  index.knnSearch(rq_flann, indices_flann, dists_flann, k, search_params);
  Rcpp::Rcout << "before knnSearch.." << std::endl;
  // arma::imat indices(indices_flann.ptr(), k, n_ref_query, true);
  arma::mat dists(dists_flann.ptr(), k, n_ref_query, true);
  Rcpp::Rcout << "before deleting stuff.." << std::endl;
  delete[] indices_flann.ptr();
  delete[] dists_flann.ptr();

  Rcpp::Rcout << "before assigning distances value.." << std::endl;
  // return only the k-th distance: 
  Rcpp::NumericVector distances(n_ref_query); 
  for(int i = 0; i < n_ref_query; i ++)
  {
    distances[i] = sqrt(dists[(i + 1) * k - 1]); 
  }

  Rcpp::Rcout << "return distances.." << std::endl;

  return distances; 
}

// // Export to R
// RcppExport SEXP rflann_Neighbour(SEXP query_SEXP,
//                                  SEXP ref_SEXP,
//                                  SEXP kSEXP,
//                                  SEXP buildSEXP,
//                                  SEXP coresSEXP,
//                                  SEXP checksSEXP) {
// BEGIN_RCPP
//     Rcpp::RObject __result;
//     Rcpp::RNGScope __rngScope;
//     Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type
//         query_(query_SEXP);
//     Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type
//         ref_(ref_SEXP);
//     Rcpp::traits::input_parameter< int >::type
//         k(kSEXP);
//     Rcpp::traits::input_parameter< std::string >::type
//         build(buildSEXP);
//     Rcpp::traits::input_parameter< int >::type
//         cores(coresSEXP);
//     Rcpp::traits::input_parameter< int >::type
//         checks(checksSEXP);
//     __result = Rcpp::wrap(Neighbour(query_, ref_, k, build, cores, checks));
//     return __result;
// END_RCPP
// }
