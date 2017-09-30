#include <Rcpp.h>
using namespace Rcpp;

NumericVector get_NN_2Set_cpp(const NumericMatrix& data, const NumericMatrix& query_data, int& D, int& ND, int& NQ, int& K, double& EPS,
                       int& SEARCHTYPE, int& USEBDTREE, double& SQRAD, NumericVector& distances);
Rcpp::List get_NN_2Set(SEXP R_data, SEXP R_query_data, SEXP R_d, SEXP R_nd, SEXP R_nq, SEXP R_k, SEXP R_error_bound,
                       SEXP R_searchtype, SEXP R_usebdtree, SEXP R_sqRad, SEXP R_distances,
                       SEXP R_verbose);
NumericVector get_points_in_radius_cpp(const NumericMatrix& data, const NumericMatrix& query_data, int& D, int& ND, int& NQ, int& K, double& EPS,
                       int& USEBDTREE, const NumericVector& SQRAD);