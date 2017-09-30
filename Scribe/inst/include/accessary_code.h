#include <Rcpp.h>
using namespace Rcpp;

double di_single_run_cpp(NumericMatrix& x, NumericMatrix& y, int n = 5); 
double di_single_run_conditioned_cpp(NumericMatrix x, NumericMatrix y, NumericMatrix& z, int n = 5);
double rdi_many_runs_cpp(NumericMatrix& x, NumericMatrix& y);
double rdi_single_run_cpp(NumericMatrix& x, NumericMatrix& y, int d = 1);
double rdi_single_run_conditioned_cpp(NumericMatrix& x, NumericMatrix& y, NumericMatrix& z, NumericVector& z_delays, int d = 1);
double lmi_multiple_run(NumericMatrix& x, NumericMatrix& y, int d = 1, IntegerVector run_vec = 0); 
List calculate_rdi_cpp(NumericMatrix& expr_data, IntegerVector delays = 1, IntegerMatrix& super_graph, IntegerVector& turning_points = 0, int method = 1); // NULL
NumericMatrix calculate_conditioned_rdi_cpp(NumericMatrix& expr_data, IntegerMatrix& super_graph, 
                                            NumericMatrix& max_rdi_value, IntegerMatrix& max_rdi_delays, int k = 5);
NumericMatrix smooth_gene(NumericMatrix expr_data, const int window_size  = 40);
double rdi_single_run_multiple_run_cpp(NumericMatrix& x, NumericMatrix& y, int d = 1, IntegerVector run_vec = 0); 
List calculate_rdi_multiple_run_cpp(NumericMatrix& expr_data, IntegerVector delays = 1, IntegerVector run_vec = 0, IntegerMatrix& super_graph, IntegerVector turning_points = 0, int method = 1); // = NULL 
NumericMatrix rdi_multiple_runs_conditioned_cpp(NumericMatrix& x, NumericMatrix& y, NumericMatrix& z, IntegerVector& z_delays, int d = 1, IntegerVector run_vec = 0);
NumericMatrix calculate_multiple_run_conditioned_rdi_cpp(NumericMatrix& expr_data, IntegerMatrix& super_graph, 
                                            NumericMatrix& max_rdi_value, IntegerMatrix& max_rdi_delays, IntegerVector run_vec, int k = 5)