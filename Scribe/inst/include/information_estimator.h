#include <Rcpp.h>
using namespace Rcpp;

List knn_density(NumericMatrix x, NumericMatrix y, int k = 5);
NumericMatrix knn_density_2d(NumericVector x, NumericVector y, IntegerVector nGrids, int k = 5);

double vd_cpp(const int d);
double entropy_cpp(const NumericMatrix& x, int k = 5);
double mi_cpp(const NumericMatrix& x, const NumericMatrix& y, int k = 5, int normalize = 0);
List cmi_cpp(const NumericMatrix& x, const NumericMatrix& y, NumericMatrix z, int k = 5, int normalize = 0); //const NumericMatrix& z
List ucmi_cpp(const NumericMatrix& x, const NumericMatrix& y, NumericMatrix z, int k = 5, int method = 1, int k_density = 0, double bw = 0); // const NumericMatrix& z,
double umi_cpp(const NumericMatrix& x, const NumericMatrix& y, int k = 5, int method = 1, int k_density = 0, double bw = 0);
// List ucmi_cpp(const NumericMatrix& x, const NumericMatrix& y, NumericMatrix z, int k, int method, int k_density, double bw, NumericVector weight); // const NumericMatrix& z, 
// double umi_cpp(const NumericMatrix& x, const NumericMatrix& y, int k = 5, int method = 1, int k_density = 0, double bw = 0, NumericVector weight = 0);
