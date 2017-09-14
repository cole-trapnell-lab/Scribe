#include <Rcpp.h>
using namespace Rcpp;

double vd_cpp(const int d);
double entropy_cpp(const NumericMatrix& x, int k = 5);
double mi_cpp(const NumericMatrix& x, const NumericMatrix& y, int k = 5);
List cmi_cpp(const NumericMatrix& x, const NumericMatrix& y, NumericMatrix z, int k = 5); //const NumericMatrix& z
