// ver 0.1; code review at Oct 2, 2017
// ver 0.2; update the documentation at Oct 2, 2017

#include <cmath>
#include <RcppArmadillo.h>
#include "../inst/include/information_estimator.h"

// [[Rcpp::depends(RcppArmadillo)]]
// #include <Rcpp.h>
// #include <omp.h>
// // [[Rcpp::plugins(openmp)]]

using namespace Rcpp;
using namespace arma;

//-------------------------------------------------------------------------------------------------------------------
/*
calculate direction information for single run
*/
//-------------------------------------------------------------------------------------------------------------------

// This function simulates the DIRECTED mutual information from X to Y when you have a SINGLE run of the processes
// Parameter n determines the the number previous time samples upon which the mi is conditioned
// double mi_cpp(double x, double y, double xy, int k = 5, int N)
// double cmi_cpp(double x, double y, double xyz, double xz, double yz, int k = 5, int N)

double di_single_run_cpp(NumericMatrix& x, NumericMatrix& y, int n = 5, bool uniformalize = false)
{
  
  int x_cols = x.cols(); int y_cols = y.cols();
  if(x_cols != y_cols)
  {
    stop("dimension of samples should be the same");
    return -1;
  }
  
  int tau = n; 
  int tot_len = x.rows() - tau; 
  
  int col_id = 0;
  int i = 0, j = 0; 
  NumericMatrix x_past(tot_len, n * y_cols);
  NumericMatrix y_past(tot_len, n * y_cols);
  
  // int max_thread = omp_get_max_threads();
  // omp_set_num_threads(max_thread);
  // #pragma omp parallel for shared(n, x_cols, y_cols, tot_len, x, y, x_past, y_past) private(i, j, col_id) //schedule(dynamic) default(none) //collapse(2) , _
  
  for(i = 0; i < tot_len; i ++)
  {
    for(j = 0; j < n * y_cols; j ++)
    { // identify the index for the column
      if(j == 0)
      {
        col_id = 0;
      }
      else
      {
        col_id = (j) % x_cols;
      }
      x_past(i, j) = x(i + (int) j / x_cols, col_id); // identify the index for the row
      y_past(i, j) = y(i + (int) j / x_cols, col_id);
    }
  }
  
  NumericMatrix y_past0 = y(Range(tau, tau + tot_len - 1), _);
  
  if(uniformalize == true) {
    // int k = 5, method = 1, k_density = 0; double bw = 0;
    List ucmi_res = ucmi_cpp(x_past, y_past0, y_past, 5, 1, 0, 0); // k, method, k_density, bw
    return ucmi_res["ucmi_res"]; 
  } else {
    List cmi_res = cmi_cpp(x_past, y_past0, y_past);
    return cmi_res["cmi_res"]; 
  } 
}

//' @title
//' di_single_run
//' @description
//' This function estimates the DIRECTED mutual information from X to Y when you have a SINGLE run of the process
//' 
//' @param x one random variable from the time-series data
//' 
//' @param y another random variable from the time-series data
//'
//' @param n Parameter n determines the the number of previous time samples upon which the mi is conditioned (delay in the formula I(x-->y)=I(x_{t-d};y_t|y_{t-1}) default: d=1)
//'
//' @param uniformalize Whether or not you want to use ucmi to calculate rdi. Default to be false. 
//' @details
//' \code{di_single_run} takes two random variables x and y as well as a delay n to estimate the direct information between variable x and y. 
//' @return a matrix for the condition mutual information estimators between all pairwise variables (x, y) in the data matrix x, y
//' @export
// [[Rcpp::export]]
double di_single_run(SEXP x, SEXP y, SEXP n, SEXP uniformalize) 
{ 
  NumericMatrix x_cpp(x); 
  NumericMatrix y_cpp(y); 
  
  int n_cpp = as<int>(n); 
  bool uniformalize_cpp = as<bool>(uniformalize); 
  
  double di_res = di_single_run_cpp(x_cpp, y_cpp, n_cpp, uniformalize_cpp);
  return di_res;
}

//-------------------------------------------------------------------------------------------------------------------
/*
calculate conditional direction information for single run
*/
//-------------------------------------------------------------------------------------------------------------------

double di_single_run_conditioned_cpp(NumericMatrix x, NumericMatrix y, NumericMatrix& z, int n = 5, bool uniformalize = false)
{
  uniformalize = false; // we cannot support kernel density estimator > 2 dim 
  int x_cols = x.cols(), y_cols = y.cols(); // , z_rows = z.rows(); // z_row is number of conditioning incoming nodes
  if(x_cols > y_cols)
  {
    stop("dimension of samples should be the same");
    return -1;
  }
  
  int tau = n; 
  int tot_len = x.rows() - tau; 
  
  NumericMatrix x_past = x(Range(tau - 1, tau - 2 + tot_len), _); 
  NumericMatrix yz_past = y(Range(tau - 1, tau - 2 + tot_len), _), tmp; 

  // need to redesign the cbind function in order to use parallel    
  // int i, j, k_x, k_yz; 
  // NumericMatrix x_past(tot_len, n * y_cols);
  // NumericMatrix yz_past(tot_len, n * y_cols * (1 + z_rows));
  
  // int max_thread = omp_get_max_threads();
  // omp_set_num_threads(max_thread);
  // #pragma omp parallel for shared(tot_len, n, y_cols, x, y, z, x_past, yz_past, z_rows) private(i, j, k_x, k_yz) //schedule(dynamic) default(none) //collapse(2) , _
  // for(i = 0; i < tot_len; i ++)
  // {
  //   for(j = 0; j < n; j ++)
  //   { // identify the index for the column
  //     for(k_x = 0; k_x < y_cols; k_x ++)
  //     {
  //       x_past(i, k_x + j * y_cols) = x(i + j, k_x); // identify the index for the row
  //     }
      
  //     for(k_yz = 0; k_yz < y_cols * (1 + z_rows); k_yz ++)
  //     {
  //       if(k_yz > y_cols * z_rows - 1) // if(k_yz < y_cols) /* k_yz is the last values */
  //       {
  //         yz_past(i, k_yz + j * y_cols * (1 + z_rows)) = y(i + j, k_yz);      
  //       }
  //       else
  //       { // i + j: column ind; k_yz % (1 + z_rows): row ind // k_yz / y_cols - 2
  //         yz_past(i, k_yz + j * y_cols * (1 + z_rows)) = z( (int) k_yz % (1 + z_rows), i + j ); // (int) k_yz / y_cols - 2       
  //       }
  //     }
  //   }
  // }

  // use cbind
  for(int i = 1; i < n + 1; i ++)
  {
    if(i > 1)
    {
      tmp = x(Range(tau - i, tau - i - 1 + tot_len), _);
      x_past = cbind(tmp, x_past);
      tmp = y(Range(tau - i, tau - i - 1 + tot_len), _);
      yz_past = cbind(tmp, yz_past);
    }
    
    for(int j = 0; j < z.cols(); j ++ )
    {
      tmp = z(Range(tau - i, tau - i - 1 + tot_len), _);
      yz_past = cbind(tmp, yz_past);
    } 
  }
  
  NumericMatrix y_past0 = y(Range(tau, tau + tot_len - 1), _);

  if(uniformalize == true) {
    List ucmi_res = ucmi_cpp(x_past, y_past0, yz_past, 5, 1, 0, 0); // k, method, k_density, bw
    return ucmi_res["ucmi_res"]; 
  } else {
    List cmi_res = cmi_cpp(x_past, y_past0, yz_past);
    return cmi_res["cmi_res"]; 
  }
}

//' @title
//' di_single_run_conditioned
//' @description
//' This function estimates the CONDITIONED DIRECTED mutual information from X to Y, conditioning on a third variable, z, when you have a SINGLE run of the processes
//' 
//' @param x one random variable from the time-series data
//' 
//' @param y another random variable from the time-series data
//'
//' @param z z is a dataframe (or matrix) containing the data of other processes upon the past of which the mi is conditioned
//' 
//' @param n Parameter n determines the the number of previous time samples upon which the mi is conditioned
//' 
//' @param uniformalize Whether or not you want to use ucmi to calculate rdi. Default to be false. 
//' @details
//' \code{di_single_run_conditioned} takes two random variables x, y and z as well as the parameter n to calculate the direct information conditioned on variable z. 
//' @return a numeric value for the condition mutual information estimator variables (x, y), conditioning on a third variable, z. 
//' @export
// [[Rcpp::export]]
double di_single_run_conditioned(SEXP x, SEXP y, SEXP z, SEXP n, SEXP uniformalize) 
{ 
  NumericMatrix x_cpp(x); 
  NumericMatrix y_cpp(y); 
  NumericMatrix z_cpp(z); 
  
  int n_cpp = as<int>(n); 
  bool uniformalize_cpp = as<bool>(uniformalize); 
  
  double di_res = di_single_run_conditioned_cpp(x_cpp, y_cpp, z_cpp, n_cpp, uniformalize_cpp);
  return di_res;
}

//-------------------------------------------------------------------------------------------------------------------
/*
calculate restricted direction information for many runs
*/
//-------------------------------------------------------------------------------------------------------------------

double rdi_many_runs_cpp(NumericMatrix& x, NumericMatrix& y, bool uniformalize = false)
{
  int dx = x.cols(), dy = y.cols();
  if(dx != dy)
  {
    stop("The dimension of time samples has to be the same for X and Y");
  }
  NumericVector ans(dx); // maybe mean or so? 
  NumericMatrix x_0, y_0, y_1;
  
  // OpenMP here too 
  int t; 
  // int max_thread = omp_get_max_threads();
  // omp_set_num_threads(max_thread);
  // #pragma omp parallel for shared(dx, x, y, ans) private(t, _, x_0, y_0, y_1) //schedule(dynamic) default(none) //collapse(2) , _
  for(t = 1; t < dx; t ++)
  {
    x_0 = x(_, Range(t - 1, t - 1));
    y_0 = y(_, Range(t, t));
    y_1 = y(_, Range(t - 1, t - 1)); 

    if(uniformalize == true) {
      List ucmi_res = ucmi_cpp(x_0, y_0, y_1, 5, 1, 0, 0); // k, method, k_density, bw
      ans[t - 1] = ucmi_res["ucmi_res"]; 
    } else {
      List cmi_res = cmi_cpp(x_0, y_0, y_1);
      ans[t - 1] = cmi_res["cmi_res"]; 
    }
  }
  
  return sum(ans); // return sum of ans? 
}

//' @title
//' rdi_many_runs
//' @description
//' This function estimates the DIRECTED mutual information from X to Y when you have multiple run of the processes.
//' 
//' @param x a random variable with multiple run of the same process.
//' 
//' @param y another random variable with multiple run of another process.
//'
//' @param uniformalize Whether or not you want to use ucmi to calculate rdi. Default to be false. 
//' @details
//' \code{rdi_many_runs} takes two random variables (with the same multiple realization of the two processes) and estimate 
//' the direct information between them at each time point. It then sums up those information estimators. This function can
//' only be used when you have hundreds runs of the same time-series experiment. 
//' @return a numeric value storing the DI from two multiple run variables
//' @export
// [[Rcpp::export]]
double rdi_many_runs(SEXP x, SEXP y, SEXP uniformalize) 
{ 
  NumericMatrix x_cpp(x); 
  NumericMatrix y_cpp(y); 
  bool uniformalize_cpp = as<bool>(uniformalize);
  
  double rdi_many_runs_res = rdi_many_runs_cpp(x_cpp, y_cpp, uniformalize_cpp);
  return rdi_many_runs_res;
}

//-------------------------------------------------------------------------------------------------------------------
/*
calculate restricted direction information for a single run
*/
//-------------------------------------------------------------------------------------------------------------------

double rdi_single_run_cpp(NumericMatrix& x, NumericMatrix& y, int d = 1, bool uniformalize = false) 
{
  int Nx = x.rows(); int Ny = y.rows(); //int dx = x.cols(); int dy = y.cols(); 
  
  if(Nx != Ny)
  {
    stop("The number of time samples has to be the same for X and Y");
    return -1;
  }
  
  NumericMatrix x_0 = x(Range(0, Nx - d - 1), _);
  NumericMatrix y_0 = y(Range(d, Ny - 1), _);
  NumericMatrix y_1 = y(Range(d - 1, Nx - 2), _);
  
  if(uniformalize == true) {
    List ucmi_res = ucmi_cpp(x_0, y_0, y_1, 5, 1, 0, 0); // k, method, k_density, bw
    return ucmi_res["ucmi_res"]; 
  } else {
    List cmi_res = cmi_cpp(x_0, y_0, y_1);
    return cmi_res["cmi_res"]; 
  }
}

//' @title
//' rdi_single_run
//' @description
//' This function estimates the RESTRICTED DIRECTED mutual information from X to Y when you have a SINGLE run of the processes
//' 
//' @param x one random variable from the time-series data
//' 
//' @param y another random variable from the time-series data
//' 
//' @param d delay in the formula I(x-->y)=I(x_{t-d};y_t|y_{t-1}) default: d=1
//'
//' @param uniformalize Whether or not you want to use ucmi to calculate rdi. Default to be false. 
//' @details
//' \code{rdi_single_run} takes two random variables x and y as well as a time delay d to estimate the restricted direct infomation between them.
//' @return a numeric value for the restricted direct information between x and y with a time delay d = 1. 
//' @export
// [[Rcpp::export]]
double rdi_single_run(SEXP x, SEXP y, SEXP d, SEXP uniformalize) 
{ 
  NumericMatrix x_cpp(x); 
  NumericMatrix y_cpp(y); 
  
  int d_cpp = as<int>(d); 
  bool uniformalize_cpp = as<bool>(uniformalize); 
  
  double rdi_single_run_res = rdi_single_run_cpp(x_cpp, y_cpp, d_cpp, uniformalize_cpp);
  return rdi_single_run_res;
}


//-------------------------------------------------------------------------------------------------------------------
/*
  calculate the lagged mutual information 
*/
//-------------------------------------------------------------------------------------------------------------------

double lmi_single_run_cpp(NumericMatrix& x, NumericMatrix& y, int delay = 1, bool uniformalize = false) 
{
  int Nx = x.rows(); int Ny = y.rows(); 
  
  if(Nx != Ny)
  {
    stop("The number of time samples has to be the same for X and Y");
    return -1;
  }
  
  NumericMatrix x_0 = x(Range(0, Nx - delay - 1), _);
  NumericMatrix y_0 = y(Range(delay, Ny - 1), _);
  
  // Rcout << "\n before running CMI " << std::endl;
  return mi_cpp(x_0, y_0); 
}

//' @title
//' lmi_single_run
//' @description
//' This subroutine calculates the lagged mutual information 
//' 
//' @param x one random variable from the time-series data
//' 
//' @param y another random variable from the time-series data
//'
//' @param delay Time lags used to estimate the RDI values  
//'
//' @param uniformalize Whether or not you want to use ucmi to calculate rdi. Default to be false. 
//' @details
//' \code{lmi_single_run} takes two random variable x and y and estimated their mutual information with a time lag d. 
//' using the KSG estimator. 
//' It relies on the ANN package to query the kNN with KDTree algorithm. 
//' @return a numeric value for the mutual information estimator between two variables (x, y) with a time lag d.
//' @export
// [[Rcpp::export]]
double lmi_single_run(SEXP x, SEXP y, SEXP delay, SEXP umi) 
{ 
  NumericMatrix x_cpp(x); 
  NumericMatrix y_cpp(y); 

  int delay_cpp = as<int>(delay); 
  bool umi_cpp = as<bool>(umi); 
  
  double lmi_res = lmi_single_run_cpp(x_cpp, y_cpp, delay_cpp, umi_cpp);
  return lmi_res;
}

// lmi for multiple runs
// [[Rcpp::export]]
double lmi_multiple_run_cpp(NumericMatrix& x, NumericMatrix& y, int d = 1, IntegerVector run_vec = 0, bool uniformalize = false) 
{
  int Nx = x.rows(); int Ny = y.rows(); int run_len; 
  
  if(Nx != Ny)
  {
    stop("The number of time samples has to be the same for X and Y");
    return -1;
  }
 
  int run_max = max(run_vec); // run_vec: vector to represent the run, need to be just numeric numbers and starts from 0. 
  int current_ind = 0; // current initial index for storing new slice of data 

  NumericMatrix x_0_res, y_0_res, tmp; // matrix passed to cmi function
  IntegerVector index_all = seq_len(run_vec.size()) - 1, current_run_ind; // cell index from 0 to maximal number of cells; run index corresponding to current run during iteration 

  if(run_max == 0) 
  {
    x_0_res = x(Range(0, Nx - d - 1), _); // x should only has 1 colum 
    y_0_res = y(Range(d, Ny - 1), _);
  } else {

    // correctly set the dimensionality for the x_0, y_0 and y_1 matrix 
    int tot_len = Nx - d * (run_max + 1); // number of samples minus the delay times total number of runs 

    mat x_0(tot_len, 1), y_0(tot_len, 1); // arma matrix 
    
    uvec pos; // arma uvec 
    vec vals; // arma vec 

    IntegerVector rng; 

    for(int i = 0; i <= run_max; i++) // concatenate the data 
    {
      current_run_ind = index_all[run_vec == i]; // get the cells that belong to a particular run 
      run_len = current_run_ind.size(); 

      if(current_run_ind.size() < d)
      {
        stop("Each run has to have more samples than the designated 'delays'");
        return -1;
      }

      rng = Range(current_ind, current_ind + run_len - d - 1); 

      pos = as<uvec>( rng ); 
      
      // assuming x is only one column (Cx, Cy is not used)
      tmp = x(Range(min(current_run_ind), min(current_run_ind) + run_len - d - 1), Range(0, 0)); //x(Range(1, 2), Range(1, 1));
      vals = as<vec>( tmp ); 
      x_0.elem(pos) = vals; 

      IntegerVector rng_vec = Range(min(current_run_ind), min(current_run_ind) + run_len - d - 1); 
      // for(int test_na = 0; test_na < vals.n_elem; test_na ++) 
      // {
      //   if(arma::is_finite(vals[test_na]) == false) 
      //   {
      //     // Rcout << "identify non-finite values for x_0 here; i is " << i << " rng_vec is " << rng_vec << std::endl; 
      //   }
      // }


      tmp = y(Range(min(current_run_ind) + d, min(current_run_ind) + run_len - 1), Range(0, 0)); //x(Range(1, 2), Range(1, 1));
      vals = as<vec>( tmp ); 
      y_0.elem(pos) = vals; 

      // for(int test_na = 0; test_na < vals.n_elem; test_na ++) 
      // {
      //   if(arma::is_finite(vals[test_na]) == false) 
      //   {
      //     Rcout << "identify non-finite values for x_0 here; i is " << i << " rng_vec is " << rng_vec << std::endl; 
      //   }
      // }

      // for(int test_na = 0; test_na < vals.n_elem; test_na ++) 
      // {
      //   if(arma::is_finite(vals[test_na]) == false) 
      //   {
      //     Rcout << "identify non-finite values for x_0 here; i is " << i << " rng_vec is " << rng_vec << std::endl; 
      //   }
      // }

      current_ind = current_ind + run_len - d; // move to the the next position after filling the current run (note that there is no - 1)
    }


    x_0_res = as<NumericMatrix>( wrap(x_0) ); // this conversion have problems? 
    y_0_res = as<NumericMatrix>( wrap(y_0) );
  }  

  return mi_cpp(x_0_res, y_0_res);  
}

//' @title
//' lmi_multiple_run
//' @description
//' This subroutine calculates the lagged mutual information with multiple realization of the same processes. 
//' 
//' @param x one random variable from the time-series data
//' 
//' @param y another random variable from the time-series data
//'
//' @param d Time lags used to estimate the RDI values  
//'
//' @param run_vec A integer vector encodes the information of the run id (run id start from 0) 
//'
//' @param uniformalize Whether or not you want to use ucmi to calculate rdi. Default to be false. 
//' @details
//' \code{lmi_multiple_run} takes two random variables x and y, each has multiple realization of the same process and estimated their mutual information with a time lag d. 
//' using the KSG estimator. 
//' It relies on the ANN package to query the kNN with KDTree algorithm. 
//' @return a numeric value for the estimated mutual information between two variables (x, y) with a time lag d. 
//' @export
// [[Rcpp::export]]
double lmi_multiple_run(SEXP x, SEXP y, SEXP d, SEXP run_vec, SEXP umi) 
{ 
  NumericMatrix x_cpp(x); 
  NumericMatrix y_cpp(y); 

  int d_cpp = as<int>(d); 

  IntegerVector run_vec_cpp(run_vec); 
  bool umi_cpp = as<bool>(umi);

  double lmi_res = lmi_multiple_run_cpp(x_cpp, y_cpp, d_cpp, run_vec_cpp, umi_cpp);
  return lmi_res;
}
//-------------------------------------------------------------------------------------------------------------------
/*
calculate conditional restricted direction information for a single run
*/
//-------------------------------------------------------------------------------------------------------------------

double rdi_single_run_conditioned_cpp(NumericMatrix& x, NumericMatrix& y, NumericMatrix& z, NumericVector& z_delays, int d = 1, bool uniformalize = false) 
{
  uniformalize = false; // we cannot support kernel density estimator > 2 dim 

  int Nx = x.rows(); int Ny = y.rows(); int Nz = z.rows(); // int dx = x.cols(); int dy = y.cols(); int dz = z.cols(); 
  if(Nx != Ny | Nx != Nz)
  {
    stop("The number of time samples has to be the same for X and Y or Z");
    return -1;
  }
  
  z_delays.push_back(d);
  int tau = max(z_delays); 
  
  int tot_len = x.rows() - tau;  
  
  NumericMatrix yz = y(Range(tau - 1, tau - 2 + tot_len), _);
  
  NumericMatrix tmp; 
  
  // OpenMP here too 
  // int max_thread = omp_get_max_threads();
  // omp_set_num_threads(max_thread);
  // #pragma omp parallel for shared(z, tau, z_delays, tot_len) private(j, tmp, yz) //schedule(dynamic) default(none) //collapse(2) , _
    
  // using cbind
  for(int j = 0; j < z.cols(); j ++ )
  {
    tmp = z(Range(tau - z_delays[j], tau - z_delays[j] + tot_len - 1), Range(j, j));
    yz = cbind(tmp, yz);
  } 
  
  NumericMatrix x_0 = x(Range(tau - d, tau - d + tot_len - 1), _);
  NumericMatrix y_0 = y(Range(tau, tau + tot_len - 1), _);
  
  if(uniformalize == true) {
    List ucmi_res = ucmi_cpp(x_0, y_0, yz, 5, 1, 0, 0); // k, method, k_density, bw
    return ucmi_res["ucmi_res"]; 
  } else {
    List cmi_res = cmi_cpp(x_0, y_0, yz);
    return cmi_res["cmi_res"]; 
  }
}

//' @title
//' rdi_single_run_conditioned
//' @description
//' This function estimates the CONDITIONED DIRECTED mutual information from X to Y CONDITIONED ON Z when you have a SINGLE run of the processes
//' 
//' @param x one random variable from the time-series data
//' 
//' @param y another random variable from the time-series data
//' 
//' @param z z is a dataframe or matrix consisting of the data for different variables which will be conditioned on. 
//' 
//' @param z_delays z_delay is also a dataframe or matrix consisting of the delays to be applied to different variables
//' 
//' @param d delay in the formula I(x-->y)=I(x_{t-d};y_t|y_{t-1}) default: d=1
//'
//' @param uniformalize Whether or not you want to use ucmi to calculate rdi. Default to be false. 
//' @details
//' \code{rdi_single_run_conditioned} takes two random variables x and y as well as the parameter d to calculate the restricted direct information conditioned on variable z. 
//' @return a numeric value for the estimated condition mutual information between variable x and y conditioning on a third variable z 
//' @export
// [[Rcpp::export]]
double rdi_single_run_conditioned(SEXP x, SEXP y, SEXP z, SEXP z_delays, SEXP d, SEXP uniformalize) 
{ 
  NumericMatrix x_cpp(x); 
  NumericMatrix y_cpp(y); 
  NumericMatrix z_cpp(z); 
  
  NumericVector z_delays_cpp(z_delays); 
  
  int d_cpp = as<int>(d); 
  bool uniformalize_cpp = as<bool>(uniformalize); 
  
  double rdi_single_run_conditioned_res = rdi_single_run_conditioned_cpp(x_cpp, y_cpp, z_cpp, z_delays_cpp, d_cpp, uniformalize_cpp);
  return rdi_single_run_conditioned_res;
}

// ----------------------------------------------------------------------------------------------------------------------------------------------- //

// Convert a NumericMatrix (R object) to an Armadillo matrix (C++ object)
// //[[Rcpp::export]]
// arma::mat R2armaMat_num(NumericMatrix rMat)
// {
//   arma::mat armaMat = arma::mat(rMat.begin(), rMat.nrow(), rMat.ncol(), false);
//   return armaMat;
// }

//-------------------------------------------------------------------------------------------------------------------
/*
calculate restricted direction information for all pairs of genes from the expr_data matrix
*/
//-------------------------------------------------------------------------------------------------------------------

// [[Rcpp::export]]
List extract_max_rdi_value_delay(NumericMatrix rdi_result, IntegerVector delays)
{
  
  int i, j, k, max_delay = delays[0], n_genes = rdi_result.rows();
  int delays_len = (int) rdi_result.cols() / n_genes; 
  double tmp, max_val;
  
  // NumericVector rdi_vec(delays_len);
  NumericMatrix max_rdi_value(n_genes, n_genes);
  IntegerMatrix max_rdi_delays(n_genes, n_genes);
  
  if(delays.length() != delays_len)
  {
    stop("Length of delays doesn't match up the rdi result");
    return -1; 
  }
  for(i = 0; i < n_genes; i ++) 
  {
    for(j = 0; j < n_genes; j ++)
    {
      if(i != j) 
      {
        max_val = rdi_result(i, j); // deal with the case where there is only a single delay tested 
        for(k = 0; k < delays_len; k ++)
        { 
          tmp = rdi_result(i, j + k * n_genes);
          if(max_val < tmp) {
            max_val = tmp;
            max_delay = delays[k];
          }
        }
        max_rdi_value(i, j) = max_val;
        max_rdi_delays(i, j) = max_delay; 
      }
    }
  }
  
  return List::create(Rcpp::Named("max_rdi_value") = max_rdi_value, 
                      Rcpp::Named("max_rdi_delays") = max_rdi_delays);
}

/** 
library(InformationEstimator)
extract_max_rdi_value_delay(AT1_RDI_parallel_res, c(5, 10, 15))
*/

// [[Rcpp::export]]
List calculate_rdi_cpp(NumericMatrix& expr_data, IntegerVector delays, IntegerMatrix& super_graph, IntegerVector& turning_points, int method, bool uniformalize = false) //, method: 1 rdi; 2: lmi 
{
  const int n_genes(expr_data.cols()), n_samples(expr_data.rows());

  if(max(delays) > n_samples)
  {
    stop("Number of samples has to be larger than the delays"); 
    return -1; 
  }
  
  IntegerMatrix max_rdi_delays(n_genes, n_genes); // declare max_rdi_delays when turning_points provided 
  
  // Convert arguments to Armadillo objects
  // const arma::mat data = R2armaMat_num(expr_data); // Continuous data does not need reindexing R->C++
  // const int n_pairs = (int)(((double)n_genes/2.0) * (double)(n_genes+1));
  
  // // Number of cores to use
  // //const int cores = sysconf(_SC_NPROCESSORS_ONLN);
  // omp_set_num_threads(cores);
  
  // Compute RDI in parallel
  int delays_len = delays.length(); 
  
  int RDI_ncols; 
  if(turning_points.length() == n_genes) 
  {
    RDI_ncols = n_genes; 
  } else {
    RDI_ncols = n_genes * delays_len; 
  }

  NumericMatrix RDI(n_genes, RDI_ncols);     
  std::fill(RDI.begin(), RDI.end(), 0); //NA_REAL 
  NumericMatrix expr_1(n_samples, 1), expr_2(n_samples, 1);
  
  // this part maybe slow if we do paralleling because RDI matrix is huge 
  
  int i, j, k; 
  IntegerVector current_pair; 
  // #pragma omp parallel for shared(n_genes, expr_data, delays_len, RDI) private(i, j, k, expr_1, expr_2, _) //schedule(dynamic) default(none) //collapse(2) , _
  for(int super_graph_ind = 0; super_graph_ind < super_graph.rows(); super_graph_ind ++) 
  {
    current_pair = super_graph(Range(super_graph_ind, super_graph_ind), _);
    i = current_pair[0]; j = current_pair[1]; 
    
    expr_1 = expr_data(_, Range(i, i)); //Range(0, n_genes - 1)
    expr_2 = expr_data(_, Range(j, j)); //Range(0, n_genes - 1)

    // int max_thread = omp_get_max_threads();
    // omp_set_num_threads(max_thread);
    // #pragma omp parallel for shared(delays_len, i, j, n_genes, expr_1, expr_2, RDI) private(k) //schedule(dynamic) default(none) //collapse(2) , _
    
    // if we provide turning_points estimation, we will use that to determine the time delay
    if(turning_points.length() == n_genes) // !Rf_isNull(turning_points)
    {
      Rcout << "using user provided information about time-delay " << turning_points.length() << std::endl;
      int current_delay = turning_points[i] - turning_points[j];
      
      if(i == j) 
      {
        continue; // keep diagnoal as 0 
      }

      if(method == 1)
      {
        RDI(i, j) = rdi_single_run_cpp(expr_1, expr_2, current_delay, uniformalize); // how to deal with delays include multiple values?
      }
      else if(method == 2)
      { // + k * n_genes
        RDI(i, j) = lmi_single_run_cpp(expr_1, expr_2, current_delay, uniformalize); // how to deal with delays include multiple values?
      }
      max_rdi_delays(i, j) = current_delay;
    } else {
      for(k = 0; k < delays_len; ++ k)
      {
        if(i == j) 
        {
          continue; // keep diagnoal as 0 
        }
        
        if(method == 1)
        {
          RDI(i, j + k * n_genes) = rdi_single_run_cpp(expr_1, expr_2, delays[k], uniformalize); // how to deal with delays include multiple values?
        }
        else if(method == 2)
        {
          RDI(i, j + k * n_genes) = lmi_single_run_cpp(expr_1, expr_2, delays[k], uniformalize); // how to deal with delays include multiple values?
        }
      }
    }
  }

  // perform the multiple run test (add a new argument to represent the vector for different runs)
  if(turning_points.length() == n_genes) // !Rf_isNull(turning_points)
  {
    return List::create(Rcpp::Named("RDI") = RDI, 
                        Rcpp::Named("delays") = R_NilValue,
                        Rcpp::Named("max_rdi_value") = RDI, 
                        Rcpp::Named("max_rdi_delays") = max_rdi_delays);
  }
  else{
    // run extract_max_rdi_value_delay 
    List max_rdi_value_delay_list = extract_max_rdi_value_delay(RDI, delays); 
    
    return List::create(Rcpp::Named("RDI") = RDI, 
                        Rcpp::Named("delays") = delays,
                        Rcpp::Named("max_rdi_value") = max_rdi_value_delay_list["max_rdi_value"], 
                        Rcpp::Named("max_rdi_delays") = max_rdi_value_delay_list["max_rdi_delays"]);
  }
}

//' @title
//' calculate_rdi_cpp_wrap
//' @description
//' This function estimates the DIRECTED mutual information for all gene pairs in the expr_data matrix when you have a SINGLE run of the processes
//' 
//' @param expr_data a matrix for all variables in a time-series where each row is a time point and each column is a variable (for example, a gene). 
//' The rows are ordered according to time, from earliest to latest. 
//' 
//' @param delays An integer vector storing the time delays between pairs of variables you would like to try. 
//' 
//' @param super_graph An integer matrix where each row is the variable IDs or gene pairs you would like to estimate restricted direct information (RDI) from. 
//'  
//' @param turning_points Either 0 or a numeric vector describing the inflection point (linear trajectory) or branch point (bifurcation point) for each gene. 
//' If the turning_point for each gene is provided, the time delay will be estimated based on the turning point. 
//'
//' @param method An integer of either 1 or 2 to determine which information metric will be used to quantify the causality. 
//' If method is 1, then lagged mutual information will be used; if method is 2, then the restricted direct information will be used. 
//'
//' @param uniformalize Whether or not you want to use ucmi to calculate rdi. Default to be false. 
//' @details
//' \code{calculate_rdi_cpp_wrap} takes an expression matrix (expr_data) and possible gene pairs (encoded in super_graph) to estimate the restricted 
//' direct information based on the time delay which can be estimated from the turning_point vector. 
//' @return a numeric matrix of conditional RDI values for all possible pairs of genes from expr_data. If the gene pairs is not encoded in the super_graph, it will remain as 0 in the matrix.  
//' @export
// [[Rcpp::export]]
List calculate_rdi_cpp_wrap(SEXP expr_data, SEXP delays, SEXP super_graph, SEXP turning_points, SEXP method, SEXP uniformalize) //, SEXP R_cores
{ 
  NumericMatrix expr_data_cpp(expr_data); 
  IntegerVector delays_cpp(delays); 
  
  IntegerMatrix super_graph_cpp(super_graph); 
  IntegerVector turning_points_cpp(turning_points);
  int method_cpp = as<int>(method); 
  bool uniformalize_cpp = as<bool>(uniformalize); 
  
  List rdi_list = calculate_rdi_cpp(expr_data_cpp, delays_cpp, super_graph_cpp, turning_points_cpp, method_cpp, uniformalize_cpp); //cores
 
  return rdi_list;
}

/** 
 library(InformationEstimator)
 sample_row = 30; sample_col = 50; 
 noise = matrix(rnorm(mean = 0, sd = 1e-10, sample_row * sample_col), nrow = sample_col)

 test <- calculate_rdi(t(lung_exprs_AT1[1:30, 1:50]) + noise, c(5, 10, 15))
 */

// max_rdi_delays could be double matrix in future (current IntegerMatrix)
// [[Rcpp::export]]
List extract_top_incoming_nodes_delays(NumericMatrix max_rdi_value, IntegerMatrix max_rdi_delays, int k = 1)
{
  int n_genes = max_rdi_value.rows(), i, j, current_k_ind;

  NumericMatrix top_incoming_values(n_genes, k + 1);
  IntegerMatrix top_incoming_delays(n_genes, k + 1), top_incoming_nodes(n_genes, k + 1);
  
  IntegerVector top_k_plus_1 = seq_len(k); top_k_plus_1.push_front(0); 
  
  NumericVector x, x_ordered;
  IntegerVector x_order_ind, top_k_plus_1_x_order_ind; 


  for(i = 0; i < n_genes; i ++)
  {
    max_rdi_value(i, i) = -5; // avoid diag is calculated as the incoming node 
    x = max_rdi_value(_, i); // first dimension is source; second dimension is target 
   
    NumericVector x_ordered = clone(x).sort(true); // decreasing sort; rdi value from highest to lowest 
    x_order_ind =  match(x_ordered, x) - 1; // order starts from 1; all duplicated points get the same location 

    top_k_plus_1_x_order_ind = x_order_ind[top_k_plus_1];  // top k + 1  ordered genes index

    max_rdi_value(i, i) = 0; // go back to the original value 

    for(j = 0; j < k + 1; j ++) // for each gene get the top k + 1 input node 
    {
      current_k_ind = top_k_plus_1_x_order_ind[j]; // index x_ordered to get the top k's id (node index)

      top_incoming_nodes(i, j) = current_k_ind; //top_k_plus_1_x_order_ind[j]; //[top_k_plus_1];
      top_incoming_delays(i, j) = max_rdi_delays(current_k_ind, i); // current_k_ind: index for current maximal incoming gene
      top_incoming_values(i, j) = x_ordered[j]; // rdi value  
    }
  } 

  return List::create(Rcpp::Named("top_incoming_nodes") = top_incoming_nodes, 
              Rcpp::Named("top_incoming_delays") = top_incoming_delays,
              Rcpp::Named("top_incoming_values") = top_incoming_values); 
}

/**
library(InformationEstimator)
sample_row = 30; sample_col = 50; 
noise = matrix(rnorm(mean = 0, sd = 1e-10, sample_row * sample_col), nrow = sample_col)
test <- calculate_rdi(t(lung_exprs_AT1[1:30, 1:50]) + noise, c(5, 10, 15))

test_extract <- extract_top_incoming_nodes_delays(test$max_rdi_value, test$max_rdi_delays, k = 2)
*/

//-------------------------------------------------------------------------------------------------------------------
/*
calculate conditional restricted direction information for all genes in a expression matrix 
*/
//-------------------------------------------------------------------------------------------------------------------

NumericMatrix calculate_conditioned_rdi_cpp(NumericMatrix& expr_data, IntegerMatrix& super_graph, 
                                            NumericMatrix& max_rdi_value, IntegerMatrix& max_rdi_delays, int k = 1, bool uniformalize = false) //, const int cores, const bool verbose
{  
  uniformalize = false; // we cannot support kernel density estimator > 2 dim 
  
  if(expr_data.cols() < 3) {
    stop("You need at least 3 genes to calculate conditional RDI values");
  }
  
  NumericVector k_ncol = NumericVector::create(k, expr_data.cols() - 2); // condition on at most n_gene - 2 incoming genes
  k = min(k_ncol); // minal value from k and the column 
  
  List top_incoming_nodes_delays_list = extract_top_incoming_nodes_delays(max_rdi_value, max_rdi_delays, k); 

  IntegerMatrix top_incoming_nodes = top_incoming_nodes_delays_list["top_incoming_nodes"];
  IntegerMatrix top_incoming_delays = top_incoming_nodes_delays_list["top_incoming_delays"]; 
  NumericMatrix top_incoming_values = top_incoming_nodes_delays_list["top_incoming_values"];  

  int n_genes = expr_data.cols(), n_sample = expr_data.rows(), tau, total_length, id3, delay_id3, current_id_tmp, delay; //tmp: current target  

  // valid_top_k_in_k_ordered is used to get the valid_top_k to the original order (setdiff function annoyingly sorts the data randomly)
  IntegerVector gene_pair(2), tmp(1), valid_top_k_in_k(k), top_k_plus_1(k + 1), valid_top_k(k), valid_top_k_in_k_ordered, valid_delay_tmp(1); // optimize by setting the correct dimension   // ind_for_top_k(k),  
  
  NumericMatrix cRDI_mat(n_genes, n_genes), expr1_t, expr2_t, past_tmp; 
  std::fill(cRDI_mat.begin(), cRDI_mat.end(), 0); // initialize the cRDI matrix as 0 instead of C ++ NA values 
  
  // OpenMP here too (not thread safe) 
  for(int i = 0; i < super_graph.rows(); i ++) 
  {
    gene_pair = super_graph(i, _); 
    top_k_plus_1 = top_incoming_nodes(Range(gene_pair[1], gene_pair[1]), _); // top k incoming nodes 
    tmp = gene_pair[0]; // current incoming node 

    if(intersect(top_k_plus_1, tmp).length() > 0) // if current node is in the top_k_plus_1, remove k; else remove the smallest one 
    {
      valid_top_k = setdiff(top_k_plus_1, tmp); // note that setdiff annoying the order of valid_top_k
    }
    else
    {
      tmp = top_k_plus_1[k]; // remove the smallest one, the k + 1 incoming node 
      valid_top_k = setdiff(top_k_plus_1, tmp); // note that setdiff annoyingly changes the order of valid_top_k
    }
    
    // assign the node indices to ind_for_top_k
    valid_top_k_in_k = match(valid_top_k, top_k_plus_1) - 1; // match gives order and the index (1 to k) starts from 1 -> valid_top_k_in_k
    valid_top_k_in_k_ordered = clone(valid_top_k_in_k).sort(); // the index for top_k_plus_1 is sorted
    valid_top_k = top_k_plus_1[valid_top_k_in_k_ordered]; // now the valid_top_k's value follows the ordering in the top_k_plus_1 
    
    if(valid_top_k.size() < k) 
    { // avoid any weird rare situation 
      Rcout << "all incoming node has the same RDI value for gene pair " << gene_pair[0] << " and " << gene_pair[1] << ". The target gene may have no expression." << std::endl;
      cRDI_mat(gene_pair[0], gene_pair[1]) = max_rdi_value(gene_pair[0], gene_pair[1]); 
      continue;
    }

    IntegerVector valid_top_k_incoming_delays(k); // re-initialize this everytime because we have push_back later

    // assign the valid top_k_incoming delays based on the sorted valid_top_k_in_k_ordered index 
    for(int j = 0; j < k; j ++ )
    {
      current_id_tmp = valid_top_k_in_k_ordered[j]; // range: 0 -> k - 1
      valid_delay_tmp = top_incoming_delays(Range(gene_pair[1], gene_pair[1]), Range(current_id_tmp, current_id_tmp)); // the current_id_tmp^{th} delay value 
      valid_top_k_incoming_delays[j] = valid_delay_tmp[0]; 
    }

    delay = max_rdi_delays(gene_pair[0], gene_pair[1]); // the delay corresponds to max RDI value for the current gene pair 
    valid_top_k_incoming_delays.push_back(delay); // identify and use the biggest value as the delay for all genes (to avoid index over-flow or out of boundary); last value corresponds to the delay for the current test gene pair 
    tau = max(valid_top_k_incoming_delays); 
    total_length = n_sample - tau; 

    expr1_t = expr_data(Range(tau - delay, tau - 1 - delay + total_length), Range(gene_pair[0], gene_pair[0])); 
    expr2_t = expr_data(Range(tau, tau - 1 + total_length), Range(gene_pair[1], gene_pair[1])); 
    NumericMatrix yz = expr_data(Range(tau - 1, tau  - 2 + total_length), Range(gene_pair[1], gene_pair[1])); // initialize yz as expr2_t with one time lag 

    // OpenMP here too 
    for(int id = 0; id < k; id ++)
    {
      id3 = valid_top_k[id];
      delay_id3 = valid_top_k_incoming_delays[id];
      past_tmp = expr_data(Range(tau - delay_id3, tau  - delay_id3 - 1 + total_length),
                            Range(id3, id3));
      yz = cbind(past_tmp, yz);
    }

    // std::ofstream file("test_cond_rdi.txt");
    // if (file.is_open())
    // {
    //   file << "current super_graph ind is " << gene_pair << '\n';
    //   file << "expr1_t\n" << expr1_t << "\nexpr2_t\n" << expr2_t << "\nyz\n" << yz << '\n';
    // }

    if(uniformalize == true) {
      List uniformalize_res = ucmi_cpp(expr1_t, expr2_t, yz, 5, 1, 0, 0); // k, method, k_density, bw
      cRDI_mat(gene_pair[0], gene_pair[1]) =  uniformalize_res["cmi_res"];    
    } else {
      List cmi_res = cmi_cpp(expr1_t, expr2_t, yz);
      cRDI_mat(gene_pair[0], gene_pair[1]) =  cmi_res["cmi_res"];    
    }
  }
  
  // return results
  return cRDI_mat; 
}

//' @title
//' calculate_conditioned_rdi_cpp_wrap
//' @description
//' This function estimates the conditional DIRECTED mutual information for all gene pairs in the expr_data matrix when you have a SINGLE run of the processes.
//' 
//' @param expr_data a matrix for all variables in a time-series where each row is a time point and each column is a variable (for example, a gene). 
//' The rows are ordered according to time, from earliest to latest. 
//' 
//' @param super_graph An integer matrix where each row is the variable IDs or gene pairs you would like to estimate restricted direct information (RDI) from. 
//'  
//' @param max_rdi_value A numeric matrix where each element corresponding to the maximal (identified from a series of time lags when calculating the rdi values) rdi value between two variables.  
//' 
//' @param max_rdi_delays An integer matrix where each element corresponding to the time delay corresponding to the maximal (identified from a series of time lags when calculating the rdi values) rdi value between two variables.  
//'
//' @param k An integer for the number of incoming nodes to be conditioned on 
//'
//' @param uniformalize Whether or not you want to use ucmi to calculate rdi. Default to be false. 
//' @details
//' \code{calculate_conditioned_rdi_cpp_wrap} takes an expression matrix (expr_data) and possible gene pairs (encoded in super_graph), as well as the matrices of maximal rdi value or the delays corresponding to those values to 
//' estimate the conditional restricted direct information, conditioning on top k incoming nodes. 
//' @return a numeric matrix of conditional RDI values for all possible pairs of genes from expr_data. If the gene pairs is not encoded in the super_graph, it will remain as 0 in the matrix.  
//' @export
// [[Rcpp::export]]
NumericMatrix calculate_conditioned_rdi_cpp_wrap(SEXP expr_data, SEXP super_graph,
                               SEXP max_rdi_value, SEXP max_rdi_delays, SEXP k, SEXP uniformalize) //, SEXP R_cores
{
  NumericMatrix expr_data_cpp(expr_data);
  IntegerMatrix super_graph_cpp(super_graph);
  NumericMatrix max_rdi_value_cpp(max_rdi_value);
  IntegerMatrix max_rdi_delays_cpp(max_rdi_delays);
  int k_cpp = as<int>(k);
  bool uniformalize_cpp = as<int>(uniformalize);
  
  NumericMatrix cRDI = calculate_conditioned_rdi_cpp(expr_data_cpp, super_graph_cpp, max_rdi_value_cpp, max_rdi_delays_cpp, k_cpp, uniformalize_cpp); //cores
  return cRDI;
}

/**
library(InformationEstimator)
dim(lung_exprs_AT1)
 
sample_row = 123; sample_col = 218; 
set.seed(2017)
noise = matrix(rnorm(mean = 0, sd = 1e-10, sample_row * sample_col), nrow = sample_row)
a <- Sys.time()
rdi_list <- calculate_rdi(lung_exprs_AT1[1:sample_row, 1:sample_col] + noise, c(5, 10, 15), method = 1)
b <- Sys.time()
rdi_time <- b - a 
  
top_delays <- extract_top_incoming_nodes_delays(rdi_list$max_rdi_value, rdi_list$max_rdi_delays, k = 1)
   
tmp <- expand.grid(1:sample_col, 1:sample_col, stringsAsFactors = F)
super_graph <- tmp[tmp[, 1] != tmp[, 2], ] - 1 # convert to C++ index 
 
max_rdi_value <- as.matrix(rdi_list$max_rdi_value); max_rdi_delays <- as.matrix(rdi_list$max_rdi_delays); 
a <- Sys.time()
cRDI <- calculate_conditioned_rdi_cpp_wrap(t(lung_exprs_AT1[1:sample_row, 1:sample_col]) + noise, as.matrix(super_graph), max_rdi_value, max_rdi_delays, 2L)
b <- Sys.time()
crdi_time <- b - a 
 
# super_graph 
a <- Sys.time()
cRDI_R <- calculate_conditioned_rdi(t(lung_exprs_AT1[1:sample_row, 1:sample_col]), super_graph = NULL, rdi_list, 2L)
b <- Sys.time()
R_crdi_time <- b - a 
 
*/

//' @title
//' smooth_gene
//' @description
//' This subroutine takes a time-series data and returns a moving average for the data. 
//' 
//' @param expr_data a matrix for all variables in a time-series where each row is a time point and each column is a variable (for example, a gene). 
//' The rows are ordered according to time, from earliest to latest. 
//' 
//' @param window_size Integer value for the smoothing window used for calculating the moving average.  
//'
//' @details
//' \code{entropy} takes a integer of dimensions and then calculate the olume of a d-dimensional unit ball for Euclidean norm
//' using the formula: 0.5 * d * log(pi) - log(gamma(0.5 * d + 1))
//' It's implimented in C++, providing a (small) increase in speed over the R equivalent.
//' @return a updated matrix with gene expression smoothed with window size equal to window_size
//' @export
// [[Rcpp::export]]
NumericMatrix smooth_gene(NumericMatrix& expr_data, const int window_size = 40)
{ // columns: genes; rows: cells 
  int win_range = expr_data.rows() - window_size;
  NumericMatrix expr_data_smooth(expr_data(Range(0, win_range), _));
  
  NumericVector tmp(win_range + 1), tmp_1;
  for(int i = 0; i < expr_data.cols(); i ++) // genes 
  {
    for(int j = 0; j <= win_range; j ++) // cells
    {
      tmp_1 = expr_data(Range(j, j + window_size - 1), Range(i, i)); 
      tmp[j] = mean(tmp_1);
    }
    expr_data_smooth(_, i) = tmp;
  }
  
  return expr_data_smooth;
}

/**
run_new_dpt <- function(cds, norm_method = 'none', root = NULL, color_by = 'Cell_type'){
  message('root should be the id to the cell not the cell name ....')
  norm_data <- monocle:::normalize_expr_data(cds, norm_method = norm_method)
  norm_data <- t(norm_data)
  duplicated_genes <- which(base::duplicated.array(norm_data))
  norm_data[duplicated_genes, 1] <- norm_data[duplicated_genes, 1] + rnorm(length(duplicated_genes), 0, 1)
  dm <- DiffusionMap(norm_data)
  dpt <- DPT(dm)
  ts <- dm@transitions
  M <- destiny:::accumulated_transitions(dm)
  if(is.null(root)){
  }
  else{
    dm <- DiffusionMap(norm_data)
    dpt <- DPT(dm, tips = root)
  }
  if('Hours' %in% colnames(pData(cds)))
    pData(cds)$Hours <- pData(cds)[, color_by]
  p1 <- qplot(DM$DC1, DM$DC2, colour = pData(cds)$Hours)
  branch <- dpt@branch
  if(is.null(root))
    root <- which.min(pData(cds)$Pseudotime)
  pt <- dpt[root, ]
  dp_res <- list(dm = dm, pt = pt, ts = ts, M = M, ev = dm@eigenvectors, p1 = p1, branch = branch)
  return(dp_res)
}

message('test di_single_run_cpp')
library(InformationEstimator)
a <- Sys.time()
di_single_run(as.matrix(XYZ[, 1]), as.matrix(XYZ[, 2]), 10)
di_single_run(as.matrix(XYZ[, 1]), as.matrix(XYZ[, 3]), 10)
di_single_run(as.matrix(XYZ[, 2]), as.matrix(XYZ[, 3]), 10)
b <- Sys.time()
b - a 
# 11.56868 mins

# 0.00824910072759
# 0.00761374668366
# 0.00794952458964
 
message('test di_multiple_run_cpp')
library(InformationEstimator)
rdi_many_runs(as.matrix(XYZ[, 1:2]), as.matrix(XYZ[, 2:3]))
rdi_many_runs(as.matrix(XYZ[, 2:3]), as.matrix(XYZ[, 1:2]))
rdi_many_runs(as.matrix(XYZ[, c(1, 2)]), as.matrix(XYZ[, c(2, 3)]))

message('test di_single_run_conditioned_cpp')
library(InformationEstimator)
a <- Sys.time()
di_single_run_conditioned(as.matrix(XYZ[, 1]), as.matrix(XYZ[, 2]), as.matrix(XYZ[, 3]))
di_single_run_conditioned(as.matrix(XYZ[, 1]), as.matrix(XYZ[, 3]), as.matrix(XYZ[, 2]))
di_single_run_conditioned(as.matrix(XYZ[, 2]), as.matrix(XYZ[, 3]), as.matrix(XYZ[, 1]))
b <- Sys.time()
b - a

# Time difference of 13.76443 mins

message('test rdi_single_run_cpp')
library(InformationEstimator)
a <- Sys.time()
rdi_single_run(as.matrix(XYZ[, 1]), as.matrix(XYZ[, 2]), 1)
rdi_single_run(as.matrix(XYZ[, 1]), as.matrix(XYZ[, 3]), 1)
rdi_single_run(as.matrix(XYZ[, 2]), as.matrix(XYZ[, 3]), 1)
b <- Sys.time()
b - a
# 
# 0.0129632127358
# 0.00873069345257
# 0.0104163795216

# Time difference of 2.090866 mins

message('test rdi_single_run_conditioned_cpp')
library(InformationEstimator)
a <- Sys.time()
rdi_single_run_conditioned(as.matrix(XYZ[, 1]), as.matrix(XYZ[, 2]), as.matrix(XYZ[, 3]), 1L, 1L)
rdi_single_run_conditioned(as.matrix(XYZ[, 1]), as.matrix(XYZ[, 3]), as.matrix(XYZ[, 2]), 1L, 1L)
rdi_single_run_conditioned(as.matrix(XYZ[, 2]), as.matrix(XYZ[, 3]), as.matrix(XYZ[, 1]), 1L, 1L)
b <- Sys.time()
b - a

library(InformationEstimator)
test_res <- smooth_gene(XYZ) #  
qplot(test[, 1], test[, 2])

# compare with the pairwise distance 
library(InformationEstimator)
N <- 10000
a <- Sys.time()
calculate_rdi(XYZ[1:N, ], 1.0, 1) # 
b <- Sys.time()

library(InformationEstimator)
N <- 10000
a <- Sys.time()
calculate_rdi(XYZ[1:N, ], 1.0, 3L)
b <- Sys.time()

calculate_conditioned_rdi(XYZ) # 

# test the neuron simulation dataset (RDI + cRDI)
# test on the neuron simulation dataset (use three cells only first):
load('/Users/xqiu/Dropbox (Personal)/Projects/Monocle2_revision/RData/fig3.RData')
  
# avoid duplicated cells
exprs(nao_sim_cds)[2, 401] <- 0.001
exprs(nao_sim_cds)[3, 801] <- 0.001

# nao_sim_cds <- reduceDimension(nao_sim_cds, max_components = 3, norm_method = 'none', verbose = T, pseudo_expr = 0, scaling = F, auto_param_selection = F, initial_method = run_dpt, reduction_method = 'SimplePPT')
nao_sim_cds <- orderCells(nao_sim_cds)
plot_cell_trajectory(nao_sim_cds, color_by = 'State') + facet_wrap(~State)
plot_cell_trajectory(nao_sim_cds, x = 2, y = 3, color_by = 'State') + facet_wrap(~State)

nao_sim_cds <- orderCells(nao_sim_cds, root_state = 4)
dimnames(nao_exprs_mat) <- dimnames(nao_sim_cds)
nao_exprs_mat <- t(nao_exprs_mat)
neuron_branch <- nao_exprs_mat[pData(nao_sim_cds)$State %in% c(4, 3), ]
astrocyte_branch <- nao_exprs_mat[pData(nao_sim_cds)$State %in% c(4, 5), ]
oligodendrocyte_branch <- nao_exprs_mat[pData(nao_sim_cds)$State %in% c(4, 5), ]

save(file = '../../RData/simulation_data', neuron_branch, astrocyte_branch, oligodendrocyte_branch)

library(InformationEstimator)
library(monocle)

noise = matrix(rnorm(mean = 0, sd = 1e-10, 200 * 13), nrow = 200)
a <- Sys.time()
res <- calculate_rdi(neuron_branch[1:200, ] + noise, 1:3, 1) # run rdi 
b <- Sys.time()

a <- Sys.time()
RDI_parallel_res <- di::calculate_and_write_pairwise_dmi(neuron_branch[1:200, ], delays = 1:3, cores = detectCores() - 2, verbose = T)
b <- Sys.time()

lung <- load_lung()
a <- Sys.time()
RDI_parallel_res <- di::calculate_and_write_pairwise_dmi(t(exprs(lung)), delays = c(5, 10, 15), cores = detectCores() - 2, verbose = T)
b <- Sys.time()

noise = matrix(rnorm(mean = 0, sd = 1e-10, nrow(lung) * ncol(lung)), nrow = ncol(lung))
res <- calculate_rdi(t(exprs(lung)) + noise, c(5, 10, 15), 1) # run rdi 

a <- Sys.time()
RDI_parallel_res <- calculate_rdi(t(exprs(lung)), delays = c(5, 10, 15), cores = detectCores() - 2, verbose = T)
b <- Sys.time()

*/
// 
// // [[Rcpp::export]]
// RcppExport SEXP loglik_MP(SEXP s_beta, SEXP s_x, SEXP s_nCpu) {
//   NumericVector x(s_x);
//   NumericVector beta(s_beta);
//   int n_cpu = IntegerVector(s_nCpu)[0];
//   double mu = beta[0];
//   double sigma = beta[1];
//   double ll = 0;
//   omp_set_dynamic(0);         // Explicitly disable dynamic teams
//   omp_set_num_threads(n_cpu); // Use n_cpu threads for all
//   // consecutive parallel regions
//   #pragma omp parallel
//   {
//     double ll_thread = 0;
//   #pragma omp for 
//     for(int i = 0; i < x.length(); i++) {
//       ll_thread -= (x[i] - mu)*(x[i] - mu);
//     }
//   #pragma omp critical
//   {
//     ll += ll_thread;
//   }
//   }
//   ll *= 0.5/sigma/sigma;
//   ll -= (0.5*log(2*M_PI) + log(sigma))*x.length();
//   NumericVector result(1, ll);
//   return result;
// }

// // Convert a NumericMatrix (R object) to an Armadillo matrix (C++ object)

// // [[Rcpp::export]]
// NumericMatrix calculate_rdi_cpp(NumericMatrix expr_data, const NumericVector delays, const int n_cores)
// {
// 
//   // Convert arguments to Armadillo objects
//   const arma::mat data = R2armaMat_num(expr_data); // Continuous data does not need reindexing R->C++
//   const int n_genes(data.n_cols), n_samples(data.n_rows);
//   const int n_pairs = (int)(((double)n_genes/2.0) * (double)(n_genes+1));
// 
//   // // Number of cores to use
//   // //const int n_cores = sysconf(_SC_NPROCESSORS_ONLN);
//   // omp_set_num_threads(n_cores);
// 
//   // Compute RDI in parallel
//   int delays_len = length(delays); 
//   NumericMatrix RDI = arma::mat(n_genes, n_genes * delays_len, arma::fill::zeros);
//   
//   // #pragma omp parallel for shared(expr_data,delays,RDI) private(i,j,k) schedule(dynamic) default(none) collapse(2)
//   // this part maybe slow because RDI matrix is huge 
//   for(int i(0); i<n_genes; ++i)
//   {
//     for(int j(0); j<n_genes; ++j)
//     {
//      for(int k(0); j < delays_len; ++ k)
//      {
//        if(i == j) 
//        {
//          RDI[i, i] = 0; 
//        }
//        RDI(i, j + k * n_genes) = rdi_single_run_cpp(expr_data(, i), expr_data(, j), delays[k]); // how to deal with delays include multiple values?
//      }
//     }
//   }
//   
//   return RDI;
// }

// [[Rcpp::plugins(openmp)]]

// find the max delay for each pair

// find top k incoming node: 

// perform the conditioinal RDI calculation 
// prepare the dataset in R before pass to C++ 
// 
// // [[Rcpp::export]]
// NumericMatrix calculate_conditioned_rdi_cpp(NumericMatrix expr_data, CharacterMatrix super_graph, 
//  NumericMatrix top_k_plus_1_incoming_id, NumericMatrix rdi_res, const int k = 1, const int n_cores, 
//  const bool verbsoe = TRUE)
// {
//  if(verbsoe == TRUE) Rcout << "Calculating the conditional restricted direct mutual information for each pair of genes" << std::endl;
//  k = min(k, ncol(expr_data) - 2); // minal value k and the column 
// 
//  NumericMatrix cRDI_mat(expr_data); 
//  for(int i = 0; i < super_graph.nrow(); i ++) 
//  {
//    NumericVector index_name = super_graph(i, _);
//    NumericVector gene_pair = super_graph(i, _); 
//    NumericVector top_k_plus_1 = top_k_plus_1_incoming_id(gene_pair[2], ); 
//    NumericVector valid_top_k = setdiff(top_k_plus_1, gene_pair[1]); 
// 
//    NumericMatrix expr1_t = genes_data(_, index_name[0]); 
//    NumericMatrix expr2_t = genes_data(_, index_name[1]); 
//    NumericMatrix past; 
//    for(int id = 0; id < length(top_k_gene_delay); id ++)
//    {
//      past = rbind(past, genes_data(_, top_k_plus_1[id])); 
//    }
// 
//         cRDI_mat(i, j) = cmi(expr1_t, expr2_t, past); 
//  }
// 
//  // return results
// 
//  return cRDI_mat; 
// }
// 
// 


// the non-deterministic crash bug (address 0x0, cause 'unknown') fixed 

// [[Rcpp::export]]
double rdi_multiple_run_cpp(NumericMatrix& x, NumericMatrix& y, int d = 1, IntegerVector run_vec = 0, bool uniformalize = false) 
{
  int Nx = x.rows(); int Ny = y.rows(); int run_len; // int Cx = x.cols(); int Cy = y.cols(); 
  
  if(Nx != Ny)
  {
    stop("The number of time-series samples has to be the same for X and Y");
    return -1;
  }

  int run_max = max(run_vec); // run_vec: vector to represent the run, need to be just numeric numbers and starts from 0. 
  int current_ind = 0; // current initial index for storing new slice of data 

  NumericMatrix x_0_res, y_0_res, y_1_res, tmp; // matrix passed to cmi function
  IntegerVector index_all = seq_len(run_vec.size()) - 1, current_run_ind; // cell index from 0 to maximal number of cells; run index corresponding to current run during iteration 
  
  if(run_max == 0) 
  {
    x_0_res = x(Range(0, Nx - d - 1), _); // x should only has 1 colum 
    y_0_res = y(Range(d, Ny - 1), _);
    y_1_res = y(Range(d - 1, Nx - 2), _); 
  } else {

    // correctly set the dimensionality for the x_0, y_0 and y_1 matrix 
    int tot_len = Nx - ( d * (run_max + 1) ); // number of samples minus, the delay times total number of runs 

    mat x_0(tot_len, 1), y_0(tot_len, 1), y_1(tot_len, 1); // arma matrix or colvec
    
    uvec pos; // arma uvec 
    vec vals; // arma vec 

    IntegerVector rng; 

    for(int i = 0; i <= run_max; i++) // concatenate the data 
    {
      current_run_ind = index_all[run_vec == i]; // get the cells that belong to a particular run 
      run_len = current_run_ind.size(); 

      if(run_len < d)
      {
        stop("Each run has to have more samples than the designated 'delays'");
        return -1;
      }

      rng = Range(current_ind, current_ind + run_len - d - 1);
      
      pos = as<uvec>( rng ); 

      // assuming x is only one column (Cx, Cy is not used)
      tmp = x(Range(min(current_run_ind), min(current_run_ind) + run_len - d - 1), Range(0, 0)); //x(Range(1, 2), Range(1, 1));
      vals = as<vec>( tmp ); 
      x_0.elem(pos) = vals; 

      IntegerVector rng_vec = Range(min(current_run_ind), min(current_run_ind) + run_len - d - 1); 

      // for(int test_na = 0; test_na < vals.n_elem; test_na ++) 
      // {
      //   if(arma::is_finite(vals[test_na]) == false) 
      //   {
      //     Rcout << "identify non-finite values for x_0 here; i is " << i << " rng_vec is " << rng_vec << std::endl; 
      //   }
      // }

      tmp = y(Range(min(current_run_ind) + d, min(current_run_ind) + run_len - 1), Range(0, 0)); //x(Range(1, 2), Range(1, 1));
      vals = as<vec>( tmp ); 
      y_0.elem(pos) = vals; 

      // for(int test_na = 0; test_na < vals.n_elem; test_na ++) 
      // {
      //   if(arma::is_finite(vals[test_na]) == false) 
      //   {
      //     Rcout << "identify non-finite values for x_0 here; i is " << i << " rng_vec is " << rng_vec << std::endl; 
      //   }
      // }

      tmp = y(Range(min(current_run_ind) + d - 1, min(current_run_ind) + run_len - 2), Range(0, 0)); //x(Range(1, 2), Range(1, 1));

      vals = as<vec>( tmp ); 
      y_1.elem(pos) = vals; 

      // for(int test_na = 0; test_na < vals.n_elem; test_na ++) 
      // {
      //   if(arma::is_finite(vals[test_na]) == false) 
      //   {
      //     Rcout << "identify non-finite values for x_0 here; i is " << i << " rng_vec is " << rng_vec << std::endl; 
      //   }
      // }

      current_ind = current_ind + run_len - d; // move to the the next position after filling the current run (note that there is no - 1)
    }

    x_0_res = as<NumericMatrix>(wrap(x_0)); // this conversion have problems? 
    y_0_res = as<NumericMatrix>(wrap(y_0));
    y_1_res = as<NumericMatrix>(wrap(y_1));

  }  

  // for(int test_na = 0; test_na < y_1_res.size(); test_na ++) 
  // {
  //   if(arma::is_finite(y_1_res[test_na]) == false) 
  //   {
  //     Rcout << "identify non-finite values for y_1_res here" << std::endl; 
  //   }
  // }

  if(uniformalize == true) {
    List ucmi_res = ucmi_cpp(x_0_res, y_0_res, y_1_res, 5, 1, 0, 0); // k, method, k_density, bw
    return ucmi_res["ucmi_res"]; 
  } else {
    List cmi_res = cmi_cpp(x_0_res, y_0_res, y_1_res);
    return cmi_res["cmi_res"]; 
  } 
}

// [[Rcpp::export]]
List calculate_rdi_multiple_run_cpp(NumericMatrix& expr_data, IntegerVector delays, IntegerVector run_vec,
                                    IntegerMatrix& super_graph, IntegerVector turning_points = 0, int method = 1, bool uniformalize = false) //, method: 1 rdi; 2: lmi 
{
  const int n_genes(expr_data.cols()), n_samples(expr_data.rows());

  if(max(delays) > n_samples)
  {
    stop("Number of samples has to be larger than the delays"); 
    return -1; 
  }
  
  IntegerMatrix max_rdi_delays(n_genes, n_genes);
  
  // Convert arguments to Armadillo objects
  // const arma::mat data = R2armaMat_num(expr_data); // Continuous data does not need reindexing R->C++
  // Rcout << "n_genes is " << n_genes << " n_samples is " << n_samples << " delay length is " << delays.length() << std::endl;
  // const int n_pairs = (int)(((double)n_genes/2.0) * (double)(n_genes+1));
  
  // // Number of cores to use
  // //const int cores = sysconf(_SC_NPROCESSORS_ONLN);
  // omp_set_num_threads(cores);
  
  // Compute RDI in parallel
  int delays_len = delays.length(); 
  
  int RDI_ncols;  
  if(turning_points.length() == n_genes) 
  {
    RDI_ncols = n_genes; 
  } else {
    RDI_ncols = n_genes * delays_len; 
  }

  NumericMatrix RDI(n_genes, RDI_ncols);    
  std::fill(RDI.begin(), RDI.end(), 0); 
  NumericMatrix expr_1(n_samples, 1), expr_2(n_samples, 1);

  // this part maybe slow if we do paralleling because RDI matrix is huge 
  
  int i, j, k; 
  IntegerVector current_pair; 
  // #pragma omp parallel for shared(n_genes, expr_data, delays_len, RDI) private(i, j, k, expr_1, expr_2, _) //schedule(dynamic) default(none) //collapse(2) , _
  for(int super_graph_ind = 0; super_graph_ind < super_graph.rows(); super_graph_ind ++) 
  {
    current_pair = super_graph(Range(super_graph_ind, super_graph_ind), _);
    i = current_pair[0]; j = current_pair[1]; 
    
    expr_1 = expr_data(_, Range(i, i)); //Range(0, n_genes - 1)
    expr_2 = expr_data(_, Range(j, j)); //Range(0, n_genes - 1)

    // int max_thread = omp_get_max_threads();
    // omp_set_num_threads(max_thread);
    // #pragma omp parallel for shared(delays_len, i, j, n_genes, expr_1, expr_2, RDI) private(k) //schedule(dynamic) default(none) //collapse(2) , _
    if(turning_points.length() == n_genes) //!Rf_isNull(turning_points) 
    {
      Rcout << "using user provided information about time-delay " << turning_points.length() << std::endl;
      int current_delay = turning_points[i] - turning_points[j];

      if(i == j) 
      {
        continue; 
      }
      
      if(method == 1)
      {
        RDI(i, j) = rdi_multiple_run_cpp(expr_1, expr_2, current_delay, run_vec, uniformalize); // how to deal with delays include multiple values?
      }
      else if(method == 2)
      {
        RDI(i, j) = lmi_multiple_run_cpp(expr_1, expr_2, current_delay, run_vec, uniformalize); // how to deal with delays include multiple values?
      } 
      max_rdi_delays(i, j) = current_delay;
    } else {
      for(k = 0; k < delays_len; ++ k)
      {
        if(i == j) 
        {
          continue; // keep diagnoal as 0 
        }
        
        if(method == 1)
        {
          RDI(i, j + k * n_genes) = rdi_multiple_run_cpp(expr_1, expr_2, delays[k], run_vec, uniformalize); // how to deal with delays include multiple values?
        }
        else if(method == 2)
        {
          RDI(i, j + k * n_genes) = lmi_multiple_run_cpp(expr_1, expr_2, delays[k], run_vec, uniformalize); // how to deal with delays include multiple values?
        }
      }
    }
  }

  // perform the multiple run test (add a new argument to represent the vector for different runs)
  if(turning_points.length() == n_genes) // Rf_isNull(turning_points) 
  {  
    return List::create(Rcpp::Named("RDI") = RDI, 
                        Rcpp::Named("delays") = R_NilValue,
                        Rcpp::Named("max_rdi_value") = RDI, 
                        Rcpp::Named("max_rdi_delays") = max_rdi_delays);
  }
  else{
    // run extract_max_rdi_value_delay 
    List max_rdi_value_delay_list = extract_max_rdi_value_delay(RDI, delays); 
    
    return List::create(Rcpp::Named("RDI") = RDI, 
                        Rcpp::Named("delays") = delays,
                        Rcpp::Named("max_rdi_value") = max_rdi_value_delay_list["max_rdi_value"], 
                        Rcpp::Named("max_rdi_delays") = max_rdi_value_delay_list["max_rdi_delays"]);
  }
}

//' @title
//' calculate_rdi_multiple_run_cpp_wrap
//' @description
//' This function estimates the DIRECTED mutual information for all gene pairs in the expr_data matrix when you have a SINGLE run of the processes
//' 
//' @param expr_data a matrix for all variables in a time-series where each row is a time point and each column is a variable (for example, a gene). 
//' The rows are ordered according to time, from earliest to latest. 
//' 
//' @param delays An integer vector storing the time delays between pairs of variables you would like to try. 
//' 
//' @param super_graph An integer matrix where each row is the variable IDs or gene pairs you would like to estimate restricted direct information (RDI) from. 
//'  
//' @param turning_points Either 0 or a numeric vector describing the inflection point (linear trajectory) or branch point (bifurcation point) for each gene. 
//' If the turning_point for each gene is provided, the time delay will be estimated based on the turning point. 
//'
//' @param method An integer of either 1 or 2 to determine which information metric will be used to quantify the causality. 
//' If method is 1, then lagged mutual information will be used; if method is 2, then the restricted direct information will be used. 
//'
//' @param uniformalize Whether or not you want to use ucmi to calculate rdi. Default to be false. 
//' @details
//' \code{calculate_rdi_multiple_run_cpp_wrap}, similar to calculate_rdi_cpp_wrap, takes an expression matrix (expr_data) and possible gene pairs (encoded in super_graph) to estimate the restricted 
//' direct information based on the time delay which can be estimated from the turning_point vector. It, however, differs to calculate_rdi_cpp_wrap, in that it can concatenate different experiments (runs) into a single run of the data for causality estimation. 
//' @return a numeric matrix of conditional RDI values for all possible pairs of genes from expr_data. If the gene pairs is not encoded in the super_graph, it will remain as 0 in the matrix.  
//' @export
// [[Rcpp::export]]
List calculate_rdi_multiple_run_cpp_wrap(SEXP expr_data, SEXP delays, SEXP run_vec, SEXP super_graph, SEXP turning_points, SEXP method, SEXP uniformalize) //, SEXP R_cores
{ 
  NumericMatrix expr_data_cpp(expr_data); 
  IntegerVector delays_cpp(delays); 
  IntegerVector run_vec_cpp(run_vec); 
  
  IntegerMatrix super_graph_cpp(super_graph); 
  IntegerVector turning_points_cpp(turning_points);
  int method_cpp = as<int>(method); 
  bool uniformalize_cpp = as<bool>(uniformalize); 
  
  List rdi_list = calculate_rdi_multiple_run_cpp(expr_data_cpp, delays_cpp, run_vec_cpp, super_graph_cpp, turning_points_cpp, method_cpp, uniformalize_cpp); //cores

  return rdi_list;
}

/** 
 library(InformationEstimator)
 sample_row = 30; sample_col = 50; 
 noise = matrix(rnorm(mean = 0, sd = 1e-10, sample_row * sample_col), nrow = sample_col)

 test <- calculate_rdi(t(lung_exprs_AT1[1:30, 1:50]) + noise, c(5, 10, 15))
 */

/** 
 library(InformationEstimator)
 sample_row = 30; sample_col = 50; 
 noise = matrix(rnorm(mean = 0, sd = 1e-10, sample_row * sample_col), nrow = sample_col)

  sample_row = 123; sample_col = 218; 
  set.seed(2017)
  noise = matrix(rnorm(mean = 0, sd = 1e-10, sample_row * sample_col), nrow = sample_row)
  a <- Sys.time()
  rdi_list <- calculate_rdi(lung_exprs_AT1[1:sample_row, 1:sample_col] + noise, c(5, 10, 15), method = 1)
  b <- Sys.time()
  rdi_time <- b - a 
    
  top_delays <- extract_top_incoming_nodes_delays(rdi_list$max_rdi_value, rdi_list$max_rdi_delays, k = 1)
     
  tmp <- expand.grid(1:sample_col, 1:sample_col, stringsAsFactors = F)
  super_graph <- tmp[tmp[, 1] != tmp[, 2], ] - 1 # convert to C++ index 
   
 // SEXP expr_data, SEXP delays, SEXP run_vec, SEXP super_graph, SEXP method
  test <- calculate_rdi_multiple_run_cpp_wrap(t(lung_exprs_AT1[1:30, 1:50]) + noise, 
            c(5, 10, 15), rep(c(1, 2), each = 25), super_graph, method)
*/

//-------------------------------------------------------------------------------------------------------------------
/*
calculate conditional restricted direction information for a single run
*/
//-------------------------------------------------------------------------------------------------------------------

// [[Rcpp::export]]
double rdi_multiple_runs_conditioned_cpp(NumericMatrix& x, NumericMatrix& y, 
                                         NumericMatrix& z, IntegerVector& z_delays, int d = 1, IntegerVector run_vec = 0, bool uniformalize = false) 
{
  uniformalize = false; // we cannot support kernel density estimator > 2 dim 
  
  int Nx = x.rows(); int Ny = y.rows(); int Nz = z.rows(); int run_len; // int Cx = x.cols(); int Cy = y.cols(); int Cz = z.cols(); 
  if(Nx != Ny | Nx != Nz)
  {
    stop("The number of time samples has to be the same for X and Y or Z");
    return -1;
  }
  
  z_delays.push_back(d);
  int tau = max(z_delays); 

  int run_max = max(run_vec); // run_vec: vector to represent the run, need to be just numeric numbers and starts from 0. 
  int current_ind = 0; // current initial index for storing new slice of data (--current_incoming_k)

  NumericMatrix x_0_res, y_0_res, yz_res, tmp; // matrix passed to cmi function
  IntegerVector index_all = seq_len(run_vec.size()) - 1, current_run_ind; // cell index from 1 to maximal number of cells; cell index corresponding to current run during iteration 

  if(run_max == 0) // note that run id starts from 0 
  {
    int tot_len = x.rows() - tau;  
    
    yz_res = y(Range(tau - 1, tau - 2 + tot_len), _);
    
    NumericMatrix tmp; 
    
    // OpenMP here too 
    // int max_thread = omp_get_max_threads();
    // omp_set_num_threads(max_thread);
    // #pragma omp parallel for shared(z, tau, z_delays, tot_len) private(j, tmp, yz) //schedule(dynamic) default(none) //collapse(2) , _
      
    // using cbind
    for(int j = 0; j < z.cols(); j ++ )
    {
      tmp = z(Range(tau - z_delays[j], tau - z_delays[j] + tot_len - 1), Range(j, j));
      yz_res = cbind(tmp, yz_res);
    } 
    
    x_0_res = x(Range(tau - d, tau - d + tot_len - 1), _);
    y_0_res = y(Range(tau, tau + tot_len - 1), _);

  } else {
    mat x_0, y_0; // pre-define the size of the matrix: (1000, 2) 
    
    uvec pos; //unsigned vector 
    vec vals; 

    int tot_len = Nx - tau * (run_max + 1); 

    mat yz(tot_len, z.cols() + 1); // tot_len: number of cells in total across all runs; z.col() + 1: number of top k incoming nodes and y vector 

    for(int i = 0; i <= run_max; i++) // concatenate the data 
    {
      current_run_ind = index_all[run_vec == i]; // get the cells that belong to a particular run 
      run_len = current_run_ind.size();

      if(current_run_ind.size() < d)
      {
        stop("Each run has to have more samples than the designated 'delays'");
        return -1;
      }
      tot_len = run_len - tau; // x.rows() 
      
      // prepare the value for yz (only the y values)
      NumericMatrix locs_tmp(2, tot_len); // first row: row index; second row: column index; each column is a coordinate
      locs_tmp(0, _) = Range(current_ind, current_ind + tot_len - 1);
      // we could change the index below to the maximal value 
      locs_tmp(1, _) = rep(0, tot_len); 
      umat locs = as<umat>( locs_tmp );  
      uvec eids = sub2ind( size(yz), locs ); // Obtain Element IDs
            
      // pos = as<uvec>(Range(current_ind, current_ind + tot_len - 1)); 
      tmp = y(Range(min(current_run_ind) + tau - 1, min(current_run_ind) + tau - 2 + tot_len), Range(0, 0)); // x(Range(min(current_run_ind), min(current_run_ind) + Nx - d - 1), Range(1, 1)); //x(Range(1, 2), Range(1, 1));

      vals = as<vec>( tmp ); 
      yz.elem( eids ) = vals;
           
      // OpenMP here too 
      // int max_thread = omp_get_max_threads();
      // omp_set_num_threads(max_thread);
      // #pragma omp parallel for shared(z, tau, z_delays, tot_len) private(j, tmp, yz) //schedule(dynamic) default(none) //collapse(2) , _
        
      // we could change the index below 
      for(int j = 0; j < z.cols(); j ++ )
      {
        locs_tmp(1, _) = rep(j + 1, tot_len); // first column belong to y values so starting from the second column 
        umat locs = as<umat>( locs_tmp );  

        uvec eids = sub2ind( size(yz), locs ); // Obtain Element IDs
        tmp = z(Range(min(current_run_ind) + tau - z_delays[j], min(current_run_ind) + tau - z_delays[j] + tot_len - 1), Range(j, j));
        vals = as<vec>( tmp ); 
        
        yz.elem( eids ) = vals;        
      } 

      tmp = x(Range(min(current_run_ind) + tau - d, min(current_run_ind) + tau - d + tot_len - 1), Range(0, 0));
      vals = as<vec>( tmp );
      x_0 = join_cols(x_0, vals); //.elem(pos) = vals;
      
      tmp = y(Range(min(current_run_ind) + tau, min(current_run_ind) + tau + tot_len - 1), Range(0, 0));      
      vals = as<vec>( tmp );
      y_0 = join_cols(y_0, vals); //.elem(pos) = vals;    

      current_ind += tot_len; // update the index (note that there is no - 1)
    }

    x_0_res = as<NumericMatrix>(wrap(x_0));
    y_0_res = as<NumericMatrix>(wrap(y_0));
    yz_res = as<NumericMatrix>(wrap(yz)); 

  }
  
  if(uniformalize == true) {
    List ucmi_res = ucmi_cpp(x_0_res, y_0_res, yz_res, 5, 1, 0, 0); // k, method, k_density, bw
    return ucmi_res["ucmi_res"]; 
  } else {
    List cmi_res = cmi_cpp(x_0_res, y_0_res, yz_res);
    return cmi_res["cmi_res"]; 
  } 
}

//' @title
//' This function simulates the CONDITIONED DIRECTED mutual information from X to Y CONDITIONED ON Z when you have a SINGLE run of the processes
//' @description
//' This subroutine calculates the volume of a d-dimensional unit ball for Euclidean norm
//' 
//' @param x one random variable from the time-series data
//' 
//' @param y another random variable from the time-series data
//' 
//' @param z z is a dataframe or matrix consisting of the data for different variables
//' 
//' @param z_delays z_delay is also a dataframe or matrix consisting of the delays to be applied to different variables
//' 
//' @param d delay in the formula I(x-->y)=I(x_{t-d};y_t|y_{t-1}) default: d=1
//'
//' @param uniformalize Whether or not you want to use ucmi to calculate rdi. Default to be false. 
//' @details
//' \code{rdi_single_run_conditioned} takes two random variables x and y as well as the parameter n to calculate the restricted direct information conditioned on variable z. 
//' @return a matrix for the condition mutual information estimators between all pairwise variables (x, y) in the data matrix x, y
//' @export
// [[Rcpp::export]]
double rdi_multiple_runs_conditioned(SEXP x, SEXP y, SEXP z, SEXP z_delays, SEXP d, SEXP run_vec, SEXP uniformalize) 
{ 
  NumericMatrix x_cpp(x); 
  NumericMatrix y_cpp(y); 
  NumericMatrix z_cpp(z); 
  
  IntegerVector z_delays_cpp(z_delays); 
  IntegerVector run_vec_cpp(run_vec); 
  
  int d_cpp = as<int>(d); 
  bool uniformalize_cpp = as<bool>(uniformalize); 
  double rdi_single_run_conditioned_res = rdi_multiple_runs_conditioned_cpp(x_cpp, y_cpp, z_cpp, z_delays_cpp, d_cpp, run_vec_cpp, uniformalize_cpp);

  return rdi_single_run_conditioned_res;
}

// test all the above functions in the following 

/** 
 library(InformationEstimator)
 sample_row = 30; sample_col = 50; 
 noise = matrix(rnorm(mean = 0, sd = 1e-10, sample_row * sample_col), nrow = sample_col)

  sample_row = 123; sample_col = 218; 
  set.seed(2017)
  noise = matrix(rnorm(mean = 0, sd = 1e-10, sample_row * sample_col), nrow = sample_row)
  a <- Sys.time()
  rdi_list <- calculate_rdi(lung_exprs_AT1[1:sample_row, 1:sample_col] + noise, c(5, 10, 15), method = 1)
  b <- Sys.time()
  rdi_time <- b - a 
    
  top_delays <- extract_top_incoming_nodes_delays(rdi_list$max_rdi_value, rdi_list$max_rdi_delays, k = 1)
     
  tmp <- expand.grid(1:sample_col, 1:sample_col, stringsAsFactors = F)
  super_graph <- tmp[tmp[, 1] != tmp[, 2], ] - 1 # convert to C++ index 
   
 // SEXP expr_data, SEXP delays, SEXP run_vec, SEXP super_graph, SEXP method
  test <- rdi_multiple_runs_conditioned_cpp(t(lung_exprs_AT1[1:30, 1:50]) + noise, 
            c(5, 10, 15), rep(c(1, 2), each = 25), super_graph, method)
*/

//-------------------------------------------------------------------------------------------------------------------
/*
calculate conditional restricted direction information for all genes in a expression matrix 
*/
//-------------------------------------------------------------------------------------------------------------------

// [[Rcpp::export]]
NumericMatrix calculate_conditioned_rdi_multiple_run_cpp(NumericMatrix& expr_data, IntegerMatrix& super_graph, 
                                            NumericMatrix& max_rdi_value, IntegerMatrix& max_rdi_delays, IntegerVector run_vec = 0, int k = 1, bool uniformalize = false) //, const int cores, const bool verbose
{
  uniformalize = false; // we cannot support kernel density estimator > 2 dim 
  
  if(expr_data.cols() < 3) {
    stop("You need at least 3 genes to calculate conditional RDI values");
  }
  
  NumericVector k_ncol = NumericVector::create(k, expr_data.cols() - 2); // conditioned at most n - 2 genes (2: two current testing genes)
  k = min(k_ncol); // minal value from k and the column 
  
  List top_incoming_nodes_delays_list = extract_top_incoming_nodes_delays(max_rdi_value, max_rdi_delays, k); // get the top k incoming node and the corresponding delays 

  IntegerMatrix top_incoming_nodes = top_incoming_nodes_delays_list["top_incoming_nodes"];
  IntegerMatrix top_incoming_delays = top_incoming_nodes_delays_list["top_incoming_delays"]; 
  NumericMatrix top_incoming_values = top_incoming_nodes_delays_list["top_incoming_values"];  

  int n_genes = expr_data.cols(), n_sample = expr_data.rows(), current_id_tmp, delay; // tau, total_length, id3, delay_id3, tmp: current target  

  // valid_top_k_in_k_ordered is used to get the valid_top_k to the original order (setdiff function annoyingly sorts the data randomly)
  IntegerVector gene_pair(2), tmp(1), valid_top_k_in_k(k), top_k_plus_1(k + 1), valid_top_k(k), valid_top_k_in_k_ordered, valid_delay_tmp(1); // ind_for_top_k(k), little computational gain by setting the correct dimension  

  NumericMatrix cRDI_mat(n_genes, n_genes), expr1_t, expr2_t, past_tmp; 
  std::fill(cRDI_mat.begin(), cRDI_mat.end(), 0); // initialize the cRDI matrix as 0 instead of the C ++ NA values (so the diagnoal will be 0 values)
  
  // OpenMP here too (not thread safe) 
  for(int i = 0; i < super_graph.rows(); i ++) 
  {
    gene_pair = super_graph(i, _); 
    top_k_plus_1 = top_incoming_nodes(Range(gene_pair[1], gene_pair[1]), _); // top k incoming nodes 
    tmp = gene_pair[0]; // get the incoming node 

    if(intersect(top_k_plus_1, tmp).length() > 0) // if current node tmp is in the top_k_plus_1, remove tmp; else remove the smallest one 
    {
      valid_top_k = setdiff(top_k_plus_1, tmp);
    }
    else
    {
      tmp = top_k_plus_1[k]; // remove the smallest one (k + 1 points in total)
      valid_top_k = setdiff(top_k_plus_1, tmp); 
    }
    
    // assign the node indices to ind_for_top_k
    valid_top_k_in_k = match(valid_top_k, top_k_plus_1) - 1; // match gives order and the index (1 to k) starts from 1 -> valid_top_k_in_k
    valid_top_k_in_k_ordered = clone(valid_top_k_in_k).sort(); // the index for top_k_plus_1 is sorted
    valid_top_k = top_k_plus_1[valid_top_k_in_k_ordered]; // now the valid_top_k's value follows the ordering in the top_k_plus_1 

    if(valid_top_k.size() < k) 
    { // avoid any weird rare situation (all incoming RDI value are 0)
      cRDI_mat(gene_pair[0], gene_pair[1]) = max_rdi_value(gene_pair[0], gene_pair[1]); 
      continue;
    }
    
    IntegerVector valid_top_k_incoming_delays(k); // re-initialize this everytime because we have push_back later

    // assign the valid top_k_incoming delays based on the sorted valid_top_k_in_k_ordered index 
    for(int j = 0; j < k; j ++ )
    {
      current_id_tmp = valid_top_k_in_k_ordered[j]; // range: 0 -> k - 1
      valid_delay_tmp = top_incoming_delays(Range(gene_pair[1], gene_pair[1]), Range(current_id_tmp, current_id_tmp)); // the current_id_tmp^{th} delay value 
      valid_top_k_incoming_delays[j] = valid_delay_tmp[0]; 
    }

    delay = max_rdi_delays(gene_pair[0], gene_pair[1]); // the delay corresponds to max RDI value for the current gene pair 
 
    // use rdi_conditioned to calculate the values: 
    NumericMatrix x = expr_data(_, Range(gene_pair[0], gene_pair[0]));
    NumericMatrix y = expr_data(_, Range(gene_pair[1], gene_pair[1]));

    // prepare the data for all top-k incoming nodes' expression data 
    NumericMatrix z(n_sample, k);
    for (int i = 0; i < valid_top_k.size(); i++) {
      z(_,i) = expr_data(_, valid_top_k(i));
    }

    cRDI_mat(gene_pair[0], gene_pair[1])  = rdi_multiple_runs_conditioned_cpp(x, y, z, valid_top_k_incoming_delays, delay, run_vec, uniformalize); 
  }
  
  // return results
  return cRDI_mat; 
}

//' @title
//' calculate_conditioned_rdi_multiple_run_wrap
//' @description
//' This function estimates the conditional DIRECTED mutual information for all gene pairs in the expr_data matrix when you have a SINGLE run of the processes.
//' 
//' @param expr_data a matrix for all variables in a time-series where each row is a time point and each column is a variable (for example, a gene). 
//' The rows are ordered according to time, from earliest to latest. 
//' 
//' @param super_graph An integer matrix where each row is the variable IDs or gene pairs you would like to estimate restricted direct information (RDI) from. 
//'  
//' @param max_rdi_value A numeric matrix where each element corresponding to the maximal (identified from a series of time lags when calculating the rdi values) rdi value between two variables.  
//' 
//' @param max_rdi_delays An integer matrix where each element corresponding to the time delay corresponding to the maximal (identified from a series of time lags when calculating the rdi values) rdi value between two variables.  
//'
//' @param run_vec An integer vector keeping the run id for each sample. 
//'
//' @param k An integer for the number of incoming nodes to be conditioned on.
//'
//' @param uniformalize Whether or not you want to use ucmi to calculate rdi. Default to be false. 
//' @details
//' \code{calculate_conditioned_rdi_multiple_run_wrap}, similar to calculate_conditioned_rdi_cpp_wrap, takes an expression matrix (expr_data) and possible gene pairs (encoded in super_graph), as well as the matrices of maximal rdi value or the delays corresponding to those values to 
//' estimate the conditional restricted direct information, conditioning on top k incoming nodes. It, however, differs to calculate_conditioned_rdi_cpp_wrap, in that it can concatenate different experiments (runs) into a single run of the data for causality estimation. 
//' @return a numeric matrix of conditional RDI values for all possible pairs of genes from expr_data. If the gene pairs is not encoded in the super_graph, it will remain as 0 in the matrix.  
//' @export
// [[Rcpp::export]]
NumericMatrix calculate_conditioned_rdi_multiple_run_wrap(SEXP expr_data, SEXP super_graph, SEXP max_rdi_value, SEXP max_rdi_delays, SEXP run_vec, SEXP k, SEXP uniformalize) //, SEXP R_cores
{ 
  NumericMatrix expr_data_cpp(expr_data), max_rdi_value_cpp(max_rdi_value); 
  IntegerMatrix super_graph_cpp(super_graph), max_rdi_delays_cpp(max_rdi_delays); 

  IntegerVector run_vec_cpp(run_vec); 
  int k_cpp = as<int>(k); 
  bool uniformalize_cpp = as<bool>(uniformalize); 
  
  NumericMatrix crdi_res = calculate_conditioned_rdi_multiple_run_cpp(expr_data_cpp, super_graph_cpp, max_rdi_value_cpp, max_rdi_delays_cpp, run_vec_cpp, k_cpp, uniformalize_cpp); //cores

  return crdi_res;
}

// implement UMI here 


/*** R
rm(list = ls())

 library(devtools)
 load_all()
 library(monocle)
 library(plyr)
 lung <- load_lung()
 lung <- buildBranchCellDataSet(lung, progenitor_method = 'duplicate')
 data <- t(exprs(lung))
 run_vec <- as.numeric(revalue(pData(lung)$Branch, c("Y_1" = 1, "Y_30" = 2)))
 noise = matrix(rnorm(mean = 0, sd = 1e-1, nrow(data) * ncol(data)), nrow = nrow(data))
 sample_gene <- 10; sample_cell <- 200
 data_sampled <- data[1:sample_cell, 1:sample_gene] + noise[1:sample_cell, 1:sample_gene]
 a <- Sys.time()
 rdi_list <- calculate_rdi(data_sampled, c(1, 2, 3), method = 1)
 b <- Sys.time()
 rdi_time <- b - a 
 top_delays <- extract_top_incoming_nodes_delays(rdi_list$max_rdi_value, rdi_list$max_rdi_delays, k = 1)
 tmp <- expand.grid(1:sample_gene, 1:sample_gene, stringsAsFactors = F)
 super_graph <- tmp[tmp[, 1] != tmp[, 2], ] - 1 # convert to C++ index 
 library(devtools)
 load_all()
 library(monocle)
 library(plyr)
 con_rid_res_test <- calculate_multiple_run_conditioned_rdi_wrap(data_sampled, as.matrix(super_graph), as.matrix(rdi_list$max_rdi_value), as.matrix(rdi_list$max_rdi_delays), run_vec[1:sample_cell] - 1, 1)
 test_extract <- extract_top_incoming_nodes_delays(rdi_list$max_rdi_value, rdi_list$max_rdi_delays, k = 1)
 calculate_rdi_multiple_run_cpp_wrap(data_sampled, c(1, 2, 3), run_vec[1:sample_cell] - 1, as.matrix(super_graph), 1)
 duplicated()
 calculate_rdi_multiple_run_cpp_wrap(data_sampled, 1:3, rep(0, length(run_vec)), as.matrix(super_graph), 1)

library(Scribe)
library(devtools)
load_all('~/Dropbox (Personal)/Projects/Causal_network/causal_network/Cpp/di-grn/Scribe')
library(monocle)

load("/Users/xqiu/Dropbox (Cole Trapnell's Lab)/Monocle 2/first_revision/Monocle2_revision/RData/fig3.RData") # cds data
load('/Users/xqiu/Dropbox (Personal)/Projects/Causal_network/causal_network/RData/neuron_network') # network data

# neuron_network not exist
neuron_network$Type <- c('Neuron', 'Oligo', 'Astro', 'Neuron', 'AO',
                        'Neuron', 'Neuron', 'Neuron', 'Neuron', "Neuron",
                        'AO', 'AO', 'Astro', 'Oligo', 'Olig', 'Astro',
                        'Astro', 'Astro', 'Olig', 'Astro', 'Oligo')

#
fData(neuron_sim_cds)$gene_short_name <- fData(neuron_sim_cds)$gene_short_names
fData(na_sim_cds)$gene_short_name <- fData(na_sim_cds)$gene_short_names

gene_name_vec <- c('Pax6', 'Mash1', 'Brn2', 'Zic1', 'Tuj1', 'Hes5', 'Scl', 'Olig2', 'Stat3', 'Myt1L', 'Aldh1L', 'Sox8', 'Mature')

debug(rdi_crdi_pseudotime)
rdi_crdi_pseudotime_res_list <- rdi_crdi_pseudotime(t(exprs(na_sim_cds)[1:12, 1:200]), window_size = 50) #13 mature gives Na values
*/

/*** R

# run multiple run on the real simulation datasets: 
 
*/

// [[Rcpp::export]]
NumericMatrix calculate_umi_cpp(NumericMatrix& expr_data, IntegerMatrix& super_graph, int k, int method, int k_density, double bw) //, method: 1 rdi; 2: lmi 
{
  const int n_genes(expr_data.cols()), n_samples(expr_data.rows());
  
  NumericMatrix uMI(n_genes, n_genes);     
  std::fill(uMI.begin(), uMI.end(), 0); //NA_REAL 
  NumericMatrix expr_1(n_samples, 1), expr_2(n_samples, 1);
  
  // this part maybe slow if we do paralleling because RDI matrix is huge 
  
  int i, j; 
  IntegerVector current_pair; 
  // #pragma omp parallel for shared(n_genes, expr_data, delays_len, RDI) private(i, j, k, expr_1, expr_2, _) //schedule(dynamic) default(none) //collapse(2) , _
  for(int super_graph_ind = 0; super_graph_ind < super_graph.rows(); super_graph_ind ++) 
  {
    current_pair = super_graph(Range(super_graph_ind, super_graph_ind), _);
    i = current_pair[0]; j = current_pair[1]; 
    
    expr_1 = expr_data(_, Range(i, i)); //Range(0, n_genes - 1)
    expr_2 = expr_data(_, Range(j, j)); //Range(0, n_genes - 1)
    
    uMI(i, j) = umi_cpp(expr_1, expr_2, k = k, method = method, k_density = k_density, bw = bw); // how to deal with delays include multiple values?
  }
  
  return uMI;
}
 
 //' @title
 //' calculate_umi_cpp_wrap
 //' @description
 //' This function estimates the uniformed mutual information for all gene pairs in the expr_data matrix
 //' 
 //' @param expr_data a matrix for all variables in a time-series where each row is a time point and each column is a variable (for example, a gene). 
 //' The rows are ordered according to time, from earliest to latest. 
 //' 
 //' @param super_graph An integer matrix where each row is the variable IDs or gene pairs you would like to estimate restricted direct information (RDI) from. 
 //'  
 //' @param k Number for nearest neighbors used in entropy calculation.
 //' 
 //' @param methodWhich 2D density estimator you would like to use. 1 is kde estimator and 2 is knn based estimator. Default to be 1. 
 //'
 //' @param k_density The number of k nearest neighbors you would like to use when calculating the density (only applicable when method == 2 or using knn based density estimation).
 //'
 //' @param bw Bindwidth used for the kernel density estimator. Currently it is not used. The bindwidth in the kde function is automatically estimated. 
 //' @details
 //' \code{calculate_umi_cpp_wrap} takes an expression matrix (expr_data) and possible gene pairs (encoded in super_graph) to estimate the uniformed 
 //' mutual information. 
 //' @return a numeric matrix of uniform mutual information values for all possible pairs of genes from expr_data. If the gene pairs is not encoded in the super_graph, it will remain as 0 in the matrix.  
 //' @export
 // [[Rcpp::export]]
 NumericMatrix calculate_umi_cpp_wrap(SEXP expr_data, SEXP super_graph, SEXP k, SEXP method, SEXP k_density, SEXP bw) 
 { 
   NumericMatrix expr_data_cpp(expr_data); 
   IntegerMatrix super_graph_cpp(super_graph); 

   int k_cpp = as<int>(k);
   int method_cpp = as<int>(method); 
   int k_density_cpp = as<int>(k_density); 

   double bw_cpp = as<double>(bw); 
   
   NumericMatrix uMI = calculate_umi_cpp(expr_data_cpp, super_graph_cpp, k_cpp, method_cpp, k_density_cpp, bw_cpp); //cores
   
   return uMI;
 }
 
