#include <cmath>
// #include <Rcpp.h>
#include <omp.h>
// // [[Rcpp::plugins(openmp)]]

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
#include "../inst/include/information_estimator.h"

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

double di_single_run_cpp(NumericMatrix& x, NumericMatrix& y, int n = 10)
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
  
  int max_thread = omp_get_max_threads();
  omp_set_num_threads(max_thread);
  #pragma omp parallel for shared(n, x_cols, y_cols, tot_len, x, y, x_past, y_past) private(i, j, col_id) //schedule(dynamic) default(none) //collapse(2) , _
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
  
  // NumericMatrix x_past = x(Range(tau - 1, tau - 2 + tot_len), _);
  // NumericMatrix y_past = y(Range(tau - 1, tau - 2 + tot_len), _);
  // NumericMatrix tmp; 
  // 
  // int i; 
  // omp_set_num_threads(8);
  // #pragma omp parallel for shared(n, tau, tot_len, x, y, x_past, y_past) private(i, _, tmp) //schedule(dynamic) default(none) //collapse(2) , _
  // for(i = 2; i < n + 1; i ++)
  // {
  //   tmp = x(Range(tau - i, tau - i - 1 + tot_len), _);
  //   x_past = cbind(tmp, x_past);
  //   tmp = y(Range(tau - i, tau - i - 1 + tot_len), _);
  //   y_past = cbind(tmp, y_past);
  // }
  // 
  NumericMatrix y_past0 = y(Range(tau, tau + tot_len - 1), _);
  
  List cmi_res = cmi_cpp(x_past, y_past0, y_past);
  return cmi_res["cmi_res"]; 
  // return cmi_cpp(x_past, y_past0, y_past);
}

//' @title
//' di_single_run
//' @description
//' This function simulates the DIRECTED mutual information from X to Y when you have a SINGLE run of the process
//' 
//' @param x one random variable from the time-series data
//' 
//' @param y another random variable from the time-series data
//'
//' @param n Parameter n determines the the number of previous time samples upon which the mi is conditioned (delay in the formula I(x-->y)=I(x_{t-d};y_t|y_{t-1}) default: d=1)
//'
//' @details
//' \code{di_single_run} takes two random variables x and y as well as a delay n to estimate the direct information between variable x and y. 
//' @return a matrix for the condition mutual information estimators between all pairwise variables (x, y) in the data matrix x, y
//' @export
// [[Rcpp::export]]
double di_single_run(SEXP x, SEXP y, SEXP n) 
{ 
  NumericMatrix x_cpp(x); 
  NumericMatrix y_cpp(y); 
  
  int n_cpp = as<int>(n); //why clone?
  
  double di_res = di_single_run_cpp(x_cpp, y_cpp, n_cpp);
  return di_res;
}

//-------------------------------------------------------------------------------------------------------------------
/*
calculate conditional direction information for single run
*/
//-------------------------------------------------------------------------------------------------------------------

double di_single_run_conditioned_cpp(NumericMatrix x, NumericMatrix y, NumericMatrix& z, int n)
{
 
  int x_cols = x.cols(), y_cols = y.cols(), z_rows = z.rows();  
  if(x_cols > y_cols)
  {
    stop("dimension of samples should be the same");
    return -1;
  }
  
  int tau = n; 
  int tot_len = x.rows() - tau; 
  
  // NumericMatrix x_past = x(Range(tau - 1, tau - 2 + tot_len), _); 
  // NumericMatrix yz_past = y(Range(tau - 1, tau - 2 + tot_len), _), tmp; 

  // need to redesign the cbind function in order to use parallel    
  int i, j, k_x, k_yz; 
  NumericMatrix x_past(tot_len, n * y_cols);
  NumericMatrix yz_past(tot_len, n * y_cols * (1 + z_rows));
  
  int max_thread = omp_get_max_threads();
  omp_set_num_threads(max_thread);
  #pragma omp parallel for shared(tot_len, n, y_cols, x, y, z, x_past, yz_past, z_rows) private(i, j, k_x, k_yz) //schedule(dynamic) default(none) //collapse(2) , _
  for(i = 0; i < tot_len; i ++)
  {
    for(j = 0; j < n; j ++)
    { // identify the index for the column
      for(k_x = 0; k_x < y_cols; k_x ++)
      {
        x_past(i, k_x + j * y_cols) = x(i + j, k_x); // identify the index for the row
    }
      
      for(k_yz = 0; k_yz <= y_cols * (1 + z_rows); k_yz ++)
      {
        if(k_yz < y_cols)
        {
          yz_past(i, k_yz + j * y_cols * (1 + z_rows)) = y(i + j, k_yz);      
        }
        else
        { // k_yz / y_cols - 2: row ind; k_yz % (1 + z_rows): column ind 
          yz_past(i, k_yz + j * y_cols * (1 + z_rows)) = z((int) k_yz / y_cols - 2, (int) k_yz % (1 + z_rows));       
        }
      }
    }
  }

  // // use cbind
  // for(int i = 1; i < n + 1; i ++)
  // {
  //   if(i > 1)
  //   {
  //     tmp = x(Range(tau - i, tau - i - 1 + tot_len), _);
  //     x_past = cbind(tmp, x_past);
  //     tmp = y(Range(tau - i, tau - i - 1 + tot_len), _);
  //     yz_past = cbind(tmp, yz_past);
  //   }
    
  //   for(int j = 0; j < z.cols(); j ++ )
  //   {
  //     tmp = z(Range(tau - i, tau - i - 1 + tot_len), _);
  //     yz_past = cbind(tmp, yz_past);
  //   } 
  // }
  
  NumericMatrix y_past0 = y(Range(tau, tau + tot_len - 1), _);
  List cmi_res = cmi_cpp(x_past, y_past0, yz_past);
  return cmi_res["cmi_res"]; 
  // return cmi_cpp(x_past, y_past0, yz_past); 
}

//' @title
//' di_single_run_conditioned
//' @description
//' This function simulates the CONDITIONED DIRECTED mutual information from X to Y when you have a SINGLE run of the processes
//' 
//' @param x one random variable from the time-series data
//' 
//' @param y another random variable from the time-series data
//'
//' @param z z is a dataframe (or matrix) containing the data of other processes upon the past of which the mi is conditioned
//' 
//' @param n Parameter n determines the the number of previous time samples upon which the mi is conditioned
//' 
//' @details
//' \code{di_single_run_conditioned} takes two random variables x and y as well as the parameter n to calculate the direct information conditioned on variable z. 
//' @return a matrix for the condition mutual information estimators between all pairwise variables (x, y) in the data matrix x, y
//' @export
// [[Rcpp::export]]
double di_single_run_conditioned(SEXP x, SEXP y, SEXP z, SEXP n) 
{ 
  NumericMatrix x_cpp(x); 
  NumericMatrix y_cpp(y); 
  NumericMatrix z_cpp(z); 
  
  int n_cpp = as<int>(n); 
  
  double di_res = di_single_run_conditioned_cpp(x_cpp, y_cpp, z_cpp, n_cpp);
  return di_res;
}

//-------------------------------------------------------------------------------------------------------------------
/*
calculate restricted direction information for many runs
*/
//-------------------------------------------------------------------------------------------------------------------

double rdi_many_runs_cpp(NumericMatrix& x, NumericMatrix& y)
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
    List cmi_res = cmi_cpp(x_0, y_0, y_1);
    ans[t - 1] = cmi_res["cmi_res"]; 
    // ans[t - 1] = cmi_cpp(x_0, y_0, y_1);
  }
  
  return sum(ans); // return sum of ans? 
}

//' @title
//' rdi_many_runs
//' @description
//' This function simulates the DIRECTED mutual information from X to Y when you have multiple run of the processes
//' 
//' @param x a random variable with multiple run of the same process
//' 
//' @param y another random variable with multiple run of another process
//'
//' @details
//' \code{rdi_many_runs} takes two random variables (with the same multiple realization of the two processes) and estimate the direct information between them. 
//' using the formula: 0.5 * d * log(pi) - log(gamma(0.5 * d + 1))
//' It's implimented in C++, providing a (small) increase in speed over the R equivalent.
//' @return a numeric vector storing the DI from two multiple run variables
//' @export
// [[Rcpp::export]]
double rdi_many_runs(SEXP x, SEXP y) 
{ 
  NumericMatrix x_cpp(x); 
  NumericMatrix y_cpp(y); 
  
  double rdi_many_runs_res = rdi_many_runs_cpp(x_cpp, y_cpp);
  return rdi_many_runs_res;
}

//-------------------------------------------------------------------------------------------------------------------
/*
calculate restricted direction information for a single run
*/
//-------------------------------------------------------------------------------------------------------------------

double rdi_single_run_cpp(NumericMatrix& x, NumericMatrix& y, int d) 
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
  
  // Rcout << "\n before running CMI " << std::endl;
  List cmi_res = cmi_cpp(x_0, y_0, y_1);
  return cmi_res["cmi_res"];
  // return cmi_cpp(x_0, y_0, y_1); 
}

//' @title
//' rdi_single_run
//' @description
//' This function simulates the RESTRICTED DIRECTED mutual information from X to Y when you have a SINGLE run of the processes
//' 
//' @param x one random variable from the time-series data
//' 
//' @param y another random variable from the time-series data
//' 
//' @param d delay in the formula I(x-->y)=I(x_{t-d};y_t|y_{t-1}) default: d=1
//'
//' @details
//' \code{rdi_single_run} takes two random variables x and y as well as a time delay d to estimate the restricted direct infomation between them.
//' @return a matrix for the condition mutual information estimators between all pairwise variables (x, y) in the data matrix x, y
//' @export
// [[Rcpp::export]]
double rdi_single_run(SEXP x, SEXP y, SEXP d) 
{ 
  NumericMatrix x_cpp(x); 
  NumericMatrix y_cpp(y); 
  
  int d_cpp = as<int>(d); 
  
  double rdi_single_run_res = rdi_single_run_cpp(x_cpp, y_cpp, d_cpp);
  return rdi_single_run_res;
}


//-------------------------------------------------------------------------------------------------------------------
/*
  calculate the lagged mutual information 
*/
//-------------------------------------------------------------------------------------------------------------------
// implement multiple run version here 
double lmi_single_run_cpp(NumericMatrix& x, NumericMatrix& y, int delay) 
{
  int Nx = x.rows(); int Ny = y.rows(); //int dx = x.cols(); int dy = y.cols(); 
  
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
//' @details
//' \code{cmi} takes two random variable x and y and estimated their mutual information conditioned on the third random variable z
//' using the KSG estimator. 
//' It relies on the ANN package to query the kNN with KDTree algorithm. 
//' @return a numeric value for the condition mutual information estimator between two variables (x, y) conditioned on variable z
//' @export
// [[Rcpp::export]]
double lmi_single_run(SEXP x, SEXP y, SEXP delay) 
{ 
  NumericMatrix x_cpp(x); 
  NumericMatrix y_cpp(y); 

  int delay_cpp = as<int>(delay); 

  double lmi_res = lmi_single_run_cpp(x_cpp, y_cpp, delay_cpp);
  return lmi_res;
}

//-------------------------------------------------------------------------------------------------------------------
/*
calculate conditional restricted direction information for a single run
*/
//-------------------------------------------------------------------------------------------------------------------

double rdi_single_run_conditioned_cpp(NumericMatrix& x, NumericMatrix& y, NumericMatrix& z, NumericVector& z_delays, int d) 
{
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
  
  List cmi_res = cmi_cpp(x_0, y_0, yz);
  return cmi_res["cmi_res"];  
  // return cmi_cpp(x_0, y_0, yz); 
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
//' @details
//' \code{rdi_single_run_conditioned} takes two random variables x and y as well as the parameter n to calculate the restricted direct information conditioned on variable z. 
//' @return a matrix for the condition mutual information estimators between all pairwise variables (x, y) in the data matrix x, y
//' @export
// [[Rcpp::export]]
double rdi_single_run_conditioned(SEXP x, SEXP y, SEXP z, SEXP z_delays, SEXP d) 
{ 
  NumericMatrix x_cpp(x); 
  NumericMatrix y_cpp(y); 
  NumericMatrix z_cpp(z); 
  
  NumericVector z_delays_cpp(z_delays); 
  
  int d_cpp = as<int>(d); 
  
  double rdi_single_run_conditioned_res = rdi_single_run_conditioned_cpp(x_cpp, y_cpp, z_cpp, z_delays_cpp, d_cpp);
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

// add the super-graph here too 
// [[Rcpp::export]]
List calculate_rdi_cpp(NumericMatrix& expr_data, IntegerVector delays, IntegerMatrix& super_graph, IntegerVector& turning_points, int method) //, method: 1 rdi; 2: lmi 
{
  const int n_genes(expr_data.cols()), n_samples(expr_data.rows());
  // Rcout << "n_genes is " << n_genes << "n_samples is " << n_samples << std::endl;
  if(max(delays) > n_samples)
  {
    stop("Number of samples has to be larger than the delays"); 
    return -1; 
  }
  
  IntegerMatrix max_rdi_delays(n_genes, n_genes); // declare max_rdi_delays when turning_points provided 
  
  // Convert arguments to Armadillo objects
  // const arma::mat data = R2armaMat_num(expr_data); // Continuous data does not need reindexing R->C++
  // Rcout << "n_genes is " << n_genes << " n_samples is " << n_samples << " delay length is " << delays.length() << std::endl;
  // const int n_pairs = (int)(((double)n_genes/2.0) * (double)(n_genes+1));
  
  // // Number of cores to use
  // //const int cores = sysconf(_SC_NPROCESSORS_ONLN);
  // omp_set_num_threads(cores);
  
  // Compute RDI in parallel
  int delays_len = delays.length(); 
  if(turning_points.length() == n_genes) 
  {
    NumericMatrix RDI(n_genes, n_genes); 
  } else {
    NumericMatrix RDI(n_genes, n_genes * delays_len);     
  }

  NumericMatrix expr_1(n_samples, 1), expr_2(n_samples, 1);
  std::fill(RDI.begin(), RDI.end(), 0); //NA_REAL 
  
  // this part maybe slow if we do paralleling because RDI matrix is huge 
  
  int i, j, k; 
  IntegerVector current_pair; 
  // #pragma omp parallel for shared(n_genes, expr_data, delays_len, RDI) private(i, j, k, expr_1, expr_2, _) //schedule(dynamic) default(none) //collapse(2) , _
  for(int super_graph_ind = 0; super_graph_ind < super_graph.rows(); super_graph_ind ++) 
  {
    current_pair = super_graph(Range(super_graph_ind, super_graph_ind), _);
    // Rcout << "current_pair is " << current_pair << std::endl;
    i = current_pair[0]; j = current_pair[1]; 
    // Rcout << "i is " << i << " j is " << j << std::endl;
    
    expr_1 = expr_data(_, Range(i, i)); //Range(0, n_genes - 1)
    expr_2 = expr_data(_, Range(j, j)); //Range(0, n_genes - 1)

    // int max_thread = omp_get_max_threads();
    // omp_set_num_threads(max_thread);
    // #pragma omp parallel for shared(delays_len, i, j, n_genes, expr_1, expr_2, RDI) private(k) //schedule(dynamic) default(none) //collapse(2) , _
    
    // if we provide turning_points estimation, we will use that to determine the time delay
    if(turning_points.length() == n_genes) // !Rf_isNull(turning_points)
    {
      Rcout << "using user provided information about time-delay " << turning_points.length() << std::endl;
      double current_delay = turning_points[i] - turning_points[j];
      if(method == 1)
      {
        RDI(i, j) = rdi_single_run_cpp(expr_1, expr_2, current_delay); // how to deal with delays include multiple values?
      }
      else if(method == 2)
      { // + k * n_genes
        RDI(i, j) = lmi_single_run_cpp(expr_1, expr_2, current_delay); // how to deal with delays include multiple values?
      }
      max_rdi_delays(i, j) = current_delay;
    } else {
      for(k = 0; k < delays_len; ++ k)
      {
        if(i == j) 
        {
          RDI(i, j + k * n_genes) = 0; 
        }
        
        // Rcout << "number of expr_1's row is " << expr_1.rows() << "number of expr_2's row  is " << expr_2.rows() << std::endl;
        
        // RDI(i, j + k * n_genes) = 0; // how to deal with delays include multiple values?
        // Rcout << "the current RDI value is " << rdi_single_run_cpp(expr_1, expr_2, delays[k]) << std::endl;
        if(method == 1)
        {
          RDI(i, j + k * n_genes) = rdi_single_run_cpp(expr_1, expr_2, delays[k]); // how to deal with delays include multiple values?
        }
        else if(method == 2)
        {
          RDI(i, j + k * n_genes) = lmi_single_run_cpp(expr_1, expr_2, delays[k]); // how to deal with delays include multiple values?
        }
      // Rcout << "The RDI value  is " << RDI(i, j + k * n_genes) << std::endl;
      }
    }
  }

  // perform the multiple run test (add a new argument to represent the vector for different runs)
  if(turning_points.length() == n_genes) // !Rf_isNull(turning_points)
  {
    // // run extract_max_rdi_value_delay 
    // List max_rdi_value_delay_list = extract_max_rdi_value_delay(RDI, delays); 
  
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
  // return max_rdi_value_delay_list;  
}

// [[Rcpp::export]]
List calculate_rdi_cpp_wrap(SEXP expr_data, SEXP delays, SEXP super_graph, SEXP turning_points, SEXP method) //, SEXP R_cores
{ 
  // Rcout << "before running calculate_rdi_cpp_wrap" << std::endl;
  NumericMatrix expr_data_cpp(expr_data); 
  IntegerVector delays_cpp(delays); 
  
  IntegerMatrix super_graph_cpp(super_graph); 
  IntegerVector turning_points_cpp(turning_points);
  int method_cpp = as<int>(method); 
  // Rcout << "before running calculate_rdi_cpp" << std::endl;
  List rdi_list = calculate_rdi_cpp(expr_data_cpp, delays_cpp, super_graph_cpp, turning_points_cpp, method_cpp); //cores
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
List extract_top_incoming_nodes_delays(NumericMatrix max_rdi_value, IntegerMatrix max_rdi_delays, int k)
{
  int n_genes = max_rdi_value.rows(), i, j, current_k_ind;
  NumericMatrix top_incoming_values(n_genes, k + 1);
  IntegerMatrix top_incoming_delays(n_genes, k + 1), top_incoming_nodes(n_genes, k + 1);
  IntegerVector top_k_plus_1 = seq_len(k); top_k_plus_1.push_front(0); 
  NumericVector x, x_ordered;
  IntegerVector x_order_ind, top_k_plus_1_x_order_ind; 

  // Rcout << top_k_plus_1 << std::endl; 

  for(i = 0; i < n_genes; i ++)
  {
    max_rdi_value(i, i) = -5; // avoid diag is calculated as the incoming node 
    x = max_rdi_value(_, i); // first dimension is source; second dimension is target 
    NumericVector x_ordered = clone(x).sort(true); // decreasing sort; rdi value from highest to lowest 
    x_order_ind =  match(x_ordered, x) - 1; // order starts from 1; all duplicated points get the same location 

    // Rcout << "here x is " << x << "x_ordered is " << x_ordered << " x_order_ind is " << x_order_ind << std::endl; 

    top_k_plus_1_x_order_ind = x_order_ind[top_k_plus_1];  
    max_rdi_value(i, i) = 0; // go back to the original value 
    // Rcout << "x is" << x << " top_k_plus_1_x_order_ind is " << top_k_plus_1_x_order_ind << std::endl; 
    for(j = 0; j < k + 1; j ++) // for each gene get the top k + 1 input node 
    {
      current_k_ind = top_k_plus_1_x_order_ind[j]; // index x_ordered to get the top k's id (node index)
      // Rcout << "top_k_plus_1_x_order_ind[j] " << top_k_plus_1_x_order_ind[j] << "x_ordered is" << x_ordered << " current_k_ind is " << current_k_ind << std::endl; 
      top_incoming_nodes(i, j) = current_k_ind; //top_k_plus_1_x_order_ind[j]; //[top_k_plus_1];
      top_incoming_delays(i, j) = max_rdi_delays(current_k_ind, i); //
      top_incoming_values(i, j) = x_ordered[j]; //;  
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
                                            NumericMatrix& max_rdi_value, IntegerMatrix& max_rdi_delays, int k) //, const int cores, const bool verbose
{
  // if(verbose == TRUE) Rcout << "Calculating the conditional restricted direct mutual information for each pair of genes" << std::endl;
  NumericVector k_ncol = NumericVector::create(k, expr_data.cols() - 2); 
  k = min(k_ncol); // minal value from k and the column 
  List top_incoming_nodes_delays_list = extract_top_incoming_nodes_delays(max_rdi_value, max_rdi_delays, k); 

  IntegerMatrix top_incoming_nodes = top_incoming_nodes_delays_list["top_incoming_nodes"];
  IntegerMatrix top_incoming_delays = top_incoming_nodes_delays_list["top_incoming_delays"]; 
  NumericMatrix top_incoming_values = top_incoming_nodes_delays_list["top_incoming_values"];  

  int n_genes = expr_data.cols(), n_sample = expr_data.rows(), tau, total_length, id3, delay_id3, current_id_tmp, delay; //tmp: current target  

  IntegerVector gene_pair(2), tmp(1), valid_top_k_in_k(k), valid_top_k_in_k_ordered, top_k_plus_1(k + 1), valid_top_k(k), ind_for_top_k(k), valid_top_k_incoming_delays(k), valid_delay_tmp(1); // optimize by setting the correct dimension 
  NumericMatrix cRDI_mat(n_genes, n_genes), expr1_t, expr2_t, past_tmp; 
  std::fill(cRDI_mat.begin(), cRDI_mat.end(), NA_REAL); // initialize the cRDI matrix as C ++ NA values 
  
  // OpenMP here too (not thread safe) 
  for(int i = 0; i < super_graph.rows(); i ++) 
  {
    // index_name = super_graph(i, _); // Range(i, i)
    gene_pair = super_graph(i, _); 
    top_k_plus_1 = top_incoming_nodes(Range(gene_pair[1], gene_pair[1]), _); // top k incoming nodes 
    tmp = gene_pair[1];

    if(intersect(top_k_plus_1, tmp).length() > 0) // if current node is in the top_k_plus_1, remove k; else remove the smallest one 
    {
      valid_top_k = setdiff(top_k_plus_1, tmp);
    }
    else
    {
      tmp = top_k_plus_1[k + 1]; // remove the smallest one 
      valid_top_k = setdiff(top_k_plus_1, tmp); 
    }
    
    // assign the node indices to ind_for_top_k
    // Rcout << "valid_top_k " << valid_top_k << "top_k_plus_1 " << top_k_plus_1 << std::endl;
    valid_top_k_in_k = match(valid_top_k, top_k_plus_1) - 1; // match gives order and the index (1 to k) starts from 1 -> valid_top_k_in_k
    // Rcout << "valid_top_k " << valid_top_k << "valid_top_k_in_k " << valid_top_k_in_k << std::endl;
    valid_top_k = top_k_plus_1[valid_top_k_in_k]; // set valid_top_k to the same order as top_k_plus_1 by valid_top_k_in_k for downstream analysis
    valid_top_k_in_k_ordered = clone(valid_top_k_in_k).sort(); 
    // Rcout << "after sort " << valid_top_k << "valid_top_k_in_k_ordered " << valid_top_k_in_k_ordered << std::endl; 
    valid_top_k = top_k_plus_1[valid_top_k_in_k_ordered]; 
    // Rcout << "New sort " << valid_top_k << std::endl;
    
    if(valid_top_k.size() < k) 
    { // avoid any weird rare situation 
      // Rcout << "all incoming node has the same RDI value for gene pair " << gene_pair[0] << " and " << gene_pair[1] << ". Either of the gene may have no expression." << std::endl;
      cRDI_mat(gene_pair[0], gene_pair[1]) = max_rdi_value(gene_pair[0], gene_pair[1]); 
      continue;
    }

    for(int j = 0; j < k; j ++ )
    {
      ind_for_top_k[j] = top_k_plus_1[valid_top_k_in_k[j]]; // for example 10, 150, 3
    }
    
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

    // Rcout << "current i is " << i << "max for top_incoming_delays is " << max(top_incoming_delays) << std::endl;
    // Rcout << "top_k is " << valid_top_k << "delay is " << delay << " valid_top_k_incoming_delays are " << valid_top_k_incoming_delays << std::endl;
    // Rcout << "n_sample is " << n_sample << " tau - delay " << tau << ", " << - delay << " tau - 1 - delay + total_length " << total_length << std::endl;
    expr1_t = expr_data(Range(tau - delay, tau - 1 - delay+total_length), Range(gene_pair[0], gene_pair[0])); 
    // Rcout << "tau" << tau << "tau - 1 + total_length" << tau - 1 + total_length << std::endl;
    expr2_t = expr_data(Range(tau, tau - 1 + total_length), Range(gene_pair[1], gene_pair[1])); 
    // Rcout << "tau - 1" << tau - 1 << "tau - 2 + total_length" << tau - 2 + total_length << std::endl;
    NumericMatrix yz = expr_data(Range(tau - 1, tau  - 2 + total_length), Range(gene_pair[1], gene_pair[1])); // initialize yz as expr2_t with one time lag 

    // OpenMP here too 
    for(int id = 0; id < k; id ++)
    {
      id3 = valid_top_k[id];
      delay_id3 = valid_top_k_incoming_delays[id];
      past_tmp = expr_data(Range(tau - delay_id3, tau  - delay_id3 - 1 + total_length),
                            Range(id3, id3));
      yz = cbind(yz, past_tmp);
    }

    // NumericMatrix past0 = transpose(yz);
    List cmi_res = cmi_cpp(expr1_t, expr2_t, yz);
    cRDI_mat(gene_pair[0], gene_pair[1]) =  cmi_res["cmi_res"];    
    // cRDI_mat(gene_pair[0], gene_pair[1]) = cmi_cpp(expr1_t, expr2_t, yz);

  }
  
  // maybe we should add multiple run analysis here 

  // return results
  return cRDI_mat; 
}

// [[Rcpp::export]]
NumericMatrix calculate_conditioned_rdi_cpp_wrap(SEXP expr_data, SEXP super_graph,
                               SEXP max_rdi_value, SEXP max_rdi_delays, SEXP k) //, SEXP R_cores
{
  NumericMatrix expr_data_cpp(expr_data);
  IntegerMatrix super_graph_cpp(super_graph);
  NumericMatrix max_rdi_value_cpp(max_rdi_value);
  IntegerMatrix max_rdi_delays_cpp(max_rdi_delays);
  int k_cpp = as<int>(k);

  NumericMatrix cRDI = calculate_conditioned_rdi_cpp(expr_data_cpp, super_graph_cpp, max_rdi_value_cpp, max_rdi_delays_cpp, k_cpp); //cores
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
//' entropy
//' @description
//' This subroutine calculates the volume of a d-dimensional unit ball for Euclidean norm
//' 
//' @param R_x number of dimension
//' 
//' @param R_k number of dimension
//'
//' @details
//' \code{entropy} takes a integer of dimensions and then calculate the olume of a d-dimensional unit ball for Euclidean norm
//' using the formula: 0.5 * d * log(pi) - log(gamma(0.5 * d + 1))
//' It's implimented in C++, providing a (small) increase in speed over the R equivalent.
//' @return a updated matrix with gene expression smoothed with window size equal to window_size
//' @export
// [[Rcpp::export]]
NumericMatrix smooth_gene(NumericMatrix& expr_data, const int window_size)
{ // columns: genes; rows: cells 
  int win_range = expr_data.rows() - window_size;
  NumericMatrix expr_data_smooth(expr_data(Range(0, win_range), _));
  
  // Rcout << "expr_data_smooth's dimensionality are " << expr_data_smooth.cols() << " " << expr_data_smooth.rows() << std::endl;
  NumericVector tmp(win_range + 1), tmp_1;
  for(int i = 0; i < expr_data.cols(); i ++) // genes 
  {
    for(int j = 0; j <= win_range; j ++) // cells
    {
      tmp_1 = expr_data(Range(j, j + window_size - 1), Range(i, i)); 
      tmp[j] = mean(tmp_1);
    }
    // Rcout << "tmp is " << tmp.length() << std::endl;
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


// add the super-graph here too 


// [[Rcpp::export]]
double rdi_single_run_multiple_run_cpp(NumericMatrix& x, NumericMatrix& y, int d, IntegerVector run_vec) 
{
  int Nx = x.rows(); int Ny = y.rows(); int Cx = x.cols(); int Cy = y.cols(); int run_len; 
  
  if(Nx != Ny)
  {
    stop("The number of time samples has to be the same for X and Y");
    return -1;
  }

  int run_max = max(run_vec); // run_vec: vector to represent the run, need to be just numeric numbers and starts from 0. 
  NumericMatrix x_0_res, y_0_res, y_1_res, tmp; // matrix passed to cmi function
  IntegerVector index_all = seq_len(run_vec.size()) - 1, current_run_ind; // cell index from 0 to maximal number of cells; cell index corresponding to current run during iteration 
  int current_ind = 0; // current initial index for storing new slice of data 

  if(run_max == 0) 
  {
    x_0_res = x(Range(0, Nx - d - 1), _); // x should only has 1 colum 
    y_0_res = y(Range(d, Ny - 1), _);
    y_1_res = y(Range(d - 1, Nx - 2), _); 
  } else {

    // correctly set the dimensionality for the x_0, y_0 and y_1 matrix 
    int tot_len = Nx - d * (run_max + 1); 
    // Rcout << "\n tot_len is" << tot_len << std::endl; 
   // for(int i = 0; i <= run_max; i++) // concatenate the data 
    // {
    //   current_run_ind = index_all[run_vec == i]; // get the cells that belong to a particular run 
    //   run_len = current_run_ind.size();

    //   if(current_run_ind.size() < d)
    //   {
    //     stop("Each run has to have more samples than the designated 'delays'");
    //     return -1;
    //   }

    //   tot_len += run_len - d; // easier way: tot_len =  Nx - d * (run_max + 1)
    // }

    mat x_0(tot_len, 1), y_0(tot_len, 1), y_1(tot_len, 1); // arma matrix 
    
    uvec pos; // arma uvec 
    vec vals; // arma vec 

    for(int i = 0; i <= run_max; i++) // concatenate the data 
    {
      current_run_ind = index_all[run_vec == i]; // get the cells that belong to a particular run 
      run_len = current_run_ind.size(); 
      // Rcout << "current_run_ind is " << current_run_ind << " run_vec is " << run_vec << "run_len is " << run_len << std::endl;
      if(current_run_ind.size() < d)
      {
        stop("Each run has to have more samples than the designated 'delays'");
        return -1;
      }

      pos = as<uvec>(wrap(Range(current_ind, current_ind + run_len - d - 1))); 
      IntegerVector test_range = Range(0, 3);
      // Rcout << "\n pos is" << pos << "test_range is " << test_range  << "current_ind is " << current_ind << std::endl; 
      
      // assuming x is only one column (Cx, Cy is not used)
      tmp = x(Range(min(current_run_ind), min(current_run_ind) + run_len - d - 1), Range(0, 0)); //x(Range(1, 2), Range(1, 1));
      vals = as<vec>(wrap(tmp)); 
      x_0.elem(pos) = vals; 

      IntegerVector rng_vec = Range(min(current_run_ind), min(current_run_ind) + run_len - d - 1); 
      for(int test_na = 0; test_na < vals.n_elem; test_na ++) 
      {
        if(arma::is_finite(vals[test_na]) == false) 
        {
          // Rcout << "identify non-finite values for x_0 here; i is " << i << " rng_vec is " << rng_vec << std::endl; 
        }
      }

      // Rcout << "x_0.elem(pos) is" << vals.t() << std::endl; 

      tmp = y(Range(min(current_run_ind) + d, min(current_run_ind) + run_len - 1), Range(0, 0)); //x(Range(1, 2), Range(1, 1));
      vals = as<vec>(wrap(tmp)); 
      y_0.elem(pos) = vals; 

      for(int test_na = 0; test_na < vals.n_elem; test_na ++) 
      {
        if(arma::is_finite(vals[test_na]) == false) 
        {
          Rcout << "identify non-finite values for x_0 here; i is " << i << " rng_vec is " << rng_vec << std::endl; 
        }
      }

      // Rcout << "y_0.elem(pos) is" << vals.t() << std::endl; 

      tmp = y(Range(min(current_run_ind) + d - 1, min(current_run_ind) + run_len - 2), Range(0, 0)); //x(Range(1, 2), Range(1, 1));
      vals = as<vec>(wrap(tmp)); 
      y_1.elem(pos) = vals; 

      for(int test_na = 0; test_na < vals.n_elem; test_na ++) 
      {
        if(arma::is_finite(vals[test_na]) == false) 
        {
          Rcout << "identify non-finite values for x_0 here; i is " << i << " rng_vec is " << rng_vec << std::endl; 
        }
      }

      // Rcout << "y_1.elem(pos) is" << vals.t() << std::endl; 

      current_ind = current_ind + run_len - d; // move to the the next position after filling the current run (note that there is no - 1)
    }

    // Rcout << "\n before converting to NumericMatrix " << "x_0 is " << x_0.t() << "y_0" << y_0.t() << "y_1" << y_1.t() << std::endl; 

    x_0_res = as<NumericMatrix>(wrap(x_0)); // this conversion have problems? 
    y_0_res = as<NumericMatrix>(wrap(y_0));
    y_1_res = as<NumericMatrix>(wrap(y_1));

  }  
  // Rcout << "\n before running CMI " << "x_0_res is " << transpose(x_0_res) << "y_0_res" << transpose(y_0_res) << "y_1_res" << transpose(y_1_res) << std::endl; 
  for(int test_na = 0; test_na < y_1_res.size(); test_na ++) 
  {
    if(arma::is_finite(y_1_res[test_na]) == false) 
    {
      // Rcout << "identify non-finite values for y_1_res here" << std::endl; 
    }
  }

  List cmi_res = cmi_cpp(x_0_res, y_0_res, y_1_res); 
  return cmi_res["cmi_res"]; 
  // return cmi_cpp(x_0_res, y_0_res, y_1_res); 
}

// [[Rcpp::export]]
List calculate_rdi_multiple_run_cpp(NumericMatrix& expr_data, IntegerVector delays, IntegerVector run_vec, IntegerMatrix& super_graph, IntegerVector turning_points, int method) //, method: 1 rdi; 2: lmi 
{
  const int n_genes(expr_data.cols()), n_samples(expr_data.rows());
  // Rcout << "n_genes is " << n_genes << "n_samples is " << n_samples << std::endl;
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
  NumericMatrix RDI(n_genes, n_genes * delays_len), expr_1(n_samples, 1), expr_2(n_samples, 1);
  std::fill(RDI.begin(), RDI.end(), 0); 
  
  // this part maybe slow if we do paralleling because RDI matrix is huge 
  
  int i, j, k; 
  IntegerVector current_pair; 
  // #pragma omp parallel for shared(n_genes, expr_data, delays_len, RDI) private(i, j, k, expr_1, expr_2, _) //schedule(dynamic) default(none) //collapse(2) , _
  for(int super_graph_ind = 0; super_graph_ind < super_graph.rows(); super_graph_ind ++) 
  {
    current_pair = super_graph(Range(super_graph_ind, super_graph_ind), _);
    // Rcout << "current_pair is " << current_pair << std::endl;
    i = current_pair[0]; j = current_pair[1]; 
    // Rcout << "i is " << i << " j is " << j << std::endl;
    
    expr_1 = expr_data(_, Range(i, i)); //Range(0, n_genes - 1)
    expr_2 = expr_data(_, Range(j, j)); //Range(0, n_genes - 1)

    // int max_thread = omp_get_max_threads();
    // omp_set_num_threads(max_thread);
    // #pragma omp parallel for shared(delays_len, i, j, n_genes, expr_1, expr_2, RDI) private(k) //schedule(dynamic) default(none) //collapse(2) , _
    if(turning_points.length() == n_genes) //!Rf_isNull(turning_points) 
    {
      double current_delay = turning_points[i] - turning_points[j];
      RDI(i, j) = rdi_single_run_multiple_run_cpp(expr_1, expr_2, current_delay, run_vec); // how to deal with delays include multiple values?
      max_rdi_delays(i, j) = current_delay;
    } else {
      for(k = 0; k < delays_len; ++ k)
      {
        if(i == j) 
        {
          continue; // keep diagnoal as NA 
          // RDI(i, j + k * n_genes) = 0; 
        }
        
        // Rcout << "number of expr_1's row is " << expr_1.rows() << "number of expr_2's row  is " << expr_2.rows() << std::endl;
        
        // RDI(i, j + k * n_genes) = 0; // how to deal with delays include multiple values?
        // Rcout << "the current RDI value is " << rdi_single_run_cpp(expr_1, expr_2, delays[k]) << std::endl;
        if(method == 1)
        {
          RDI(i, j + k * n_genes) = rdi_single_run_multiple_run_cpp(expr_1, expr_2, delays[k], run_vec); // how to deal with delays include multiple values?
        }
        else if(method == 2)
        {
          // RDI(i, j + k * n_genes) = lmi_single_run_cpp(expr_1, expr_2, delays[k], run_vec); // how to deal with delays include multiple values?
        }
        // Rcout << "The RDI value  is " << RDI(i, j + k * n_genes) << std::endl;
      }
    }
  }

  // perform the multiple run test (add a new argument to represent the vector for different runs)
  if(turning_points.length() == n_genes) // Rf_isNull(turning_points) 
  {
    // // run extract_max_rdi_value_delay 
    // List max_rdi_value_delay_list = extract_max_rdi_value_delay(RDI, delays); 
  
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
  // return max_rdi_value_delay_list;  
}

// [[Rcpp::export]]
List calculate_rdi_multiple_run_cpp_wrap(SEXP expr_data, SEXP delays, SEXP run_vec, SEXP super_graph, SEXP turning_points, SEXP method) //, SEXP R_cores
{ 
  // Rcout << "before running calculate_rdi_cpp_wrap" << std::endl;
  NumericMatrix expr_data_cpp(expr_data); 
  IntegerVector delays_cpp(delays); 
  IntegerVector run_vec_cpp(run_vec); 
  
  IntegerMatrix super_graph_cpp(super_graph); 
  IntegerVector turning_points_cpp(turning_points);
  int method_cpp = as<int>(method); 
  // Rcout << "before running calculate_rdi_cpp" << std::endl;
  List rdi_list = calculate_rdi_multiple_run_cpp(expr_data_cpp, delays_cpp, run_vec_cpp, super_graph_cpp, turning_points_cpp, method_cpp); //cores
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
double rdi_multiple_runs_conditioned_cpp(NumericMatrix& x, NumericMatrix& y, NumericMatrix& z, IntegerVector& z_delays, int d, IntegerVector run_vec) 
{
  int Nx = x.rows(); int Ny = y.rows(); int Nz = z.rows(); int Cx = x.cols(); int Cy = y.cols(); int Cz = z.cols(); int run_len; 
  if(Nx != Ny | Nx != Nz)
  {
    stop("The number of time samples has to be the same for X and Y or Z");
    return -1;
  }
  
  z_delays.push_back(d);
  int tau = max(z_delays); 

  int run_max = max(run_vec); // run_vec: vector to represent the run, need to be just numeric numbers and starts from 0. 
  NumericMatrix x_0_res, y_0_res, yz_res, tmp; // matrix passed to cmi function
  IntegerVector index_all = seq_len(run_vec.size()) - 1, current_run_ind; // cell index from 1 to maximal number of cells; cell index corresponding to current run during iteration 
  int current_ind = 0; // current initial index for storing new slice of data (--current_incoming_k)

  // Rcout << "index_all is: " << index_all << std::endl;
  // Rcout << "run_max is: " << run_max << std::endl;

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
    
    uvec pos; //, test;
    vec vals; 

    int tot_len = Nx - tau * (run_max + 1); 
    // Rcout << "Pass here 2; run_max is "  << run_max << std::endl;

    // // correctly set the dimensionality for the yz matrix 
    // for(int i = 0; i <= run_max; i++) // concatenate the data 
    // {
    //   current_run_ind = index_all[run_vec == i]; // get the cells that belong to a particular run 
    //   run_len = current_run_ind.size();

    //   if(current_run_ind.size() < d)
    //   {
    //     stop("Each run has to have more samples than the designated 'delays'");
    //     return -1;
    //   }

    //   tot_len += run_len - tau; // easier way: Nx - tau * (run_max + 1)
    // }
    
    // Rcout << "yz dimensionality is "  << tot_len << ", " << z.cols() + 1 << " Nx is " << Nx << " tau is " << tau << std::endl;
    
    mat yz(tot_len, z.cols() + 1); // tot_len: number of cells in total across all runs; z.col() + 1: number of top k incoming nodes

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
      

      NumericMatrix locs_tmp(2, tot_len); // first row: row index; second row: column index; each column is a coordinate
      locs_tmp(0, _) = Range(current_ind, current_ind + tot_len - 1);
      locs_tmp(1, _) = rep(0, tot_len); 
      umat locs = as<umat>(wrap(locs_tmp));  
      // Rcout << "locs is (first) "  << locs << " size(yz) is " << size(yz) << std::endl;
      uvec eids = sub2ind( size(yz), locs ); // Obtain Element IDs
      
      // Rcout << "pass eids"  << eids << std::endl;
      
      // pos = as<uvec>(wrap(Range(current_ind, current_ind + tot_len - 1))); 
      tmp = y(Range(min(current_run_ind) + tau - 1, min(current_run_ind) + tau - 2 + tot_len), Range(0, 0)); // x(Range(min(current_run_ind), min(current_run_ind) + Nx - d - 1), Range(1, 1)); //x(Range(1, 2), Range(1, 1));

      vals = as<vec>(wrap(tmp)); 
      yz.elem( eids ) = vals;
           
      // Rcout << "pass yz join_rows "  << yz << std::endl;
      // OpenMP here too 
      // int max_thread = omp_get_max_threads();
      // omp_set_num_threads(max_thread);
      // #pragma omp parallel for shared(z, tau, z_delays, tot_len) private(j, tmp, yz) //schedule(dynamic) default(none) //collapse(2) , _
        
      // using cbind
      for(int j = 0; j < z.cols(); j ++ )
      {
        locs_tmp(1, _) = rep(j + 1, tot_len); // first column belong to y values so starting from the second column 
        umat locs = as<umat>(wrap(locs_tmp));  
        // Rcout << "locs is (second) "  << locs << " size(yz) is " << size(yz) << std::endl;

        uvec eids = sub2ind( size(yz), locs ); // Obtain Element IDs
        // current_incoming_k = valid_top_k[j];
        tmp = z(Range(min(current_run_ind) + tau - z_delays[j], min(current_run_ind) + tau - z_delays[j] + tot_len - 1), Range(j, j));
        vals = as<vec>(wrap(tmp)); 
        // yz = join_cols(vals, yz);
        
        yz.elem( eids ) = vals;
        
      } 
      // Rcout << "pass yz join_rows iteration, yz (rows and columns) "  << yz.n_rows << ", " << yz.n_cols << std::endl;
     
      // Rcout << "pos is"  << pos << std::endl;
      // test = as<uvec>(wrap(d - 1));
      // x_0.elem(test) = as<vec>(wrap(1));

      tmp = x(Range(min(current_run_ind) + tau - d, min(current_run_ind) + tau - d + tot_len - 1), Range(0, 0));
      // Rcout << "locs is (third) "  << locs << " size(yz) is " << size(yz) << std::endl;
      vals = as<vec>(wrap(tmp));
      x_0 = join_cols(x_0, vals); //.elem(pos) = vals;
      
      tmp = y(Range(min(current_run_ind) + tau, min(current_run_ind) + tau + tot_len - 1), Range(0, 0));      
      vals = as<vec>(wrap(tmp));
      y_0 = join_cols(y_0, vals); //.elem(pos) = vals;    

      current_ind += tot_len; // update the index (note that there is no - 1)
    }

    x_0_res = as<NumericMatrix>(wrap(x_0));
    y_0_res = as<NumericMatrix>(wrap(y_0));
    // yz.reshape(y_0.n_rows, 2); // reshape yz matrix
    yz_res = as<NumericMatrix>(wrap(yz)); 

  }
  
  // Rcout << "before cmi_cpp " << std::endl;
  // Rcout << "dimension of each matrix " << "x_0_res dimension are " << x_0_res.rows()  << " ; " << x_0_res.cols() << " y_0_res dimension are " << y_0_res.rows() << " ; " << y_0_res.cols() << std::endl;
  // Rcout << "dimension of each matrix " << "yz_res dimension are " << yz_res.rows() << "; " << yz_res.cols()<< std::endl;
  List cmi_res = cmi_cpp(x_0_res, y_0_res, yz_res);
  return cmi_res["cmi_res"]; 
  // return cmi_cpp(x_0_res, y_0_res, yz_res)
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
//' @details
//' \code{rdi_single_run_conditioned} takes two random variables x and y as well as the parameter n to calculate the restricted direct information conditioned on variable z. 
//' @return a matrix for the condition mutual information estimators between all pairwise variables (x, y) in the data matrix x, y
//' @export
// [[Rcpp::export]]
double rdi_multiple_runs_conditioned(SEXP x, SEXP y, SEXP z, SEXP z_delays, SEXP d, SEXP run_vec) 
{ 
  NumericMatrix x_cpp(x); 
  NumericMatrix y_cpp(y); 
  NumericMatrix z_cpp(z); 
  
  IntegerVector z_delays_cpp(z_delays); 
  IntegerVector run_vec_cpp(run_vec); 
  
  int d_cpp = as<int>(d); 
  // (NumericMatrix& x, NumericMatrix& y, NumericMatrix& z, IntegerVector& z_delays, int d, IntegerVector run_vec) 
  double rdi_single_run_conditioned_res = rdi_multiple_runs_conditioned_cpp(x_cpp, y_cpp, z_cpp, z_delays_cpp, d_cpp, run_vec_cpp);
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
NumericMatrix calculate_multiple_run_conditioned_rdi_cpp(NumericMatrix& expr_data, IntegerMatrix& super_graph, 
                                            NumericMatrix& max_rdi_value, IntegerMatrix& max_rdi_delays, IntegerVector run_vec, int k) //, const int cores, const bool verbose
{
  Rcout << "max_rdi_value(gene_pair[0], gene_pair[1]) is " << max_rdi_value << std::endl;
  // if(verbose == TRUE) Rcout << "Calculating the conditional restricted direct mutual information for each pair of genes" << std::endl;
  NumericVector k_ncol = NumericVector::create(k, expr_data.cols() - 2); // conditioned at most n - 2 genes (2: two current testing genes)
  k = min(k_ncol); // minal value from k and the column 
  List top_incoming_nodes_delays_list = extract_top_incoming_nodes_delays(max_rdi_value, max_rdi_delays, k); // get the top k incoming node and the corresponding delays 

  IntegerMatrix top_incoming_nodes = top_incoming_nodes_delays_list["top_incoming_nodes"];
  IntegerMatrix top_incoming_delays = top_incoming_nodes_delays_list["top_incoming_delays"]; 
  NumericMatrix top_incoming_values = top_incoming_nodes_delays_list["top_incoming_values"];  

  Rcout << "top_incoming_nodes is " << top_incoming_nodes << std::endl;

  int n_genes = expr_data.cols(), n_sample = expr_data.rows(), tau, total_length, id3, delay_id3, current_id_tmp, delay; //tmp: current target  

  IntegerVector gene_pair(2), tmp(1), valid_top_k_in_k(k), valid_top_k_in_k_ordered, top_k_plus_1(k + 1), valid_top_k(k), valid_top_k_incoming_delays(k), valid_delay_tmp(1); // ind_for_top_k(k), little computational gain by setting the correct dimension 
  NumericMatrix cRDI_mat(n_genes, n_genes), expr1_t, expr2_t, past_tmp; 
  std::fill(cRDI_mat.begin(), cRDI_mat.end(), 0); // NA_REAL initialize the cRDI matrix as C ++ NA values (so the diagnoal will be NA values)
  
  // OpenMP here too (not thread safe) 
  for(int i = 0; i < super_graph.rows(); i ++) 
  {
   // index_name = super_graph(i, _); // Range(i, i)
    gene_pair = super_graph(i, _); 
    // Rcout << "before running rdi with gene_pair as " << gene_pair << std::endl;
    top_k_plus_1 = top_incoming_nodes(Range(gene_pair[1], gene_pair[1]), _); // top k incoming nodes 
    // Rcout << "top_k_plus_1 is " << top_k_plus_1 << std::endl;
    tmp = gene_pair[0]; // get the incoming node 

    if(intersect(top_k_plus_1, tmp).length() > 0) // if current node tmp is in the top_k_plus_1, remove tmp; else remove the smallest one 
    {
      valid_top_k = setdiff(top_k_plus_1, tmp);
    }
    else
    {
      Rcout << "length of top_k_plus_1 is " << top_k_plus_1.size() << "; k is " << k << std::endl;
      tmp = top_k_plus_1[k]; // remove the smallest one (k + 1 points in total)
      valid_top_k = setdiff(top_k_plus_1, tmp); 
    }
    
    // assign the node indices to ind_for_top_k
    // Rcout << "valid_top_k " << valid_top_k << "top_k_plus_1 " << top_k_plus_1 << std::endl;
    valid_top_k_in_k = match(valid_top_k, top_k_plus_1) - 1; // match gives order and the index (1 to k) starts from 1 -> valid_top_k_in_k
    // Rcout << "valid_top_k " << valid_top_k << "valid_top_k_in_k " << valid_top_k_in_k << std::endl;
    // valid_top_k = top_k_plus_1[valid_top_k_in_k]; // set valid_top_k to the same order as top_k_plus_1 by valid_top_k_in_k for downstream analysis
    valid_top_k_in_k_ordered = clone(valid_top_k_in_k).sort();  // sort the index valid_top_k_in_k 
    // Rcout << "after sort " << valid_top_k << "valid_top_k_in_k_ordered " << valid_top_k_in_k_ordered << std::endl; 
    valid_top_k = top_k_plus_1[valid_top_k_in_k_ordered]; 
    Rcout << "New sort " << valid_top_k << std::endl;

    if(valid_top_k.size() < k) 
    { // avoid any weird rare situation (all incoming RDI value are 0)
      // Rcout << "all incoming node has the same RDI value for gene pair " << gene_pair[0] << " and " << gene_pair[1] << ". Either of the gene may have no expression." << std::endl;
      cRDI_mat(gene_pair[0], gene_pair[1]) = max_rdi_value(gene_pair[0], gene_pair[1]); 
      // Rcout << "max_rdi_value(gene_pair[0], gene_pair[1]) is " << max_rdi_value << std::endl;

      continue;
    }

    // for(int j = 0; j < k; j ++ )
    // {
    //   ind_for_top_k[j] = top_k_plus_1[valid_top_k_in_k[j]]; // for example 10, 150, 3
    // }
    
    for(int j = 0; j < k; j ++ )
    {
      current_id_tmp = valid_top_k_in_k_ordered[j]; // range: 0 -> k - 1
      valid_delay_tmp = top_incoming_delays(Range(gene_pair[1], gene_pair[1]), Range(current_id_tmp, current_id_tmp)); // the current_id_tmp^{th} delay value 
      valid_top_k_incoming_delays[j] = valid_delay_tmp[0]; 
    }

    delay = max_rdi_delays(gene_pair[0], gene_pair[1]); // the delay corresponds to max RDI value for the current gene pair 
    // valid_top_k_incoming_delays.push_back(delay); // identify and use the biggest value for all delays (to avoid index over-flow); last value is not used 
    // tau = max(valid_top_k_incoming_delays); 
    // total_length = n_sample - tau; 

    // // Rcout << "current i is " << i << "max for top_incoming_delays is " << max(top_incoming_delays) << std::endl;
    // expr1_t = expr_data(Range(tau - delay, tau - 1 - delay+total_length), Range(gene_pair[0], gene_pair[0])); 
    // // Rcout << "tau" << tau << "tau - 1 + total_length" << tau - 1 + total_length << std::endl;
    // expr2_t = expr_data(Range(tau, tau - 1 + total_length), Range(gene_pair[1], gene_pair[1])); 
    // // Rcout << "tau - 1" << tau - 1 << "tau - 2 + total_length" << tau - 2 + total_length << std::endl;
    // NumericMatrix yz = expr_data(Range(tau - 1, tau  - 2 + total_length), Range(gene_pair[1], gene_pair[1])); // initialize yz as expr2_t with one time lag 

    // // OpenMP here too 
    // for(int id = 0; id < k; id ++)
    // {
    //   id3 = valid_top_k[id];
    //   delay_id3 = valid_top_k_incoming_delays[id];
    //   past_tmp = expr_data(Range(tau - delay_id3, tau  - delay_id3 - 1 + total_length),
    //                         Range(id3, id3));
    //   yz = cbind(yz, past_tmp);
    // }

    // // NumericMatrix past0 = transpose(yz);
    // cRDI_mat(gene_pair[0], gene_pair[1]) = cmi_cpp(expr1_t, expr2_t, yz);

    // do the multiple-run analysis here 

    // Rcout << "before running rdi" << std::endl;

    // use rdi_conditioned to calculate the values: 
    NumericMatrix x = expr_data(_, Range(gene_pair[0], gene_pair[0]));
    NumericMatrix y = expr_data(_, Range(gene_pair[1], gene_pair[1]));

    // prepare the data for all top-k incoming nodes' expression data 
    NumericMatrix z(n_sample, k);
    for (int i = 0; i < valid_top_k.size(); i++) {
      z(_,i) = expr_data(_, valid_top_k(i));
    }

    // // // x, y, z, valid_top_k_incoming_delays, delay, run_vec
    Rcout << "x is (before rdi_multiple_runs_conditioned_cpp) " << x << ", " << x.cols(); 
    // Rcout << "y is (before rdi_multiple_runs_conditioned_cpp) " << y << ", " << x.cols(); 
    // Rcout << "z is (before rdi_multiple_runs_conditioned_cpp) " << z << ", " << x.cols(); 
    // Rcout << "\n y is (row, column) " << y.rows() << ", " << y.cols(); 
    // Rcout << "\n z is (row, column) " << z.rows() << ", " << z.cols() ; 
    // Rcout << "\n valid_top_k_incoming_delays is " << valid_top_k_incoming_delays; 
    // Rcout << "\n delay is " << delay; 
    // Rcout << "\n run_vec is " << run_vec << endl; 

    cRDI_mat(gene_pair[0], gene_pair[1])  = rdi_multiple_runs_conditioned_cpp(x, y, z, valid_top_k_incoming_delays, delay, run_vec); 
    Rcout << "after running rdi with gene_pair as " << gene_pair << "the current crDI is" << cRDI_mat(gene_pair[0], gene_pair[1]) << "results is " << std::endl;
  }
  
  // return results
  return cRDI_mat; 
}

// [[Rcpp::export]]
NumericMatrix calculate_multiple_run_conditioned_rdi_wrap(SEXP expr_data, SEXP super_graph, SEXP max_rdi_value, SEXP max_rdi_delays, SEXP run_vec, SEXP k) //, SEXP R_cores
{ 
  // Rcout << "before running calculate_rdi_cpp_wrap" << std::endl;
  NumericMatrix expr_data_cpp(expr_data), max_rdi_value_cpp(max_rdi_value); 
  IntegerMatrix super_graph_cpp(super_graph), max_rdi_delays_cpp(max_rdi_delays); 

  IntegerVector run_vec_cpp(run_vec); 
  int k_cpp = as<int>(k); 
  // Rcout << "before running calculate_rdi_cpp" << std::endl;
  NumericMatrix crdi_res = calculate_multiple_run_conditioned_rdi_cpp(expr_data_cpp, super_graph_cpp, max_rdi_value_cpp, max_rdi_delays_cpp, run_vec_cpp, k_cpp); //cores
  return crdi_res;
}

/** 
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
 
Var2 Var1
2      0    1
3      0    2
4      0    3
5      0    4
6      0    5
7      0    6
8      0    7
9      0    8
10     0    9
11     0   10
12     0   11
13     1    0
15     1    2
16     1    3
17     1    4
18     1    5
19     1    6
20     1    7
21     1    8
22     1    9
23     1   10
24     1   11
25     2    0
26     2    1
28     2    3
29     2    4
30     2    5
31     2    6
32     2    7
33     2    8
34     2    9
35     2   10
36     2   11
37     3    0
38     3    1
39     3    2
41     3    4
42     3    5
43     3    6
44     3    7
45     3    8
46     3    9
47     3   10
48     3   11
49     4    0
50     4    1
51     4    2
52     4    3
54     4    5
55     4    6
56     4    7
57     4    8
58     4    9
59     4   10
60     4   11
61     5    0
62     5    1
63     5    2
64     5    3
65     5    4
67     5    6
68     5    7
69     5    8
70     5    9
71     5   10
72     5   11
73     6    0
74     6    1
75     6    2
76     6    3
77     6    4
78     6    5
80     6    7
81     6    8
82     6    9
83     6   10
84     6   11
85     7    0
86     7    1
87     7    2
88     7    3
89     7    4
90     7    5
91     7    6
93     7    8
94     7    9
95     7   10
96     7   11
97     8    0
98     8    1
99     8    2
100    8    3
101    8    4
102    8    5
103    8    6
104    8    7
106    8    9
107    8   10
108    8   11
109    9    0
110    9    1
111    9    2
112    9    3
113    9    4
114    9    5
115    9    6
116    9    7
117    9    8
119    9   10
120    9   11
121   10    0
122   10    1
123   10    2
124   10    3
125   10    4
126   10    5
127   10    6
128   10    7
129   10    8
130   10    9
132   10   11
133   11    0
134   11    1
135   11    2
136   11    3
137   11    4
138   11    5
139   11    6
140   11    7
141   11    8
142   11    9
143   11   10 
*/

/**
 // run multiple run on the real simulation datasets: 
 

 */











