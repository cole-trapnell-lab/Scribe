// ver 0.1; code review at Oct 2, 2017
// ver 0.2; update the documentation at Oct 2, 2017

#include <RcppArmadillo.h>
#include <boost/math/special_functions/digamma.hpp>
#include <RANNinf.h>

// flann package
#include <string>

#include "../inst/include/information_estimator.h"
#include "../inst/include/density.h"

// #include <ANN/ANN.h> // ANN declarations
// #include <omp.h>
// #include "../inst/include/ann_neighbor_radius.h"
// #include "../inst/include/radius_search.h"
// #include "../inst/include/neighbour.h"

using namespace Rcpp;
using namespace arma;
using namespace boost::math;
using namespace RANNinf;
// using namespace Rcpp::stat

// [[Rcpp::plugins(openmp)]]

// [[Rcpp::depends(RcppArmadillo)]]

//-------------------------------------------------------------------------------------------------------------------
/*
calculate local density  
*/

// knn_dis = [tree_x.query(point, k_density + 1, p=np.inf)[0][k_density] for point in x]
// density_estimate = np.array([float(k_density) / N / knn_dis[i] ** dx for i in range(len(knn_dis))])
//   weight = (1 / density_estimate) / np.mean(1 / density_estimate)

List knn_density_cpp(NumericMatrix x, NumericMatrix y, int k = 5)
{
  int dx = x.cols();
  int N = y.rows(); // number of query points
  NumericVector weight(N), density_estimate(N); 

  int dimension = N;  // k *
  Rcpp::NumericVector data_y_distances(dimension); 

  double error_bound = 0.0; int searchtype = 1; int usebdtree = 0; double sqRad = 0.0; // don't use bd tree (usebdtree = 0)

  get_NN_2Set_cpp(x, y, dx, N, N, k, error_bound, searchtype, usebdtree, sqRad, data_y_distances); //data_xyz_nn_index,
  
  // NumericVector vd = 0.5 * d * log(M_PI) - log(gamma(0.5 * d + 1)) + log(data_xyz_distances); // 
  // density_estimate = k / (N * exp(vd)); // this is for the L2 norm (the one below is the infinity norm) 
  density_estimate = k / (N * pow(data_y_distances, dx)); 
  
  weight = (1 / density_estimate / mean(1 / density_estimate));
  
  return List::create(Rcpp::Named("density_estimate") = density_estimate, 
                      Rcpp::Named("weight") = weight);
}

//' @title
//' knn_density
//' @description
//' This subroutine calculates the 1d density of x at positions y 
//' 
//' @param x 1d vector of the data 
//' @param y 1d vector of querying points, positions used to estimate the density 
//' @param k number of nearest neighbors used to estimate the 1d density, default to be 5 
//' 
//' @details
//' \code{knn_density} takes a vector of original data points x and a vector of querying points y
//' to calculate the density at each point y using the k-nearest neighbors. 
//' @return a List where the element is the density estimate (name: density_estimate), 
//' the second one is the weight calculated based on density_estimate.  
//' @export
// [[Rcpp::export]]
List knn_density(SEXP x, SEXP y, SEXP k) //&
{ 
  NumericMatrix x_cpp = as<int>(x); 
  NumericMatrix y_cpp = as<int>(y); 
  int k_cpp = as<int>(k); 

  List knn_density_res = knn_density_cpp(x_cpp, y_cpp, k_cpp);
  return knn_density_res;
}

NumericMatrix knn_density_2d_cpp(NumericVector x, NumericVector y, IntegerVector nGrids, int k = 5)
{
  double x_min = min(x), x_max = max(x);
  double y_min = min(y), y_max = max(y);
  double x_step = (x_max - x_min) / (nGrids[0]), y_step = (y_max - y_min) / (nGrids[1]);

  NumericMatrix x_y_grid(nGrids[0] * nGrids[1], 2);
  NumericMatrix x_y_density(nGrids[0], nGrids[1]);

  for(int i = 0; i < nGrids[0]; i++) 
  {
    for(int j = 0; j < nGrids[1]; j++)
    {
      // make sure we use the middle of each grid 
      x_y_grid(i * nGrids[0] + j, _) = NumericVector::create(x_min + x_step * ((double) i + 0.5), 
                                          y_min + y_step * ((double) j + 0.5) ); 
    }
  }

  NumericMatrix xy = cbind(x, y);
  int dxy = 2;
  int N = xy.rows();
  int NQ = nGrids[0] * nGrids[1]; 

  Rcpp::NumericVector data_xy_distances(nGrids[0] * nGrids[1]), density_estimate(data_xy_distances); 

  double error_bound = 0.0; int searchtype = 1; int usebdtree = 0; double sqRad = 0.0; // don't use bd tree (usebdtree = 0)

  get_NN_2Set_cpp(xy, x_y_grid, dxy, N, NQ, k, error_bound, searchtype, usebdtree, sqRad, data_xy_distances); //data_xyz_nn_index,
  
  // NumericVector vd = 0.5 * d * log(M_PI) - log(gamma(0.5 * d + 1)) + log(data_xyz_distances); // 
  // density_estimate = k / (N * exp(vd)); // this is for the L2 norm (the one below is the infinity norm) 
  density_estimate = k / (N * pow(data_xy_distances, dxy));  

  for(int i = 0; i < nGrids[0]; i++) 
  {
    for(int j = 0; j < nGrids[1]; j++)
    {
      x_y_density(i, j) = density_estimate[i * nGrids[0] + j]; 
    }
  }

  return(x_y_density);
}

// enable Gaussian smoothing 
//' @title
//' knn_density_2d
//' @description
//' This subroutine calculates the density for a 2d space. 
//' 
//' @param x A vector for the values of the data on the first dimension 
//' @param y A vector for the values of the data on the second dimension 
//' @param nGrids A vector of two for the grid numbers on the first and second dimension 
//' @param k number of nearest neighbors used to calculate the 2d density 
//' 
//' @details
//' \code{knn_density_2d} 
//' @return a numeric value for the d-dimensional unit ball for Euclidean norm
//' @export a matrix of density estimate, calculated on the center of each grid from the data x and y. 
// [[Rcpp::export]]
NumericMatrix knn_density_2d(SEXP x, SEXP y, SEXP nGrids, SEXP k) //&
{ 
  NumericVector x_cpp = as<int>(x); 
  NumericVector y_cpp = as<int>(y); 
  IntegerVector nGrids_cpp = as<int>(nGrids); 
  int k_cpp = as<int>(k); 

  NumericMatrix knn_density_2d_res = knn_density_2d_cpp(x_cpp, y_cpp, nGrids_cpp, k_cpp);
  return knn_density_2d_res;
}

/*
 knn_density(1:10, 5)
 knn_density(matrix(1:100, ncol = 1), 5)
 
 density(1:100, bw = 1, n = 100, from = 1, to = 100)$y
 
 cpp <- knn_density(matrix(1:100, ncol = 1), 5)
 Rversion <- 1 / density(1:100, bw = 0.2, n = 100, from = 1, to = 100)$y / mean(1 / density(1:100, bw = 0.2, n = 100, from = 1, to = 100)$y)
 
 qplot(cpp, Rversion)
 
// NumericVector Cquantile(NumericVector x, NumericVector q) 
// {
//   NumericVector y = clone(x);
//   std::sort(y.begin(), y.end());
//   return y[x.size()*(q - 0.000000001)];
// }

// NumericVector Cquantile(NumericVector x, NumericVector probs) {
//   Environment stats("package:stats");
//   Function quantile = stats["quantile"];
//   int npr = probs.size();
//   NumericVector ans(npr);
//   for(int i=0; i<npr; i++){
//     ans[i] = as<double>(quantile(x, probs[i]));
//   }
//   return ans;
// }

// double bandwidth_nrd_cppfunction(NumericVector x)
// {
//     NumericVector r = Cquantile(x, NumericVector::create(0.25, 0.75));
//     double h = (r[2] - r[1])/1.34;
//     double res = 4 * 1.06 * min(NumericVector::create(sqrt(var(x)), h)) * pow(x.size(), (-1/5));

//     return(res);
// }
å
 */

// [[Rcpp::export]]
double digamma_0(double x) // this modified digamma function avoids the undefined case for 0 
{
  if(x == 0)
  {
   return digamma(1); // arma::datum::inf
  }
  else {
    return digamma(x);
  }
}
//-------------------------------------------------------------------------------------------------------------------
/*
  calculate vd 
*/
//-------------------------------------------------------------------------------------------------------------------

//'  This subroutine calculates the volume of a d-dimensional unit ball for Euclidean norm
double vd_cpp(const int d)
{
  double vd_res = 0.5 * d * log(M_PI) - log(gamma(0.5 * d + 1));
  return vd_res;
}

//' @title
//' vd
//' @description
//' This subroutine calculates the volume of a d-dimensional unit ball for Euclidean norm
//' 
//' @param d number of dimension
//' 
//' @details
//' \code{vd} takes a integer of dimensions and then calculate the volume of a d-dimensional unit ball for Euclidean norm
//' using the formula: 0.5 * d * log(pi) - log(gamma(0.5 * d + 1))
//' It's implimented in C++, providing a (small) increase in speed over the R equivalent.
//' @return a numeric value for the d-dimensional unit ball for Euclidean norm
//' @export
// [[Rcpp::export]]
double vd(SEXP d) //&
{ 
  int d_cpp = as<int>(d); 
  double vd_res = vd_cpp(d_cpp);
  return vd_res;
}

//-------------------------------------------------------------------------------------------------------------------
/*
  calculate entropy 
*/
//-------------------------------------------------------------------------------------------------------------------

//' This function estimates the entropy of a continuous random variable
double entropy_cpp(const NumericMatrix& x, int k) // = 5 NumericVector
{
  int N = x.rows();
  if(k > N)
  {
    stop("k is larger than total number of samples");
    return -1;
  }
  
  int d = x.cols();
  
  int dimension = N;  //k *
  Rcpp::NumericVector distances(dimension);

  // Rcpp::IntegerVector nn_index; 
  // nn_index = IntegerVector(dimension);

  // use ANN: 
  double error_bound = 0.0; int searchtype = 1; int usebdtree = 0; double sqRad = 0.0; 
  
  get_NN_2Set_cpp(x, x, d, N, N, k, error_bound, searchtype, usebdtree, sqRad, distances); //, *nn_index, *distances
  
  // // use flann: 
  // int cores = 1, checks = 1; std::string build("kdtree");  
  // NumericVector distances = Neighbour(x, k, build, cores, checks); 

  NumericVector log_knn_dist(N);
  for(int i = 0; i < N; i++)
  {
    log_knn_dist[i] = log(distances[i]); //(i + 1) * k - 1
  }

  double entropy = - digamma_0(k) + digamma_0(N) + d * mean(log_knn_dist); 
  return entropy;
}

//' @title
//' entropy
//' @description
//' This subroutine estimates the entropy of a continuous random variable
//' 
//' @param x data matrix used for calculating the entropy
//' 
//' @param k number for nearest neighbors used in entropy calculation
//'
//' @details
//' \code{entropy} takes a continuous random variable and then estimates
//' entropy using the KSG estimator. 
//' It relies on the ANN package to query the kNN with KDTree algorithm.  
//' @return a numeric value of entropy estimate
//' @export
// [[Rcpp::export]]
double entropy(SEXP x, SEXP k) //
{ 
  NumericMatrix x_cpp(x);
  int k_cpp = as<int>(k);

  double entropy = entropy_cpp(x_cpp, k_cpp); 
  return(entropy);
}

//-------------------------------------------------------------------------------------------------------------------
/*
  calculate mutual information 
*/
//-------------------------------------------------------------------------------------------------------------------

// This function estimates the mutual information of two random variables based on their observed values
double mi_cpp(const NumericMatrix& x, const NumericMatrix& y, int k, int normalize) 
{

  k = k + 1; // the first nearest point is itself in ANN

  int N = x.rows();
  NumericMatrix data = cbind(x, y); int d_data = data.cols();
  
	int dx = x.cols();	int dy = y.cols(); 

	if(k > N)
  {
    stop("k is larger than total number of samples");
    return -1;
  }
  if(y.rows() != N)
  {
    stop("Number of samples should be the same");
    return -1;
  }

  // use ANN 
  double error_bound = 0.0; int searchtype = 1; int usebdtree = 0; double sqRad = 0.0; int x_y_k = 0;

  int dimension =  N;  //k *
  Rcpp::NumericVector data_distances(dimension); 
  
  get_NN_2Set_cpp(data, data, d_data, N, N, k, error_bound, searchtype, usebdtree, sqRad,  data_distances); //data_nn_index,
  
  // // use flann: 
  // int cores = 1, checks = 1; std::string build("kdtree"); 
  // NumericVector data_distances = Neighbour(data, k, build, cores, checks); 
 
	NumericVector information_samples(N, digamma_0(N));

  NumericVector k_xy = get_points_in_radius_cpp(data, data, d_data, N, N, x_y_k, error_bound, usebdtree, data_distances);
  NumericVector cnt_x = get_points_in_radius_cpp(x, x, dx, N, N, x_y_k, error_bound, usebdtree, data_distances);
  NumericVector cnt_y = get_points_in_radius_cpp(y, y, dy, N, N, x_y_k, error_bound, usebdtree, data_distances); 
  
  // remove it self 
  k_xy = k_xy - 1;
  cnt_x = cnt_x - 1;
  cnt_y = cnt_y - 1; 
  // // use flann: 
  // NumericVector k_xy = RadiusSearch(data, data_distances, build, cores, checks); 
  // NumericVector cnt_x = RadiusSearch(x, data_distances, build, cores, checks); 
  // NumericVector cnt_y = RadiusSearch(y, data_distances, build, cores, checks); 

  int i; 

  // int max_thread = omp_get_max_threads();
  // omp_set_num_threads(max_thread);
  // #pragma omp parallel for shared(information_samples, N, cnt_x, cnt_y) private(i) //schedule(dynamic) default(none) //collapse(2) , _
  for(i = 0; i < N; i++)
  {
    information_samples[i] += digamma_0(k_xy[i]) - digamma_0(cnt_x[i]) - digamma_0(cnt_y[i]);
  }

  NumericVector weight(N);
  weight = weight + 1.0;

  List knn_density_res; 

  if(normalize != 0) {
    knn_density_res = knn_density_cpp(x, x, k); 
    weight = knn_density_res["weight"];
  }

	double mi_res = mean(information_samples * weight); //

	return mi_res;
}

//' @title
//' mi
//' @description
//' This function estimates the mutual information of two random variables, x, y, based on their observed values
//' 
//' @param x one random variable from the time-series data
//' 
//' @param y another random variable from the time-series data
//'
//' @param k number for nearest neighbors used in entropy calculation
//'
//' @param normalize A logic flag to determine whether or not you want to normalize the MI value by the density of x. 
//'
//' @details
//' \code{mi} takes two random variables x and y to estimate the mutual information between them 
//' using the KSG estimator
//' It relies on the ANN package to query the kNN with KDTree algorithm. 
//' @return a estimated mutual information value between two variables (x, y)
//' @export
// [[Rcpp::export]]
double mi(SEXP x, SEXP y, SEXP k, SEXP normalize) //&
{ 
  NumericMatrix x_cpp(x); 
  NumericMatrix y_cpp(y); 
  int k_cpp = as<int>(k);
  int normalize_cpp = as<int>(normalize);

  double mi_res = mi_cpp(x_cpp, y_cpp, k_cpp, normalize_cpp);

  return mi_res;
}

//-------------------------------------------------------------------------------------------------------------------
/*
  calculate conditional mutual information 
*/
//-------------------------------------------------------------------------------------------------------------------

//' This function estimates the CONDITIONAL mutual information of X and Y given Z
//' Multiply a number by two
//'
List cmi_cpp(const NumericMatrix& x, const NumericMatrix& y, NumericMatrix z, int k, int normalize){ // const NumericMatrix& z, 
  
  k = k + 1; 	

  int N = x.rows();
  if(y.rows() != N | z.rows() != N)
  {
    stop("all Matrix should include the same number of samples");
  }

  NumericMatrix data_xz = cbind(x, z);
  NumericMatrix data_yz = cbind(y, z);
  NumericMatrix data_xyz = cbind(cbind(x, y), z);

  int dxz = data_xz.cols();  int dyz = data_yz.cols(); int dz = z.cols(); int d_data_xyz = data_xyz.cols();

  double error_bound = 0.0; int searchtype = 1; int usebdtree = 0; double sqRad = 0.0; int cmbn_k = 0; // don't use bd tree (usebdtree = 0)

  int dimension = N;  // number of cells 
  Rcpp::NumericVector data_xyz_distances(dimension); 

  get_NN_2Set_cpp(data_xyz, data_xyz, d_data_xyz, N, N, k, error_bound, searchtype, usebdtree, sqRad, data_xyz_distances); //data_xyz_nn_index,

  NumericVector information_samples(N);

  NumericVector k_xyz = get_points_in_radius_cpp(data_xyz, data_xyz, d_data_xyz, N, N, cmbn_k, error_bound, usebdtree, data_xyz_distances);
  NumericVector cnt_xz = get_points_in_radius_cpp(data_xz, data_xz, dxz, N, N, cmbn_k, error_bound, usebdtree, data_xyz_distances);
  NumericVector cnt_yz = get_points_in_radius_cpp(data_yz, data_yz, dyz, N, N, cmbn_k, error_bound, usebdtree, data_xyz_distances); 

  NumericVector cnt_z = get_points_in_radius_cpp(z, z, dz, N, N, cmbn_k, error_bound, usebdtree, data_xyz_distances); 

  k_xyz = k_xyz - 1;
  cnt_xz = cnt_xz - 1; 
  cnt_yz = cnt_yz - 1; 
  cnt_z = cnt_z - 1; 

  int i;
  // int max_thread = omp_get_max_threads();
  // omp_set_num_threads(max_thread);
  // #pragma omp parallel for shared(information_samples, N, cnt_xz, cnt_yz, cnt_z) private(i) //schedule(dynamic) default(none) //collapse(2) , _
  for(i = 0; i < N; i++)
  {
    if(k_xyz[i] == 0 | cnt_xz[i] == 0 | cnt_yz[i] == 0 | cnt_z[i] == 0)
      continue; 
    
    information_samples[i] += digamma_0(k_xyz[i]) - digamma_0(cnt_xz[i]) - digamma_0(cnt_yz[i]) + digamma_0(cnt_z[i]); 
  }

  NumericVector weight(N);
  weight = weight + 1.0;

  List knn_density_res; 

  if(normalize != 0) {
    knn_density_res = knn_density_cpp(x, x, k);  // try data_xz
    weight = knn_density_res["weight"]; 
  }
  
  double cmi_res = mean(information_samples * weight);//mean(information_samples);

  // if(arma::is_finite(cmi_res) == false) 
  // {
  //   Rcout << "weight is " << weight << std::endl;
  //   Rcout << "cmi_res is " << cmi_res << std::endl;
  //   stop("cmi_res is NaN or NA values!");
  // }

  return List::create(Rcpp::Named("cmi_res") = cmi_res, 
            Rcpp::Named("information_samples") = information_samples); 
}

//' @title
//' cmi
//' @description
//' This subroutine calculates the volume of a d-dimensional unit ball for Euclidean norm
//' 
//' @param x one random variable from the time-series data
//' 
//' @param y another random variable from the time-series data
//'
//' @param z condition random variable for variables (x, y) from the time-series data
//'
//' @param k number for nearest neighbors used in entropy calculation
//'
//' @param normalize A logic flag to determine whether or not you want to normalize the MI value by the density of x. 
//'
//' @details
//' \code{cmi} takes two random variable x and y and estimated their mutual information conditioned on the third random variable z
//' using the KSG estimator. 
//' It relies on the ANN package to query the kNN with KDTree algorithm. 
//' @return a estimated conditional mutual information value between two variables (x, y), conditioning on a third variable z. 
//' @export
// [[Rcpp::export]]
List cmi(SEXP x, SEXP y, SEXP z, SEXP k, SEXP normalize) 
{ 
  NumericMatrix x_cpp(x); 
  NumericMatrix y_cpp(y); 
  NumericMatrix z_cpp(z); 

  int k_cpp = as<int>(k); 
  int normalize_cpp = as<int>(normalize); 

  List cmi_res = cmi_cpp(x_cpp, y_cpp, z_cpp, k_cpp, normalize_cpp);
  return cmi_res;
}

// List ucmi_res = ucmi_cpp(x_cpp, y_cpp, z_cpp, k_cpp, method_cpp, k_density_cpp, bw_cpp);
// [[Rcpp::interfaces(r, cpp)]]
// [[Rcpp::export]]
List ucmi_cpp(const NumericMatrix& x, const NumericMatrix& y, NumericMatrix z, int k, int method, int k_density, double bw){ // , NumericVector weight const NumericMatrix& z, 
  
  k = k + 1;  

  int N = x.rows();
  if(y.rows() != N | z.rows() != N)
  {
    stop("all Matrix should include the same number of samples");
  }

  NumericMatrix data_xz = cbind(x, z);
  NumericMatrix data_yz = cbind(y, z);
  NumericMatrix data_xyz = cbind(cbind(x, y), z);

  int dxz = data_xz.cols();  int dyz = data_yz.cols(); int dz = z.cols(); int d_data_xyz = data_xyz.cols();

  double error_bound = 0.0; int searchtype = 1; int usebdtree = 0; double sqRad = 0.0; int cmbn_k = 0; // don't use bd tree (usebdtree = 0)

  // kde_cpp(NumericMatrix data, int k = 1, int b = 1, int pdf = 1, int density_sample_type = 1)
  NumericVector weight; //(N)
  
  if(method == 1) {
    int k = 1, b = 1, pdf = 1, density_sample_type = 1;
    NumericMatrix kde_mat; // NumericMatrix kde_cpp
    
    kde_mat = kde_cpp(data_xz, k, b, pdf, density_sample_type); // get the density
    weight = kde_mat(_, 0); // get the density estimator
    // weight = exp(weight);

    weight = (1 / weight) / mean(1 / weight); // convert the density estimator into weight estimate
  } else if(method == 2) {
    List knn_density_res; 
    knn_density_res = knn_density_cpp(x, x, k);  // try data_xz
    weight = knn_density_res["weight"]; 
  }

  int dimension = N;  // number of cells 
  Rcpp::NumericVector data_xyz_distances(dimension); 

  get_NN_2Set_cpp(data_xyz, data_xyz, d_data_xyz, N, N, k, error_bound, searchtype, usebdtree, sqRad, data_xyz_distances); //data_xyz_nn_index,

  NumericVector information_samples(N);

  NumericVector k_xyz = get_points_in_radius_cpp(data_xyz, data_xyz, d_data_xyz, N, N, cmbn_k, error_bound, usebdtree, data_xyz_distances);
  NumericVector cnt_xz = get_points_in_radius_cpp(data_xz, data_xz, dxz, N, N, cmbn_k, error_bound, usebdtree, data_xyz_distances);
  
  List yz_list = get_points_indices_in_radius_cpp(data_yz, data_yz, dyz, N, N, cmbn_k, error_bound, usebdtree, data_xyz_distances); 
  IntegerMatrix yz_neighbors = yz_list["nn_index"];
  NumericVector yz_num_cnts = yz_list["num_cnts"]; 
  
  List z_list = get_points_indices_in_radius_cpp(z, z, dz, N, N, cmbn_k, error_bound, usebdtree, data_xyz_distances); 
  IntegerMatrix z_neighbors = z_list["nn_index"]; 
  NumericVector z_num_cnts = z_list["num_cnts"]; 
  
  k_xyz = k_xyz - 1;
  cnt_xz = cnt_xz - 1; 
  // cnt_yz = cnt_yz - 1; 
  // cnt_z = cnt_z - 1; 

  int i, current_neighbor_cnt;
  IntegerVector yz_neighbors_ind, z_neighbors_ind; 
  weight[is_na(weight)] = 1; // convert the na values into 1  

  Rcout << "weight vector is " << weight << std::endl;
  // int max_thread = omp_get_max_threads();
  // omp_set_num_threads(max_thread);
  // #pragma omp parallel for shared(information_samples, N, cnt_xz, cnt_yz, cnt_z) private(i) //schedule(dynamic) default(none) //collapse(2) , _
  // std::ofstream file("test_ucmi_cpp.txt");
  
  for(i = 0; i < N; i++)
  {
    if(k_xyz[i] == 0 | cnt_xz[i] == 0) // | cnt_yz[i] == 0 | cnt_z[i] == 0
      continue; 
    
    // information_samples[i] += weight[i]* digamma(len(tree_xyz.query_ball_point(data_xyz[i], knn_dis[i], p=np.inf )) -1)
    // information_samples[i] += weight[i]* -digamma(len(tree_xz.query_ball_point(data_xz[i], knn_dis[i], p=np.inf )) - 1)
    // information_samples[i] += weight[i]* -digamma( np.sum( weight[j] for j in tree_yz.query_ball_point(data_yz[i], knn_dis[i], p=np.inf )) - weight[i])
    // information_samples[i] += weight[i]* digamma( np.sum( weight[j] for j in tree_z.query_ball_point(z[i], knn_dis[i], p=np.inf)) - weight[i])
    current_neighbor_cnt = yz_num_cnts[i] - 1;
    yz_neighbors_ind = yz_neighbors(i, _); 
    yz_neighbors_ind = yz_neighbors_ind[Range(0, current_neighbor_cnt)]; 
    
    current_neighbor_cnt = z_num_cnts[i] - 1;
    z_neighbors_ind = z_neighbors(i, _); 
    z_neighbors_ind = z_neighbors_ind[Range(0, current_neighbor_cnt)];

    NumericVector sub_weight_yz = weight[yz_neighbors_ind]; double digamma_0_yz = sum(sub_weight_yz) - weight[i];
    NumericVector sub_weight_z = weight[z_neighbors_ind]; double digamma_0_z = sum(sub_weight_z) - weight[i];
    
    // if (file.is_open() && i < 5)
    // {
    //   file << "current yz_neighbors_ind ind is " << yz_neighbors_ind << '\n';
    //   file << "current z_neighbors_ind ind is " << z_neighbors_ind << '\n';
    //   file << "weight is " << weight << "\n"; 
    //   file << "sub_weight_yz is " << sub_weight_yz << "\n"; 
    //   file << "sub_weight_z is " << sub_weight_z << "\n"; 
    //   file << "digamma_0_yz is " << digamma_0_yz << "\n"; 
    //   file << "digamma_0_z is " << digamma_0_z << "\n"; 
    // }
    
    information_samples[i] += weight[i] * digamma_0(k_xyz[i]); 
    information_samples[i] -= weight[i] * digamma_0(cnt_xz[i]);
    information_samples[i] -= weight[i] * digamma_0(digamma_0_yz);
    information_samples[i] += weight[i] * digamma_0(digamma_0_z); 
  }

  double ucmi_res = mean(information_samples); //mean(information_samples);

  // if(arma::is_finite(cmi_res) == false) 
  // {
  //   Rcout << "weight is " << weight << std::endl;
  //   Rcout << "cmi_res is " << cmi_res << std::endl;
  //   stop("cmi_res is NaN or NA values!");
  // }

  return List::create(Rcpp::Named("ucmi_res") = ucmi_res, 
            Rcpp::Named("information_samples") = information_samples,
            Rcpp::Named("weight") = weight); 
}

//' @title
//' ucmi
//' @description
//' This subroutine calculates the uniformed conditional mutual information where 
//' the distribution for x and z is replaced by a uniform distribution.  
//' 
//' @param x one random variable from the time-series data
//' @param y another random variable from the time-series data
//' @param z condition random variable for variables (x, y) from the time-series data
//' @param k number for nearest neighbors used in entropy calculation
//' @param method Which 2D density estimator you would like to use. 1 is kde estimator and 2 is knn based estimator. 
//' @param k_density The number of k nearest neighbors you would like to use when calculating the density 
//' (only applicable when method == 2 or using knn based density estimation). 
//' @param bw Bindwidth used for the kernel density estimator. Currently it is not used. The bindwidth in the kde function is automatically estimated. 
//' 
//' @details
//' \code{ucmi} takes two random variable x and y and estimated their mutual information conditioned on the third random variable z
//' using the KSG estimator while x, y is replaced by a uniform distribution. It relies on a C++ implmentation of kde estimator 
//' (https://github.com/timnugent/kernel-density) and the ANN package to query the kNN with KDTree algorithm. 
//' @return  a estimated conditional mutual information value between two variables (x, y), conditioning on a third variable z where 
//' the distribution for the x, z is replaced by a uniform distribution. 
//'  
//' @export
// [[Rcpp::export]]
List ucmi(SEXP x, SEXP y, SEXP z, SEXP k, SEXP method, SEXP k_density, SEXP bw) //, SEXP weight
{ 
  NumericMatrix x_cpp(x); 
  NumericMatrix y_cpp(y); 
  NumericMatrix z_cpp(z); 

  int k_cpp = as<int>(k); 

  int method_cpp = as<int>(method); 
  int k_density_cpp = as<int>(k_density);
  double bw_cpp = as<double>(bw); 

  // NumericVector weight_cpp = as<NumericVector>(weight); 
  
  List ucmi_res = ucmi_cpp(x_cpp, y_cpp, z_cpp, k_cpp, method_cpp, k_density_cpp, bw_cpp); //, weight_cpp
  return ucmi_res;
}

// [[Rcpp::interfaces(r, cpp)]]
// [[Rcpp::export]]
double umi_cpp(const NumericMatrix& x, const NumericMatrix& y, int k, int method, int k_density, double bw) // , NumericVector weight
{ 
  
  k = k + 1; // the first nearest point is itself in ANN

  int N = x.rows();
  NumericMatrix data = cbind(x, y); int d_data = data.cols();
  
  int dx = x.cols();  int dy = y.cols(); 

  if(k > N)
  {
    stop("k is larger than total number of samples");
    return -1;
  }
  if(y.rows() != N)
  {
    stop("Number of samples should be the same");
    return -1;
  }

  // kde_cpp(NumericMatrix data, int k = 1, int b = 1, int pdf = 1, int density_sample_type = 1)
  NumericVector weight; //(N)

  if(method == 1) {
    int k = 1, b = 1, pdf = 1, density_sample_type = 1;
    NumericMatrix kde_mat; // NumericMatrix kde_cpp

    kde_mat = kde_cpp(x, k, b, pdf, density_sample_type); // get the density
    weight = kde_mat(_, 0); // get the density estimator
    // weight = exp(weight);

    weight = (1 / weight) / mean(1 / weight); // convert the density estimator into weight estimate
  } else if(method == 2) {
    List knn_density_res;
    knn_density_res = knn_density_cpp(x, x, k);  // try data_xz
    weight = knn_density_res["weight"];
  }

  // use ANN 
  double error_bound = 0.0; int searchtype = 1; int usebdtree = 0; double sqRad = 0.0; int x_y_k = 0;

  int dimension =  N;  
  Rcpp::NumericVector data_distances(dimension); 
  
  get_NN_2Set_cpp(data, data, d_data, N, N, k, error_bound, searchtype, usebdtree, sqRad,  data_distances); //data_nn_index,
  
  NumericVector k_xy = get_points_in_radius_cpp(data, data, d_data, N, N, x_y_k, error_bound, usebdtree, data_distances);
  NumericVector cnt_x = get_points_in_radius_cpp(x, x, dx, N, N, x_y_k, error_bound, usebdtree, data_distances);
  // NumericVector cnt_y = get_points_in_radius_cpp(y, y, dy, N, N, x_y_k, error_bound, usebdtree, data_distances); 

  List y_list = get_points_indices_in_radius_cpp(y, y, dy, N, N, x_y_k, error_bound, usebdtree, data_distances); 
  IntegerMatrix y_neighbors = y_list["nn_index"];
  NumericVector y_num_cnts = y_list["num_cnts"]; 
    
  // remove it self 
  k_xy = k_xy - 1;
  cnt_x = cnt_x - 1;

  NumericVector weight_y(N); //(N)

  double ans = digamma_0(k - 1) + 2 * log(N - 1) - digamma_0(N); // + vd(dx) + vd(dy) - vd(dx + dy) 

  int i, current_neighbor_cnt, nx;
  double ny;
  IntegerVector y_neighbors_ind;
  NumericVector tmp; 
  
  // int max_thread = omp_get_max_threads();
  // omp_set_num_threads(max_thread);
  // #pragma omp parallel for shared(information_samples, N, cnt_x, cnt_y) private(i) //schedule(dynamic) default(none) //collapse(2) , _
  for(i = 0; i < N; i++)
  {
    nx = cnt_x[i];

    current_neighbor_cnt = y_num_cnts[i] - 1;
    y_neighbors_ind = y_neighbors(i, _); 
    y_neighbors_ind = y_neighbors_ind[Range(0, current_neighbor_cnt)]; 
    tmp = weight[y_neighbors_ind]; 
    ny = sum(tmp) - weight[i];

    ans += -weight[i] * log(nx) / N;
    ans += -weight[i] * log(ny) / N; 
  }

  return ans;
}

//' @title
//' umi
//' @description
//' This subroutine calculates the uniformed  mutual information where 
//' the distribution for x is replaced by a uniform distribution.  
//' 
//' @param x one random variable from the time-series data
//' @param y another random variable from the time-series data
//' @param k number for nearest neighbors used in entropy calculation
//' @param method Which 2D density estimator you would like to use. 1 is kde estimator and 2 is knn based estimator. 
//' @param k_density The number of k nearest neighbors you would like to use when calculating the density 
//' (only applicable when method == 2 or using knn based density estimation). 
//' @param bw Bindwidth used for the kernel density estimator. Currently it is not used. The bindwidth in the kde function is automatically estimated. 
//' 
//' @details
//' \code{umi} takes two random variable x and y and estimated their mutual using the KSG estimator while x is replaced by a uniform distribution. 
//' It relies on a C++ implmentation of kde estimator (https://github.com/timnugent/kernel-density) and the ANN package to query the kNN with KDTree algorithm. 
//' @return  A estimated uniform mutual information value between two variables (x, y) where the distribution for the x is replaced by a uniform distribution. 
//'  
//' @export
// [[Rcpp::export]]
double umi(SEXP x, SEXP y, SEXP k, SEXP method, SEXP k_density, SEXP bw)
{ 
  NumericMatrix x_cpp(x); 
  NumericMatrix y_cpp(y); 

  int k_cpp = as<int>(k); 

  int method_cpp = as<int>(method); 
  int k_density_cpp = as<int>(k_density);
  double bw_cpp = as<double>(bw); 
  
  double ucmi_res = umi_cpp(x_cpp, y_cpp, k_cpp, method_cpp, k_density_cpp, bw_cpp);
  return ucmi_res;
}

/*** R
library(MASS)
library(Scribe)
lung <- load_lung()
x <- exprs(lung)[, 1]
y <- exprs(lung)[, 2]
x <- log(x + 1)
y <- log(y + 1)
data_xy <- matrix(c(x, y), ncol = 2, byrow = T)
res <- kde2d(x, y)
contour(res, xlab = "previous duration",
        ylab = "duration", levels  =  c(0.05, 0.1, 0.2, 0.4) )
kde_cpp_res <- kde_cpp(data_xy)
x_d = x[-length(y)]
y_d = y[-length(y)]
y_t = y[-1]
ucmi(x_d, y_t, y_d, k = 5, method == 1, k_density = 0, bw = 0)

attach(geyser)
plot(duration, waiting, xlim = c(0.5,6), ylim = c(40,100))
f1 <- kde2d(duration, waiting, n = 50, lims = c(0.5, 6, 40, 100))

data_xy <- matrix(c(duration, duration), ncol = 2, byrow = F)
kde_cpp_res <- kde_cpp(data_xy)

##############################################################################################################################
# test the result vs the python implmentation 
##############################################################################################################################
x <- read.csv('/Users/xqiu/Dropbox (Personal)/Projects/Causal_network/causal_network/Python_code/x_ucmi.txt', header = F)
y <- read.csv('/Users/xqiu/Dropbox (Personal)/Projects/Causal_network/causal_network/Python_code/y_ucmi.txt', header = F)
z <- read.csv('/Users/xqiu/Dropbox (Personal)/Projects/Causal_network/causal_network/Python_code/z_ucmi.txt', header = F)

cmi_res <- cmi(as.matrix(x, ncol = 1), as.matrix(y, ncol = 1), as.matrix(z, ncol = 1), k = 5L, 0L)
ucmi_res <- ucmi(as.matrix(x, ncol = 1), as.matrix(y, ncol = 1), as.matrix(z, ncol = 1), k = 5, method = 1, k_density = 0, bw = 0)
cmi_res$cmi_res
ucmi_res$ucmi_res

cmi_res <-cmi(as.matrix(x^2, ncol = 1), as.matrix(y^2, ncol = 1), as.matrix(z^2, ncol = 1), k = 5L, 0L)
ucmi_res <-ucmi(as.matrix(x^2, ncol = 1), as.matrix(y^2, ncol = 1), as.matrix(z^2, ncol = 1), k = 5, method = 1, k_density = 0, bw = 0)
cmi_res$cmi_res
ucmi_res$ucmi_res

cmi_res <-cmi(as.matrix(x^3, ncol = 1), as.matrix(y^3, ncol = 1), as.matrix(z^3, ncol = 1), k = 5L, 0L)
ucmi_res <-ucmi(as.matrix(x^3, ncol = 1), as.matrix(y^3, ncol = 1), as.matrix(z^3, ncol = 1), k = 5, method = 1, k_density = 0, bw = 0)
cmi_res$cmi_res
ucmi_res$ucmi_res

cmi_res <-cmi(as.matrix(x^4, ncol = 1), as.matrix(y^4, ncol = 1), as.matrix(z^4, ncol = 1), k = 5L, 0L)
ucmi_res <- ucmi(as.matrix(x^4, ncol = 1), as.matrix(y^4, ncol = 1), as.matrix(z^4, ncol = 1), k = 5, method = 1, k_density = 0, bw = 0)
cmi_res$cmi_res
ucmi_res$ucmi_res

x5 <- read.csv('/Users/xqiu/Dropbox (Personal)/Projects/Causal_network/causal_network/Python_code/x10_ucmi.txt', header = F)
y5 <- read.csv('/Users/xqiu/Dropbox (Personal)/Projects/Causal_network/causal_network/Python_code/y10_ucmi.txt', header = F)
z5 <- read.csv('/Users/xqiu/Dropbox (Personal)/Projects/Causal_network/causal_network/Python_code/z10_ucmi.txt', header = F)
weight <- read.csv('/Users/xqiu/Dropbox (Personal)/Projects/Causal_network/causal_network/Python_code/weight.txt', header = F)

data_xy <- matrix(c(x5[, 1], z5[, 1]), ncol = 2, byrow = F)
kde_cpp_res <- kde_cpp(as.matrix(data_xy)) # this is very slow 
kde_res <- kde2d(data_xy)

a <- Sys.time()
cmi_res <-cmi(as.matrix(x5, ncol = 1), as.matrix(y5, ncol = 1), as.matrix(z5, ncol = 1), k = 5L, 0L)
b <- Sys.time()

k = 1
b = 1
pdf = 1
density_sample_type = 1;

data_xz <- matrix(x5, z5, ncol = 2)

a <- Sys.time()
ucmi_res <-ucmi(as.matrix(x5, ncol = 1), as.matrix(y5, ncol = 1), as.matrix(z5, ncol = 1), k = 5, method = 1, k_density = 0, bw = 0)
b <- Sys.time()

ucmi_res <-ucmi(as.matrix(x5, ncol = 1), as.matrix(y5, ncol = 1), as.matrix(z5, ncol = 1), k = 5, method = 1, k_density = 0, bw = 0)


a <- Sys.time()
kde_cpp_res <- kde_cpp(as.matrix(data_xz), k = 1, b = 1, pdf = 1, density_sample_type = 2) # this is very slow 
b <- Sys.time()
b - a 

# compare between python and cpp for the neuron dataset: 
data_xy <- matrix(c(x5[, 1], z5[, 1]), ncol = 2, byrow = F)
a <- Sys.time()
kde_cpp_res <- kde_cpp(as.matrix(data_xz), k = 1, b = 1, pdf = 1, density_sample_type = 2) # this is very slow 
  b <- Sys.time()
  b - a 

umi(as.matrix(x5, ncol = 1), as.matrix(y5, ncol = 1), k = 5, method = 1, k_density = 5, bw = 0.01)
weight_umi
weight_umi <- read.csv('/Users/xqiu/Dropbox (Personal)/Projects/Causal_network/causal_network/Python_code/weight_umi.txt', header = F)
kde_cpp_res <- kde_cpp(as.matrix(x5)) # this is very slow 
qplot(t((1 / kde_cpp_res) / mean(1 / kde_cpp_res)), weight_umi)

#!/usr/bin/Rscript
filename <- "univariate_pdf.png"
data <- read.table("/Users/xqiu/Dropbox (Personal)/Projects/Causal_network/causal_network/Cpp/Real_deal/Scribe/kernel-density/data/multivariate.csv", header=FALSE, sep="," ,comment.char="#")
  png(filename)
  plot(data$V1,data$V2,xlab="x",ylab="density",main="Univariate PDF")
  x <-dev.off()
  
  filename <- "univariate_cdf.png"
data <- read.table("/Users/xqiu/Dropbox (Personal)/Projects/Causal_network/causal_network/Cpp/Real_deal/Scribe/kernel-density/data/multivariate.csv", header=FALSE, sep="," ,comment.char="#")
  png(filename)
  plot(data$V1,data$V2,xlab="x",ylab="density",main="Univariate CDF")
  x <-dev.off()
  
  filename <- "pdf_default_gaussian.png"
png(filename, width = 1000, height = 1000)
  tab = read.csv("matrix_default.csv", header=F ,comment.char="#") # read in data in a dataframe
  x = as.numeric(tab[1,-1]) # assign values of first axis
  y = as.numeric(tab[-1,1]) # assign values of second axis
  z = as.matrix(tab[-1,-1]) # cast the remaining data into a  matrix
  nrz <- nrow(z)
  ncz <- ncol(z)
  jet.colors <- colorRampPalette( c("red", "yellow") )
  nbcol <- 20
color <- jet.colors(nbcol)
  zfacet <- z[-1, -1] + z[-1, -ncz] + z[-nrz, -1] + z[-nrz, -ncz]
facetcol <- cut(zfacet, nbcol)
  persp(x,y,z,col = color[facetcol],zlab="density",main="PDF - default bandwidth, gaussian kernel",phi = 30, theta = -30)
  t <-dev.off()
  
  filename <- "pdf_optimal_secant_gaussian.png"
png(filename, width = 1000, height = 1000)
  tab = read.csv("matrix_optimal.csv", header=F ,comment.char="#") # read in data in a dataframe
  x = as.numeric(tab[1,-1]) # assign values of first axis
  y = as.numeric(tab[-1,1]) # assign values of second axis
  z = as.matrix(tab[-1,-1]) # cast the remaining data into a  matrix
  nrz <- nrow(z)
  ncz <- ncol(z)
  jet.colors <- colorRampPalette( c("red", "yellow") )
  nbcol <- 20
color <- jet.colors(nbcol)
  zfacet <- z[-1, -1] + z[-1, -ncz] + z[-nrz, -1] + z[-nrz, -ncz]
facetcol <- cut(zfacet, nbcol)
  persp(x,y,z,col = color[facetcol],zlab="density",main="PDF - optimal bandwidth (secant), gaussian kernel",phi = 30, theta = -30)
  t <-dev.off()
  
  filename <- "pdf_optimal_bisection_guassian.png"
png(filename, width = 1000, height = 1000)
  tab = read.csv("matrix_optimal_safe.csv", header=F ,comment.char="#") # read in data in a dataframe
  x = as.numeric(tab[1,-1]) # assign values of first axis
  y = as.numeric(tab[-1,1]) # assign values of second axis
  z = as.matrix(tab[-1,-1]) # cast the remaining data into a  matrix
  nrz <- nrow(z)
  ncz <- ncol(z)
  jet.colors <- colorRampPalette( c("red", "yellow") )
  nbcol <- 20
color <- jet.colors(nbcol)
  zfacet <- z[-1, -1] + z[-1, -ncz] + z[-nrz, -1] + z[-nrz, -ncz]
facetcol <- cut(zfacet, nbcol)
  persp(x,y,z,col = color[facetcol],zlab="density",main="PDF - optimal bandwidth (bisection), gaussian kernel",phi = 30, theta = -30)
  t <-dev.off()
  
  filename <- "pdf_default_box.png"
png(filename, width = 1000, height = 1000)
  tab = read.csv("matrix_default_box.csv", header=F ,comment.char="#") # read in data in a dataframe
  x = as.numeric(tab[1,-1]) # assign values of first axis
  y = as.numeric(tab[-1,1]) # assign values of second axis
  z = as.matrix(tab[-1,-1]) # cast the remaining data into a  matrix
  nrz <- nrow(z)
  ncz <- ncol(z)
  jet.colors <- colorRampPalette( c("red", "yellow") )
  nbcol <- 20
color <- jet.colors(nbcol)
  zfacet <- z[-1, -1] + z[-1, -ncz] + z[-nrz, -1] + z[-nrz, -ncz]
facetcol <- cut(zfacet, nbcol)
  persp(x,y,z,col = color[facetcol],zlab="density",main="PDF - default bandwidth, box kernel",phi = 30, theta = -30)
  t <-dev.off()
  
  filename <- "pdf_default_epanechnikov.png"
png(filename, width = 1000, height = 1000)
  tab = read.csv("matrix_default_epa.csv", header=F ,comment.char="#") # read in data in a dataframe
  x = as.numeric(tab[1,-1]) # assign values of first axis
  y = as.numeric(tab[-1,1]) # assign values of second axis
  z = as.matrix(tab[-1,-1]) # cast the remaining data into a  matrix
  nrz <- nrow(z)
  ncz <- ncol(z)
  jet.colors <- colorRampPalette( c("red", "yellow") )
  nbcol <- 20
color <- jet.colors(nbcol)
  zfacet <- z[-1, -1] + z[-1, -ncz] + z[-nrz, -1] + z[-nrz, -ncz]
facetcol <- cut(zfacet, nbcol)
  persp(x,y,z,col = color[facetcol],zlab="density",main="PDF - default bandwidth, epanechnikov kernel",phi = 30, theta = -30)
  t <-dev.off()
*/




