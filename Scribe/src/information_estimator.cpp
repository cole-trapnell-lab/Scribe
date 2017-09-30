#include <RcppArmadillo.h>
#include <ANN/ANN.h>	// ANN declarations
#include <boost/math/special_functions/digamma.hpp>
#include <omp.h>

// flann package
#include <string>

#include "../inst/include/ann_neighbor_radius.h"
#include "../inst/include/radius_search.h"
#include "../inst/include/neighbour.h"
#include "../inst/include/information_estimator.h"

using namespace Rcpp;
using namespace arma;
using namespace boost::math;
// using namespace Rcpp::stat


// [[Rcpp::plugins(openmp)]]

// [[Rcpp::depends(RcppArmadillo)]]

//-------------------------------------------------------------------------------------------------------------------
/*
calculate local density  
*/
//-----
// knn_dis = [tree_x.query(point, k_density + 1, p=np.inf)[0][k_density] for point in x]
// density_estimate = np.array([float(k_density) / N / knn_dis[i] ** dx for i in range(len(knn_dis))])
//   weight = (1 / density_estimate) / np.mean(1 / density_estimate)

// [[Rcpp::export]]
NumericVector knn_density(NumericMatrix x, NumericMatrix y, int k = 5)
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
  
  return(weight);
}

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

// [[Rcpp::export]]
NumericMatrix knn_density_2d(NumericVector x, NumericVector y, IntegerVector nGrid, int k = 5)
{
  double x_min = min(x), x_max = max(x);
  double y_min = min(y), y_max = max(y);
  double x_step = (x_max - x_min) / (nGrid[0] - 1), y_step = (y_max - y_min) / (nGrid[1] - 1);

  NumericMatrix x_y_grid(nGrid[0] * nGrid[1], 2);
  NumericMatrix x_y_density(nGrid[0], nGrid[1]);

  for(int i = 0; i < nGrid[0]; i++) 
  {
    for(int j = 0; j < nGrid[1]; j++)
    {
      x_y_grid(i * nGrid[0] + j, _) = NumericVector::create(x_min + x_step * i, y_min + y_step * j); 
    }
  }

  NumericMatrix xy = cbind(x, y);
  int dxy = 2;
  int N = xy.rows();
  int NQ = nGrid[0] * nGrid[1]; 

  Rcpp::NumericVector data_xy_distances(nGrid[0] * nGrid[1]), density_estimate(data_xy_distances); 

  double error_bound = 0.0; int searchtype = 1; int usebdtree = 0; double sqRad = 0.0; // don't use bd tree (usebdtree = 0)

  get_NN_2Set_cpp(xy, x_y_grid, dxy, N, NQ, k, error_bound, searchtype, usebdtree, sqRad, data_xy_distances); //data_xyz_nn_index,
  
  // NumericVector vd = 0.5 * d * log(M_PI) - log(gamma(0.5 * d + 1)) + log(data_xyz_distances); // 
  // density_estimate = k / (N * exp(vd)); // this is for the L2 norm (the one below is the infinity norm) 
  density_estimate = k / (N * pow(data_xy_distances, dxy));  

  for(int i = 0; i < nGrid[0]; i++) 
  {
    for(int j = 0; j < nGrid[1]; j++)
    {
      x_y_density(i, j) = density_estimate[i * nGrid[0] + j]; 
    }
  }

  // return List::create(Rcpp::Named("cmi_res") = cmi_res, 
  //           Rcpp::Named("information_samples") = information_samples); 
  return(x_y_density);
}

/*
 knn_density(1:10, 5)
 knn_density(matrix(1:100, ncol = 1), 5)
 
 density(1:100, bw = 1, n = 100, from = 1, to = 100)$y
 
 cpp <- knn_density(matrix(1:100, ncol = 1), 5)
 Rversion <- 1 / density(1:100, bw = 0.2, n = 100, from = 1, to = 100)$y / mean(1 / density(1:100, bw = 0.2, n = 100, from = 1, to = 100)$y)
 
 qplot(cpp, Rversion)
 
 */

// [[Rcpp::export]]
double digamma_0(double x) // this modified digamma function avoids the undefined case for 0 
{
  if(x == 0)
  {
   // Rcout << "count is 0" << std::endl;
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
  Rcpp::IntegerVector nn_index; 
  Rcpp::NumericVector distances(dimension);
  nn_index = IntegerVector(dimension);

  // use ANN: 
  double error_bound = 0.0; int searchtype = 1; int usebdtree = 0; double sqRad = 0.0; 

  // // test that there is no na or nan values before running knn-graph search 
  // int is_na_tmp = 0;
  // for(int test_na = 0; test_na < x.size(); test_na ++) 
  // {
  //   if(arma::is_finite(x[test_na]) == false) 
  //   {
  //     is_na_tmp = 1;
  //   }
  // }
  // if(is_na_tmp == 1)
  // {
  //   stop("your data have some NaN or NA values!");
  // }
  
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
//' @return a vector of entropy values
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
double mi_cpp(const NumericMatrix& x, const NumericMatrix& y, int k) //, int normalize
{

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
  
  // // test that there is no na or nan values before running knn-graph search 
  // int is_na_tmp = 0;
  // for(int test_na = 0; test_na < data.size(); test_na ++) 
  // {
  //   if(arma::is_finite(data[test_na]) == false) 
  //   {
  //     is_na_tmp = 1;
  //   }
  // }
  // if(is_na_tmp == 1)
  // {
  //   stop("your data have some NaN or NA values!");
  // }
  
  get_NN_2Set_cpp(data, data, d_data, N, N, k, error_bound, searchtype, usebdtree, sqRad,  data_distances); //data_nn_index,
  
  // // use flann: 
  // int cores = 1, checks = 1; std::string build("kdtree"); 
  // NumericVector data_distances = Neighbour(data, k, build, cores, checks); 
 
	NumericVector information_samples(N, digamma_0(N));

	// use ANN 		
  NumericVector all_kth_dist(N);
  for(int i = 0; i < N; i++)
  {
    all_kth_dist[i] = data_distances[i]; 
  }	

  NumericVector k_xy = get_points_in_radius_cpp(data, data, d_data, N, N, x_y_k, error_bound, usebdtree, data_distances);
  NumericVector cnt_x = get_points_in_radius_cpp(x, x, dx, N, N, x_y_k, error_bound, usebdtree, data_distances);
  NumericVector cnt_y = get_points_in_radius_cpp(y, y, dy, N, N, x_y_k, error_bound, usebdtree, data_distances); 
  
  // // use flann: 
  // NumericVector k_xy = RadiusSearch(data, data_distances, build, cores, checks); 
  // NumericVector cnt_x = RadiusSearch(x, data_distances, build, cores, checks); 
  // NumericVector cnt_y = RadiusSearch(y, data_distances, build, cores, checks); 
 
  int i; 
  int max_thread = omp_get_max_threads();
  // Rcout << "Max number of threads available " << max_thread << std::endl; 
  omp_set_num_threads(max_thread);
  #pragma omp parallel for shared(information_samples, N, cnt_x, cnt_y) private(i) //schedule(dynamic) default(none) //collapse(2) , _
  for(i = 0; i < N; i++)
  {
    information_samples[i] += digamma_0(k_xy[i]) - digamma_0(cnt_x[i]) - digamma_0(cnt_y[i]);
  }

  NumericVector weight = knn_density(x, x, k); 
  
	double mi_res = mean(information_samples * weight);//mean(information_samples);
	return mi_res;
}

//' @title
//' mi
//' @description
//' This function estimates the mutual information of two random variables based on their observed values
//' 
//' @param x one random variable from the time-series data
//' 
//' @param y another random variable from the time-series data
//'
//' @param k number for nearest neighbors used in entropy calculation
//'
//' @details
//' \code{mi} takes two random variables x and y to estimate the mutual information between them 
//' using the KSG estimator
//' It relies on the ANN package to query the kNN with KDTree algorithm. 
//' @return a numeric value for the mutual information estimator between two variables (x, y)
//' @export
// [[Rcpp::export]]
double mi(SEXP x, SEXP y, SEXP k) //&
{ 
  NumericMatrix x_cpp(x); 
  NumericMatrix y_cpp(y); 
  int k_cpp = as<int>(k);

  double mi_res = mi_cpp(x_cpp, y_cpp, k_cpp);

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
List cmi_cpp(const NumericMatrix& x, const NumericMatrix& y, NumericMatrix z, int k){ // const NumericMatrix& z, 
	
  // Rcout << " dimension of x is: " << x.rows() << "," << x.cols() << std::endl;
  // Rcout << " dimension of y is: " << y.rows() << "," << y.cols() << std::endl;
  // Rcout << " dimension of z is: " << z.rows() << "," << z.cols() << std::endl;

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

  int dimension = N;  // k *
  Rcpp::NumericVector data_xyz_distances(dimension); 

  // Rcout << "before get_NN_2Set_cpp" << std::endl;
  // Rcout << " x is: " << x << std::endl;
  // Rcout << " y is: " << y << std::endl;
  // Rcout << " z is: " << z << std::endl;
  
  // // test whether or not there is na or nan values before running knn-graph search 
  // int is_na_tmp = 0, is_duplicated_tmp = 0;
  // bool sum_diff;
  // for(int test_na = 0; test_na < data_xyz.size(); test_na ++) 
  // {
  //   if(arma::is_finite(data_xyz[test_na]) == false) 
  //   {
  //     is_na_tmp = 1;
  //   }
  // }
  

  // // std::ofstream file2("k_xyz.txt");
  // // if (file2.is_open())
  // // {
  // //   // file2 << "Here is the matrix k_xyz:\n" << k_xyz << '\n';
  // //   file2 << "data_xyz :\n" << data_xyz << '\n';
  // //   file2 << "data_xz :\n" << data_xz << '\n';
  // //   file2 << "data_yz :\n" << data_yz << '\n';
  // //   file2 << "z :\n" << z << '\n';
  // // }

  // // the following loop tests whether or not we have duplicated points in the data 
  // // how to make this step quicker? duplicated -> index -> interate the index to count -> add noise for xyz, xy, yz
  // NumericVector testa, testb, cnt_xyz(data_xyz.rows() - 2);  
  // for(int test_dup_i = 0; test_dup_i < data_xyz.rows() - 2; test_dup_i ++)
  // {
  //   cnt_xyz[test_dup_i] = 0;
  //   for(int test_dup_j = test_dup_i + 1; test_dup_j < data_xyz.rows() - 1; test_dup_j ++)
  //   { 
  //     testa = data_xyz(test_dup_i, _) ; 
  //     testb = data_xyz(test_dup_j, _) ; 

  //     // data_xyz
  //     sum_diff = is_true(all(testa == testb)); 
  //     if(sum_diff)
  //     {
  //       data_xyz(test_dup_j, _) = data_xyz(test_dup_j, _) + runif(d_data_xyz, 1e-6, 1e-4); // high dimension has smaller variance 
  //       is_duplicated_tmp = 1; 
  //       cnt_xyz[test_dup_i] += 1; // update the number of duplicated point 
  //     }

  //     // data_xy, data_yz and data_z
  //     testa = data_xz(test_dup_i, _) ; 
  //     testb = data_xz(test_dup_j, _) ; 

  //     sum_diff = is_true(all(testa == testb)); 
  //     if(sum_diff)
  //     {
  //       data_xz(test_dup_j, _) = data_xz(test_dup_j, _) + runif(dxz, 1e-10, 1e-7); // rnorm(dxz, 0, 1e-10); 
  //       is_duplicated_tmp = 1; 
  //     }

  //     testa = data_yz(test_dup_i, _) ; 
  //     testb = data_yz(test_dup_j, _) ; 

  //     sum_diff = is_true(all(testa == testb)); 
  //     if(sum_diff)
  //     {
  //       data_yz(test_dup_j, _) = data_yz(test_dup_j, _) + runif(dyz, 1e-10, 1e-6); // rnorm(dyz, 0, 1e-10); 
  //       is_duplicated_tmp = 1; 
  //     }

  //     testa = z(test_dup_i, _); 
  //     testb = z(test_dup_j, _); 

  //     sum_diff = is_true(all(testa == testb)); 
  //     if(sum_diff)
  //     {
  //       z(test_dup_j, _) = z(test_dup_j, _) + runif(dz, 1e-10, 1e-6); // rnorm(dz, 0, 1e-10); 
  //       is_duplicated_tmp = 1; 
  //     }

  //   }
    
  // }

  // if(is_na_tmp == 1 | is_duplicated_tmp == 1)
  // {
  //   Rcout << "your data have some NaN, NA or duplicated values!" << std::endl;
  //   // stop("your data have some NaN, NA or duplicated values!");
  // }

  // Rcout << "the duplicated counts for cnt_xyz are: " << cnt_xyz << std::endl;
  // Rcout << "duplicated index are: " << cnt_xyz == k << std::endl;
  
  // check whether or not we have duplicated points (duplicated points will lead to error)
  // Rcout << "before get_NN_2Set_cpp" << std::endl;

  // NumericVector get_NN_2Set_cpp(const NumericMatrix& data, int& D, int& ND, int& NQ, int& K, double& EPS,
  //                      int& SEARCHTYPE, int& USEBDTREE, double& SQRAD, NumericVector& distances) // const NumericVector& query,  IntegerVector& nn_index, 

  get_NN_2Set_cpp(data_xyz, data_xyz, d_data_xyz, N, N, k, error_bound, searchtype, usebdtree, sqRad, data_xyz_distances); //data_xyz_nn_index,

  // Rcout << "run get_NN_2Set_cpp to get data_xyz_distances \n" << data_xyz_distances << std::endl;
  // Rcout << "after get_NN_2Set_cpp" << std::endl;
  NumericVector information_samples(N);

  // std::ofstream file("test.txt");
  // if (file.is_open())
  // {
  //   file << "length of data_xyz_distances is:\n" << data_xyz_distances.size() << '\n';
  //   file << "Here is the matrix data_xyz_distances:\n" << data_xyz_distances << '\n';
  // }

  // NumericVector get_points_in_radius_cpp(const NumericMatrix& data, int& D, int& ND, int& NQ, int& K, double& EPS,
  //                        int& USEBDTREE, const NumericVector& SQRAD) 
  // avoid ties of cells 
  // Rcout << "get_points_in_radius_cpp 1 data_xyz " << data_xyz  << d_data_xyz << " N "  << N << " cmbn_k " << cmbn_k << " error_bound " << error_bound << " usebdtree " << usebdtree << " data_xyz_distances " << data_xyz_distances << std::endl;
  NumericVector k_xyz = get_points_in_radius_cpp(data_xyz, data_xyz, d_data_xyz, N, N, cmbn_k, error_bound, usebdtree, data_xyz_distances);
  // Rcout << "get_points_in_radius_cpp 2 data_xz " << data_xz << " dxz " << dxz << " N "  << N << " cmbn_k " << cmbn_k << " error_bound " << error_bound << " usebdtree " << usebdtree << " data_xyz_distances " << data_xyz_distances << std::endl;
  NumericVector cnt_xz = get_points_in_radius_cpp(data_xz, data_xz, dxz, N, N, cmbn_k, error_bound, usebdtree, data_xyz_distances);
  // Rcout << "get_points_in_radius_cpp 3 data_yz " << data_yz << " dyz " << dyz << std::endl;
  NumericVector cnt_yz = get_points_in_radius_cpp(data_yz, data_yz, dyz, N, N, cmbn_k, error_bound, usebdtree, data_xyz_distances); 
  // Rcout << "get_points_in_radius_cpp 4" << std::endl;

  // try {
  //   NumericVector cnt_z = get_points_in_radius_cpp(c_z, c_z, dz, N, N, cmbn_k, error_bound, usebdtree, data_xyz_distances); 
  //   throw 42;x
  // } catch (std::exception& e) {
  //     Rcout << " the c_z values are: " << c_z << '\n';
  //     Rcout << " the data_xyz_distances values are: " << data_xyz_distances << '\n';
  //     cerr << "Exception catched : " << e.what() << std::endl;
  // }
  NumericVector cnt_z = get_points_in_radius_cpp(z, z, dz, N, N, cmbn_k, error_bound, usebdtree, data_xyz_distances); 
  // Rcout << "after running get_points_in_radius_cpp z" << z << " dz " << dz << std::endl;

  // NumericVector cnt_z = cnt_xz; 
  // Rcout << "k_xyz" << k_xyz << " cnt_xz " << cnt_xz << " cnt_yz " << cnt_yz << " cnt_z " << cnt_z << std::endl;
  // Rcout << "cnt_z" << std::endl;
  // k_xyz[k_xyz==0] = 1; cnt_xz[cnt_xz == 0] = 1; cnt_yz[cnt_yz == 0] = 1; cnt_z[cnt_z == 0] = 1; // avoid the case where digamma_0 cannot deal with 0
  // Rcout << " dimension of x is: " << x.rows() << "," << x.cols() << std::endl;
  // Rcout << "\n before paralleling " << std::endl;

  // std::ofstream file3("cnt_xz.txt");
  // if (file3.is_open())
  // {
  //   file3 << "Here is the matrix m:\n" << cnt_xz << '\n';
  // }

  // std::ofstream file4("cnt_yz.txt");
  // if (file4.is_open())
  // {
  //   file4 << "Here is the matrix m:\n" << cnt_yz << '\n';
  // }

  // std::ofstream file5("cnt_z.txt");
  // if (file5.is_open())
  // {
  //   file5 << "Here is the matrix m:\n" << cnt_z << '\n';
  // }

  int i;
  // int max_thread = omp_get_max_threads();
  // omp_set_num_threads(max_thread);
  // #pragma omp parallel for shared(information_samples, N, cnt_xz, cnt_yz, cnt_z) private(i) //schedule(dynamic) default(none) //collapse(2) , _
  for(i = 0; i < N; i++)
  {
    if(k_xyz[i] == 0 | cnt_xz[i] == 0 | cnt_yz[i] == 0 | cnt_z[i] == 0)
      continue; 
    // if(k_xyz[i] < cnt_xyz[i]){
    //   k_xyz[i] = cnt_xyz[i];
    //   // Rcout << "k_xyz[i]" << k_xyz[i] << " cnt_xz[i] " << cnt_xz[i] << " cnt_z[i] " << cnt_z[i] << std::endl;
    // }
    // if(i == 0)
    //   Rcout << "k_xyz[i]" << k_xyz[i] << " cnt_xyz[i] " << cnt_xyz[i] << std::endl;

    // Rcout << "current i is " << i << std::endl;
    information_samples[i] += digamma_0(k_xyz[i]) - digamma_0(cnt_xz[i]) - digamma_0(cnt_yz[i]) + digamma_0(cnt_z[i]); 
  }

  // // std::ofstream file2("k_xyz.txt");
  // if (file2.is_open())
  // {
  //   file2 << "Here is the matrix k_xyz:\n" << k_xyz << '\n';
  //   // file2 << "data_xyz :\n" << data_xyz << '\n';
  //   // file2 << "data_xz :\n" << data_xz << '\n';
  //   // file2 << "data_yz :\n" << data_yz << '\n';
  //   // file2 << "z :\n" << z << '\n';
  // }

  // return the information_samples and study the relationship between the sample and the pseudotime 
  // NumericVector weight = knn_density(data_xz, k);  // try data_xz
  NumericVector weight = knn_density(x, x, k);  // try data_xz
  
  // double cmi_res = mean(information_samples * weight);//mean(information_samples);
  double cmi_res = mean(information_samples);
  
  // test that there is no na or nan values before running knn-graph search 
  int is_na_tmp = 0;
  // if(rcpp::is.finite(weight) == false) 
  // {
  //   Rcout << "weight is " << weight << std::endl;
  //   stop("weight is NaN or NA values!");
  // }
  if(arma::is_finite(cmi_res) == false) 
  {
    Rcout << "weight is " << weight << std::endl;
    Rcout << "cmi_res is " << cmi_res << std::endl;
    stop("cmi_res is NaN or NA values!");
  }

	// double cmi_res = mean(information_samples);
	// return cmi_res;
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
//' @details
//' \code{cmi} takes two random variable x and y and estimated their mutual information conditioned on the third random variable z
//' using the KSG estimator. 
//' It relies on the ANN package to query the kNN with KDTree algorithm. 
//' @return a numeric value for the condition mutual information estimator between two variables (x, y) conditioned on variable z
//' @export
// [[Rcpp::export]]
List cmi(SEXP x, SEXP y, SEXP z, SEXP k) 
{ 
  NumericMatrix x_cpp(x); 
  NumericMatrix y_cpp(y); 
  NumericMatrix z_cpp(z); 

  int k_cpp = as<int>(k); 

  List cmi_res = cmi_cpp(x_cpp, y_cpp, z_cpp, k_cpp);
  return cmi_res;
}

/*** R
message('test vd')
library(InformationEstimator)
vd(10)

# 0.9361576864649548

message('test get_NN_2Set')
library(InformationEstimator)
load('/Users/xqiu/Dropbox (Personal)/Projects/Causal_network/causal_network/csv_data/XYZ.txt')
XYZ <- t(XYZ)
dimensioon <- 3; ND <- 100000; NQ <- 100000; k <- 5; eps <- 0; searchtypeInt <- 1; tree_type <- 1; radius <- 0; 
nn.idx <- 1L; dists <- 0.1;
res <- get_NN_2Set(XYZ[, ], XYZ[, ], 3L, 100000L, 100000L, 5L, as.double(0.0), 1L, 0L, as.double(0.0), dists, T)

data_nmslib <- read.table("/Users/xqiu/Dropbox (Personal)/Projects/Genes_Inference_in_Cell_Differentiation_Process/notebook/data_nmslib.txt", sep = '\t')
res <- get_NN_2Set(XYZ[, ], XYZ[, ], 3L, 100000L, 100000L, 5L, as.double(0.0), 1L, 0L, as.double(0.0), dists, T)
a <- Sys.time(); res <- Neighbour(data_nmslib, data_nmslib, 1000, cores = 1); b <- Sys.time()
# 

message('test entropy')
library(InformationEstimator)
entropy(as.matrix(XYZ[, 1]), 5)

# 1.1298735949169973

message('test mi')
library(InformationEstimator)
a <- Sys.time()
mi(as.matrix(XYZ[, 1]), as.matrix(XYZ[, 2]), 5)
mi(as.matrix(XYZ[, 1]), as.matrix(XYZ[, 3]), 5)
mi(as.matrix(XYZ[, 2]), as.matrix(XYZ[, 3]), 5)
mi(as.matrix(XYZ[, 1]), as.matrix(XYZ[, 1]), 5)
mi(as.matrix(XYZ[, 3]), as.matrix(XYZ[, 3]), 5)
mi(as.matrix(XYZ[, 2]), as.matrix(XYZ[, 2]), 5)
b <- Sys.time()
b - a 

# 1.88421442831
# 1.4081763245
# 1.82672752049
# 10.0068027965
# 10.0068027965
# 10.0068027965

# Time difference of 16.8602 secs

message('test cmi')
library(InformationEstimator)
a <- Sys.time()
cmi(as.matrix(XYZ[, 1]), as.matrix(XYZ[, 2]), as.matrix(XYZ[, 3]), 5);
cmi(as.matrix(XYZ[, 1]), as.matrix(XYZ[, 3]), as.matrix(XYZ[, 2]), 5);
cmi(as.matrix(XYZ[, 2]), as.matrix(XYZ[, 3]), as.matrix(XYZ[, 1]), 5);
b <- Sys.time()
b - a 

# 0.47189332607
# 0.00569163896183
# 0.411085051308

# Time difference of 35.11695 secs

*/

/*** R
# test the relationship between cmi value and the pseudotime 

################################################################################################################################################################################
# real data 
################################################################################################################################################################################
lung <- load_lung()
lung_AT1 <- lung[, pData(lung)$State %in% c(2, 3)]
lung_AT1 <- lung_AT1[, order(pData(lung_AT1)$Pseudotime)]
delay <- 5
x <- exprs(lung_AT1[2, 5:(ncol(lung_AT1) - 1)])
y <- exprs(lung_AT1[45, 6:(ncol(lung_AT1))])
z <- exprs(lung_AT1[45, 5:(ncol(lung_AT1) - 1)])
cmi_res <- cmi(t(as.matrix(x)), t(as.matrix(y)), t(as.matrix(z)), k = 5)
qplot(1:length(cmi_res$information_samples), cmi_res$information_samples)

gene_num <- nrow(lung)
pseudotime_cmi_res_df <- matrix(nrow = gene_num * (gene_num - 1), ncol = length(cmi_res$information_samples))
cnt <- 1

################################################################################################################################################################################
# make a heatmap for all pairwise relationship
################################################################################################################################################################################
for(i in 1:nrow(lung)) {
  message('current i is ', i)
  for(j in 1:nrow(lung)) {
    if(i == j)
      next; 
    x <- exprs(lung_AT1[i, delay:(ncol(lung_AT1) - 1)])
    y <- exprs(lung_AT1[j, (delay + 1):(ncol(lung_AT1))])
    z <- exprs(lung_AT1[j, delay:(ncol(lung_AT1) - 1)])
    cmi_res <- cmi(t(as.matrix(x)), t(as.matrix(y)), t(as.matrix(z)), k = 5)
    pseudotime_cmi_res_df[cnt, ] <- cmi_res$information_samples
    cnt <- cnt + 1
  }
}

# make smooth curves for the temporal variation of the information samples 
# only check for the positive gene-pairs
mean_cmi <- rowMeans(pseudotime_cmi_res_df)
test <- di::smooth_genes(t(pseudotime_cmi_res_df[which(mean_cmi > 0), ]), window_size = 20)
pheatmap::pheatmap(t(test), cluster_cols = F)
  
################################################################################################################################################################################
# simulation data 
################################################################################################################################################################################
load('/Users/xqiu/Dropbox (Personal)/Projects/Causal_network/causal_network/RData/neuron_network')
gene_name_vec <- c('Pax6', 'Mash1', 'Brn2', 'Zic1', 'Tuj1', 'Hes5', 'Scl', 'Olig2', 'Stat3', 'Myt1L', 'Aldh1L', 'Sox8', 'Mature')

cell_simulate <- readMat('/Users/xqiu/Dropbox (Personal)/Projects/DDRTree_fstree/DDRTree_fstree/mat_data/cell_simulate.mat')
all_cell_simulation <- cell_simulate$cell.simulate[, 1:400, ] #time 0-20 are the period two branches appear
example_data <- all_cell_simulation[, , 1]

tmp <- expand.grid(1:ncol(t(example_data)), 1:ncol(t(example_data)), stringsAsFactors = F)
super_graph <- tmp[tmp[, 1] != tmp[, 2], ] - 1 # convert to C++ index
super_graph <- super_graph[, c(2, 1)]

delay <- 1
x <- example_data[1, delay:(ncol(example_data) - 1)]
y <- example_data[2, (delay + 1):(ncol(example_data))]
z <- example_data[2, (delay):(ncol(example_data) - 1)]
cmi_res <- cmi(matrix(x, ncol = 1), matrix(y, ncol = 1), matrix(z, ncol = 1), k = 5)
qplot(1:length(cmi_res$information_samples), cmi_res$information_samples)
cmi_res$information_samples

# simulation data (when n = 1)
data_ori <- read.table('/Users/xqiu/Dropbox (Personal)/Projects/Causal_network/causal_network/Python_code/simulation_expr_mat_n_equal_1.txt', sep = '\t', skip = 1)
x3d <- array(dim = c(13, 101, 20))
    
data_ori_gene_unique <- unique(data_ori$V1)
    
data_ori$V2
for(i in 1:dim(x3d)[1]) {
  for(j in 1:dim(x3d)[3]) {
    x3d[i, , j] <- as.numeric(data_ori[data_ori$V1 == data_ori_gene_unique[i] & data_ori$V2 == paste0("R", j), -c(1:2)]) # make sure data_ori$V2 has R#
  }
}
    
dim(x3d)  
example_data <- x3d[, , 1]
################################################################################################################################################################################
# show the gene expression 
################################################################################################################################################################################

# make a heatmap for all pairwise relationship
# real data 
gene_num <- nrow(example_data)
pseudotime_cmi_res_df <- matrix(nrow = gene_num * (gene_num - 1), ncol = length(cmi_res$information_samples))
cnt <- 1

cmi_res_df <- matrix(0, nrow = 13, ncol = 13)
for(i in 1:nrow(example_data)) {
  message("current i is ", i)
  for(j in 1:nrow(example_data)) {
    if(i == j)
      next; 
    x <- example_data[i, delay:(ncol(example_data) - 1)]
    y <- example_data[j, (delay + 1):(ncol(example_data))]
    z <- example_data[j, (delay):(ncol(example_data) - 1)]
    cmi_res <- cmi(matrix(x, ncol = 1), matrix(y, ncol = 1), matrix(z, ncol = 1), k = 5)
    pseudotime_cmi_res_df[cnt, ] <- cmi_res$information_samples
    cnt <- cnt + 1
    
    cmi_res_df[i, j] <- cmi_res$cmi_res
  }
}

# only check for the positive gene-pairs
mean_cmi <- rowMeans(pseudotime_cmi_res_df)
test <- di::smooth_genes(t(pseudotime_cmi_res_df[which(mean_cmi > 0), ]), window_size = 20)
# set the name for the gene pair below: 
gene_uniq <- unique(c(as.character(neuron_network[, 1]), as.character(neuron_network[, 2])))
all_cmbns <- expand.grid(gene_name_vec, gene_name_vec)
valid_all_cmbns <- all_cmbns[all_cmbns$Var1 != all_cmbns$Var2, ]
valid_all_cmbns_df <- data.frame(pair = paste(tolower(valid_all_cmbns$Var2), tolower(valid_all_cmbns$Var1), sep = '_'), pval = 0)
row.names(valid_all_cmbns_df) <- valid_all_cmbns_df$pair
valid_all_cmbns_df[paste(tolower(neuron_network$V1), tolower(neuron_network$V2), sep = '_'), 2] <- 1

test <- t(test)
row.names(test) <- valid_all_cmbns_df$pair[which(mean_cmi > 0)]

pheatmap::pheatmap(test[, ], cluster_cols = F, cluster_rows = F, annotation_names_row = T)

row.names(pseudotime_cmi_res_df) <- valid_all_cmbns_df$pair
pheatmap::pheatmap(pseudotime_cmi_res_df, cluster_cols = F, cluster_rows = F, annotation_names_row = T)

################################################################################################################################################################################
# show the result for cRDI 
################################################################################################################################################################################

# calculate rdi values
run_vec <- rep(1, ncol(example_data))
a <- Sys.time()
rdi_list <- calculate_rdi_multiple_run_cpp(t(example_data), delay = 1, run_vec - 1, as.matrix(super_graph), method = 1) #calculate_rdi(data_noise, delay, method = 1)
b <- Sys.time()
rdi_time <- b - a

top_k_list <- extract_top_incoming_nodes_delays(rdi_list$max_rdi_value, rdi_list$max_rdi_delays, 1)
con_pseudotime_cmi_res_df <- matrix(nrow = gene_num * (gene_num - 1), ncol = length(cmi_res$information_samples))

# di::di_single_run_conditioned
# function (x, y, z, n = 10) 
# {
#   if (is.numeric(x)) 
#     x <- as.matrix(x)
#     if (is.numeric(y)) 
#       y <- as.matrix(y)
#       if (is.numeric(z)) 
#         z <- as.matrix(z)
#         if (ncol(x) != ncol(y)) 
#           stop("The number of time samples has to be the same for X and Y")
#           if (nrow(x) != nrow(z)) 
#             stop("The number of time samples has to be the same for X and all Zs")
#             tau <- n
#             tot_len <- nrow(x) - tau
#             x_past <- x[(tau):(tau - 1 + tot_len), ]
#           yz_past <- y[(tau):(tau - 1 + tot_len), ]
#           for (i in 1:n) {
#             if (i > 1) {
#               x_past <- cbind(x[(tau - i + 1):(tau - i + tot_len), 
#               ], x_past)
#               yz_past <- cbind(y[(tau - i + 1):(tau - i + tot_len), 
#               ], yz_past)
#             }
#             for (j in 1:ncol(z)) {
#               yz_past <- cbind(z[(tau - i + 1):(tau - i + tot_len), 
#                                j], yz_past)
#             }
#           }
#           return(cmi(x_past, y[(tau + 1):(tau + tot_len), ], yz_past))
# }


cnt <- 1
cmi_res_df <- matrix(0, nrow = 13, ncol = 13)
for(i in 1:nrow(example_data)) {
  message("current i is ", i)
  for(j in 1:nrow(example_data)) {
    if(i == j)
      next; 
    
    z_top_k_ind <- top_k_list$top_incoming_nodes[j, 1] + 1
    if(z_top_k_ind == i)
      z_top_k_ind <- top_k_list$top_incoming_nodes[j, 2] + 1
    x_past <- example_data[i, delay:(ncol(example_data) - 1)]
    yz_past <- example_data[j, (delay):(ncol(example_data) - 1)]
    yz_past <- rbind(example_data[z_top_k_ind, delay:(ncol(example_data) - 1)], yz_past)
    y_past <- example_data[j, (delay + 1):(ncol(example_data))]
    
    cmi_res <- cmi(matrix(x_past, ncol = 1), matrix(y_past, ncol = 1), t(yz_past), k = 5)
    con_pseudotime_cmi_res_df[cnt, ] <- cmi_res$information_samples
    cnt <- cnt + 1
    cmi_res_df[i, j] <- cmi_res$cmi_res
  }
}

con_rdi_res_test <- calculate_multiple_run_conditioned_rdi_wrap(t(example_data), as.matrix(super_graph), as.matrix(rdi_list$max_rdi_value), as.matrix(rdi_list$max_rdi_delays), run_vec - 1, 1)

cmi_res_df - con_rdi_res_test
a <- Sys.time()
rdi_list <- calculate_rdi_multiple_run_cpp(t(example_data), delay = c(1), run_vec - 1, as.matrix(super_graph), method = 1) #* 100 + noise
# rdi_list <- calculate_rdi_multiple_run_cpp(data_rem_dup, delay = c(1), run_vec[!duplicated(data)] - 1, as.matrix(super_graph_remove_dup), method = 1) # + noise_remove_dup
b <- Sys.time()
rdi_time <- b - a

dimnames(rdi_list$max_rdi_value) <- list(uniq_gene, uniq_gene)
con_rdi_res_test <- calculate_multiple_run_conditioned_rdi_wrap(t(example_data), as.matrix(super_graph), as.matrix(rdi_list$max_rdi_value), as.matrix(rdi_list$max_rdi_delays), run_vec - 1, 1)
dimnames(con_rdi_res_test) <- list(uniq_gene, uniq_gene)

row.names(con_pseudotime_cmi_res_df) <- valid_all_cmbns_df$pair
pheatmap::pheatmap(con_pseudotime_cmi_res_df, cluster_cols = F, cluster_rows = F, annotation_names_row = T)


 */




