#include <RcppArmadillo.h>

#include <kde/kde.hpp>  // ANN declarations
using namespace Rcpp;
using namespace arma;

// [[Rcpp::interfaces(r, cpp)]]
// [[Rcpp::export]]
NumericMatrix kde_cpp(NumericMatrix data, int k = 1, int b = 1, int pdf = 1, int density_sample_type = 1)
{
  NumericMatrix density; 
  int N = data.rows(), d = data.cols();

  // set up the size of the density matrix 
	if(density_sample_type == 1) {
		 density = NumericMatrix(N, 1);
	} else if(density_sample_type == 2) {
		if(d == 1) {
		  density = NumericMatrix(1001, 1);
		} else if(d == 2) {
		  density = NumericMatrix(201, 201);
		}
	}

  NumericVector tmp; 
  
	KDE* kde = new KDE();
	vector<double> data_row;

	kde->set_kernel_type(k);
	kde->set_bandwidth_opt_type(b);

	for(int i = 0; i < N; i ++){
		tmp = data(i, _); 
	  data_row = as<vector <double> > (tmp);
	  // Rcout << "tmp is " << tmp << std::endl; 
	  kde->add_data(data_row);
	}

	int vars_count = kde->get_vars_count(); 

	if(vars_count == 1){
    
		double min_x = kde->get_min(0);
		double max_x = kde->get_max(0);

		if(density_sample_type == 1) {
		  for(int i = 0; i < N; i ++) {
				double x = data(i, 0); 

				if(pdf == 1){
					density(i, 0) = kde->pdf(x);
				}else if (pdf == 2){
					density(i, 0) = kde->cdf(x);
				}
			}
		} else if(density_sample_type == 2){
			double x_increment = (max_x-min_x)/1000.0;
	    
	    int i = 0;
			for(double x = min_x; x <= max_x; x += x_increment){
				if(pdf == 1){
					density(i, 0) = kde->pdf(x);
				}else{
					density(i, 0) = kde->cdf(x);
				}
				i ++; 
			}
		}

	}else if(vars_count == 2){

		double min_x = kde->get_min(0);
		double max_x = kde->get_max(0);
		double min_y = kde->get_min(1);
		double max_y = kde->get_max(1);

		if(density_sample_type == 1) {
		  
		  for(int i = 0; i < N; i++){
		    double x = data(i, 0);
		    double y = data(i, 1); 
		    // Rcout << "x is " << x << " y is " << y << std::endl; 
		    
		    if(pdf == 1){
		      density(i, 0) = kde->pdf(x, y);
		    } else if(pdf == 2){
		      density(i, 0) = kde->cdf(x, y);
		    }
		  }
		} else if(density_sample_type == 2){
		  int nd = 200;
		  double x_increment = (max_x-min_x)/nd;
		  double y_increment = (max_y-min_y)/nd;
		  double x = min_x;
		  double y = min_y;
		  
		  for(int i = 0; i < nd; i++, x += x_increment){
		    y = min_y;
		    
		    for(int j = 0; j < nd; j++, y += y_increment){
		      if(pdf == 1){
		        density(i, j) = kde->pdf(x, y);
		      }else{
		        density(i, j) = kde->cdf(x, y);
		      }
		    }
		  }
		}
	}
	delete kde;

	return(density);
}

//' @title
//' kde
//' @description
//' This subroutine calculates the 2d density for a two dimensional matrix using kernel density estimator. 
//' 
//' @param data A two dimensional matrix, data, where the first row correspondds to the sample and column corresponds to the dimensions.
//' @param k Kernel type. 1 = Gaussian (default); 2 = Box; 3 = Epanechnikov
//' @param b Bandwidth optimisation (Gaussian only). 1 = Default; 2 = AMISE optimal, secant method; 3 = AMISE optimal, bisection method
//' @param pdf Calculate PDF or CDF. 1 = PDF (default); 2 = CDF.
//' @param density_sample_type Where do you want to calculate the kde density. 1 = Original data points (default); 2 = on a pre-defined grid.  
//' 
//' @details
//' \code{kde} takes a 2D matrix and uses kernel density estimator to calculate the density 
//' at location of the original data points (default) or a grid of points. For a one dimension matrix, 
//' the grid is calculated on 1001 evenly spaced points while for a two dimensional matrix, it will be 
//' calculated on 201 evenly spaced points on each dimension (201 x 201 points in total).  
//' @return a NumericMatrix where the element is the density estimate (name: density_estimate), 
//' the second one is the weight calculated based on density_estimate. This function is based on from Tim Nugent
//' (https://github.com/timnugent/kernel-density).  
//' @export
// [[Rcpp::export]]
NumericMatrix kde(SEXP data, SEXP k, SEXP b, SEXP pdf, SEXP density_sample_type) 
{ 
  NumericMatrix data_cpp(data);
   
  int k_cpp = as<int>(k); 
  int b_cpp = as<int>(b); 
  int pdf_cpp = as<int>(pdf); 
  int density_sample_type_cpp = as<int>(density_sample_type); 

  NumericMatrix kde_res = kde_cpp(data_cpp, k_cpp, b_cpp, pdf_cpp, density_sample_type_cpp);
  return kde_res;
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

data_xy <- matrix(c(x5[, 1], z5[, 1]), ncol = 2, byrow = F)
kde_cpp_res <- kde_cpp(as.matrix(data_xy[, ]), density_sample_type = 2) # this is very slow 
kde_res <- kde2d(data_xy[, 1], data_xy[, 2], 200)

test_data <- read.csv("/Users/xqiu/Dropbox (Personal)/Projects/Causal_network/causal_network/Cpp/Real_deal/Scribe/kernel-density/data/multivariate.csv", sep = ',', header = F)
kde_res <- kde2d(test_data[, 1], test_data[, 2], n = 200)
kde_cpp_res <- kde_cpp(as.matrix(test_data), k = 1, b = 1, pdf = 1, density_sample_type = 2) # this is very slow 
kde_cpp_res <- kde_cpp(as.matrix(test_data), k = 1, b = 1, pdf = 1, density_sample_type = 1) # this is very slow 
t((1 / kde_cpp_res) / mean(1 / kde_cpp_res))[1:10]

kde(test_data, eval.points = test_data)
ks_res <- kde(test_data, eval.points = test_data)
1 / ks_res$estimate / (mean(1 / ks_res$estimate))

kde_sci_learn <- read.csv('/Users/xqiu/Dropbox (Personal)/Projects/Causal_network/causal_network/Python_code/kde_sci_learn.csv', header = F)

qplot(ks_res$estimate, kde_sci_learn) + geom_abline()
qplot(kde_cpp_res, kde_sci_learn) + geom_abline()
qplot(1 / ks_res$estimate / (mean(1 / ks_res$estimate)), 1 / kde_sci_learn$V1 / (mean(1 / kde_sci_learn$V1)))

*/

// // [[Rcpp::interfaces(r, cpp)]]
// // [[Rcpp::export]]
// List kde2d_mass(NumericVector x, NumericVector y, NumericVector h, IntegerVector n = 25, NumericVector lims)
// {
// 	int nx = x.length(); 
// 	if(y.length() != nx) 
// 	{
// 		stop("data vectors must be the same length");
// 	}
	
// 	for(int test_na = 0; test_na < nx; test_na ++)
// 	{
// 		if(arma::is_finite(x[test_na]) || arma::is_finite(y[test_na])) 
// 		{
// 			stop("Missing or infinite values in the data are not allowed");
// 			return -1; 
// 		}
// 	}

// 	for(int test_na = 0; test_na < lims.length(); test_na ++)
// 	{
// 		if(arma::is_finite(lims[test_na])) 
// 		{
// 			stop("only finite values are allowed in 'lims'"); 
// 		}
// 	}
	
// 	NumericVector n = rep(n, 2); // length.out 2 
// 	IntegerVector gx = seq_len(lims[0], lims[1], n[0]);
// 	IntegerVector gy = seq_len(lims[2], lims[3], n[1]);

// 	for(int test_na = 0; test_na < h.length(); test_na ++)
// 	{
// 		if(h[test_na] <= 0) 
// 		{
// 			stop("Bandwidth must be strictly positive"); 
// 			return -1; 
// 		}

// 		if(arma:is_finite(h[test_na]))
// 		{
// 			h[0] = bandwidth_nrd(x);
// 			h[1] = bandwidth_nrd(y);		
// 			break; 
// 		}
// 	}

// 	h = h / 4; 

// 	vec gx_arma = as<vec>(gx);
// 	vec gx_arma = as<vec>(gy);
// 	vec ax = (gx_arma - x.t()) / h[0];
// 	vec ax = (gx_arma - x.t()) / h[1];

// 	mat z_a = matrix(dnorm(ax), , nx); 
// 	mat z_b = matrix(dnorm(ay), , nx);

//     mat z = z_a * trans(zb)/(nx * h[0] * h[1]); 
    
//     // use inteperlation to predict the density at each original data point

//     return List::create(Named("x") = gx, Named("y") = gy, Named("z") = z); 

// }

// /*** R  
	
// If I understand what you want to do, it could be achieved by fitting a smoothing model to the grid density estimate and then using that to predict the density at each point you are interested in. For example:

// # Simulate some data and put in data frame DF
// n <- 100
// x <- rnorm(n)
// y <- 3 + 2* x * rexp(n) + rnorm(n)
// # add some outliers
// y[sample(1:n,20)] <- rnorm(20,20,20)
// DF <- data.frame(x,y)

// # Calculate 2d density over a grid
// library(MASS)
// dens <- kde2d(x,y)

// # create a new data frame of that 2d density grid
// # (needs checking that I haven't stuffed up the order here of z?)
// gr <- data.frame(with(dens, expand.grid(x,y)), as.vector(dens$z))
// names(gr) <- c("xgr", "ygr", "zgr")

// # Fit a model
// mod <- loess(zgr~xgr*ygr, data=gr)

// # Apply the model to the original data to estimate density at that point
// DF$pointdens <- predict(mod, newdata=data.frame(xgr=x, ygr=y))

// # Draw plot
// library(ggplot2)
// ggplot(DF, aes(x=x,y=y, color=pointdens)) + geom_point()
// */






