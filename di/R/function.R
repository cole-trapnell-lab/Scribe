#' This subroutine calculates the volume of a d-dimensional unit ball for Euclidean norm
#'
#' @param d number of deminsion
#' @return a numeric value for the d-dimensional unit ball for Euclidean norm
#' @export
vd <- function(d){
	0.5 * d * log(pi) - log(gamma(0.5 * d + 1))
}

#' This function estimates the entropy of a continuous random variable
#'
#' @param x data matrix used for calculating the entropy
#' @param k number for nearest neighbors used in entropy calculation
#' @importFrom RANN nn2
#' @return a vector of entropy values
#' @export
entropy <- function(x, k = 5) {
  if(is.numeric(x))
    x <- as.matrix(x)
 	N <- nrow(x) # The number of observed samples
 	d <- ncol(x) # The number of the dimensions of the data

 	knn_dis <- RANN::nn2(x, k = k + 1)$nn.dists[, k + 1]
 	ans <- -digamma(k) + digamma(N)

 	return(ans + d * mean(log(knn_dis)))
}

#' This function estimates the mutual information of two random variables based on their observed values
#'
#' @param x one random variable from the time-series data
#' @param y another random variable from the time-series data
#' @param k number for nearest neighbors used in entropy calculation
#' @importFrom RANN nn2
#' @return a numeric value for the mutual information estimator between two variables (x, y)
#' @export
mi <- function(x, y, k = 5) {
  if(is.numeric(x))
    x <- as.matrix(x)
  if(is.numeric(y))
      y <- as.matrix(y)

  N <- nrow(x)
	if(nrow(y) != N)
		warning('Matrix should have the same length')
	dx <- ncol(x)
	dy <- ncol(y)
	data <- cbind(x, y)

	x_knn_res <- RANN::nn2(x, k = nrow(x))
	y_knn_res <- RANN::nn2(y, k = nrow(y))
	data_knn_res <- RANN::nn2(data, k = nrow(x))

	knn_dis <- data_knn_res$nn.dists[, k + 1]
	information_samples <- rep(digamma(N), N)

	RANN::nn2(data, radius = knn_dis[1])$nn.idx
	#run the nn2 only once for all data in a vector
	for(i in 1:N){
	  # print(c(sum(data_knn_res$nn.dists[i, ] <= knn_dis[i]), sum(x_knn_res$nn.dists[i, ] <= knn_dis[i]), sum(y_knn_res$nn.dists[i, ] <= knn_dis[i])))
	  information_samples[i] <- information_samples[i] + digamma(sum(data_knn_res$nn.dists[i, ] <= knn_dis[i]) - 1) #always k
	  information_samples[i] <- information_samples[i] - digamma(sum(x_knn_res$nn.dists[i, ] <= knn_dis[i]) - 1)
	  information_samples[i] <- information_samples[i] - digamma(sum(y_knn_res$nn.dists[i, ] <= knn_dis[i]) - 1)
	}

   return(mean(information_samples))
}

#' This function estimates the CONDITIONAL mutual information of X and Y given Z
#'
#' @param x one random variable from the time-series data
#' @param y another random variable from the time-series data
#' @param z condition random variable for variables (x, y) from the time-series data
#' @param k number for nearest neighbors used in entropy calculation
#' @importFrom RANN nn2
#' @return a numeric value for the condition mutual information estimator between two variables (x, y) conditioned on variable z
#' @export
cmi <- function(x, y, z, k=5) {
  if(is.numeric(x))
    x <- as.matrix(x)
  if(is.numeric(y))
    y <- as.matrix(y)
  if(is.numeric(z))
    z <- as.matrix(z)

  N <- nrow(x)
	if(nrow(y) != N | nrow(z) != N)
		warning('Matrix should have the same length')

	dx <- ncol(x)
	dy <- ncol(y)
	dz <- ncol(z)

	data_xyz <- Reduce(cbind, list(x, y, z))
	data_xz <- cbind(x, z)
	data_yz <- cbind(y, z)

	xz_knn_res <- RANN::nn2(data_xz, k = nrow(data_xz))
	yz_knn_res <- RANN::nn2(data_yz, k = nrow(data_yz))
	z_res <- RANN::nn2(z, k = nrow(z))
	data_knn_res <- RANN::nn2(data_xyz, k = nrow(data_yz))

	knn_dis <- data_knn_res$nn.dists[, k + 1]
	information_samples <- rep(0, N)

	#run the nn2 only once for all data in a vector
	for(i in 1:N){
	  # print(c(sum(data_knn_res$nn.dists[i, ] <= knn_dis[i]), sum(xz_knn_res$nn.dists[i, ] <= knn_dis[i]), sum(yz_knn_res$nn.dists[i, ] <= knn_dis[i]), sum(z_res$nn.dists[i, ] <= knn_dis[i])))
	  # print(c(digamma(sum(data_knn_res$nn.dists[i, ] <= knn_dis[i]) - 1), digamma(sum(xz_knn_res$nn.dists[i, ] <= knn_dis[i]) - 1), digamma(sum(yz_knn_res$nn.dists[i, ] <= knn_dis[i]) - 1), digamma(sum(z_res$nn.dists[i, ] <= knn_dis[i]) - 1)))

	  information_samples[i] <- information_samples[i] + digamma(sum(data_knn_res$nn.dists[i, ] <= knn_dis[i]) - 1)
	  information_samples[i] <- information_samples[i] - digamma(sum(xz_knn_res$nn.dists[i, ] <= knn_dis[i]) - 1)
	  information_samples[i] <- information_samples[i] - digamma(sum(yz_knn_res$nn.dists[i, ] <= knn_dis[i]) - 1)
	  information_samples[i] <- information_samples[i] + digamma(sum(z_res$nn.dists[i, ] <= knn_dis[i]) - 1)
	}

	# print(information_samples)
  return(mean(information_samples))
}

#' This function simulates the DIRECTED mutual information from X to Y when you have a SINGLE run of the process
#'
#' @param x one random variable from the time-series data
#' @param y another random variable from the time-series data
#' @param d delay in the formula I(x-->y)=I(x_{t-d};y_t|y_{t-1}) default: d=1
#' @return a matrix for the condition mutual information estimators between all pairwise variables (x, y) in the data matrix x, y
#' @export
di_single_run <- function(x, y, n=10) {
  if(is.numeric(x))
    x <- as.matrix(x)
  if(is.numeric(y))
    y <- as.matrix(y)

  if(ncol(x) != ncol(y))
 		stop('The number of time samples has to be the same for X and Y')
    tau <- n
    tot_len <- nrow(x) - tau
    x_past <- x[(tau):(tau - 1 + tot_len), ]
    y_past <- y[(tau):(tau - 1 + tot_len), ]
    for(i in 2:n) {
        x_past = cbind(x[(tau - i + 1):(tau - i + tot_len), ], x_past)
        y_past = cbind(y[(tau - i + 1):(tau - i + tot_len), ], y_past)
    }

    return(cmi(x_past, as.matrix(y[(tau + 1):(tau+tot_len), ]), y_past))
}

#' This function simulates the CONDITIONED DIRECTED mutual information from X to Y when you have a SINGLE run of the processes
#'
#' @param x one random variable from the time-series data
#' @param y another random variable from the time-series data
#' @param z z is a dataframe (or matrix) containing the data of other processes upon the past of which the mi is conditioned
#' @param n Parameter n determines the the number of previous time samples upon which the mi is conditioned
#' @return a matrix for the condition mutual information estimators between all pairwise variables (x, y) in the data matrix x, y
#' @export
di_single_run_conditioned <- function(x, y, z, n = 10) {
  if(is.numeric(x))
    x <- as.matrix(x)
  if(is.numeric(y))
    y <- as.matrix(y)
  if(is.numeric(z))
    z <- as.matrix(z)

  if(ncol(x) != ncol(y))
 		stop('The number of time samples has to be the same for X and Y')
 	if(nrow(x) != nrow(z))
 		stop('The number of time samples has to be the same for X and all Zs')

    tau <- n
    tot_len <- nrow(x) - tau
    x_past <- x[(tau):(tau - 1 + tot_len), ]
    yz_past <- y[(tau):(tau - 1 + tot_len), ]
    for(i in 1:n){
		if(i > 1) {
	        x_past <- cbind(x[(tau - i + 1):(tau - i + tot_len), ], x_past)
	        yz_past <- cbind(y[(tau - i + 1):(tau - i + tot_len), ], yz_past)
        }
      	for(j in 1:ncol(z)){
      		yz_past <- cbind(z[(tau - i + 1):(tau - i + tot_len), j], yz_past)
      	}
    }

    return(cmi(x_past, y[(tau + 1):(tau+tot_len), ], yz_past))
}

#' This function simulates the RESTRICTED DIRECTED mutual information from X to Y when you have a SINGLE run of the processes
#'
#' @param x one random variable from the time-series data
#' @param y another random variable from the time-series data
#' @param d delay in the formula I(x-->y)=I(x_{t-d};y_t|y_{t-1}) default: d=1
#' @return a matrix for the condition mutual information estimators between all pairwise variables (x, y) in the data matrix x, y
#' @export
rdi_single_run <- function(x, y, d=1) {
  if(is.numeric(x))
    x <- as.matrix(x)
  if(is.numeric(y))
    y <- as.matrix(y)

  if(nrow(x) != nrow(y))
 		stop('The number of time samples has to be the same for X and Y')

    return(cmi(x[1:(nrow(x) - d), ], y[-(1:d), ], y[d:(nrow(y) - 1), ]))
}

#' This function simulates the CONDITIONED DIRECTED mutual information from X to Y CONDITIONED ON Z when you have a SINGLE run of the processes
#'
#' @param x one random variable from the time-series data
#' @param y another random variable from the time-series data
#' @param z z is a dataframe or matrix consisting of the data for different variables
#' @param z_delay z_delay is also a dataframe or matrix consisting of the delays to be applied to different variables
#' @param d delay in the formula I(x-->y)=I(x_{t-d};y_t|y_{t-1}) default: d=1
#' @return a matrix for the condition mutual information estimators between all pairwise variables (x, y) in the data matrix x, y
#' @export
rdi_single_run_conditioned <- function(x, y, z, z_delays, d = 1) {
  if(is.numeric(x))
    x <- as.matrix(x)
  if(is.numeric(y))
    y <- as.matrix(y)
  if(is.numeric(z))
    z <- as.matrix(z)

  if(nrow(x) != nrow(y))
 		stop('The number of time samples has to be the same for X and Y')
 	if(nrow(x) != nrow(z))
 		stop('The number of time samples has to be the same for X and all Zs')

    tau <- max(c(z_delays, d))
    tot_len <- nrow(x) - tau
    yz <- y[(tau):(tau - 1 + tot_len), ]
    for(i in 1:ncol(z)){
        # yz = np.concatenate( (z[key][tau-z_delays[key]:tau-z_delays[key]+tot_len], yz), axis=1)
        yz <- cbind(z[(tau - z_delays[i] + 1):(tau-z_delays[i]+tot_len), i], yz)
    }

    return(cmi(x[(tau - d + 1):(tau - d + tot_len), ], y[(tau + 1):(tau + tot_len), ], yz))
}

#' other functions need to implement:
#' sc (Shannon capacity);
#' csc (Conditioned Shannon capacity);
#' causal_sc (Causal Shannon capacity);
#' causal_sc_conditioned
#' d_regularizer
#' get_obj
#' get_grad
#' projection


#' #' This function simulates the DIRECTED mutual information from X to Y when you have MANY runs of the processes
#' #'
#' #' @param x one random variable from the time-series data
#' #' @param y another random variable from the time-series data
#' #' @return a matrix for the condition mutual information estimators between all pairwise variables (x, y) in the data matrix x, y
#' #' @export
#' # rdi_many_runs <- function(x, y) {
#' #  	if(ncol(x) != ncol(y))
#' #  		stop('The number of time samples has to be the same for X and Y')
#'
#' #     T <- ncol(x)
#' #     ans <- 0
#'
#' #     for(t in 1:(T - 1)) {
#' #         ans <- ans + cmi(x[, (t):t], y[, (t + 1):(t + 1)], y[, (t):t])
#' #         print(t, '\n')
#' #     }
#'
#' #     return(ans)
#' # }

