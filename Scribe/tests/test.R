library(InformationEstimator)
library(destiny)
library(monocle)

load('/Users/xqiu/Dropbox (Personal)/Projects/Causal_network/causal_network/csv_data/XYZ.txt')
XYZ <- t(XYZ)
dimensioon <- 3; ND <- 100000; NQ <- 100000; k <- 5; eps <- 0; searchtypeInt <- 1; tree_type <- 1; radius <- 0; 

library(InformationEstimator)
N <- 10000
a <- Sys.time()
calculate_rdi_cpp(XYZ[1:N, ], 1.0, 3L)
b <- Sys.time()
# 
# # entropy_cpp(XYZ[, 1], k = 5, N = 100000)
# # vd_cpp(10)
# nn.idx = 1L; dists = 1.0;
# 
# # res <- get_NN_2Set(as.double(XYZ[, 1]), as.double(XYZ[, 1]), 1L, 100000L, 100000L, 5L, as.double(0.0), 1L, 0L, as.double(0.0), nn.idx, dists, T)
# 
# 
# res <- get_NN_2Set(XYZ[, ], XYZ[, ], 3L, 100000L, 100000L, 5L, as.double(0.0), 1L, 0L, as.double(0.0), nn.idx, dists, T)

a <- Sys.time(); res <- get_NN_2Set(XYZ[, ], XYZ[, ], 3L, 100000L, 100000L, 1000L, as.double(0.0), 1L, 0L, as.double(0.0), nn.idx, dists); b <- Sys.time()
a <- Sys.time(); res <- Neighbour(XYZ, XYZ, 1000, cores = 8); b <- Sys.time()

nr <- 1000
nc <- 10
p <- LshParameterSetter$new(nr, nc)
X <- matrix(rnorm(nr * nc), nr, nc)
tab <- LshNnTable$new(t(X), p)
tab$find_nearest_neighbor(as.vector(X[1, ]))
