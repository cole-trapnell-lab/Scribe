#include <Rcpp.h>
#include <ANN/ANN.h>  // ANN declarations

#include "../inst/include/ann_neighbor_radius.h"
using namespace Rcpp;

/**
Things to update: 
1. avoid the conversion of NumericMatrix to NumericVector 
2. replace ANN with FLANN 
*/

//------------------------------------------------------------------------------------------------
//				 Near Neighbours Program
//------------------------------------------------------------------------------------------------
// 
// //[[Rcpp::export]]
// https://stackoverflow.com/questions/30129968/how-to-initialize-numericvector-to-a-specific-size-after-it-has-been-declared
// https://stackoverflow.com/questions/13782943/how-to-resize-a-numericvector

//-------------------------------------------------------------------------------------------------------------------
/*
  get the kth neighbor's distance 
*/
//-------------------------------------------------------------------------------------------------------------------

NumericVector get_NN_2Set_cpp(const NumericMatrix& data, const NumericMatrix& query_data, int& D, int& ND, int& NQ, int& K, double& EPS,
                       int& SEARCHTYPE, int& USEBDTREE, double& SQRAD, NumericVector& distances) // const NumericVector& query,  IntegerVector& nn_index, 
{
  const int d = D;		// Number of Dimensions for points
  const int nd = ND;	// Number of Data points
  const int nq= NQ;		// Number of Query points
  const int k = K;		// Maximum number of Nearest Neighbours
  
  const int searchtype = SEARCHTYPE;
  const bool usebdtree = USEBDTREE?true:false;
  
  const double error_bound = EPS;	// enough said!
  const double sqRad = SQRAD;		// Squared Radius for rad search
  
  ANNkd_tree	*the_tree;	// Search structure
  
  ANNpointArray data_pts 	= annAllocPts(nd,d);		// Allocate data points
  // ANNidxArray nn_idx 		= new ANNidx[k];		// Allocate near neigh indices
  // ANNdistArray dists 		= new ANNdist[k];		// Allocate near neighbor dists
  ANNidxArray nn_idx;    // Allocate near neigh indices
  ANNdistArray dists;   // Allocate near neighbor dists

  int *d_ptr = new int[d]; // pointer for updating ind
  int ptr = 0;

  NumericVector ANNkdFRPtsInRange(nq);
  
  // set up column offsets for query point matrix (to convert Row/Col major)
  for(int i = 0; i < d; i++)
  {
    d_ptr[i] = i*nd;
  }
  
  for(int i = 0; i < nd; i++) // now construct the points (to convert Row / Col major)
  {
    for(int j = 0; j < d; j++)
    {    
      data_pts[i][j]=data[ d_ptr[j]++ ];
    }
  }
  
  if(usebdtree){
      the_tree = new ANNbd_tree(	// Build search structure
      data_pts,			// The data points
      nd,					// Number of data points
      d);					// Dimension of space
  } else {
    the_tree = new ANNkd_tree( data_pts, nd, d);
  }
  
  // set up offsets for query point matrix (to convert Row / Col major)
  for(int i = 0; i < d; i++)
  {
    d_ptr[i] = i*nq;
  }
  
  ANNpoint pq;
  pq = annAllocPt(d);
  int i, j;

  // omp_set_num_threads(8); // thread-safe not satisfied 
  // #pragma omp parallel for shared(nq, d, data, d_ptr, searchtype, the_tree, k, error_bound, distances, ANNkdFRPtsInRange, sqRad) private(i, j, nn_idx, dists, pq) //schedule(dynamic) default(none) //collapse(2) , _
  for(i = 0; i < nq; i++)	// Run all query points against tree
  {
    pq = annAllocPt(d);
    nn_idx    = new ANNidx[k];    // Allocate near neigh indices
    dists    = new ANNdist[k];   // Allocate near neighbor dists
    // read coords of current query point
    for(j = 0; j < d; j++)
    {
      pq[j]=query_data[ d_ptr[j]++ ]; // update to query point 
    }

   if(searchtype == 1) {
      the_tree->annkSearch( // search
          pq, // query point
          k,    // number of near neighbors
          nn_idx,   // nearest neighbors (returned)
          dists,    // distance (returned)
          error_bound); // error bound
    } else if(searchtype == 2) { // Priority search
      the_tree->annkPriSearch(pq, k, nn_idx, dists, error_bound);
    } else if(searchtype == 3) {  // Fixed radius search
      ANNkdFRPtsInRange[i] =  the_tree->annkFRSearch( pq, sqRad, k, nn_idx, dists,error_bound);
    }

    distances[ptr++] = ANN_ROOT(dists[k - 1]); // only return the kth distance
  }

  // Do a little bit of memory management......
  annDeallocPt(pq);
  annDeallocPts(data_pts);
  delete [] nn_idx;
  delete [] dists;
  delete [] d_ptr;
  delete the_tree;
  annClose();  

  return ANNkdFRPtsInRange;
}

// //[[Rcpp::export]]
Rcpp::List get_NN_2Set(SEXP R_data, SEXP R_query_data, SEXP R_d, SEXP R_nd, SEXP R_nq, SEXP R_k, SEXP R_error_bound,
                       SEXP R_searchtype, SEXP R_usebdtree, SEXP R_sqRad, SEXP R_distances,
                       SEXP R_verbose) //SEXP R_nn_index, SEXP R_query, 
{ 
  bool verbose = as<bool>(R_verbose);
  
  if (verbose)
    Rcpp::Rcout << "Mapping R_data, R_query_data" << std::endl;
  Rcpp::NumericMatrix data(R_data); 
  Rcpp::NumericMatrix query_data(R_query_data); 
  
  if (verbose)
    Rcpp::Rcout << "Mapping R_d, R_nd, R_k, R_error_bound, R_searchtype, R_usebdtree, R_sqRad" << std::endl;
  int d = as<int>(R_d);
  int nd = as<int>(R_nd);
  int nq = as<int>(R_nq);
  int k = as<int>(R_k);
  double error_bound = as<double>(R_error_bound);
  int searchtype = as<int>(R_searchtype);
  int usebdtree = as<int>(R_usebdtree);
  double sqRad = as<double>(R_sqRad);

  if (verbose)
    Rcpp::Rcout << "Mapping R_nn_index, R_distance" << std::endl;
  
  // Rcpp::IntegerVector nn_index(R_nn_index); 
  Rcpp::NumericVector distances(R_distances); 
  
  int dimension = nq; //only store the kth distance (not k * nq)   
  distances = NumericVector(dimension);

  NumericVector ANNkdFRPtsInRange = get_NN_2Set_cpp(data, query_data, d, nd, nq, k, error_bound, searchtype, usebdtree, sqRad,  distances); //nn_index, query,
  
  return Rcpp::List::create(Rcpp::Named("distances") = Rcpp::wrap(distances), 
                          Rcpp::Named("num_points") = Rcpp::wrap(ANNkdFRPtsInRange));    
}
  
//-------------------------------------------------------------------------------------------------------------------
/*
  get number of cells within a radius 
*/
//-------------------------------------------------------------------------------------------------------------------

NumericVector get_points_in_radius_cpp(const NumericMatrix& data, const NumericMatrix& query_data, int& D, int& ND, int& NQ, int& K, double& EPS,
                       int& USEBDTREE, const NumericVector& SQRAD) //& const NumericVector& query, 
{
  const int d = D;    // Number of Dimensions for points
  const int nd = ND;  // Number of Data points
  const int nq= NQ;   // Number of Query points
  const int k = K;    // Maximum number of Nearest Neighbours
  
  const bool usebdtree = USEBDTREE?true:false;
  
  const double error_bound = EPS; // enough said!
  
  ANNkd_tree  *the_tree;  // Search structure
  
  ANNpointArray data_pts  = annAllocPts(nd,d);    // Allocate data points
  ANNidxArray nn_idx    = new ANNidx[k];    // Allocate near neigh indices
  ANNdistArray dists    = new ANNdist[k];   // Allocate near neighbor dists
  
  int *d_ptr = new int[d]; // pointer for updating ind

  NumericVector ANNkdFRPtsInRange(nq);
  
  // set up column offsets for query point matrix (to convert Row/Col major)
  for(int i = 0; i < d; i++)
  {
    d_ptr[i] = i*nd;
  }
  
  for(int i = 0; i < nd; i++) // now construct the points (to convert Row / Col major)
  {
    for(int j = 0; j < d; j++)
    {    
      data_pts[i][j]=data[ d_ptr[j]++ ];
    }
  }
  
  if(usebdtree){
      the_tree = new ANNbd_tree(  // Build search structure
      data_pts,     // The data points
      nd,         // Number of data points
      d);         // Dimension of space
  } else {
    the_tree = new ANNkd_tree( data_pts, nd, d);
  }
  
  // set up offsets for query point matrix (to convert Row / Col major)
  for(int i = 0; i < d; i++)
  {
    d_ptr[i] = i*nq;
  }
  
  ANNpoint pq = annAllocPt(d);
  int i, j;
  // omp_set_num_threads(8);
  // #pragma omp parallel for shared(nq, d, data, d_ptr, the_tree, SQRAD, k, error_bound, ANNkdFRPtsInRange) private(i, j, nn_idx, dists, pq) //schedule(dynamic) default(none) //collapse(2) , _
  for(i = 0; i < nq; i++) // Run all query points against tree
  {
    // read coords of current query point
    for(j = 0; j < d; j++)
    {
      pq[j]=query_data[ d_ptr[j]++ ];
    }

    ANNkdFRPtsInRange[i] =  the_tree->annkFRSearch( pq, SQRAD[i], k, nn_idx, dists,error_bound);
  }
  
  // Do a little bit of memory management......
  annDeallocPt(pq);
  annDeallocPts(data_pts);
  delete [] nn_idx;
  delete [] dists;
  delete [] d_ptr;
  delete the_tree;
  annClose();  

  return ANNkdFRPtsInRange;
}
