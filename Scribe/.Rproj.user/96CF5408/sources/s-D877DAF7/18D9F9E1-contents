// this function should be updated to implement the network sparsifier algorithm 

// #include <omp.h>
#include <Rcpp.h>
using namespace Rcpp;

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//

// [[Rcpp::export]]
NumericMatrix clr(NumericMatrix x) {
  int n_gene = x.rows();
  int n_col = x.cols();
  if(n_gene != n_col)
  {
    throw std::logic_error( "the input matrix shohuld be a square matrix" );
  }
  // assume each column as a target 
  int i, j;
  double m_tmp, v_tmp; 
  NumericVector m(n_gene), v(n_gene);
  for(i = 0; i < n_gene; i++)
  {
    m_tmp = 0; v_tmp = 0;  
    for(j = 0; j < n_gene; j ++)
    {
      m_tmp += x(i, j);
    }
    m[i] = m_tmp / n_gene;
    for(j = 0; j < n_gene; j ++)
    {
      v_tmp += pow(x(i, j) - m[i], 2);
    }
    v[i] = v_tmp / (n_gene - 1); // denominator as n - 1
  }

  for(i = 0; i < n_gene; i++)
  {
    for(j = 0; j < n_gene; j ++)
    {
      double zij = x(i, j) - m[i];
      zij = (zij < 0 || v[i] == 0) ? 0 : sqrt(zij*zij/v[i]);
      x(i, j) = zij;
    }
  }

  return x; 
}

// maybe also implement Qi's sparsifer method here 

// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//

/*** R
library(Scribe)
library(monocle)
lung <- load_lung()
clr(exprs(lung)[1:5, 1:5]) 
apply(exprs(lung)[1:5, 1:5], 1, function(x) (x - mean(x)) / (sd(x))) # note that apply flip the x, y order
clr(exprs(lung)) 
*/
