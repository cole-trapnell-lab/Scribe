## **Scribe**: Towards inferring causal regulations with single cell dynamics-coupled measurements

Single-cell transcriptome sequencing now routinely samples thousands of cells, potentially providing enough data to reconstruct *causal gene regulatory networks* from observational data. Here, we developed **Scribe**, a toolkit for detecting and visualizing causal regulations, and explore the potential for single-cell experiments to power network reconstruction. **Scribe** employs *Restricted Directed Information* to determine causality by estimating the strength of information transferred from a potential regulator to its downstream target by taking advantage of time-delays. We apply **Scribe** and other leading approaches for network reconstruction to several types of single-cell measurements and show that there is a dramatic drop in performance for "pseudotime” ordered single-cell data compared to live imaging data. We demonstrate that performing causal inference requires temporal coupling between measurements. We show that methods such as “*RNA velocity*” restore some degree of coupling through an analysis of chromaffin cell fate commitment. These analyses therefore highlight an important shortcoming in experimental and computational methods for analyzing gene regulation at single-cell resolution and point the way towards overcoming it. 

## Installation

Note that this is still an alpha version of **Scribe**. Stable version of **Scribe** will be released when it is ready. Until then, please use **Scribe** with caution. We welcome any bugs reports of **Scribe** and any comments, suggestions regarding to our manuscript (See below). You can install this alpha version of **Scribe** via the following steps:

```sh
# Scribe depends on RANNinf package, first clone the RANNinf package 
git clone https://github.com/cole-trapnell-lab/RANNinf
# You can then cd to the parent folder of `RANNinf` and build and install it with the following command in the terminal: 
R CMD build RANNinf # this build the RANNinf package and produce a file RANNinf_2.5.0.99.tar.gz
R CMD install RANNinf_2.5.0.99.tar.gz # install the built RANNinf package, similar to the RANN package but the infinity norm is used 
# then clone the github repo with 
git clone git@github.com:cole-trapnell-lab/Scribe.git
# then cd to the directory where the cloned github repo located
R CMD install Scribe_0.1.tar.gz # install the built Scribe package. 
```
Package dependencies issues may incur when you try to install RANNinf or Scribe, you can easily install those packages from CRAN or bioconductor. 
On Mac OS, you may confronted the error (*clang: error: linker command failed with exit code 1 (use -v to see invocation)*) wheny you try to install **Scribe**. That is because Scribe depends on gfortran and your system should use the updated gfortran. See https://thecoatlessprofessor.com/programming/rcpp-rcpparmadillo-and-os-x-mavericks--lgfortran-and--lquadmath-error/. 

## Citation
Xiaojie Qiu, Arman Rahimzamani, Li Wang, Qi Mao, Timothy Durham, Jose L McFaline-Figueroa, Lauren Saunders, Cole Trapnell, Sreeram Kannan (2018): Towards inferring causal gene regulatory networks from single cell expression measurements. BioRxiv

biorxiv link: https://www.biorxiv.org/content/early/2018/09/25/426981

twitter link: https://twitter.com/coletrapnell/status/1044986820520435712 
