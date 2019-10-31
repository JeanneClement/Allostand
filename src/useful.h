#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]
using namespace arma;

// Prototype of useful functions

// Voxelize the crown centered around the point (0, 0, 0), given a voxel size (dvox).
arma::mat CrownToVoxels(double crown_rx, double crown_ry, double crown_rz, double dvox);

// End