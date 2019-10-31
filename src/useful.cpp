#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

// Voxelize the crown centered around the point (0, 0, 0), given a voxel size (dvox).
arma::mat CrownToVoxels(double crown_rx, double crown_ry, double crown_rz, double dvox) {
  
  // Crown radius in voxels
  int rx = std::ceil(crown_rx / dvox);
  int ry = std::ceil(crown_ry / dvox);
  int rz = std::ceil(crown_rz / dvox);
  
  // Create an eighth of the ellipsoid along (1, 1, 1)
  arma::mat crown_q111;
  // Count of rows 
  int i = 0;
  for (int ix = 0; ix < rx; ix++) {
    for (int iy = 0; iy < ry; iy++) {
      for (int iz = 0; iz < rz; iz++) {
        if ( ix*ix/(rx*rx) + iy*iy/(ry*ry) + iz*iz/(rz*rz) <= 1) {
          arma::irowvec voxels_int = {ix, iy, iz};
          arma::rowvec voxels = arma::conv_to<arma::rowvec>::from(voxels_int);
          crown_q111.insert_rows(i, voxels);
          i+=1;
        }
      }
    }
  }
  
  // Project the other 7 eighth along axes 
  arma::rowvec q000 = {-1,-1,-1};
  arma::rowvec q001 = {-1,-1,1};
  arma::rowvec q010 = {-1,1,-1};
  arma::rowvec q011 = {-1,1,1};
  arma::rowvec q100 = {1,-1,-1};
  arma::rowvec q101 = {1,-1,1};
  arma::rowvec q110 ={1,1,-1};
  
  arma::mat crown_q000 = crown_q111.each_row()%q000;
  arma::mat crown_q001 = crown_q111.each_row()%q001;
  arma::mat crown_q010 = crown_q111.each_row()%q010;
  arma::mat crown_q011 = crown_q111.each_row()%q011; 
  arma::mat crown_q100 = crown_q111.each_row()%q100;
  arma::mat crown_q101 = crown_q111.each_row()%q101;
  arma::mat crown_q110 = crown_q111.each_row()%q110;
  
  arma::mat crown_mat = arma::join_cols(crown_q111,arma::join_cols(crown_q000, crown_q001));
  crown_mat = arma::join_cols(crown_mat, arma::join_cols(crown_q010, crown_q011));
  crown_mat = arma::join_cols(crown_mat, arma::join_cols(crown_q100, crown_q101));
  crown_mat = arma::join_cols(crown_mat, crown_q110);
  
  return crown_mat;
}
