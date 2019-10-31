#include <RcppArmadillo.h>
#include <RcppArmadilloExtensions/sample.h>
#include "useful.h"

// [[Rcpp::depends(RcppArmadillo)]]

using namespace arma;

// [[Rcpp::export]]
Rcpp::List AlloStand(arma::mat data, int stand_attempt_max, int tree_attempt_max, 
                     double overlap_max, double plot_length, double quadrat_length, double voxel_size,
                     bool force, bool verbose) {
  
  int tree_index_max = data.n_rows;
  Rcpp::List stand;
  // Initialises output Rcpp::List and assign observed variables 
  // if the precedent replication has failed
  arma::vec quadrat_x(tree_index_max); quadrat_x = data.col(1); // Quadra x-coordinate
  arma::vec quadrat_y(tree_index_max); quadrat_y = data.col(2); // Quadra y-coordinate
  arma::vec stand_dbh = data.col(0); // Diametre at breast height
  arma::vec stand_height(tree_index_max); stand_height = data.col(3); // height
  arma::vec stand_trunk_height(tree_index_max); stand_trunk_height = data.col(4); // Trunk height
  arma::vec stand_crown_radius(tree_index_max); stand_crown_radius = data.col(5); // Crown radius
  arma::vec stand_x(tree_index_max); // X coordinate of the tree in the plot
  arma::vec stand_y(tree_index_max); // Y coordinate of the tree in the plot
  arma::vec overlap_vec(tree_index_max); // overlap rate for each tree's crown 
  
  // index of quadrat (0,0),(0,1) to coordinates in meter of the bottom left corner
  quadrat_x = quadrat_x * quadrat_length;
  quadrat_y = quadrat_y * quadrat_length;
  
  // Define overlap as a vector
  arma::vec overlap_max_vec(tree_attempt_max); overlap_max_vec.fill(overlap_max);
  
  // Find highest tree (that will determine mockup vertical dimension)
  // Extract maximal tree height, metre
  double tree_height_max = stand_height.max()+voxel_size;
  
  // Mockup dimensions
  arma::rowvec plot_dim(3); plot_dim = {plot_length,plot_length,tree_height_max};
  arma::rowvec plot_dim_voxel = arma::ceil(plot_dim / voxel_size);
  
  // Loop while the stand is not successfully populated, up to the maximal stand attempt
  for (int stand_attempt_index=1; stand_attempt_index <= stand_attempt_max+1; stand_attempt_index++) {
    
    if(verbose){
      Rprintf("  Attempt %i starting... \n", stand_attempt_index);
    }
    
    // Initialises the mockup with zeros
    arma::Cube<double> mockup; mockup.zeros(plot_dim_voxel(0)+1,plot_dim_voxel(1)+1,plot_dim_voxel(2)+1);
    
    // Loop over the trees
    int n_fail = 0;
    for (int tree_index=0; tree_index < tree_index_max; tree_index++) {
      if(verbose){
        Rprintf("Tree %i / %i \n",tree_index+1, tree_index_max);
      }
      
      // Tree not located by default
      bool located = false;
      
      // Loop while tree is not located, up to the maximal tree attempt
      for (int tree_attempt_index=1;  tree_attempt_index <= tree_attempt_max; tree_attempt_index++) { 
        if(verbose){
          Rprintf(".");
        }
        
        // Choose a random x y coordinate within the quadrat, in metre
        double tree_x = arma::as_scalar(Rcpp::RcppArmadillo::sample(arma::regspace(std::max(quadrat_x(tree_index), voxel_size),voxel_size, std::min(quadrat_x(tree_index)+quadrat_length, plot_length-voxel_size)),1, false));
        double tree_y = arma::as_scalar(Rcpp::RcppArmadillo::sample(arma::regspace(std::max(quadrat_y(tree_index), voxel_size),voxel_size, std::min(quadrat_y(tree_index)+quadrat_length, plot_length-voxel_size)),1, false));
        
        
        // Create voxelised crown
        arma::mat crown = CrownToVoxels(stand_crown_radius(tree_index), stand_crown_radius(tree_index), 0.5*(stand_height(tree_index)-stand_trunk_height(tree_index)), voxel_size);
        
        // Position the tree in the mockup
        arma::rowvec tree_centre = {tree_x, tree_y, 0.5*(stand_height(tree_index) + stand_trunk_height(tree_index))};
        arma::rowvec tree_centre_voxel = arma::ceil(tree_centre / voxel_size) + 1;
        arma::mat tree_voxel = crown.each_row() + tree_centre_voxel;
        
        // Crown must be strictly above ground
        if (arma::sum(tree_voxel.col(2)<0)!=0) {
          Rcpp::stop("Error with crown depth for tree %i", tree_index+1);
        }
        
        // Simless plot
        int nvox_tree = tree_voxel.n_rows;
        double sum_occuped_vox = 0.0;
        for (int i=0; i<nvox_tree; i++){
          for (int j=0; j<tree_voxel.n_cols; j++){
            // Negative coordinates are re-injected in opposite edge
            if(tree_voxel(i,j) < 0){
              tree_voxel(i,j) += plot_dim_voxel(j);
            }
            
            // Coordinates exceeding plot dimension are re-injected in opposite edge
            if(tree_voxel(i,j) > plot_dim_voxel(j)){
              tree_voxel(i,j) -= plot_dim_voxel(j);
            }
          }
          sum_occuped_vox += mockup(tree_voxel(i,0),tree_voxel(i,1),tree_voxel(i,2));
        }
        
        // Computes overlap rate of current tree in the mockup
        double overlap = sum_occuped_vox/nvox_tree;
        
        if (overlap <= overlap_max_vec(tree_attempt_index-1)) {
          // Overlap is acceptable, the tree position is validated
          for (int i=0; i< nvox_tree; i++) {
            mockup(tree_voxel(i,0),tree_voxel(i,1),tree_voxel(i,2)) = 1.0;
          }
          stand_x(tree_index) = tree_x;
          stand_y(tree_index) = tree_y;
          overlap_vec(tree_index) = overlap; 
          
          located = true;
          if(verbose){
            Rprintf("[DONE] \n");
          }
          
          // Exit the tree attempt loop
          break;
        }
        continue;
        
      } // end of loop for (tree_attempt_index in 1:tree_attempt_max)
      
      // Tree could not be positioned in the current stand attempt
      if (!located) {
        if(verbose){
          Rprintf("\n [FAILED] too many unsuccessful positions \n");
        }
        
        // Exit the tree loop and start a new stand attempt
        if(stand_attempt_index < stand_attempt_max){
          break;
        }
        
        // If the last stand attempt failed try one more time 
        if( stand_attempt_index==stand_attempt_max){
          if(verbose){
            Rprintf("Simulation [FAILED] %i unsuccessful attempt to locate all trees \n", stand_attempt_max);
          }
          
          // On supplementary attempt maximal overlap for some trees is adapted to locate them if force=true 
          if(force){
            overlap_max_vec = arma::regspace(overlap_max,(1-overlap_max)/(tree_attempt_max-1),1);
          }
          break;
        }
        
        // On supplementary attempt some trees may not be located if force=false
        if( stand_attempt_index == stand_attempt_max+1 & !force ){
          stand_x(tree_index) = NA_REAL;
          stand_y(tree_index) = NA_REAL;
          n_fail += 1;
          continue;
        }
        
      }
      
      // If all trees are located before last attempt return results
      if (located & tree_index == tree_index_max-1 & stand_attempt_index <= stand_attempt_max) {
        stand_attempt_index=stand_attempt_max+2;      
      }
      
    } // end of loop for(tree_index in 0:(tree_index_max-1))
    // Move to the following tree
    
    if(verbose & n_fail!=0){
      Rprintf(" Return results with %i / %i trees not located \n", n_fail, tree_index_max);
    }
    
    if(verbose & force & (stand_attempt_index==stand_attempt_max+1)){
      Rprintf("Last attempt by adapting maximal overlap to locate all trees \n", stand_attempt_max);
    }    
  } // end of loop for (stand_attempt_index in 1:stand_attempt_max)
  
  if(verbose){
    Rprintf("Stand - Simulation terminated \n");
  }
  
  stand = Rcpp::List::create(Rcpp::Named("quadrat.x") = quadrat_x, // Quadra x-coordinate
                             Rcpp::Named("quadrat.y")= quadrat_y, // Quadra y-coordinate
                             Rcpp::Named("dbh") = stand_dbh, // Diametre at breast height
                             Rcpp::Named("height") = stand_height, // height
                             Rcpp::Named("trunk.height") = stand_trunk_height, // Trunk height
                             Rcpp::Named("crown.radius") = stand_crown_radius, // Crown radius
                             Rcpp::Named("overlap") = overlap_vec, // overlap of each tree's crown
                             Rcpp::Named("x")= stand_x , Rcpp::Named("y") = stand_y);// X and Y coordinate of the tree in the plot
  
  // Rcpp::List converted in R list
  return stand;
}

/*** R

# Allometric function's parameters 
## DBH to tree height allometry parameters
H_a <- 35
H_b <- 0.03
H_c <- 0.8
## DBH to crown radius allometry parameters
Cr_x2 <- 3.2489
Cr_y2 <- 1.6941
Cr_inter1 <- 0.3939
Cr_coef1 <- 0.4002
Cr_coef2 <- 0.4102
Cr_RSE <- 0.494
## Tree height to trunk height allometry parameters
Ht_inter <- -0.70087
Ht_coef <- 1.04896
Ht_RSE <- 0.3885
# create dataframe
allom.param <- data.frame(H_a, H_b, H_c, Cr_x2, Cr_y2, Cr_inter1, Cr_coef1, Cr_coef2, Cr_RSE, Ht_inter, Ht_coef, Ht_RSE)

source('/home/beauclair/Documents/Allostand/R/randomTreeInventory.R')
# Data simulation 
data <- randomTreeInventory(ntree=200, nquadrat=5, measurement.rate=1, allom.param = allom.param)
# Sort dbh to position the highest trees first 
data <- data[sort(data$dbh,decreasing=T,index.return=T)$ix,]

# Function to locate trees 
stand <- AlloStand(as.matrix(data), stand_attempt_max=10, tree_attempt_max =350,
                   overlap_max=0.2, plot_length=100.0, quadrat_length=20.0, voxel_size=0.5, force=T, verbose=T)
stand <- as.data.frame(stand)

## Representation of results 
source('/home/beauclair/Documents/Allostand/R/tree3d.R')
plot.stand(stand = stand[!is.na(stand$x),], res=20)
**/
  