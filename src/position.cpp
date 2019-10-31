#include <RcppArmadillo.h>
#include <RcppArmadilloExtensions/sample.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppGSL)]]

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

// [[Rcpp::export]]
Rcpp::List AlloStand(arma::mat data, int stand_attempt_max, int tree_attempt_max, 
                     double overlap_max, double plot_length, double quadrat_length, double voxel_size) {
  
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
  
  // index of quadrat (0,0),(0,1) to coordinates in meter of the bottom left corner
  quadrat_x = quadrat_x * quadrat_length;
  quadrat_y = quadrat_y * quadrat_length;
  
  // Find highest tree (that will determine mockup vertical dimension)
  // Extract maximal tree height, metre
  double tree_height_max = stand_height.max()+voxel_size;
  
  // Mockup dimensions
  arma::rowvec plot_dim(3); plot_dim = {plot_length,plot_length,tree_height_max};
  arma::rowvec plot_dim_voxel = arma::ceil(plot_dim / voxel_size);
  
  // Loop while the stand is not successfully populated, up to the maximal stand attempt
  for (int stand_attempt_index=1; stand_attempt_index <= stand_attempt_max; stand_attempt_index++) {
    
    Rprintf("  Attempt %i starting... \n", stand_attempt_index);
    
    // Initialises the mockup with zeros
    arma::Cube<double> mockup; mockup.zeros(plot_dim_voxel(0)+1,plot_dim_voxel(1)+1,plot_dim_voxel(2)+1);
    
    // Loop over the trees
    for (int tree_index=0; tree_index < tree_index_max; tree_index++) {
      Rprintf("Tree %i / %i \n",tree_index+1, tree_index_max);
      
      // Tree not located by default
      bool located = false;
      
      // Loop while tree is not located, up to the maximal tree attempt
      for (int tree_attempt_index=1;  tree_attempt_index <= tree_attempt_max; tree_attempt_index++) { 
      Rprintf(".");
        
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
        if (overlap <= overlap_max) {
          // Overlap is acceptable, the tree position is validated
          for (int i=0; i< nvox_tree; i++) {
            mockup(tree_voxel(i,0),tree_voxel(i,1),tree_voxel(i,2)) = 1.0;
          }
          stand_x(tree_index) = tree_x;
          stand_y(tree_index) = tree_y;
          
          located = true;
          Rprintf("[DONE] \n");
          
          // Exit the tree attempt loop
          break;
        }
      continue;
      } // end of loop for (tree_attempt_index in 1:tree_attempt_max)
      
      // Tree could not be positioned in the current stand attempt
      if (!located) {
        Rprintf("[FAILED] too many unsuccessful positions \n");
        // Exit the tree loop and start a new stand attempt
        if( stand_attempt_index==stand_attempt_max)
          Rprintf(" Simulation [FAILED] %i unsuccessful attempt to locate trees \n", stand_attempt_max);
        break; 
        }      
      
    } // end of loop for(tree_index in 0:(tree_index_max-1))
    // Move to the following tree

  } // end of loop for (stand_attempt_index in 1:stand_attempt_max)
  Rprintf("Stand - Simulation terminated \n");
  
  stand = Rcpp::List::create(Rcpp::Named("quadrat.x") = quadrat_x, // Quadra x-coordinate
                             Rcpp::Named("quadrat.y")= quadrat_y, // Quadra y-coordinate
                             Rcpp::Named("dbh") = stand_dbh, // Diamtre at breast height
                             Rcpp::Named("height") = stand_height, // height
                             Rcpp::Named("trunk.height") = stand_trunk_height, // Trunk height
                             Rcpp::Named("crown.radius") = stand_crown_radius, // Crown radius
                             Rcpp::Named("x")= stand_x , Rcpp::Named("y") = stand_y);// X and Y coordinate of the tree in the plot
  
  // Rcpp::List converted in R list
  return stand;
}

/*** R

# Data simulation 

# check whether the numeric is integer
.is.int <- function(n) { return (n%%1==0) }

randomTreeInventory <- function(ntree=150, dbh.min=0.1, dbh.max=1.5, nquadrat=0, measurement.rate=0, shape=1, allom.param) {
  # Creates a random tree inventory as a data frame. 
  #
  # Args:
  #   ntree: strictly positive integer, the number of trees in the plot
  #   dbh.min, dbh.max, strictly positive float numbers, the minimum and maximum diameter at breast height, in meter
  #   nquadrat: a positive integer, the number of quadrats along one side of the square plot (i.e nquadrat * nquadrat quadrats in the plot)
  #   measurement.rate: fraction of trees that may provide additional measurements such as height, trunk height and crown radius
  #   shape: a positive float number, the shape parameter of the truncated pareto distribution that is used to generate the DBH distribution
  # 
  # Returns: 
  #   A data frame with the following vectors: 
  #   dbh, diameter at breast height in meter
  #   quadrat.x, quadrat x-coordinate
  #   quadrat.y, quadrat y-coordinate
  #   height, tree height in meter
  #   trunk.height, trunk height in meter
  #   crown.radius, crown radius in meter
  ###
  
  # require VGAM and stats packages
  stopifnot(require(VGAM))
  stopifnot(require(stats))
  
  # check input parameters
  stopifnot(.is.int(ntree) && ntree > 0)
  stopifnot(is.numeric(dbh.min) && dbh.min > 0)
  stopifnot(is.numeric(dbh.max) && dbh.max > 0)
  stopifnot(.is.int(nquadrat) && nquadrat >= 0)
  stopifnot(is.numeric(measurement.rate) && measurement.rate >= 0 && measurement.rate <= 1)
  stopifnot(is.numeric(shape) && shape > 0)
  
  # creates NA dataframe
  trees <- as.data.frame(matrix(NA, ntree, 6))
  names(trees) <- c("dbh", "quadrat.x", "quadrat.y", "height", "trunk.height", "crown.radius")
  
  # DBH distribution from truncated pareto law
  trees$dbh <- VGAM::rtruncpareto(ntree, dbh.min, dbh.max, shape)
  
  # random scattering in the quadrats
  if (nquadrat > 0) {
    quadrat.coord <- array(round(runif(2 * ntree, 0, nquadrat - 1)), dim=c(ntree, 2))
    trees$quadrat.x <- quadrat.coord[,1]
    trees$quadrat.y <- quadrat.coord[,2]
  }
  
  # provide additional field measurements (tree height, trunk height, crown radius) 
  if (measurement.rate > 0) {
    # random mask to select which trees will provide field measurmement
    measurement.mask <- array(round(runif(3 * ntree) <= measurement.rate), dim=c(ntree, 3))
    # tree height
    trees$height[which(measurement.mask[,1]>0)] <- dbhToTreeHeight(trees$dbh[which(measurement.mask[,1]>0)],
                                                                   allom.param$H_a,  allom.param$H_b,  allom.param$H_c, err=0)
    # crown radius
    trees$crown.radius[which(measurement.mask[,2]>0)] <- dbhToCrownRadius(trees$dbh[which(measurement.mask[,2]>0)],
                                                                          allom.param$Cr_x2,  allom.param$Cr_y2,  allom.param$Cr_inter1,
                                                                          allom.param$Cr_coef1,  allom.param$Cr_coef2,  allom.param$Cr_RSE, err=0)
    # trunk height
    trees$trunk.height[which(measurement.mask[,3]>0)] <- treeHeightToTrunkHeight(dbhToTreeHeight(trees$dbh[which(measurement.mask[,3]>0)],
                                                                                 allom.param$H_a,  allom.param$H_b,  allom.param$H_c, err=0),
                                                                                 allom.param$Ht_inter,  allom.param$Ht_coef,  allom.param$Ht_RSE, err=0)
  }
  
  # return tree inventory
  return(trees)
}
# DBH to tree height allometry parameters
H_a <- 35
H_b <- 0.03
H_c <- 0.8
# DBH to crown radius allometry parameters
Cr_x2 <- 3.2489
Cr_y2 <- 1.6941
Cr_inter1 <- 0.3939
Cr_coef1 <- 0.4002
Cr_coef2 <- 0.4102
Cr_RSE <- 0.494
# Tree height to trunk height allometry parameters
Ht_inter <- -0.70087
Ht_coef <- 1.04896
Ht_RSE <- 0.3885
# create dataframe
allom.param <- data.frame(H_a, H_b, H_c, Cr_x2, Cr_y2, Cr_inter1, Cr_coef1, Cr_coef2, Cr_RSE, Ht_inter, Ht_coef, Ht_RSE)

data <- randomTreeInventory(ntree=200, nquadrat=5, measurement.rate=1, allom.param = allom.param)
# Sort dbh to position the highest trees first 
data <- data[sort(data$dbh,decreasing=T,index.return=T)$ix,]

# Function to locate trees 
stand <- AlloStand(as.matrix(data), stand_attempt_max=10, tree_attempt_max =350,
                   overlap_max=0.1, plot_length=100.0, quadrat_length=20.0, voxel_size=0.5)
stand <- as.data.frame(stand)

## Representation of results 
source('/home/beauclair/Documents/allostand/v1/tree3d.R')
source('/home/beauclair/Documents/allostand/v1/allostand_utils.R')
plot.stand(stand = stand, res=20)
**/
  