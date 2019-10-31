source('/home/beauclair/Documents/Allostand/R/tree_allometry.R')
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
    trees$trunk.height[which(measurement.mask[,3]>0)] <- treeHeightToTrunkHeight(dbhToTreeHeight(trees$dbh[which(measurement.mask[,3]>0)], allom.param$H_a,  allom.param$H_b,  allom.param$H_c, err=0),
                                                                                 allom.param$Ht_inter,  allom.param$Ht_coef,  allom.param$Ht_RSE, err=0)
  }
  
  # return tree inventory
  return(trees)
}
