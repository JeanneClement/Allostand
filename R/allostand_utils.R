Crown_allometry_error <- function(err.type="none", sd=0) {
  # Computes the error for the crown allometry.
  #
  # Args:
  #   err.type: a character, the type of error, among "none", "negative" or "maximal"
  #     err.type = none, computes allometry without any error
  #     err.type = negative, computes allometry with negative error (a smaller crown)
  #     err.type = maximal, computes allometry with maximized error (the smallest crown)
  #
  # Returns:
  #   The error for the crown allometry, as a double
  switch(err.type,
         none={err<-0},
         negative={err <- -abs(sample(rnorm(1000, mean = 0, sd = sd), 1))},
         maximal={err <- min(rnorm(1000, mean = 0, sd = sd))},
         stop(paste("Function Crown_allometry_error, argument err.type=\"",err.type, "\", it must be either \"none\", \"negative\" or \"maximal\"", sep=""))
  )
  return(err)
}

CrownToVoxels <- function(crown.rx, crown.ry, crown.rz, dvox) {
  
  # Crown radius in voxels
  crown.vx.rx <- ceiling(crown.rx / dvox) + 1
  crown.vx.ry <- ceiling(crown.ry / dvox) + 1
  crown.vx.rz <- ceiling(crown.rz / dvox) + 1
  
  x <- array(rep(seq(0,crown.vx.rx),times=crown.vx.ry*crown.vx.rz), c(crown.vx.rx+1,crown.vx.ry+1,crown.vx.rz+1))
  y <- array(rep(seq(0,crown.vx.ry),each=crown.vx.rx+1,times=crown.vx.rz), c(crown.vx.rx+1,crown.vx.ry+1,crown.vx.rz+1))
  z <- array(rep(seq(0,crown.vx.rz),each=(crown.vx.rx+1)*(crown.vx.ry+1),times=1), c(crown.vx.rx+1,crown.vx.ry+1,crown.vx.rz+1))
  
  crown.box <- (x^2/(crown.vx.rx)^2+y^2/(crown.vx.ry)^2+z^2/(crown.vx.rz)^2 <= 1)
  vx.in <- which(crown.box)
  
  crown.q111 <- array(c(x[vx.in], y[vx.in], z[vx.in]), c(length(vx.in),3))
  q000 <- c(-1,-1,-1)
  q001 <- c(-1,-1, 1)
  q010 <- c(-1, 1,-1)
  q011 <- c(-1, 1, 1)
  q100 <- c( 1,-1,-1)
  q101 <- c( 1,-1, 1)
  q110 <- c( 1, 1,-1)
  
  crown <- rbind(crown.q111, t(apply(crown.q111,1,function(n) n*q000)), t(apply(crown.q111,1,function(n) n*q001)), t(apply(crown.q111,1,function(n) n*q010)), t(apply(crown.q111,1,function(n) n*q011)), t(apply(crown.q111,1,function(n) n*q100)), t(apply(crown.q111,1,function(n) n*q101)), t(apply(crown.q111,1,function(n) n*q110)))
  colnames(crown) <- c('x', 'y', 'z')
  
  return(unique(crown))
}

CrownToPixels <- function(crown.rx, crown.ry, pixel.size) {
  
  # Crown radius in voxels
  crown.px.rx <- ceiling(crown.rx / pixel.size) + 1
  crown.px.ry <- ceiling(crown.ry / pixel.size) + 1
  
  x <- array(rep(seq(0,crown.px.rx),times=crown.px.ry), c(crown.px.rx+1,crown.px.ry+1))
  y <- array(rep(seq(0,crown.px.ry),each=crown.px.rx+1), c(crown.px.rx+1,crown.px.ry+1))
  
  crown.box <- (x^2/(crown.px.rx)^2+y^2/(crown.px.ry)^2 <= 1)
  px.in <- which(crown.box)
  
  crown.q11 <- array(c(x[px.in], y[px.in]), c(length(px.in),2))
  q00 <- c(-1,-1)
  q01 <- c(-1, 1)
  q10 <- c( 1,-1)
  
  crown <- rbind(crown.q11, t(apply(crown.q11,1,function(n) n*q00)), t(apply(crown.q11,1,function(n) n*q01)), t(apply(crown.q11,1,function(n) n*q10)))
  colnames(crown) <- c('x', 'y')
  
  return(unique(crown))
}

read.stand <- function(file) {
  # Read a stand CSV file and return it as a data frame. The CSV file must be tab-separated and dot as decimal separator.
  #
  # Args:
  #   file: a character, the file name of the CSV file
  #
  # Returns:
  #   The stand as a data frame, with the trees in rows and at least the following columns:
  #   stand$x the x coordinate
  #   stand$y the y coordinate
  #   stand$dbh the diametre at breast height (in centimetre)
  #   stand$height the tree height (in metre)
  #   stand$trunk.height the trunk height (in metre)
  #   stand$crown.radius the crown radius (in metre)
  stand.file <- file.path(file)
  if (!file.exists(stand.file)) {
    stop("Function read.stand, could not find stand file ", stand.file)
  }
  stand <- read.csv(stand.file, header=TRUE, sep="\t", na.strings="NA", dec=".", strip.white=TRUE)
  # Update column names if old ALLOSTAND version
  colnames(stand)[colnames(stand)=="ID_arbre"] <- "id"
  colnames(stand)[colnames(stand)=="D"] <- "dbh"
  colnames(stand)[colnames(stand)=="SP"] <- "species"
  colnames(stand)[colnames(stand)=="H"] <- "height"
  colnames(stand)[colnames(stand)=="Ht"] <- "trunk.height"
  colnames(stand)[colnames(stand)=="Cr"] <- "crown.radius"
  colnames(stand)[colnames(stand)=="X_temp"] <- "x"
  colnames(stand)[colnames(stand)=="Y_temp"] <- "y"
  stand$RAS <- NULL
  
  # return stand
  return(stand)
}

TreeToVoxels <- function(tree, voxel.size, plot.dim) {
  # Transform a tree into a vector of voxels, in a given plot and a given spatial resolution.
  #
  # Args:
  #   tree: a row of the stand data frame, with at least the following columns
  #     tree$crown.radius
  #     tree$height
  #     tree$trunk.height
  #     tree$x
  #     tree$y
  #   voxel.size: an integer, spatial resolution of the voxelised plot, in metre
  #   plot.dim (optional: a 3d integer vector, the three dimensions, x, y, z, in metre, of the plot
  #
  # Returns:
  #   The tree as a vector of voxels
    
  # Voxelised crown
  crown <- CrownToVoxels(tree$crown.radius, tree$crown.radius, 0.5*(tree$height-tree$trunk.height), voxel.size)
  # Tree centre
  tree.centre = c(tree$x, tree$y, 0.5*(tree$height+tree$trunk.height))
  tree.centre.vx <- ceiling(tree.centre / voxel.size) + 1
  # Voxelised tree
  tree.vx <- t(apply(crown,1,function(r)r+tree.centre.vx))
  
  # Trim voxels outside plot
  if (!missing(plot.dim)) {
    # mockup dimension in voxels
    plot.dim.voxel = ceiling(plot.dim / voxel.size)
    # indexes of voxels outside plot
    voxel.out <- t(apply(tree.vx, 1, function(x) any(x<1 | x>plot.dim.voxel)))
    if (length(which(voxel.out))) {
      tree.vx <- tree.vx[!voxel.out,]
    }
  }
  
  # Return list of voxels
  return(tree.vx)
}

TreeToPixels <- function(tree, pixel.size, plot.dim) {
  # Pixelised crown
  crown <- CrownToPixels(tree$crown.radius, tree$crown.radius, pixel.size)
  # Tree centre
  tree.centre = c(tree$x, tree$y)
  tree.centre.px <- ceiling(tree.centre / pixel.size) + 1
  # Pixelised tree
  tree.px <- t(apply(crown,1,function(r)r+tree.centre.px))
  
  # Trim pixels outside plot
  if (!missing(plot.dim)) {
    # mockup dimension in pixels
    if (length(plot.dim)==1) plot.dim <- rbind(plot.dim, plot.dim, deparse.level = 0)
    plot.dim.pixel = ceiling(plot.dim[1:2] / pixel.size)
    # indexes of pixels outside plot
    pixel.out <- t(apply(tree.px, 1, function(x) any(x<1 | x>plot.dim.pixel)))
    if (length(which(pixel.out))) {
      tree.px <- tree.px[!pixel.out,]
    }
  }
  
  # Return list of pixels
  return(tree.px)
}

TreeHMax <- function(tree, pixel.size, plot.dim) {
  # Computes the height in metre of every pixel of the tree in the 2D projection
  #
  # Args:
  #   tree: a row of the stand data frame, with at least the following columns
  #     tree$crown.radius
  #     tree$height
  #     tree$trunk.height
  #     tree$x
  #     tree$y
  #   pixel.size: an integer, spatial resolution of the voxelised plot, in metre
  #   plot.dim (optional: a 3d integer vector, the three dimensions, x, y, z, in metre, of the plot
  #
  # Return:
  #   A vector of the same length than the TreeToPixels function that contains the height in metre of every pixel of the tree crown (2D projection).
  #   If plot.dim is specified, every pixel outside the plot boundaries are discarded.
  
  # Pixelised crown
  crown <- CrownToPixels(tree$crown.radius, tree$crown.radius, pixel.size)
  # Crown radius in pixel
  radius <- max(crown)
  
  # Trim pixels outside plot
  if (!missing(plot.dim)) {
    # mockup dimension in pixels
    if (length(plot.dim)==1) plot.dim <- rbind(plot.dim, plot.dim, deparse.level = 0)
    plot.dim.pixel = ceiling(plot.dim[1:2] / pixel.size)
    # Tree centre
    tree.centre.px <- ceiling(c(tree$x, tree$y) / pixel.size) + 1
    # Pixelised tree
    tree.px <- t(apply(crown,1,function(r)r+tree.centre.px))
    # indexes of pixels outside plot
    pixel.out <- t(apply(tree.px, 1, function(x) any(x<1 | x>plot.dim.pixel)))
    if (length(which(pixel.out))) {
      crown <- crown[!pixel.out,]
    }
  }
  
  # Half crown height
  rz <- 0.5*(tree$height - tree$trunk.height)
  # Crown centre height
  cz <- 0.5*(tree$height+tree$trunk.height)
  # Crown maximal height in metre
  return(cz + rz * sqrt(abs(1 - (crown[,"x"]/radius)^2 - (crown[,"y"]/radius)^2)))
}

TreeInMockup <- function(tree, plot.dim, voxel.size) {
  
  # tree as voxels
  tree.voxel <- TreeToVoxels(tree, voxel.size)
  # mockup dimension in voxels
  plot.dim.voxel = ceiling(plot.dim / voxel.size)
  # initialises mockup
  mockup <- array(rep(0,times=prod(plot.dim.voxel)), plot.dim.voxel)
  # locate tree in mockup
  mockup[tree.voxel] <- 1
  
  # return mockup
  return(mockup)
}

Canopy <- function(stand, covering.max, plot.dim, pixel.size) {
  
  # vector of tree indexes forming the canopy
  canopy <- NULL
  # stand sorted by tree height
  standh <- stand[order(stand$height, decreasing=TRUE),]
  # mockup dimension in pixels
  if (length(plot.dim)==1) plot.dim <- rbind(plot.dim, plot.dim, deparse.level = 0)
  plot.dim.pixel = ceiling(plot.dim[1:2] / pixel.size)
  # initialises mockup
  mockup <- array(rep(0,times=prod(plot.dim.pixel)), plot.dim.pixel)
  # loop over trees
  for (tree.index in seq(1, nrow(stand))) {
    # Pixelised tree
    tree.pixel <- TreeToPixels(standh[tree.index,], pixel.size, plot.dim)
    covering <- sum(mockup[tree.pixel]>0)/nrow(tree.pixel)
    if (covering <= covering.max) {
      # Assign tree index into mockup
      mockup[tree.pixel] <- as.numeric(rownames(standh[tree.index,]))
      canopy <- rbind(canopy, as.numeric(rownames(standh[tree.index,])))
    }
  }
  
# Debug purpose: this plot and stand.plot2d(stand[Canopy(stand, plot.dim, px), ], plot.dim, px) must match  
#   ncolor <- 7
#   tree.color <- rainbow(ncolor)[mockup[which(mockup>0,1:2)]%%ncolor+1]
#   plot(which(mockup>0, c(1,2)), pch='.', col=tree.color)
  
  return(canopy)
}

plot.mockup2d <- function(stand, plot.dim, pixel.size) {
  
  # get mockup2d
  mockup <- StandToMockup2d(stand, plot.dim, pixel.size)
  # colorbar
  ncolor <- 7
  tree.color <- topo.colors(ncolor)[mockup[which(mockup>0,1:2)]%%ncolor+1]
  # scatter plot
  plot(which(mockup>0, c(1,2)), pch=16, col=tree.color)
}

plot.mockup3d <- function(stand, plot.dim, voxel.size) {
  
  # get mockup3d
  mockup <- StandToMockup3d(stand, plot.dim, voxel.size)
  # colorbar
  ncolor <- 7
  tree.color <- rainbow(ncolor)[mockup[which(mockup>0,1:2)]%%ncolor+1]
  # rgl plot3d
  library('rgl')
  plot3d(which(mockup > 0, 1:3), col=tree.color)
}

plot.stand <- function(stand, res=48) {
  # Draw 3D stand with RGL library. Must source tree3d.R script beforehand.
  # For an animation type play3d(spin3d(axis = c(0, 0, 1), rpm = 1)) after calling this function
  
  library(rgl)
  #open3d()
  greens <- colours()[grep("green", colours())]
  icol <- round(runif(nrow(stand), 1, length(greens)))
  greens <- greens[icol]
  for (i in seq(nrow(stand))) {
    tree <- stand[i,]
    draw.tree(x=tree$x, y=tree$y, dbh=tree$dbh, tree.height=tree$height, crown.width = 2*tree$crown.radius, crown.base.height = tree$trunk.height, col=greens[i], res=res)
  }  
}

StandToMockup3d <- function(stand, plot.dim, voxel.size) {
  
  # mockup dimension in voxels
  plot.dim.voxel = ceiling(plot.dim / voxel.size)
  # initialises mockup
  mockup <- array(rep(0,times=prod(plot.dim.voxel)), plot.dim.voxel)
  # loop over trees
  by(stand, nrow(stand):1, function(tree) mockup[TreeToVoxels(tree, voxel.size, plot.dim)] <<- as.numeric(rownames(tree)))
  # return mockup
  return(mockup)  
}

StandToMockup2d <- function(stand, plot.dim, pixel.size) {
  
  # mockup dimension in pixels
  if (length(plot.dim)==1) plot.dim <- rbind(plot.dim, plot.dim, deparse.level = 0)
  plot.dim.pixel = ceiling(plot.dim[1:2] / pixel.size)
  # initialises mockup
  mockup <- array(rep(0,times=prod(plot.dim.pixel)), plot.dim.pixel)
  # loop over trees
  by(stand, nrow(stand):1, function(tree) mockup[TreeToPixels(tree, pixel.size, plot.dim)] <<- as.numeric(rownames(tree)))
  # return mockup
  return(mockup)  
}

StandToDart <- function(stand, fin, fout, version = "5.6.1") {
  
  if (missing(stand)) {
    # read stand data from CSV file
    cat("Reading stand from file", fin, "\n")
    stand <- read.stand(fin)
  }
  
  if (missing(fout)) {
    fout <- paste(sub("\\.[[:alnum:]]+$", "", fin), ".dart.csv", sep="")
  }
  
  # initilise tree data frame
  tree <- as.data.frame(matrix(NA, nrow(stand), 11))
  
  # fill up tree data frame
  # tree species (index)
  #tree[,1] <- sapply(stand$species,function(x)which(levels(stand$species)==x))
  tree[,1] <- 0
  # Tree x coordinate
  tree[,2] <- stand$x
  # Tree Y coordinate
  tree[,3] <- stand$y
  # Trunk height below crown
  tree[,4] <- stand$trunk.height
  # Trunk height within crown
  tree[,5] <- 1
  # Trunk diameter below crown
  tree[,6] <- stand$dbh / 100.
  # Crown type (0 ellipsoid)
  tree[,7] <- 0
  # Crown height
  tree[,8] <- stand$height - stand$trunk.height
  # Crown azimuth
  tree[,9] <- 0
  # Crown 1st axis
  tree[,10] <- stand$crown.radius*2
  # Crown 2nd axis
  tree[,11] <- stand$crown.radius*2
  
  # write dart tree file
  cat("Writing DART tree_position.txt file (version", version, ") in", fout, "\n")
  if (version < "5.6.1") {
    write.table(tree, fout, row.names=F, col.names=F, sep="\t", dec=".")
  } else {
    # since DART 5.6.1, tree file contains column headers
    tree <- tree[, c(-5, -9)]
    #colnames(tree) <- c("SPECIES_ID", "POS_X", "POS_Y", "T_HEI_BELOW", "T_HEI_WITHIN", "T_DIA_BELOW", "C_TYPE", "C_HEI", "C_GEO_1", "C_GEO_2")
    colnames(tree) <- c("SPECIES_ID", "POS_X", "POS_Y", "T_HEI_BELOW", "T_DIA_BELOW", "C_TYPE", "C_HEI", "C_GEO_1", "C_GEO_2")
    write.table(tree, fout, row.names=F, col.names=T, sep="\t", quote = F, dec=".")
  }
}
  
read.dart <- function(file, as.stand = F, version = "5.6.1") {
  #
  # Read DART tree file (handles both DART version <= 5.6.0 and new format since 5.6.1)
  # 
  ## Args:
  #   file: a character, the file name of the CSV file
  #   as.stand: whether to convert the DART tree file as a stand data frame (refer to read.stand for details of the stand data frame)
  #
  # Returns:
  #   The stand as a data frame, with the trees in rows and at least the following columns:
  #   stand$x the x coordinate
  #   stand$y the y coordinate
  #   stand$dbh the diametre at breast height (in centimetre)
  #   stand$height the tree height (in metre)
  #   stand$trunk.height the trunk height (in metre)
  #   stand$crown.radius the crown radius (in metre)
 
  #
  dart.file <- file.path(file)
  if (!file.exists(dart.file)) {
    stop("Function read.stand, could not find DART tree file ", dart.file)
  }
  tree <- read.csv(dart.file, header=TRUE, sep="\t", na.strings="NA", dec=".", strip.white=TRUE, comment.char = "*")
  if (version < "5.6.1") {
    # prior to DART 5.6.1 tree file does not contain mandatory headers
  colnames(tree) <- c("SPECIES_ID",
                      "POS_X", "POS_Y",
                      "T_HEI_BELOW", "T_HEI_WITHIN", "T_DIA_BELOW",
                      "C_TYPE", "C_HEI", "C_ROT_INT", "C_GEO_1", "C_GEO_2")
  }
  if (as.stand) {
    stand <- as.data.frame(matrix(NA, nrow(tree), 6))
    stand[,1] <- tree$POS_X
    stand[,2] <- tree$POS_Y
    stand[,3] <- 100. * tree$T_DIA_BELOW
    stand[,4] <- tree$T_HEI_BELOW + tree$C_HEI
    stand[,5] <- tree$T_HEI_BELOW
    stand[,6] <- 0.25 * (tree$C_GEO_1 + tree$C_GEO_2)
    colnames(stand) <- c("x", "y", "dbh", "height", "trunk.height", "crown.radius")
    return(stand)
  } else
    return(tree)  
}

HMaxMockup <- function(stand, plot.dim, pixel.size) {
  # Returns a 2D mockup of the plot with the highest point (in metre) occupied by the trees in every pixel. Similar to StandToMockup2d function but returns heights instead of tree indexes.
  
  # mockup dimension in pixels
  if (length(plot.dim)==1) plot.dim <- rbind(plot.dim, plot.dim, deparse.level = 0)
  plot.dim.pixel = ceiling(plot.dim[1:2] / pixel.size)
  # initialises mockup
  mockup <- array(rep(0,times=prod(plot.dim.pixel)), plot.dim.pixel)
  # loop over trees
  for (itree in seq(nrow(stand))) {
    tree <- stand[itree,]
    hmax <- TreeHMax(tree, pixel.size, plot.dim)
    tree.px <- TreeToPixels(tree, pixel.size, plot.dim)
    mockup[tree.px] <- apply(cbind(mockup[tree.px], hmax), 1, max)
  }
  # return mockup
  return(mockup)  
}

GapFraction <- function(mockup.hmax, hmin = 2, pixel.size, plot=FALSE) {
  # Computes the gap fraction of the plot given a HMax 2D mockup (result of the HMaxMockup function) and the minimum vegetation height to be taken into account.
  # 
  
  if (plot) {
    plot(which(mockup.hmax >= hmin, c(1,2)), col="red",pch=".", xlab="X axis (m)", ylab="Y axis (m)", xaxt="n", yaxt="n",las=1)
    x.ticks <- pretty(seq(dim(mockup.hmax)[1]))
    x.ticklabs <- as.character(x.ticks * pixel.size)
    axis(1,at=x.ticks, labels=x.ticklabs, cex.axis=0.8)
    y.ticks <- pretty(seq(dim(mockup.hmax)[2]))
    y.ticklabs <- as.character(y.ticks * pixel.size)
    axis(2,at=y.ticks, labels=y.ticklabs, cex.axis=0.8)
  }
  return(length(which(mockup.hmax < hmin)) / length(mockup.hmax))  
}

Filename <- function(file) {
  # Extract the filename without extension of a full filename
  
  if (missing(file))
    stop("Argument file must be provided")
  
  # Extract the filename witout the extension
  return(sub("^([^.]*).*", "\\1", basename(file)))
}



