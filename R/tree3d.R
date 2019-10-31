ellipsoid3d <- function(rx=1,ry=1,rz=1,n=30,ctr=c(0,0,0),...) {
  # Adapted from demo("shapes3d") package rgl
  # Returns a 3D ellipsoid from rgl package
  
  trans <- diag(4)
  degvec <- seq(0,pi,length=n)
  ecoord2 <- function(p) {
    c(rx*cos(p[1])*sin(p[2]),ry*sin(p[1])*sin(p[2]),rz*cos(p[2])) }
  v <- apply(expand.grid(2*degvec,degvec),1,ecoord2)
  v <- rbind(v,rep(1,ncol(v))) ## homogeneous
  e <- expand.grid(1:(n-1),1:n)
  i1 <- apply(e,1,function(z)z[1]+n*(z[2]-1))
  i2 <- i1+1
  i3 <- (i1+n-1) %% n^2 + 1
  i4 <- (i2+n-1) %% n^2 + 1
  i <- rbind(i1,i2,i4,i3)
  return(rotate3d(qmesh3d(v,i,material=list(...)),matrix=trans))
}

trunk <- function(x, y, trunk.height, dbh) {
  # Returns a tree trunk as a 3D cylinder from rgl package
  #
  # Args:
  #   x: x-coordinate in metre from origin (0, 0) of the plot
  #   y: x-coordinate in metre from origin (0, 0) of the plot
  #   dbh: diametre at breast height in metre 
  #   trunk.height: trunk height in metre
  #
  # Returns:
  #   A cylinder3d object from rgl package
  
  by <- trunk.height / 3
  z <- seq(0, to=trunk.height, by=by)
  nz <- length(z)
  radius = 0.5*dbh
  return(cylinder3d(cbind(rep(x, nz), rep(y, nz), z), radius = radius, sides=96))
}

draw.tree <- function(x=0, y=0, dbh=50, tree.height=33, crown.width=10, crown.base.height=25, col="green", trunk.draw = TRUE, res=30) {
  # Draw a 3D tree
  #
  # Args:
  #   x: x-coordinate in metre from origin (0, 0) of the plot
  #   y: x-coordinate in metre from origin (0, 0) of the plot
  #   dbh: diametre at breast height in metre 
  #   tree.height: tree height in metre
  #   crown.width: crown width in metre
  #   crown.base.height: crown base height in metre
  #   col: colour of the crown
  #   trunk.draw: whether to draw the trunk
  #   res: the resolution for drawing the crown. res = 30 for instance means that the ellipsoid will be drawn with 30 crescents of 30 facets each.
  
  # draw trunk
  if (trunk.draw)
    shade3d(trunk(x, y, crown.base.height, dbh), col="chocolate4")
  # draw crown
  crown.height <- tree.height - crown.base.height
  crown.centre <- crown.base.height + 0.4 * crown.height
  shade3d(translate3d(rotate3d(scale3d(ellipsoid3d(n=res), 0.5*crown.width, 0.5*crown.width, 0.5*crown.height), angle=pi, x=0, y=1, z=0), x=x, y=y, z=crown.centre), col=col)
}