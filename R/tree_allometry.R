
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

dbhToCrownRadius <- function(DBH, x2, y2, inter1, coef1, coef2, RSE, err=0) {
  # DBH in meter
  # returns Cr meter
  
  DBHcm <- 100 * DBH
  Cd <- ifelse(log(DBHcm) <= x2,
               exp(inter1 + log(DBHcm) * coef1 + err) * exp(RSE^2 / 2),
               exp(y2 + (log(DBHcm) - x2) * (coef1 + coef2) + err)*exp(RSE^2 / 2))
  
  # return crown radius
  return(0.5 * Cd)
}

dbhToTreeHeight <- function(DBH, a, b, c, err=0) {
  # DBH in meter
  # returns H in meter
  
  DBHcm <- 100 * DBH
  H <- a*(1-exp(-b*(DBHcm^c))) + err
  return(H)
}

treeHeightToTrunkHeight <- function(H, inter, coef, RSE, err=0) {
  # H in meter
  # returns Ht in meter
  
  Ht <- exp(inter + RSE^2/2 + coef*log(H)) + err
  return(Ht)
}
