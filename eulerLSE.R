
#' Least square estimates in Euler scheme
#' 
#' @param ID Vector of track identifiers
#' @param xy Matrix of locations
#' @param time Vector of observation times
#' @param gradarray Array of gradients of covariates, evaluated at
#' the locations of the observations
#' 
#' @return A list of: Bhat, the vector of estimates, and var, the
#' covariance matrix of the estimates.
eulerLSE <- function(ID=NULL, xy, time, gradarray) 
{
    n <- nrow(xy)
    if(is.null(ID))
        ID <- rep(1,n)
    
    i0 <- which(ID[-1]!=ID[-n])
    i1 <- c(1, i0+1) # first obs of each track
    i2 <- c(i0, n) # last obs of each track
    
    # matrix of covariate derivatives
    D <- 0.5*rbind(gradarray[-i2,1,], gradarray[-i2,2,])
    # diagonal matrix of time intervals
    dt <- time[-i1] - time[-i2]
    TT <- diag(rep(dt,2))
    # vector of 2-d steps
    dxy <- xy[-i1,] - xy[-i2,]
    TY <- matrix(dxy, ncol=1)
    
    # estimation
    DTTD <- t(D) %*% TT %*% D
    DTTDinv <- solve(DTTD)
    Bhat <- DTTDinv %*% t(D) %*% TY
    
    return(list(Bhat=as.vector(Bhat), var=DTTDinv))
}