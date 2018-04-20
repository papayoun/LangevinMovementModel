
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
    i2 <- c(i0, n) # last obs of each track
    
    # matrix of covariate derivatives
    X <- 0.5*rbind(gradarray[-i2,1,], gradarray[-i2,2,])
    # diagonal matrix of time intervals
    TT <- diag(rep(diff(time)[-i0],2))
    # vector of 2-d steps
    TZ <- matrix(apply(xy, 2, function(c) diff(c)[-i0]), ncol=1)
    
    # estimation
    XTTX <- t(X) %*% TT %*% X
    XTTXinv <- solve(XTTX)
    Bhat <- XTTXinv %*% t(X) %*% TZ
    
    return(list(Bhat=as.vector(Bhat), var=XTTXinv))
}