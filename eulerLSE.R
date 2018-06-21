
#' Least square estimates in Euler scheme
#' 
#' @param ID Vector of track identifiers
#' @param xy Matrix of locations
#' @param time Vector of observation times
#' @param gradarray Array of gradients of covariates, evaluated at
#' the locations of the observations
#' @param withspeed Logical. If TRUE, the speed parameter is estimated.
#' 
#' @return A list of: est, the vector of estimates, and var, the
#' covariance matrix of the estimates.
eulerLSE <- function(ID=NULL, xy, time, gradarray, withspeed=TRUE) 
{
    n <- nrow(xy)
    if(is.null(ID))
        ID <- rep(1,n)
    
    i0 <- which(ID[-1]!=ID[-n])
    i1 <- c(1, i0+1) # first obs of each track
    i2 <- c(i0, n) # last obs of each track
    
    # matrix of covariate derivatives
    D <- 0.5*rbind(gradarray[-i2,1,], gradarray[-i2,2,])
    # vector of time intervals (requires less memory than diagonal matrix)
    dt <- time[-i1] - time[-i2]
    T <- rep(sqrt(dt),2)
    TT <- rep(dt,2)
    
    # vector of 2-d steps
    dxy <- xy[-i1,] - xy[-i2,]
    Y <- matrix(dxy/sqrt(dt), ncol=1)
    TY <- matrix(dxy, ncol=1)
    
    # estimation
    DTTD <- t(D*TT) %*% D
    DTTDinv <- solve(DTTD)
    Bhat <- DTTDinv %*% t(D) %*% TY
    
    if(withspeed) {
        # estimate speed
        fitted <- t(T*t(D)) %*% Bhat
        SSE <- sum((Y-fitted)^2)
        speed <- SSE/(2*(n-1)-length(Bhat))        

        est <- c(Bhat/speed,speed)
        var <- DTTDinv/speed
    } else {
        est <- Bhat
        var <- DTTDinv
    }

    return(list(est=est, var=var))
}
