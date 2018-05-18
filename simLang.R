
#' Simulate from the Langevin movement model
#' 
#' This function is based on the Euler approximation.
#' 
#' @param beta Vector of resource selection coefficients
#' @param speed Scalar speed parameter
#' @param time Vector of times of observations
#' @param xy0 Vector of coordinates (x,y) of initial location
#' @param xgrid Grid on which the covariates are known
#' @param ygrid Grid on which the covariates are known
#' @param covarray Array of values of the covariates at the points given by
#' xgrid and ygrid, of dimensions (length(xgrid),length(ygrid),length(beta)).
simLang <- function(beta, speed, time, xy0, xgrid, ygrid, covarray) {
    nbObs <- length(time)
    xy <- matrix(0,nbObs,2)
    xy[1,] <- xy0
    dt <- diff(time)
    g <- matrix(NA,nbObs,2)
    for(t in 2:nbObs) {
        cat("\rSimulating Langevin process...",round(100*t/nbObs),"%")
        g[t-1,] <- gradLogUD(beta=beta,xy=xy[t-1,],xgrid=xgrid,ygrid=ygrid,covarray=covarray)
        xy[t,] <- xy[t-1,] + 0.5 * g[t-1,] * speed * dt[t-1] + rnorm(2, 0, sqrt(speed*dt[t-1]))
    }
    cat("\n")
    return(xy)
}