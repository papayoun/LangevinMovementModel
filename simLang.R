
##' Simulate from the Langevin movement model
##' 
##' @param nbObs Number of locations to simulate
##' @param beta Vector of resource selection coefficients
##' @param h Discretisation step, for Euler approximation
##' @param xgrid Grid on which the covariates are known
##' @param ygrid Grid on which the covariates are known
##' @param covarray Array of values of the covariates at the points given by
##' xgrid and ygrid, of dimensions (length(xgrid),length(ygrid),length(beta)).
simLang <- function(nbObs, beta, h, xgrid, ygrid, covarray) {
    xy <- matrix(0,nbObs,2)
    for(t in 2:nbObs) {
        xy[t,] <- xy[t-1,] + 
            gradLogUD(beta=beta,xy=xy[t-1,],xgrid=xgrid,ygrid=ygrid,covarray=covarray) * h + 
            rnorm(2,0,sqrt(h))
    }
    return(xy)
}