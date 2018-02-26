
##' Resource selection function
##' 
##' @param xy Matrix of points where the RSF should be evaluated
##' @param beta Vector of resource selection coefficients
##' @param covlist List of covariate rasters
RSF <- function(xy,beta,covlist)
{
    xy <- matrix(xy,ncol=2)
    lim <- extent(covlist[[1]])
    cov <- sapply(covlist,extract,xy)
    reg <- as.vector(cov%*%beta)
    reg[which(xy[,1]<lim[1] | xy[,1]>lim[2] | xy[,2])<lim[3] | xy[,2]>lim[4]] <- -Inf
    return(exp(reg))
}

##' Smooth resource selection function
##' 
##' This function uses interpCov to define a smooth interpolation of the
##' resource selection function.
##' 
##' @param beta Vector of resource selection coefficients
##' @param xy Point at which the RSF should be evaluated
##' @param xgrid Grid on which the covariates are known
##' @param ygrid Grid on which the covariates are known
##' @param covarray Array of values of the covariates at the points given by
##' xgrid and ygrid, of dimensions (length(xgrid),length(ygrid),length(beta)).
##' 
##' @return Value of the RSF at the point xy.
smoothRSF <- function(xy, beta, xgrid, ygrid, covarray) {
    covvals <- apply(covarray, 3, interpCov, xy=xy, xgrid=xgrid, ygrid=ygrid)
    return(exp(sum(beta*covvals)))
}