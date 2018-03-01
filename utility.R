
#' Interpolate 2D covariate
#' 
#' At the moment, this is based on the function \code{interp.surface} from
#' the package fields (bilinear interpolation).
#' 
#' @param xy Point where the covariate should be interpolated
#' @param xgrid Grid on which the covariate is known
#' @param ygrid Grid on which the covariate is known
#' @param covmat Matrix of values of the covariate at the points given by
#' xgrid and ygrid.
#' 
#' @return Interpolated value of the covariate at the point xy.
#' 
#' @details Note that covmat needs to rotated (as e.g. with "image"), so you 
#' might need to use something like
#'   covmat <- t(apply(as.matrix(covraster),2,rev))
#' before passing it to this function
interpCov <- function(xy,xgrid,ygrid,covmat) {
    xy <- matrix(xy,ncol=2)
    covxyz <- list(x=xgrid,y=ygrid,z=covmat)
    return(interp.surface(covxyz,xy))
}

#' Gradient of the log of the utilisation distribution
#' 
#' @param beta Vector of resource selection coefficients
#' @param xy Point at which the gradient should be evaluated
#' @param xgrid Grid on which the covariates are known
#' @param ygrid Grid on which the covariates are known
#' @param covarray Array of values of the covariates at the points given by
#' xgrid and ygrid, of dimensions (length(xgrid),length(ygrid),length(beta)).
#' 
#' @return Gradient of the log-UD in xy.
gradLogUD <- function(beta, xy, xgrid, ygrid, covarray) {
    gradvals <- apply(covarray, 3, 
                      function(covmat) grad(interpCov, x=xy, xgrid=xgrid, ygrid=ygrid,covmat=covmat))
    return(gradvals%*%beta)
}