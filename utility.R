
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
    gradvals <- apply(covarray, 3, function(covmat) 
        grad(interpCov, x=xy, xgrid=xgrid, ygrid=ygrid,covmat=covmat))
    return(gradvals%*%beta)
}

#' Gradient of covariate field
#' 
#' @param xy Matrix of locations where the gradient should be evaluated
#' @param xgrid Grid on which the covariates are known
#' @param ygrid Grid on which the covariates are known
#' @param covarray Array of values of the covariates at the points given by
#' xgrid and ygrid, of dimensions (length(xgrid),length(ygrid),length(beta)).
#' 
#' @return Three-dimensional array of gradients of covariate fields. The rows index time, 
#' the columns are the dimensions (x and y), and the layers index the covariates.
covGrad <- function(xy, xgrid, ygrid, covarray) {
    gradarray <- array(NA,dim=c(nrow(xy),2,dim(covarray)[3]))
    for(i in 1:dim(covarray)[3]) {
        gradarray[,,i] <- t(apply(xy,1,function(x) 
            grad(interpCov,x,xgrid=xgrid,ygrid=ygrid,covmat=covarray[,,i])))
    }
    return(gradarray)
}

#' Hessian of covariate field
#' 
#' @param xy Matrix of locations where the gradient should be evaluated
#' @param xgrid Grid on which the covariates are known
#' @param ygrid Grid on which the covariates are known
#' @param covarray Array of values of the covariates at the points given by
#' xgrid and ygrid, of dimensions (length(xgrid),length(ygrid),length(beta)).
#' 
#' @return Three-dimensional array of second derivatives of covariate fields. The rows
#' index time, the four columns correspond to the elements of the Hessian matrix 
#' (d/dx^2, d/dydx,d/dxdy, d/dy^2), and the layers index the covariates.
covHessian <- function(xy, xgrid, ygrid, covarray) {
    hessarray <- array(NA,dim=c(nrow(xy),4,dim(covarray)[3]))
    for(i in 1:dim(covarray)[3]) {
        hessarray[,,i] <- t(apply(xy,1, function(x)
            c(hessian(interpCov,x,xgrid=xgrid,ygrid=ygrid,covmat=covarray[,,i]))))
    }
    return(hessarray)
}