source("OzakiFunctions.R")
#' Negative log-likelihood function for Langevin movement model
#' 
#' @param par Vector of (working) model parameters: parameters of the 
#' utilisation distribution, and speed parameter.
#' @param xy Matrix of observed locations (two columns: x, y)
#' @param time Vector of times of observations
#' @param gradarray Three-dimensional array of gradients of covariate fields. 
#' The rows index time, the columns are the dimensions (x and y), and the 
#' layers index the covariates.
#' @param hessarray Optional three-dimensional array of hessians of covariate fields. 
#' The rows index time, the columns are (d/dx^2, d/dxdy, d/dydx, d/dy^2), and the 
#' layers index the covariates. Only needed if \code{method="ozaki"}.
#' @param method Method used for the approximation of the likelihood. Either
#' "euler" (default) or "ozaki".
#' 
#' @return Negative log-likelihood
nllkLang <- function(par, xy, time, ID = NULL, gradarray, 
                     hessarray = NULL, method = "euler") {
    n <- nrow(xy)
    
    # unpack parameters
    ncov <- dim(gradarray)[3]
    if(length(par)==ncov+1) {
        beta <- par[1:ncov]
        speed <- exp(par[ncov+1])
    } else if(length(par)==ncov) {
        beta <- par
        speed <- 1
    }
    
    # multiply gradients by beta coefficients
    gradmat <- 0.5 * apply(gradarray, 2, function(mat) mat %*% beta)
    
    # only one track
    if(is.null(ID))
        ID <- rep(1,nrow(xy))
    
    # indices of first and last obs in tracks
    i0 <- which(ID[-1]!=ID[-n])
    i1 <- c(1, i0+1) # first obs
    i2 <- c(i0, n) # last obs
    
    dt <- time[-i1] - time[-i2]
    
    if(method == "euler") {
        llk <- sum(dnorm(xy[-i1,], xy[-i2,] + speed * dt * gradmat[-i2,], sqrt(speed * dt), log=TRUE))
    }
    if(method == "ozaki"){
        hessmat <- 0.5 * apply(hessarray, 2, function(mat) mat %*% beta)
        Means <- getOzakiMean(xy, dt, gradmat, hessmat) # matrix of size (n-1) * 2
        CovsInv <- getOzakiCovariance(xy, dt, gradmat, hessmat, Inv = T) # matrix of size (n-1) * 4
        DetsInv <- apply(CovsInv, 1, vecDetM) # compute determinant to spot potential problems
        ToEuler <- NumProblemDet(DetsInv) # use Euler approximation for potential problematic segments
        if(any(ToEuler)){
            Nr <- sum(ToEuler)
            Means[ToEuler, ] <- xy[-n,][ToEuler,] + dt[ToEuler] * (gradmat[-n,][ToEuler, ])
            CovsInv[ToEuler, ] <- cbind(1/ dt[ToEuler], rep(0, Nr), 
                                        rep(0, Nr), 1 / dt[ToEuler]) 
            DetsInv[ToEuler] <- 1 / (dt[ToEuler] * dt[ToEuler])
        }
        llk <- sum(logdnorm2d(xs = xy[-1, ], vecMeans = Means, vecCovInvs = CovsInv, DetInv = DetsInv))
    }
    return(-llk)
}
