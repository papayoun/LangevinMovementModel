source("OzakiFunctions.R")
#' Negative log-likelihood function for Langevin movement model
#' 
#' @param beta Parameters of the utilisation distribution
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
nllkLang <- function(beta, xy, time, gradarray, hessarray = NULL, method = "euler") {
  dt <- diff(time)
  LangevinFactor <- 0.5 # 0.5 to fit with the Langevin 0.5 factor
  gradmat <- LangevinFactor * apply(gradarray, 2, function(mat) mat %*% beta)
  n <- nrow(xy)
  if(method == "euler"){
    llk <- sum(dnorm(xy[-1, ], xy[-n,] + dt * gradmat[-n, ], sqrt(dt), log=TRUE))
  }
  if(method == "ozaki"){
    hessmat <- LangevinFactor * apply(hessarray, 2, function(mat) mat %*% beta)
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