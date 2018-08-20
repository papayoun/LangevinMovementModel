source("OzakiFunctions.R")
#' Negative log-likelihood function for Langevin movement model
#' 
#' @param par Vector of (working) model parameters: parameters of the 
#' utilisation distribution, and diffusion parameter gamma.
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
        gamma <- exp(par[ncov+1])
    } else if(length(par)==ncov) {
        beta <- par
        gamma <- 1
    }
    
    # multiply gradients by beta coefficients
    gradmat <- 0.5 * gamma * apply(gradarray, 2, function(mat) mat %*% beta)
    
    # only one track
    if(is.null(ID))
        ID <- rep(1,nrow(xy))
    
    # indices of first and last obs in tracks
    i0 <- which(ID[-1]!=ID[-n])
    i1 <- c(1, i0+1) # first obs
    i2 <- c(i0, n) # last obs
    
    dt <- time[-i1] - time[-i2]
    
    if(method == "euler") {
        llk <- sum(dnorm(xy[-i1,], xy[-i2,] + dt * gradmat[-i2,], sqrt(gamma * dt), log=TRUE))
    }
    if(method == "ozaki"){
        hessmat <- 0.5 * gamma * apply(hessarray, 2, function(mat) mat %*% beta)
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
nllkLangForOptim <- function(vecPars, xy, time, ID = NULL, gradarray, 
                             hessarray = NULL, method = "euler"){
  driftParams <- vecPars[-length(vecPars)]
  diffusionParam <- vecPars[length(vecPars)]
  nllkLang(beta = driftParams, gamma = diffusionParam, xy = xy, time =  time, ID =  ID,
           gradarray = gradarray, hessarray = hessarray, method = method)
}
eulerEstimate <- function(xy, time, gradarray){
  n <- nrow(xy)
  d <- dim(gradarray)[3] # Number of covariates
  timeIncrementMatrix <- diag(rep(diff(time), 2))# 2 * (n - 1) * 2 * (n - 1) matrix
  designMatrix <- 0.5 * rbind(gradarray[-n , 1 , 1 : d ], gradarray[-n , 2 , 1:d])# 2 (n - 1) * d matrix
  incrementVector <- matrix(apply(xy, 2, diff)  , ncol=1)# 2 * (n - 1) * d matrix
  estimatedCovariance <- solve(t(designMatrix) %*% timeIncrementMatrix %*% designMatrix)
  betaEstimation <- estimatedCovariance %*% t(designMatrix) %*% incrementVector;
  return(list(coefficients = as.numeric(betaEstimation), covariance = estimatedCovariance))
}

eulerEstimateUnknownGamma <- function(xy, time, gradarray){
  n <- nrow(xy)
  J <- dim(gradarray)[3] # Number of covariates
  timeIncrementMatrix <- rep(diff(time), 2)# 2 * (n - 1) * 2 * (n - 1) matrix
  designMatrix <- 0.5 * rbind(gradarray[-n , 1 , 1:J], gradarray[-n , 2 , 1:J])# 2 (n - 1) * d matrix
  TtimesD <- timeIncrementMatrix * designMatrix
  incrementVector <- matrix(apply(xy, 2, diff)  , ncol = 1)# 2 * (n - 1) * d matrix
  observations <- matrix(apply(xy, 2, function(z) diff(z) / sqrt(diff(time)))  , ncol = 1)# 2 * (n - 1) * J matrix
  estimatedCovariance <- solve(t(designMatrix) %*% TtimesD)
  estimatedNu <- estimatedCovariance %*% t(designMatrix) %*% incrementVector;
  degreeFreedom <- 2 * (n - 1) - J
  estimatedGamma <- colSums((observations -  sqrt(timeIncrementMatrix) * designMatrix %*% estimatedNu)^2) / degreeFreedom
  estimatedBeta <- estimatedNu / estimatedGamma * (degreeFreedom - 2) / degreeFreedom
  estimatedBetaCovariance <- (2 * estimatedBeta %*% t(estimatedBeta) / (degreeFreedom - 4) 
                              + estimatedCovariance / estimatedGamma * (1 + 2 / (degreeFreedom - 4)))
  confidenceIntervals <- t(sapply(1:length(estimatedBeta), function(j){
    estimatedBeta[j] + c(1, -1) * qnorm(0.025) *  sqrt(estimatedBetaCovariance[j, j])
  }))
  rownames(confidenceIntervals) <- rownames(estimatedBetaCovariance) <- colnames(estimatedBetaCovariance) <- paste0("beta", 1:J)
  return(list(betaHat = as.numeric(estimatedBeta), gammaHat  = estimatedGamma,
              betaHatCovariance = estimatedBetaCovariance, betaHat95CI = confidenceIntervals))
}
