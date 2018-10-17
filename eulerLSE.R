
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
#' 
#' 
eulerLSE <- function(ID=NULL, xy, time, gradarray, withspeed=TRUE) 
{
  n <- nrow(xy)
  J <- dim(gradarray)[3]
  # if(any(gradarray == 0)){
  #   valRef <- min(abs(gradarray[gradarray != 0])) / 1000
  #   newGradArray <- do.call(function(...) abind(..., along = 3), lapply(1:J, function(j){
  #     x <- gradarray[ ,, j]
  #     nChange <- sum(x == 0)
  #     if(nChange > 0){
  #       x[x == 0] <- runif(nChange,  - valRef, valRef)
  #     }
  #     return(x)
  #   }))
  #   gradarray <- newGradArray
  # }
  if(is.null(ID))
    ID <- rep(1,n)
  
  i0 <- which(ID[-1] != ID[-n])
  i1 <- c(1, i0 + 1) # first obs of each track
  i2 <- c(i0, n) # last obs of each track
  
  # matrix of covariate derivatives
  designMatrix <- 0.5 * rbind(gradarray[-i2, 1,], gradarray[-i2, 2,])
  
  # vector of time intervals (requires less memory than diagonal matrix)
  timeIncrements <- rep(time[-i1] - time[-i2], 2)
  rootTimeIncrements <- sqrt(timeIncrements)
  
  # vector of 2-d steps
  
  positionIncrements <- as.numeric(xy[-i1,] - xy[-i2,])
  Y <- positionIncrements / rootTimeIncrements
  # duplicatedLines <- (duplicated(t(designMatrix)) | apply(designMatrix, 2, function(x) all(x == 0)))
  # if(any(duplicatedLines)){
  #   designMatrix <- designMatrix[, !duplicatedLines]
  #   timeIncrements <- timeIncrements[!duplicatedLines]
  #   rootTimeIncrements <- rootTimeIncrements[!duplicatedLines]
  #   positionIncrements <- positionIncrements[!duplicatedLines]
  #   Y <- Y[!duplicatedLines]
  # }
  # estimation
  estimatedBetaCovariance <- estimatedCovariance <- solve (t(designMatrix * timeIncrements) %*% designMatrix)
  degreeFreedom <- length(Y) - J
  estimatedBeta <- estimatedNu <- estimatedCovariance %*% t(designMatrix) %*% positionIncrements;
  estimatedGamma <- NULL
  CIgamma <- NULL
  if(withspeed){
    estimatedGamma <- sum((Y -  sqrt(timeIncrements) * designMatrix %*% estimatedNu)^2) / degreeFreedom
    estimatedBeta <- estimatedNu / estimatedGamma * (degreeFreedom - 2) / degreeFreedom
    estimatedBetaCovariance <- (2 * estimatedBeta %*% t(estimatedBeta) / (degreeFreedom - 4) 
                                + estimatedCovariance / estimatedGamma * (1 + 2 / (degreeFreedom - 4)))
    
    CIgamma <- estimatedGamma * degreeFreedom/qchisq(c(0.975,0.025), degreeFreedom)
  }
  
  # Derive R squared
  RSS1 <- sum((Y -  sqrt(timeIncrements) * designMatrix %*% estimatedNu)^2)
  RSS0 <- sum(Y^2)
  R2 <- 1 - RSS1/RSS0
  
  confidenceIntervals <- t(sapply(1:length(estimatedBeta), function(j){
    estimatedBeta[j] + c(1, -1) * qnorm(0.025) *  sqrt(estimatedBetaCovariance[j, j])
  }))
  
  rownames(confidenceIntervals) <- rownames(estimatedBetaCovariance) <- colnames(estimatedBetaCovariance) <- paste0("beta", 1:J)
  return(list(betaHat = as.numeric(estimatedBeta), gammaHat  = estimatedGamma,
              betaHatCovariance = estimatedBetaCovariance, betaHat95CI = confidenceIntervals,
              gammaHat95CI = CIgamma, R2=R2))
}
# eulerLSE <- function(ID=NULL, xy, time, gradarray, withspeed=TRUE) 
# {
#     n <- nrow(xy)
#     if(is.null(ID))
#         ID <- rep(1,n)
#     
#     i0 <- which(ID[-1]!=ID[-n])
#     i1 <- c(1, i0+1) # first obs of each track
#     i2 <- c(i0, n) # last obs of each track
#     
#     # matrix of covariate derivatives
#     D <- 0.5 * rbind(gradarray[-i2,1,], gradarray[-i2,2,])
#     # vector of time intervals (requires less memory than diagonal matrix)
#     dt <- time[-i1] - time[-i2]
#     T <- rep(sqrt(dt),2)
#     TT <- rep(dt,2)
#     
#     # vector of 2-d steps
#     dxy <- xy[-i1,] - xy[-i2,]
#     Y <- matrix(dxy/sqrt(dt), ncol=1)
#     TY <- matrix(dxy, ncol=1)
#     
#     # estimation
#     DTTD <- t(D*TT) %*% D
#     DTTDinv <- solve(DTTD)
#     Bhat <- DTTDinv %*% t(D) %*% TY
#     
#     if(withspeed) {
#         # estimate speed
#         fitted <- t(T*t(D)) %*% Bhat
#         SSE <- sum((Y-fitted)^2)
#         speed <- SSE/(2*(n-1)-length(Bhat))        
# 
#         est <- c(Bhat/speed,speed)
#         var <- DTTDinv/speed
#     } else {
#         est <- Bhat
#         var <- DTTDinv
#     }
# 
#     return(list(est=est, var=var))
# }
