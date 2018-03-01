source("LinearAlgebraFunctions.R")
#' Obtaining Ozaki required mean
#' 
#' At the moment, this function is based on the Euler approximation.
#' 
#' @param xy Matrix of observed locations (two columns: x, y)
#' @param time Vector of times of observations
#' @param gradmat Computed gradient at xy positions. (In the case of the RSF, it
#' is already the sum of all gradients). It is a T by 2 matrix where T is the number 
#' of observations for the given trajectory
#' @param jacobmat Computed Jacobian at xy positions. It must be T by 4 matrix
#' in the following order (d/dxdx, d/dxdy, d/dydx, d/dydy)
#' @return A matrix having size of xy, corresponding of the expected mean 
#' of the ozaki linearization, which is the muO of the article of Gloaguen et al
#' 
getOzakiMean <- function(xy, Deltas, gradmat, jacobmat) {
  jacobmatInv <- t(apply(jacobmat, 1, vecSolveM))  #For each row, computing the inverse row
  Idvec <- c(1, 0, 0, 1)
  return(t(sapply(1:nrow(xy), function(i){
    expJiMinusId <- (vecExpM(jacobmat[i, ] * Deltas[i]) - IdVec)
    xy[i, ] + vecMtimesV(expJiMinusId, vecMtimesM(jacobmatInv[i, ], gradmat[i, ]))
  })))
}

getOzakiCov <- function(xy, Deltas, gradmat, jacobmat) {
  JksJ <- t(apply(jacobmat, 1, vecSolveM))  #For each row, computing the inverse row
  Idvec <- c(1, 0, 0, 1)
  I4 <- diag(1, 4)
  return(t(sapply(1:nrow(xy), function(i){
    JksJ <- MksumM(jacobmat[i, ])
    JksJinv <- MksumMInv(jacobmat[i, ])
    as.numeric(JksJinv %*% (expm::expm(JksJ * Deltas[i]) - I4) %*% Idvec)
  })))
}

