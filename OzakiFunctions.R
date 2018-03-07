source("LinearAlgebraFunctions.R")
#' Mean of Ozaki transition density
#' 
#' @param xy Matrix of observed locations (two columns: x, y)
#' @param Deltas Vector of time intervals
#' @param gradmat Computed gradient at xy positions. (In the case of the RSF, it
#' is already the sum of all gradients). It is a T by 2 matrix where T is the number 
#' of observations for the given trajectory
#' @param jacobmat Computed Jacobian at xy positions. It must be T by 4 matrix
#' in the following order (d/dxdx, d/dxdy, d/dydx, d/dydy)
#' 
#' @return A matrix having size of xy, corresponding to the expected mean 
#' of the Ozaki linearization, denoted muO in Gloaguen et al.
#' 
getOzakiMean <- function(xy, Deltas, gradmat, jacobmat, jacobmatInv) {
    jacobmatInv <- t(apply(jacobmat, 1, vecSolveM))  # For each row, compute the inverse row
    Idvec <- c(1, 0, 0, 1) # I2 in vectorized version
    return(t(sapply(1:(nrow(xy) - 1), function(i){
        if(any(is.na(jacobmatInv[i, ]))){ # If numerical problem, returns the euler approx
            return(xy[i, ] + Deltas[i] * gradmat[i, ])
        }
        else{
            expJiMinusId <- (vecExpM(jacobmat[i, ] * Deltas[i]) - Idvec)
            return(xy[i, ] + vecMtimesV(expJiMinusId, vecMtimesV(jacobmatInv[i, ], gradmat[i, ])))
        }
    })))
}

#' Covariance matrix of Ozaki transition density
#' 
#' @param xy Matrix of observed locations (two columns: x, y)
#' @param Deltas Vector of time intervals
#' @param gradmat Computed gradient at xy positions. (In the case of the RSF, it
#' is already the sum of all gradients). It is a T by 2 matrix where T is the number 
#' of observations for the given trajectory
#' @param jacobmat Computed Jacobian at xy positions. It must be T by 4 matrix
#' in the following order (d/dxdx, d/dxdy, d/dydx, d/dydy)
getOzakiCovariance <- function(xy, Deltas, gradmat, jacobmat, Inv = F) {
    Idvec <- c(1, 0, 0, 1)
    I4 <- diag(1, 4)
    return(t(sapply(1:(nrow(xy) - 1), function(i){
        if(all.equal(jacobmat[i, 1], 0) == T | all.equal(jacobmat[i, 4], 0) == T) # All betas equal to 0
            vecCov <- c(Deltas[i], 0, 0, Deltas[i]) # Euler approximation
        else{
            JksJ <- MksumM(jacobmat[i, ]) # Kronecker sum
            JksJinv <- MksumMInv(jacobmat[i, ]) # Inverse Kronecker sum
            if(any(is.na(JksJinv))){ # Case of numerical problem
                vecCov <- c(Deltas[i], 0, 0, Deltas[i]) # Euler approximation
            }
            else{
                vecCov <- as.numeric(JksJinv %*% (expm::expm(JksJ * Deltas[i]) - I4) %*% Idvec)
                if(mode(all.equal(vecCov[2], vecCov[3], tolerance = 10^(-4))) == "character"){
                    print("Line", i)
                    print("jacobmat[i,]")
                    print(jacobmat[i, ])
                    print("Non symmetric matrix")
                }
            }
        }
        if(Inv){
            return(vecSolveM(vecCov)) # directly to inverse
        }
        else{
            return(vecCov)
        }
    })))
}

