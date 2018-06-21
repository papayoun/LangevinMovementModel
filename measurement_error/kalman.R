
##' Kalman filter for Langevin diffusion (Euler scheme)
##' 
##' @param par Vector of parameters (betas, speed, and observation error)
##' @param xy Matrix of observed locations
##' @param time Vector of times of locations
##' @param ID Vector of track IDs
##' @param gradarray Array of gradients of covariate fields, evaluated at
##' observed locations
##' 
##' @return Negative log-likelihood
kalman <- function(par, xy, time, ID, gradarray)
{
    nbObs <- nrow(xy)
    dt <- c(diff(time),1)
    
    # unpack parameters
    nbCovs <- dim(gradarray)[3]
    beta <- par[1:nbCovs]
    speed <- exp(par[nbCovs+1])
    h <- exp(par[nbCovs+2]) # measurement error SD
    
    # multiply gradients by beta coefficients
    gradmat <- 0.5 * apply(gradarray, 2, function(mat) mat %*% beta)
    
    # covariance matrix of observation process
    H <- diag(h^2,2)
    
    # initial state mean and covariance
    aest <- matrix(xy[1,] + speed * dt[1] * gradmat[1,], ncol=1)
    Pest <- H
    
    llk <- 0
    for(i in 2:nbObs) {
        if(ID[i]!=ID[i-1]) {
            aest <- matrix(xy[i,] + speed * dt[i] * gradmat[i,], ncol=1)
            Pest <- H
        } else {
            Bc <- matrix(speed * dt[i] * gradmat[i,], ncol=1)
            Q <-  speed * dt[i] * diag(2)
            
            # measurement residual
            u <- matrix(xy[i,],ncol=1) - aest
            
            # residual covariance
            F <- Pest + H
            
            # add contribution to log-likelihood
            detF <- F[1,1]*F[2,2] - F[2,1] * F[1,2]
            llk <- llk - (log(detF) + sum(u*solve(F,u)))/2
            
            # Kalman gain
            K <- Pest %*% solve(F)
            
            # update state estimate and state covariance
            aest <- aest + K %*% u + Bc
            Pest <- Pest %*% t(diag(2)-K) + Q
        }
    }
    
    return(-llk)
}
