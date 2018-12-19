
#' Evaluate the log-RSF based on interpolated covariates
logRSFinterp <- function(xy, beta, xgrid, ygrid, covarray) {
    c <- apply(covarray, 3, function(covmat) 
        interpCov(xy=xy, xgrid=xgrid, ygrid=ygrid, covmat=covmat))
    return(sum(beta*c))
}

#' Simulate from the Langevin movement model
#' 
#' This function is based on the Euler approximation.
#' 
#' @param beta Vector of resource selection coefficients
#' @param speed Scalar speed parameter
#' @param time Vector of times of observations
#' @param xy0 Vector of coordinates (x,y) of initial location
#' @param xgrid Grid on which the covariates are known
#' @param ygrid Grid on which the covariates are known
#' @param covarray Array of values of the covariates at the points given by
#' xgrid and ygrid, of dimensions (length(xgrid),length(ygrid),length(beta)).
simLangMH <- function(beta, speed=1, time, xy0, xgrid, ygrid, covarray) {
    nbObs <- length(time)
    dt <- diff(time)
    
    # Initialise location
    xy <- matrix(0,nbObs,2)
    xy[1,] <- xy0
    
    # Log-RSF and gradient of log-UD at initial position
    g <- gradLogUD(beta=beta,xy=xy[1,],xgrid=xgrid,ygrid=ygrid,covarray=covarray)
    logRSF <- logRSFinterp(xy[1,], beta, xgrid, ygrid, covarray=covarray)
    
    # Count acceptances and rejections
    acc <- 0
    rej <- 0
    for(t in 2:nbObs) {
        if(t%%100==0) {
            cat(paste("\rSimulating Langevin process...",round(100*t/nbObs),
                      "% -- acc =",round(acc/(rej+acc)*100),"%"))  
        }
        
        # Propose new location
        xyprime <- xy[t-1,] + 0.5 * g * speed * dt[t-1] + rnorm(2, 0, sqrt(speed*dt[t-1]))
        
        # Log-RSF and gradient of log-UD at proposed location
        gprime <- gradLogUD(beta=beta,xy=xyprime,xgrid=xgrid,ygrid=ygrid,covarray=covarray)
        logRSFprime <- logRSFinterp(xyprime, beta, xgrid, ygrid, covarray=covarray)
        
        # Log-proposals in both directions, for the acceptance ratio
        logProp <- sum(dnorm(xy[t-1,], xyprime + speed*dt[t-1]*gprime/2, 
                             sqrt(speed*dt[t-1]), log=TRUE))
        logPropPrime <- sum(dnorm(xyprime, xy[t-1,] + speed*dt[t-1]*g/2, 
                                  sqrt(speed*dt[t-1]), log=TRUE))
        
        # Log acceptance ratio
        logAR <- logRSFprime + logProp - logRSF - logPropPrime
        if(log(runif(1))<logAR) {
            # If accepted, update location, gradient, and log-RSF
            xy[t,] <- xyprime
            g <- gprime
            logRSF <- logRSFprime
            acc <- acc+1
        } else {
            xy[t,] <- xy[t-1,]
            rej <- rej + 1
        } 
    }
    cat("\n")
    return(list(xy=xy, acc=acc, rej=rej))
}
