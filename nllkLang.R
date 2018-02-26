
#' Langevin movement model negative log-likelihood function
#' 
#' @param beta Parameters of the utilisation distribution
#' @param xy Matrix of observed locations (two columns: x, y)
#' @param gradarray Three-dimensional array of gradients of covariate fields. The rows index time, 
#' the columns are the dimensions (x and y), and the layers index the covariates.
#' 
#' @return Negative log-likelihood
nllkLang <- function(beta,xy,gradarray,dt) {
    gradmat <- apply(gradarray, 2, function(mat) mat%*%beta)
    llk <- sum(dnorm(xy[-1,],xy[-nrow(xy),]+dt*gradmat[-nrow(xy),],sqrt(dt),log=TRUE))
    return(-llk)
}