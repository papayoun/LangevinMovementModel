
# Functions for basic algebra when a 4-vector codes a 2 by 2 matrix

# Matrix multiplication
vecMtimesM <- function(vec1, vec2){
    return(c(vec1[1] * vec2[1] + vec1[3] * vec2[2],
             vec1[2] * vec2[1] + vec1[4] * vec2[2],
             vec1[1] * vec2[3] + vec1[3] * vec2[4],
             vec1[2] * vec2[3] + vec1[4] * vec2[4]))
}

# Determinant
vecDetM <- function(vecM){
    vecM[1] * vecM[4] - vecM[2] * vecM[3]
}

# Inverse
vecSolveM <- function(vecM){
    Det <- vecM[1] * vecM[4] - vecM[2] * vecM[3] # determinant
    if(is.na(Det) | Det == 0)
        return(rep(NA, 4))
    return(c(vecM[4], -vecM[2], -vecM[3], vecM[1]) / Det)
}

# Matrix-vector multiplication
vecMtimesV <- function(vecM, vec2D){
    return(c(vecM[1] * vec2D[1] + vecM[3] * vec2D[2],
             vecM[2] * vec2D[1] + vecM[4] * vec2D[2]))
}

# Matrix exponential
vecExpM <- function(vecM){
    return(as.numeric(expm::expm(matrix(vecM, nrow = 2))))
}

# Inner Kronecker sum (output is a matrix)
MksumM <- function(vecM){
    a <- vecM[1]; b <- vecM[2]; c <- vecM[3]; d <- vecM[4]
    VecMat <- c(2 * a, c    , c    , 0,
                b    , a + d, 0    , c,
                b    , 0    , a + d, c,
                0    , b    , b    , 2 * d )
    matrix(VecMat, nrow = 4, ncol = 4, byrow = T)
}

# Inverse of inner Kronecker sum (output is a matrix)
MksumMInv <- function(vecM){
    a <- vecM[1]; b <- vecM[2]; c <- vecM[3]; d <- vecM[4]
    Det <- 2 * (a + d) * (a*d - b*c)
    if(is.na(Det)){
      return(matrix(NA, nrow = 4, ncol = 4))
    }
    else if(Det == 0){
        return(matrix(NA, nrow = 4, ncol = 4))
    }
    VecMat <- c(d*d + a*d - b*c, -c*d       , -c*d       , c*c          ,
                -b*d           , 2*a*d - b*c, b*c        , -a*c         ,
                -b*d           , b*c        , 2*a*d - b*c, -a*c         ,
                b*b            , -a*b       , -a*b       ,a*a + d*a - b*c)
    matrix(VecMat / Det, nrow = 4, ncol = 4, byrow = T)
}

#' Fast function for computing bivariate quadratic form
#' 
#' Calculates t(v) %*% M %*% v for multiple input vectors v and matrices M.
#' 
#' @param allv Matrix with two columns. Each row corresponds to a 2-vector v.
#' @param allM Matrix with 4 columns, and as many rows as allv. Each row corresponds
#' to a 2x2 matrix M (stored as a 4-vector).
#' 
#' @return Vector of length the number of rows of the input matrices.
#' For each row v of allv and M of allM, returns t(v) times matrix(M) times v.
vecVMV <- function(allv, allM){
    return(allM[, 1] * allv[, 1] * allv[, 1] 
           + allv[, 1] * allv[, 2] * (allM[, 2] + allM[, 3])
           + allM[ ,4] * allv[, 2] * allv[, 2])
}

#' Fast function for log density of a bivariate normal distribution
#' 
#' @param xs Matrix with 2 columns and, say, n lines
#' @param vecMeans Matrix of same dimension as xs, corresponding to evaluated means
#' @param vecCovInvs Matrix of n lines and 4 columns, vectorization of the PRECISION matrix 
#' (inverse of covariance)
#' @param DetInv Number or vector of length n giving the determinant of each matrix of vecCovInvs.
#' If NULL (default), it will be computed.
#' 
#' @return For each row x xs, m of vecMeans, C of vecCovInvs, compute the pdf (evaluated in x) of a 
#' bivariate normal N(m, C^(-1))
logdnorm2d <- function(xs, vecMeans, vecCovInvs, DetInv = NULL){
    if(is.null(DetInv)){
        DetInv <- apply(vecCovInvs, 1, vecDetM)
    }
    return(- log(2 * pi) + 0.5 * log(DetInv) - 0.5 * vecVMV(xs - vecMeans, vecCovInvs))
}


#' Spot abnormal values of determinants to avoid numerical problems in the Ozaki method
#' 
#' @param Dets vector of determinants
#' 
#' @return A vector of booleans. FALSE if the determinant is acceptable, TRUE otherwise.
NumProblemDet <- function(Dets){
    sapply(Dets, function(x) isTRUE(x < 10^(-5) | x > 10^(10)))
}
