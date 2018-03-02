#' Functions for basic algebra when a 4-vector codes a 2 by 2 matrix
#' 
#' @param vec1,vec2,vecM vectors of length 4 vectors, for matrix M 2 by 2, vectors are coded
#' as expected by R. I.E, a matrix 
#' a11 a12
#' a21 a22 is coded as the vector (a11, a21, a12, a22)
#' @param vec2D a vector of length 2
#' @return vectors of length 4 corresponing to the linear algebra computation, 
#' except for the MksM and MksM that returns the matrix correspnding to the inner kronecker sum
#' and its inverse respectively
vecMtimesM <- function(vec1, vec2){
  return(c(vec1[1] * vec2[1] + vec1[3] * vec2[2],
           vec1[2] * vec2[1] + vec1[4] * vec2[2],
           vec1[1] * vec2[3] + vec1[3] * vec2[4],
           vec1[2] * vec2[3] + vec1[4] * vec2[4]))
}

vecDetM <- function(vecM){
  vecM[1] * vecM[4] - vecM[2] * vecM[3]
}

vecSolveM <- function(vecM){
  Det <- vecM[1] * vecM[4] - vecM[2] * vecM[3]#determinant
  if(is.na(Det) | Det == 0)
    return(rep(NA, 4))
  return(c(vecM[4], -vecM[2], -vecM[3], vecM[1]) / Det)
}

vecMtimesV <- function(vecM, vec2D){
  return(c(vecM[1] * vec2D[1] + vecM[3] * vec2D[2],
           vecM[2] * vec2D[1] + vecM[4] * vec2D[2]))
}
vecExpM <- function(vecM){
  return(as.numeric(expm::expm(matrix(vecM, nrow = 2))))
}
MksumM <- function(vecM){# Inner Kronecker sum
  a <- vecM[1]; b <- vecM[2]; c <- vecM[3]; d <- vecM[4]
  VecMat <- c(2 * a, c    , c    , 0,
              b    , a + d, 0    , c,
              b    , 0    , a + d, c,
              0    , b    , b    , 2 * d )
  matrix(VecMat, nrow = 4, ncol = 4, byrow = T)
}

MksumMInv <- function(vecM){# Inverse of the inner Kronecker sum
  a <- vecM[1]; b <- vecM[2]; c <- vecM[3]; d <- vecM[4]
  Det <- 2 * (a + d) * (a*d - b*c)
  if(Det == 0){
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
#' @param vec2ds matrix with 2 columns and, say, n lines
#' as 
#' @param vecMs  matrix with 4 columns and n line(same number of lines as x)
#' @return For each row v of vecDs and M of vecMs, return t(v) times matrix(M) times v
vecVMV <- function(vec2Ds, vecMs){
  return(vecMs[, 1] * vec2Ds[, 1] * vec2Ds[, 1] 
         + vec2Ds[, 1] * vec2Ds[, 2] * (vecMs[, 2] + vecMs[, 3])
         + vecMs[ ,4] * vec2Ds[, 2] * vec2Ds[, 2])
}

#' Fast function for log density of a bivariate normal distribution
#' 
#' @param xs matrix with 2 columns and, say, n lines
#' as 
#' @param vecMeans  matrix of same dimension as x, corresponding to evaluated means
#' @param vecCovInvs matrix of n lines and 4 columns, vectorization of the PRECISION matrix 
#' (inverse of covariance)
#' @param DetInv number or vector of length n giving the determinant of each matrix of vecCovInvs,
#' can be omitted, it will be computed in this case
#' @return For each row x xs, m of vecMeans, C of vecCovInvs, compute the pdf (evaluated in x) of a 
#' bivariate normal N(m, C^(-1))
logdnorm2d <- function(xs, vecMeans, vecCovInvs, DetInv = NULL){
  if(is.null(DetInv)){
    DetInv <- apply(vecCovInvs, 1, vecDetM)
  }
  - log(2 * pi) + 0.5 * log(DetInv) - 0.5 * vecVMV(xs - vecMeans, vecCovInvs)
}


#' Spot abnormal values of determinants to avoid numerical problems in the Ozaki method
#' 
#' @param Dets vector of determinants
#' @return A vector of booleans. False is the determinant is acceptable, True if is problematic
NumProblemDet <- function(Dets){
  sapply(Dets, function(x) isTRUE(x < 10^(-5) | x > 10^(10)))
}
