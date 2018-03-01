rm(list = ls())
#' Functions for basic algebra when a 4-vector codes a 2 by 2 matrix
#' 
#' @param vec,vec1,vec2,vecM vectors of length 4 vectors, for matrix M 2 by 2, vectors are coded
#' as 
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

vecSolveM <- function(vecM){
  Det <- vecM[1] * vecM[4] - vecM[2] * vecM[3]#determinant
  return(c(vecM[4], -vecM[2], -vecM[3], vecM[1]) / Det)
}

vecMtimesV <- function(vecM, vec2D){
  return(c(vecM[1] * vec2D[1] + vecM[3] * vec2D[2],
           vecM[2] * vec2D[1] + vecM[4] * vec2D[2]))
}
vecExpM <- function(vecM){
  return(as.numeric(expm::expm(matrix(vecM, nrow = 2))))
}

MksumM <- function(vecM){
  a <- vecM[1]; b <- vecM[2]; c <- vecM[3]; d <- vecM[4]
  VecMat <- c(2 * a, c    , c    , 0,
              b    , a + d, 0    , c,
              b    , 0    , a + d, c,
              0    , b    , b    , 2 * d )
  matrix(VecMat, nrow = 4, ncol = 4, byrow = T)
}

MksumMInv <- function(vecM){
  a <- vecM[1]; b <- vecM[2]; c <- vecM[3]; d <- vecM[4]
  Det <- 2 * (a + d) * (a*d - b*c)
  VecMat <- c(d*d + a*d - b*c, -c*d       , -c*d       , c*c          ,
              -b*d           , 2*a*d - b*c, b*c        , -a*c         ,
              -b*d           , b*c        , 2*a*d - b*c, -a*c         ,
              b*b            , -a*b       , -a*b       ,a*a + d*a - b*c)
  matrix(VecMat / Det, nrow = 4, ncol = 4, byrow = T)
}

