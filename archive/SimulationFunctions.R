rm(list = ls())
library(mixtools)
ScalarFun <- function(x, C, D, Om){
  exp(-(x - C)^2 / D) * sin(Om * (x -C))
}

ScalarFirstDer <- function(x, C, D, Om){
  -2 * (x - C) / D * ScalarFun(x, C, D, Om) + Om * exp(-(x - C)^2 / D) * cos(Om * (x -C))
}

ScalarSecDer <- function(x, C, D, Om){
  ((-2 / D - Om^2) * ScalarFun(x, C, D, Om) 
   -2 * (x - C) / D * (ScalarFirstDer(x, C, D, Om) +  Om * exp(-(x - C)^2 / D) * cos(Om * (x -C))))
}
# 
# Fun <- function(x, CenterParam, DilatationParam, IntensityParam, FrequenceParam){
#   MainTerms <- exp(-(x - CenterParam)^2 / DilatationParam) * sin(FrequenceParam * (x -CenterParam))
#   return(IntensityParam * prod(MainTerms))
# }

Fun <- function(x, Cs, Ds, Oms, I){
  if(!(all(c(length(x), length(Cs), length(Ds), length(Oms)) == 2))){
    stop("x, Cs, Ds, Oms must all be of length 2")
  }
  I * prod(ScalarFun(x, Cs, Ds, Oms))
}

GradFun <- function(x, Cs, Ds, Oms, I){
  if(!(all(c(length(x), length(Cs), length(Ds), length(Oms)) == 2))){
    stop("x, Cs, Ds, Oms must all be of length 2")
  }
  I * ScalarFun(x, Cs, Ds, Oms)[2:1] * ScalarFirstDer(x, Cs, Ds, Oms)
}

HessFun <- function(x, Cs, Ds, Oms, I){
  f <- ScalarFun(x, Cs, Ds, Oms)
  df <- ScalarFirstDer(x, Cs, Ds, Oms)
  ddf <- ScalarSecDer(x, Cs, Ds, Oms)
  I * c(f[2] * ddf[1], df[1] * df[2], df[2] * df[1], ddf[2] * f[1])
}
# UD parameters
# Beta1 <- -1;        Beta2 <- 0.5
# Alpha <- 0;
ScalC1 <- function(x) IP1 * prod(ScalarFun(x, C = CP1, DP1, FP1))
GradC1 <- function(x) GradFun(x, Cs = CP1, Ds = DP1, Oms = FP1, I = IP1)
ScalC2 <- function(x) IP2 * prod(ScalarFun(x, C = CP2, D = DP2, Om = FP2))
GradC2 <- function(x) GradFun(x, Cs = CP2, Ds = DP2, Oms = FP2, I = IP2)
GradDist <- function(x) - 2 * x
HessDist <- function(x) as.numeric(diag(- 2, 2))
HessC1 <- function(x) HessFun(x, Cs = CP1, Ds = DP1, Oms = FP1, I = IP1)
HessC2 <- function(x) HessFun(x, Cs = CP2, Ds = DP2, Oms = FP2, I = IP2)

UDMap <- function(x, Alpha, Beta1, Beta2, log = T){
  res <- - Alpha * sum(x^2) + Beta1 * ScalC1(x) + Beta2 * ScalC2(x)
  if(!log)
    res <- exp(res)
  return(res)
}



MeanStepEuler <- function(x, Delta, Alpha, Beta1, Beta2, gamma = 1){
  return(x + (0.5 * gamma * Delta * (Alpha * GradDist(x)
                             + Beta1 * GradC1(x) 
                             + Beta2 * GradC2(x))))
}
OzakiMoments <- function(x, Delta, Alpha, Beta1, Beta2, gamma = 1){
  HessVec <- 0.5 * gamma * (Alpha * HessDist(x) + Beta1 * HessC1(x) + Beta2 * HessC2(x))
  HessVecInv <- vecSolveM(HessVec)
  expJMinusId <- vecExpM(HessVec * Delta) - c(1, 0, 0, 1)
  GradVec <- 0.5 *  (Alpha * GradDist(x)
                     + Beta1 * GradC1(x) 
                     + Beta2 * GradC2(x))
  OzMean <- x + vecMtimesV(expJMinusId, vecMtimesV(HessVecInv, GradVec))
  JksJ <- MksumM(HessVec) # Kronecker sum
  JksJinv <- MksumMInv(HessVec) # Inverse Kronecker sum
  OzCov <- matrix(JksJinv %*% (expm::expm(JksJ * Delta) - diag(1, 4)) %*% c(gamma, 0, 0, gamma),
                  ncol = 2, nrow = 2)
  return(list(Mean = OzMean, Cov = OzCov))
} 


rNextObs <- function(x0, Delta, Alpha, Beta1, Beta2, gamma = 1, method = "Euler"){
  if(method == "Euler"){
    return(MeanStepEuler(x = x0, Delta = Delta, Alpha =  Alpha,
                         Beta1 = Beta1, Beta2 = Beta2, gamma = gamma)
           + rnorm(2, 0, sqrt(Delta * gamma))  )
  }
  if(method == "Ozaki"){
    OM <- OzakiMoments(x = x0, Delta = Delta, Alpha =  Alpha,
                       Beta1 = Beta1, Beta2 = Beta2, gamma = gamma) 
    return(rmvnorm(1, OM$Mean, OM$Cov))
  }
}
# 
# dNextObs <- function(y, x0, Delta, Alpha, Beta1, Beta2, gamma = 1, Log = F){
#   Mn <- x0 + MeanStep(x0, Delta, Alpha =  Alpha, Beta1 = Beta1, Beta2 = Beta2)
#   if(Log)
#     return(logdmvnorm(y, mu = Mn, sigma = diag(Delta, length(x0))))
#   else
#     return(dmvnorm(y, mu = Mn, sigma = diag(Delta, length(x0))))
# }

SimulationProcess <- function(Start, Npoints, Delta, A, B, C, gamma = 1, method = "Euler"){
  Times <- seq(0, (Npoints - 1) * Delta, by = Delta)
  OutputPositions = matrix(Start, ncol = 2, byrow = T, nrow = Npoints)
  for(i in 2:Npoints){
    OutputPositions[i, ] <- rNextObs(OutputPositions[i - 1, ], 
                                     Delta,Alpha =  A, Beta1 =  B, Beta2 =  C,
                                     gamma = gamma,
                                     method = method)
  }
  Output <- cbind(OutputPositions, Times)
  colnames(Output) <- c("X1", "X2", "Time")
  return(Output)
} 
