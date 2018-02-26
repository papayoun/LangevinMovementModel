rm(list = ls())
library(mixtools)
Fun <- function(x, CenterParam, DilatationParam, IntensityParam, FrequenceParam){
  MainTerms <- exp(-(x - CenterParam)^2 / DilatationParam) * sin(FrequenceParam * (x -CenterParam))
  return(IntensityParam * prod(MainTerms))
}

GradFun <- function(x, CenterParam, DilatationParam, IntensityParam, FrequenceParam){
  MainTerms <- exp(-(x - CenterParam)^2 / DilatationParam) * sin(FrequenceParam * (x -CenterParam) )
  DerTerms <- (- 2 * (x - CenterParam) / DilatationParam * MainTerms# Derivation of exp 
                + FrequenceParam * exp(-(x - CenterParam)^2 / DilatationParam) * cos(FrequenceParam * (x -CenterParam)))
  IntensityParam * sapply(1:length(MainTerms), function(i){
     prod(MainTerms[-i]) * DerTerms[i]
  })
}

# Covariates parameters
IP1 <- 6;           IP2 <- 6
DP1 <- c(40, 40);   DP2 <- c(40, 40)
CP1 <- c(0, 0);     CP2 <- c(-2, pi/2)
FP1 <- c(0.6, 0.2); FP2 <- c(0.1, 0.5)

# UD parameters
Beta1 <- -1;        Beta2 <- 0.5
Alpha <- 1 / 20
MeanStep <- function(x, Delta, Alpha, Beta1, Beta2){
  (0.5 * Delta * (-2 * x * Alpha 
                  + Beta1 * GradFun(x, CP1, DP1, IP1, FP1) 
                  + Beta2 * GradFun(x, CP2, DP2, IP2, FP2)))
}
GradC1 <- function(x) GradFun(x, CP1, DP1, IP1, FP1)
GradC2 <- function(x) GradFun(x, CP2, DP2, IP2, FP2)
rNextObs <- function(x0, Delta, Alpha, Beta1, Beta2){
  x0 + MeanStep(x0, Delta, Alpha, Beta1, Beta2) + rnorm(2, 0, sqrt(Delta))  
}

dNextObs <- function(y, x0, Delta, Alpha, Beta1, Beta2, Log = F){
  Mn <- x0 + MeanStep(x0, Delta, Alpha, Beta1, Beta2)
  if(Log)
    return(logdmvnorm(y, mu = Mn, sigma = diag(Delta, length(x0))))
  else
    return(dmvnorm(y, mu = Mn, sigma = diag(Delta, length(x0))))
}

SimulationProcess <- function(Start, Npoints, Delta){
  Times <- seq(0, (Npoints - 1) * Delta, by = Delta)
  OutputPositions = matrix(Start, ncol = 2, byrow = T, nrow = Npoints)
  for(i in 2:Npoints){
    OutputPositions[i, ] <- rNextObs(OutputPositions[i - 1, ], 
                                     Delta, Alpha, Beta1, Beta2)
  }
  Output <- cbind(OutputPositions, Times)
  colnames(Output) <- c("X1", "X2", "Time")
  return(Output)
} 
