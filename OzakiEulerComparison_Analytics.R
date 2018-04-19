rm(list = ls())
source("archive/SimulationFunctions.R")
source("archive/CovariateParameters.R")
source("nllkLang.R")
A <- 0.05; B1 <- -1; B2 <- 0.5
DataSetLength <- 500
SimulEstim <- function(seed){
  set.seed(seed)
  CompleteSample <- SimulationProcess(Start = c(0, 0), Npoints = 5 * 10^3, 
                                      Delta = 0.01, A, B1, B2, method = "Ozaki")
  Sel <- sort(sample(floor(seq(1, nrow(CompleteSample), by = 10)), 
                    replace = F, size = DataSetLength))
  xy <- as.matrix(CompleteSample[Sel, c("X1", "X2")])
  time <- CompleteSample[Sel, "Time"]
  #Array of covariates gradient
  gradarray <- array(dim = c(length(Sel), 2,  3))
  gradarray[,,1] <- t(apply(xy, 1, GradC1)) ; gradarray[,,2] <- t(apply(xy, 1, GradC2))
  gradarray[,,3] <- t(apply(xy, 1, GradDist))
  #Array of covariates hessian
  jacobarray <- array(dim = c(length(Sel), 4,  3))
  jacobarray[,,1] <- t(apply(xy, 1, HessC1)) ; jacobarray[,,2] <- t(apply(xy, 1, HessC2))
  jacobarray[,,3] <- t(apply(xy, 1, HessDist))
  CFE <- function(ParamVector){
    nllkLang(beta = ParamVector, xy = xy, time = time, 
             gradarray = gradarray, jacobarray)
  }
  CFO <- function(ParamVector){
    nllkLang(beta = ParamVector, xy = xy, time = time, 
             gradarray = gradarray, jacobarray, method = "ozaki")
  }
  fitE <- nlminb(start=c(0, 0, 0), objective = CFE,
                 control=list(trace = 0))
  fitO <- nlminb(start = c(0, 0, 0), objective = CFO, 
                 control = list(trace = 0, x.tol = 10^(-4)))
  return(list(Data = CompleteSample[Sel, ], fitE = fitE, fitO = fitO))
}
library(parallel)
AllResults <- mclapply(1:30, SimulEstim, mc.cores = 3)
save(AllResults, file = "Results_OzakiEulerComparison_Analytics.RData")

