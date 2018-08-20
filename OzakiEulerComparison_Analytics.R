rm(list = ls())
source("archive/SimulationFunctions.R")
source("archive/CovariateParameters.R")
source("nllkLang.R")
source("eulerLSE.R")
source("utility.R")
library(fields)
library(numDeriv)

DATAPATH <- "simulatedData/analyticCase/"

doubleEstimation <- function(seed, covArray, dataSetLength){
  set.seed(seed)
  CompleteSample <- read.table(paste0(DATAPATH, "simAnalyticData_seed", seed, ".txt"),
                               header = T, sep = ";")
  Sel <- sort(sample(floor(seq(1, floor(nrow(CompleteSample)), by = 10)), 
                    replace = F, size = dataSetLength))
  xy <- as.matrix(CompleteSample[Sel, c("X1", "X2")])
  time <- CompleteSample[Sel, "Time"]
  #True array of covariates gradient
  trueGradArray <- array(dim = c(length(Sel), 2,  3))
  trueGradArray[,, 1] <- t(apply(xy, 1, GradC1)) ; trueGradArray[,,2] <- t(apply(xy, 1, GradC2))
  trueGradArray[,, 3] <- t(apply(xy, 1, GradDist))
  #True array of covariates hessian
  trueJacobArray <- array(dim = c(length(Sel), 4,  3))
  trueJacobArray[,, 1] <- t(apply(xy, 1, HessC1)) ; trueJacobArray[,,2] <- t(apply(xy, 1, HessC2))
  trueJacobArray[,, 3] <- t(apply(xy, 1, HessDist))
  #Interpolated covArray
  interpGradArray <- array(dim = c(length(Sel), 2,  3))
  interpGradArray[,, 3] <- trueGradArray[,, 3]#The gradient for distance can always be computed exactly
  interpGradArray[,, 1:2] <- covGrad(xy, xGrid, yGrid, covArray)
  #Interpolated covArray
  interpJacobArray <- array(dim = c(length(Sel), 4,  3))
  interpJacobArray[,, 3] <- trueJacobArray[,, 3]#The gradient for distance can always be computed exactly
  interpJacobArray[,, 1:2] <- covHessian(xy, xGrid, yGrid, covArray)
  objFunTrueArrays <- function(ParamVector){
    # nllkLang(beta = ParamVector, xy = xy, time = time, 
    #          gradarray = trueGradArray, hessarray = trueJacobArray, method = "ozaki")
    nllkLangForOptim(vecPars = ParamVector, xy = xy, time = time, 
             gradarray = trueGradArray, hessarray = trueJacobArray, method = "ozaki")
  }
  objFunInterpArrays <- function(ParamVector){
    nllkLang(par = ParamVector, xy = xy, time = time, 
             gradarray = interpGradArray, hessarray = interpJacobArray, method = "ozaki")
  }
  fitEulerTrueGrad <- eulerLSE(xy = xy, time = time, gradarray = trueGradArray)
  fitEulerInterpGrad <- eulerLSE(xy = xy, time = time, gradarray = interpGradArray)
  fitOzTrueHess <- nlminb(start = c(0, 0, 0, 1), objective = objFunTrueArrays, 
                 control = list(trace = 0, x.tol = 10^(-4)))
  fitOzInterpHess <- nlminb(start = c(0, 0, 0), objective = objFunInterpArrays, 
                          control = list(trace = 1, x.tol = 10^(-4)))
  result <- list(fitEulTG = fitEulerTrueGrad, fitEulIG = fitEulerInterpGrad,
              fitOzTG = fitOzTrueHess, fitOzIG = fitOzInterpHess,
              selectedLines = Sel)
  save(result, file = paste0(DATAPATH, "tmp", seed, ".RData"))
  return(result)
}
library(parallel)
trueParams <- c(B1 = -1, B2 = 0.5, A = 0.05)
dataSetLength <- 300
xGrid <- yGrid <- seq(-15, 15, length = 8)
covariateSampling <- as.matrix(expand.grid(xGrid, yGrid))
covariateArray <- array(dim = c(length(xGrid), length(yGrid),  2))
covariateArray[,, 1] <- matrix(apply(covariateSampling, 1, ScalC1), nrow = length(xGrid)); 
covariateArray[,, 2] <- matrix(apply(covariateSampling, 1, ScalC2), nrow = length(xGrid))

T1 <- Sys.time()
AllResults <- mclapply(1:3, function(i) doubleEstimation(i, covArray = covariateArray, dataSetLength = dataSetLength), mc.cores = 3)
T2 <- Sys.time()
save(AllResults, file = "Results_OzakiEulerComparison_AnalyticsAndGradient.RData")
