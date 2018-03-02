source("archive/SimulationFunctions.R")
source("archive/CovariateParameters.R")
source("nllkLang.R")
CompleteSample <- as.matrix(read.table("archive/SimulatedSample.txt", header = T))
#True Beta = c(-1, 0.5)
set.seed(123)
Npoints = 2000
Sel <- seq(1, nrow(CompleteSample), length.out = 1000)
xy <- as.matrix(CompleteSample[Sel, c("X1", "X2")])
time <- CompleteSample[Sel, "Time"]
#Array of covariates gradient
gradarray <- array(dim = c(length(Sel), 2,  2))
gradarray[,,1] <- t(apply(xy, 1, GradC1))
gradarray[,,2] <- t(apply(xy, 1, GradC2))
#Array of covariates hessian
jacobarray <- array(dim = c(length(Sel), 4,  2))
jacobarray[,,1] <- t(apply(xy, 1, HessC1))
jacobarray[,,2] <- t(apply(xy, 1, HessC2))

CFE <- function(ParamVector){
  nllkLang(beta = ParamVector, xy = xy, time = time, 
           gradarray = gradarray, jacobarray)
}
CFO <- function(ParamVector){
  nllkLang(beta = ParamVector, xy = xy, time = time, 
           gradarray = gradarray, jacobarray, method = "ozaki")
}
# fit <- nlminb(start=c(-1, 0.5, 0.05),objective=ContrastFunction,
#               control=list(trace=1))
fitE <- nlminb(start=c(0,0),objective=CFE,
              control=list(trace=1))
fitO <- nlminb(start=c(0,0),objective=CFO, control=list(trace=1, x.tol = 10^(-4)))
