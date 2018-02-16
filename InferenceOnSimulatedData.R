source("ExampleSimulationFunctionsAndParameters.R")
CompleteSample <- as.matrix(read.table("SimulatedSample.txt", header = T))

OptimContrast <- function(Param0, MyData, maxiter = 10,
                          PosCols = c("X1", "X2"), TimeCol = "Time"){
  n <- nrow(MyData)
  Positions <- MyData[, PosCols]
  Times <- MyData[, TimeCol]
  Deltas <- (Times[-1] - Times[-n])
  ContrastFunction <- function(ParamVector){
    sum(sapply(1:(n-1), function(i){
      dNextObs(y = Positions[i + 1, ], x0 =  Positions[i, ], 
               Delta = Times[i + 1] - Times[i], Alpha = ParamVector[1], 
               Beta1 = ParamVector[2], Beta2 = ParamVector[3], Log = T)
    }))
  }
  GradC1s <-  t(apply(Positions[-n, ], 1, GradC1))
  GradC2s <-  t(apply(Positions[-n, ], 1, GradC2))
  Res <- matrix(NA, ncol = length(Param0), nrow = maxiter + 1)
  Contrast <- rep(NA, maxiter + 1)
  Res[1, ] <- Param0
  Contrast[1] <- ContrastFunction(Res[1, ])
  colnames(Res) <- c("Alpha", "Beta1", "Beta2")
  DenomAlpha <- sum(rowSums(Positions[-n, ]^2) * Deltas)
  Dxis <- Positions[-1, ] - Positions[-n, ]
  # Function to update Alpha
  Num2Alpha <- sum(rowSums(Positions[-n, ] * Dxis))# Do not change
  GetAlpha <- function(Beta1, Beta2){
    CovGradTerm <- Beta1 * GradC1s + Beta2 * GradC2s
    Num1 <- 0.5 * sum(Deltas * rowSums((Positions[-n, ] * CovGradTerm)))
    (Num1 - Num2Alpha) / DenomAlpha
  }
  #Function to update Beta1
  DenomBeta1 <- 0.5 * sum(rowSums(GradC1s^2) * Deltas)# Do not change
  GetBeta1 <- function(Alpha, Beta2){
    NumTerm <- (Dxis
                + 0.5 * Deltas * (2 * Alpha * Positions[-n, ] - Beta2 * GradC2s))
    Num <- sum(rowSums(GradC1s * NumTerm))
    Num / DenomBeta1
  }
  #Function to update Beta2
  DenomBeta2 <- 0.5 * sum(rowSums(GradC2s^2) * Deltas)# Do not change
  GetBeta2 <- function(Alpha, Beta1){
    NumTerm <- (Dxis
                + 0.5 * Deltas * (2 * Alpha * Positions[-n, ] - Beta1 * GradC1s))
    Num <- sum(rowSums(GradC2s * NumTerm))
    Num / DenomBeta2
  }
  for(i in 1:maxiter){
    Res[i + 1, "Alpha"] <- GetAlpha(Beta1 = Res[i, "Beta1"], 
                                    Beta2 = Res[i, "Beta2"])
    Res[i + 1, "Beta1"] <- GetBeta1(Alpha = Res[i + 1, "Alpha"], 
                                    Beta2 = Res[i, "Beta2"])
    Res[i + 1, "Beta2"] <- GetBeta2(Alpha = Res[i + 1, "Alpha"], 
                                    Beta1 =  Res[i + 1, "Beta1"])
    Contrast[i + 1] <- ContrastFunction(Res[i + 1, ])
  }
  list(ParamHistory = Res, ContrastHistory = Contrast)
}


library(parallel)
MyData <- CompleteSample[seq(1, nrow(CompleteSample), by = 200), ]
DiffRes <- mclapply(1:10, function(i){
  Param0 <- c(runif(1, 0, 5), runif(1, -5, 5), runif(1, -5, 5))
  OptimContrast(Param0, MyData, maxiter = 5)
}, mc.cores = 3)
MatContVal <- sapply(DiffRes, function(x) x$ContrastHistory)
matplot(MatContVal)


FigPath <- "/home/gloaguen/Documents/LangevinMovementModel/figures/"
FigHeight <- 600
FigWidth <- 800
png(paste0(FigPath, "ParamEstimation.png"), height = FigHeight, width = FigWidth)
par(mar = c(4, 4, 1, 1))
plot(c(1, 5), c(-5, 5), type = "n", xlab = "", ylab = "")
MyCols <- function(alpha) 
  c(rgb(1, 0, 0, alpha), rgb(0, 1, 0, alpha), rgb(0, 0, 1, alpha))
mtext(side = c(1, 2), text = c("Iteration", "Param. Value"), line = c(2, 2),
      cex = 1.5)
lapply(DiffRes, function(x) matplot(x$ParamHistory, 
                                    col = MyCols(0.5) ,
                                    lty = 1, lwd = 2, type = "l", add = T))
abline(h = c(Alpha, Beta1, Beta2), col = MyCols(1), lty = 2, lwd = 2)
legend(x = 4, y = 5, lty = 1, col = MyCols(1), legend = c(expression(alpha),
                                                          expression(beta[1]),
                                                          expression(beta[2])),
       cex = 1.5, lwd = 2)
dev.off()


# Estimation over 100 trajectories ----------------------------------------

Procedure <- function(i){
  set.seed(i)
  CompleteSample <- SimulationProcess(Start = c(0, 0), Npoints = 10^5, Delta = 0.005)
  MyData <- CompleteSample[seq(1, nrow(CompleteSample), by = 200), ]
  OptimContrast(c(Alpha, Beta1, Beta2), MyData, maxiter = 4)$ParamHistory[5,]
}
Try <- mclapply(1:300, Procedure, mc.cores = 3)  
Results <- do.call(rbind, Try)
boxplot(Results, col = MyCols(0.5))
abline(h = c(Alpha, Beta1, Beta2), col = MyCols(1), lty = 2, lwd = 2)

png(paste0(FigPath, "BoxplotParams.png"), height = FigHeight, width = FigWidth)
par(mfrow = c(1, 3), mar = c(4, 4, 1, 1), oma = c(0, 0, 0, 0))
TruePar <- c(Alpha, Beta1, Beta2)
Names <- c(expression(alpha),
           expression(beta[1]),
           expression(beta[2]))
sapply(1:3, function(i){
  boxplot(Results[, i, drop = F],
                                col = MyCols(0.5)[i])
  mtext(side = 1, text = Names[i], cex = 1.5, line = 2)
  abline(h = TruePar[i], col = MyCols(1)[i], lwd = 2, lty = 2)
  })
dev.off()
