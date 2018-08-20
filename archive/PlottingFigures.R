# source("ExampleSimulationFunctionsAndParameters.R")
CompleteSample <- read.table("archive/SimulatedSample.txt", header = T)

FigHeight <- 600
FigWidth <- 800
FigPath <- "/home/gloaguen/Documents/LangevinMovementModel/figures/"

MarDef <- c(0.1, 0.1, 0.1, 0.1)
# Computing required matrix ---------------------------------------------------------


Xs <- seq(-10, 10, length.out = 250)
Grid <- expand.grid(Xs, Xs)
Z1 <- matrix(apply(Grid, 1, Fun, DilatationParam = DP1, CenterParam = CP1,
                   IntensityParam = IP1, FrequenceParam = FP1), nrow = length(Xs))
Z2 <- matrix(apply(Grid, 1, Fun, DilatationParam = DP2, CenterParam = CP2,
                   IntensityParam = IP2, FrequenceParam = FP2), nrow = length(Xs))
NormVect <- matrix(apply(Grid, 1, function(x) sum(x^2) * Alpha), nrow = length(Xs))
GradZ1 <- apply(Grid, 1, GradFun, DilatationParam = DP1, CenterParam = CP1,
                IntensityParam = IP1, FrequenceParam = FP1)
GradZ2 <- apply(Grid, 1, GradFun, DilatationParam = DP1, CenterParam = CP2,
                IntensityParam = IP2, FrequenceParam = FP2)
NormGrad <- matrix(sqrt(colSums((t(-2 * Grid * Alpha)
                                 + (Beta1 * GradZ1)  
                                 +  (Beta2 * GradZ2))^2)), nrow = length(Xs))

MyPotFun <- exp(-NormVect + Beta1 * Z1 + Beta2 * Z2)

# Covariatex ----------------------------------------------------------
png(paste0(FigPath, "Covariate1.png"), height = FigHeight, width = FigWidth)
par(mar = MarDef)
image(Xs, Xs, Z1, col = terrain.colors(n = 100), 
      xlab = "", ylab = "", xaxt = "n",
      yaxt = "n")
box()
dev.off()
png(paste0(FigPath, "Covariate2.png"), height = FigHeight, width = FigWidth)
par(mar = MarDef)
image(Xs, Xs, Z2, col = terrain.colors(n = 100), 
      xlab = "", ylab = "", xaxt = "n",
      yaxt = "n")
box()
dev.off()
library(fields)

# Utilization Distribution ------------------------------------------------
png(paste0(FigPath, "ResultingUD.png"), height = FigHeight, width = FigWidth)
par(mar = MarDef)
image(Xs, Xs, MyPotFun, col = terrain.colors(n = 100), 
      xlab = "", ylab = "", xaxt = "n",
      yaxt = "n")
box()
dev.off()

png(paste0(FigPath, "TrajPlusUD.png"), height = FigHeight, width = FigWidth)
par(mar = MarDef)
image(Xs, Xs, MyPotFun, col = terrain.colors(n = 100), 
      xlab = "", ylab = "", xaxt = "n",
      yaxt = "n")
Sel <- seq(1, nrow(CompleteSample), by = 100)
points(CompleteSample[Sel, 1:2], col = 1, pch = 20, cex = 0.5, type = "b")
points(CompleteSample[1, 1:2, drop = F], col = "red", cex = 2, pch = 20)
box()
dev.off()
png(paste0(FigPath, "TrajPlusC1.png"), height = FigHeight, width = FigWidth)
par(mar = MarDef)
image(Xs, Xs, Z1, col = terrain.colors(n = 100), 
      xlab = "", ylab = "", xaxt = "n",
      yaxt = "n")
points(CompleteSample[Sel, 1:2], col = 1, pch = 20, cex = 0.5, type = "b")
points(CompleteSample[1, 1:2, drop = F], col = "red", cex = 2, pch = 20)
box()
dev.off()
png(paste0(FigPath, "TrajPlusC2.png"), height = FigHeight, width = FigWidth)
par(mar = MarDef)
image(Xs, Xs, Z2, col = terrain.colors(n = 100), 
      xlab = "", ylab = "", xaxt = "n",
      yaxt = "n")
points(CompleteSample[Sel, 1:2], col = 1, pch = 20, cex = 0.5, type = "b")
points(CompleteSample[1, 1:2, drop = F], col = "red", cex = 2, pch = 20)
box()
dev.off()

png(paste0(FigPath, "GradientNorm.png"), height = FigHeight, width = FigWidth)
par(mar = MarDef)
image(Xs, Xs, NormGrad, col = terrain.colors(n = 100), 
      xlab = "", ylab = "", xaxt = "n",
      yaxt = "n")
box()
dev.off()

image(Xs, Xs, MyPotFun, col = terrain.colors(n = 100), xlim = c(-10, 10),
           ylim = c(-10, 10))
image(Xs, Xs, log(MyPotFun), col = terrain.colors(n = 100), xlim = c(-10, 10),
           ylim = c(-10, 10))
contour(Xs, Xs, log(MyPotFun), col = terrain.colors(n = 100), xlim = c(-10, 10),
              ylim = c(-10, 10))
image(Xs, Xs, NormGrad, xlim = c(-10, 10),col = terrain.colors(100), ylim = c(-10, 10))
Sel <- seq(1, nrow(CompleteSample), by = 100)
points(CompleteSample[Sel, 1:2])