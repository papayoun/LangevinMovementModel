
library(ggplot2)
library(viridis)
library(parallel)
library(extrafont)
source("archive/SimulationFunctions.R")
source("archive/CovariateParameters.R")
source("nllkLang.R")
source("eulerLSE.R")
source("utility.R")

DATAPATH <- "simulatedData/analyticCase/"
AllResults <- lapply(1:600, function(seed){
  res <-  try(load(paste0(DATAPATH, "tmp", seed, ".RData")))
  if(inherits(res, "try-error"))
    return(NULL)
  else 
    return(result)
})

sampleIndex <- floor(seq(1, 600, length.out = 10))
sampleData <- do.call(rbind, mclapply(sampleIndex, function(seed){
  load(paste0(DATAPATH, "tmp", seed, ".RData"))
  CompleteSample <- read.table(paste0(DATAPATH, "simAnalyticData_seed", seed, ".txt"),
                               header = T, sep = ";")
  myData <- CompleteSample[result$selectedLines, ]
  cbind(myData, Traj = rep(seed, nrow(myData)))
}, mc.cores = 3))

myX <- myY <- seq(-20 , 20, length = 201)

getCovariates <-  function(xGrid, yGrid, params = c(A = 0.05, B1 = -1, B2 = 0.5)){
  covariateSampling <- as.matrix(expand.grid(xGrid, yGrid))
  A = params["A"]; B1 = params["B1"]; B2 = params["B2"]
  covariateArray <- array(dim = c(length(xGrid), length(yGrid),  3))
  covariateArray[,, 1] <- B1 * matrix(apply(covariateSampling, 1, ScalC1), nrow = length(xGrid)); 
  covariateArray[,, 2] <- B2 * matrix(apply(covariateSampling, 1, ScalC2), nrow = length(xGrid))
  covariateArray[,, 3] <- -A * matrix(apply(covariateSampling, 1, function(x) sum(x^2)), nrow = length(xGrid)); 
  covariateArray
}

# Plot covariates
grid <- expand.grid(myX,myY)
c1map <- data.frame(x=grid[,1], y=grid[,2], val=apply(grid, 1, ScalC1))
c2map <- data.frame(x=grid[,1], y=grid[,2], val=apply(grid, 1, ScalC2))
c3map <- data.frame(x=grid[,1], y=grid[,2], val=apply(grid, 1, function(x) sum(x^2)))

ggopts <- theme(axis.title = element_blank(), axis.text = element_blank(), 
                axis.ticks = element_blank(), legend.title = element_text(size=25), 
                legend.text = element_text(size=20), legend.key.height=unit(3,"line"))

pdf("covAnalytic1.pdf", width = 7, height = 6)
c1plot <- ggplot(c1map, aes(x,y)) + geom_raster(aes(fill=val)) +
    coord_equal() + scale_fill_viridis(name=expression(c[1])) + ggopts
plot(c1plot)
dev.off()

pdf("covAnalytic2.pdf", width = 7, height = 6)
c2plot <- ggplot(c2map, aes(x,y)) + geom_raster(aes(fill=val)) +
    coord_equal() + scale_fill_viridis(name=expression(c[2])) + ggopts
plot(c2plot)
dev.off()

pdf("covAnalytic3.pdf", width = 7, height = 6)
c3plot <- ggplot(c3map, aes(x,y)) + geom_raster(aes(fill=val)) +
    coord_equal() + scale_fill_viridis(name=expression(c[3])) + ggopts
plot(c3plot)
dev.off()


covariatesModel <- function(xGrid, yGrid, params = c(A = 0.05, B1 = -1, B2 = 0.5)){
  A = params["A"]; B1 = params["B1"]; B2 = params["B2"]
  covariateSampling <- as.matrix(expand.grid(xGrid, yGrid))
  covariateArray <- array(dim = c(length(xGrid), length(yGrid),  2))
  covariateArray[,, 1] <- matrix(apply(covariateSampling, 1, ScalC1), nrow = length(xGrid)); 
  covariateArray[,, 2] <- matrix(apply(covariateSampling, 1, ScalC2), nrow = length(xGrid))
  distance <- matrix(apply(covariateSampling, 1, function(x) sum(x^2)), nrow = length(xGrid)); 
  exp(- A * distance + B1 * covariateArray[,, 1] + B2 * covariateArray[,, 2])
}
png(paste0(FIGPATH, "UDAnalytic.png"), height = 600, width = 800)
par(mfrow = c(1, 1), mar = c(0.1, 0.1, 3, 0.1), family = "LM Roman 10")
covArrays <- getCovariates(myX, myY)
image(myX, myY, covariatesModel(myX, myY), xaxt = "n", yaxt = "n", xlab = "", ylab = "", col = terrain.colors(30),
      main = expression(pi(z)), cex.main = 3)
dev.off()
save(sampleData, file = paste0(FIGPATH, "sampleData.RData"))
png(paste0(FIGPATH, "sampleTrajectories.png"), height = 600, width = 800)
par(mfrow = c(1, 1), mar = c(0.1, 0.1, 3, 0.1), family = "LM Roman 10")
covArrays <- getCovariates(myX, myY)
image(myX, myY, covariatesModel(myX, myY), xaxt = "n", yaxt = "n", xlab = "", ylab = "", col = "white")
sapply(split(sampleData[,1:2], sampleData$Traj), function(x) {
  points(x, cex = 0.5, type = "b", pch = 20)})
sapply(split(sampleData[,1:2], sampleData$Traj), function(x) {
  points(x[1, 1], x[1, 2], col = "red", pch = 20, cex = 2)})
dev.off()
load(file = paste0(FIGPATH, "sampleData.RData"))
fullCovariates <- getCovariates(myX, myY)
image(myX, myY, getCovariates(myX, myY)[,,1], col = terrain.colors(30))
points()
points(sampleData[, 1:2], pch = 20, cex = 0.5, type = "b")
library(MASS)
kernelEst <- kde2d(sampleData[, 1], sampleData[, 2], n = 201, lims = c(-20,20,-20,20))
image(kernelEst, col = terrain.colors(30))
xGrid <- yGrid <- seq(-15, 15, length = 8)
image(myX, myY, getCovariates(myX, myY)[,,1], col = terrain.colors(30))
points(expand.grid(xGrid, yGrid), col = "red")
points(sampleData, type = "b", pch = 20, cex = 0.5)
 
results <- do.call(rbind.data.frame,
                   lapply(1:length(AllResults), function(i){
                     x <- AllResults[[i]]
                     if(is.null(x))
                       return(NULL)
                     out <- data.frame(B1 = rep(NA, 4), B2 = rep(NA, 4), A = rep(NA, 4), 
                                       method = rep(c("euler", "ozaki"), c(2,2)),
                                       gradient = rep(c("analytic", "interpolated"), 2),
                                       seed = rep(i, 4))
                     out[out$method == "euler" & out$gradient == "analytic", 1:3] <- x$fitEulTG$Bhat
                     out[out$method == "ozaki" & out$gradient == "analytic", 1:3] <- x$fitOzTG$par
                     out[out$method == "euler" & out$gradient == "interpolated", 1:3] <- x$fitEulIG$Bhat
                     out[out$method == "ozaki" & out$gradient == "interpolated", 1:3] <- x$fitOzIG$par
                     out
                   }))
trueParams <-  c(A = 0.05, B1 = -1, B2 = 0.5)
myBox <- function(x) {
  boxplot(results[, x] ~ results$method + results$gradient, xlab = "", ylab = "", pch = 20,
          main = "", names = c("ETG", "OTG", "EIG", "OIG"), cex.lab = 4)
  abline(h = trueParams[x])
}
png(paste0(FIGPATH, "estBeta1Sc1.png"), height = 600, width = 800)
par(mfrow = c(1, 1), mar = c(3, 3, 0.1, 0.1), family = "LM Roman 10", las = 2)
myBox("B1")
dev.off()
png(paste0(FIGPATH, "estBeta2Sc1.png"), height = 600, width = 800)
par(mfrow = c(1, 1), mar = c(3, 3, 0.1, 0.1), family = "LM Roman 10", las = 2)
myBox("B2")
dev.off()
png(paste0(FIGPATH, "estBeta3Sc1.png"), height = 600, width = 800)
par(mfrow = c(1, 1), mar = c(3, 3, 0.1, 0.1), family = "LM Roman 10", las = 2)
myBox("A")
dev.off()
boxplot(B2 ~ method + gradient, data = results)
boxplot(A ~ method + gradient, data = results)

medianEstimateOz <- apply(results[results$method == "ozaki" & results$gradient == "analytic",1:3], 2, quantile, prob = 0.5)
medianEstimateEu <- apply(results[results$method == "euler" & results$gradient == "analytic",1:3], 2, quantile, prob = 0.5)

image(myX, myY, covariatesModel(myX, myY), col = terrain.colors(30))
image(myX, myY, covariatesModel(myX, myY, medianEstimateOz), col = terrain.colors(30))
image(myX, myY, covariatesModel(myX, myY, medianEstimateEu), col = terrain.colors(30))

