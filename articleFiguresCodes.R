rm(list = ls())

library(raster)
library(ggplot2)
library(viridis)

source("archive/SimulationFunctions.R")
source("archive/CovariateParameters.R")
source("LinearAlgebraFunctions.R")
# source("optionsParsingScript.R")
SAVEPATH <- "simulatedData/analyticCase/"
library(extrafont)
A <- 0.05; B1 <- -1; B2 <- 0.5
library(parallel)
gammas <- c(1, 100)
samples <- mclapply(gammas, function(gam)
  SimulationProcess(Start = c(0, 0), Npoints = 2 * 10^6 + 1, 
                    Delta = 0.01, A = A, B = B1, C = B2, gamma = gam, method = "Euler"), mc.cores = 2)
xlim <- range(sapply(samples, function(x) range(x[ ,1])))
ylim <- range(sapply(samples, function(x) range(x[ ,2])))

ncells <- 401
myX <- seq(xlim[1], xlim[2], length.out = ncells + 1)
myY <- seq(ylim[1], ylim[2], length.out = ncells + 1)
grille <- as.matrix(expand.grid(myX, myY))
ud <-apply(grille, 1, function(x) UDMap(x, A, B1, B2, F))

pdf("UDFig1.pdf", width = 7, height = 6)
UDmap <- data.frame(x=grille[,1], y=grille[,2], val=ud)
UDplot <- ggplot(UDmap, aes(x,y)) + geom_raster(aes(fill=val)) +
    coord_equal() + scale_fill_viridis(name=expression(pi)) +
    theme(axis.title = element_blank(), axis.text = element_blank(), 
          axis.ticks = element_blank(), legend.title = element_text(size=15), 
          legend.text = element_text(size=12))

plot(UDplot)
dev.off()

# par(mfrow = c(1, 2))
# sapply(samples, function(x){
#   sel <- seq(1, nrow(x), by = 1000)
#   kernelEst <-   kde2d(x[sel, 1], x[sel, 2], n = ncells, lims = c(xlim, ylim))
#   image(kernelEst, col = terrain.colors(30), , xaxt = "n", yaxt = "n",
#         xlab = "", ylab = "")
#   # points(x[sel, 1:2], pch =20, cex = 0.2)
# })

pdf("ExPaths.pdf", width=8, height=8)
par(mfrow = c(2, 2), mar = rep(0, 4))
mapply(function(x, y){
  sel1 <- seq(1, 5001, by = 10)
  t1 <- x[sel1, 3][length(sel1)]
  plot(x[sel1, 1:2], type = "l", lwd = 0.5, asp = 1, xlim = xlim, ylim = ylim, 
       col="dimgrey", xlab = "", ylab = "", xaxt = "n", yaxt = "n")
  points(x[sel1, 1:2], type = "p", pch = 20, cex = 0.2)
  points(x[1, 1:2, drop = F], type = "p", pch = 20, cex = 2, col = "red")
  text(xlim[1] + 5, ylim[1] + 2, substitute(paste(Delta, " = ", b), list(b = t1)), cex = 2)
  text(xlim[1] + 5, ylim[2] - 2, substitute(paste(gamma, " = ", b), list(b = y)), cex = 2)
  
  sel2 <- seq(1, 100001, by = 10)
  t2 <- x[sel2, 3][length(sel2)]
  plot(x[sel2, 1:2], type = "l", lwd = 0.5, asp=1, xlim = xlim, ylim = ylim, 
       col="dimgrey", xlab = "", ylab = "", xaxt = "n", yaxt = "n")
  points(x[sel2, 1:2], type = "p", pch = 20, cex = 0.2)
  points(x[1, 1:2, drop = F], type = "p", pch = 20, cex = 2, col = "red")
  text(xlim[1] + 5, ylim[1] + 2, substitute(paste(Delta, " = ", b), list(b = t2)), cex = 2)
  text(xlim[1] + 5, ylim[2] - 2, substitute(paste(gamma, " = ", b), list(b = y)), cex = 2)
}, samples, gammas)
dev.off()
par(mfrow = c(1, 1), mar = c(5, 5, 4, 1))
speeds <- lapply(samples, function(x){
  n <- nrow(x)
  sqrt(rowSums((x[-1, 1:2] - x[-n, 1:2])^2)) / diff(x[, 3])
})
speedsHist <- lapply(speeds, function(x){
  res <- density(x)
  res$y[res$x < 0] <- 0
  res
})
yHistLims <- c(0, max(sapply(speedsHist, function(z) max(z$y))))
xHistLims <- c(0, max(sapply(speedsHist, function(z) max(z$x))))
plot(xHistLims, yHistLims, main = "Speed distribution", type = "n", ylab = "Density", xlab = "Speed")
sapply(1:2, function(i) lines(speedsHist[[i]], col = i, lwd = 2))
legend("topright", bty = "n", col = 1:2, lty = 1, lwd = 2, cex = 2,
       legend  = sapply(gammas, 
                                 function(gam) as.expression(substitute(gamma~"="~B, list(B = gam)))))
xlab <- bquote(.(assay) ~ AC50 ~ (mu*M))
plot(sample1[, 1:2], pch =20, cex = 0.2)
points(sample2[, 1:2], pch = 20, cex = 0.2, col = "red")
library(MASS)
?kde2d
sel <- seq(1, nrow(sample1), by = 100)
myEst <- kde2d(sample1[sel, 1], sample1[sel, 2], n = 101)
image(myEst)
myEst <- kde2d(sample2[sel, 1], sample2[sel, 2], n = 101)
image(myEst)
speeds