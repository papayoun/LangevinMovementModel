
# Load packages
library(raster)
library(numDeriv)
library(fields)
library(expm)
library(ggplot2)
library(viridis)
library(ash)
# Load other source files
source("simCov.R")
source("utility.R")
source("nllkLang.R")
source("eulerLSE.R")
source("simLangMH.R")

###############################
## Define covariates and RSF ##
###############################
# Simulate covariate fields
set.seed(1)
lim <- c(-50,50,-50,50) # limits of map
res <- 1 # grid resolution
rho <- c(10,10) # correlation coefficients
ncov <- 2
covlist <- list()
for(i in 1:ncov)
    covlist[[i]] <- simCov(rho = rho[i], lim = lim, res = res)

# Coordinates of centres of grid cells
xgrid <- seq(lim[1]+res/2,lim[2]-res/2,by=res)
ygrid <- seq(lim[3]+res/2,lim[4]-res/2,by=res)
xygrid <- expand.grid(xgrid, ygrid)

# Include squared distance to centre of map as covariate
dist2 <- ((xygrid[,1])^2+(xygrid[,2])^2)/100
covlist[[3]] <- rasterFromXYZ(cbind(xygrid,dist2))

# Store covariates in an array 
# (note that each layer is rotated, as required by interp.surface)
covarray <- array(NA, dim = c(length(xgrid),length(ygrid),length(covlist)))
for(i in 1:length(covlist))
    covarray[,,i] <- t(apply(as.matrix(covlist[[i]]),2,rev))

# Define artificial RSF
beta <- c(2, 4, -1) # true parameter values (used in simulations)
rsfRaster <- 0
for(i in 1:length(covlist))
    rsfRaster <- rsfRaster + beta[i]*covlist[[i]]
rsfRaster <- exp(rsfRaster)

###################
## Simulate data ##
###################
set.seed(1)
dt <- 0.01
nobs <- 1e6
speed <- 1
time <- dt*(1:nobs)

t0 <- Sys.time()
simdat <- simLangMH(beta=beta, speed=speed, time=time, xy0=c(0,0), 
                    xgrid=xgrid, ygrid=ygrid, covarray=covarray)
print(Sys.time()-t0)

############################
## Empirical distribution ##
############################
# Count locations in grid cells
xy <- simdat$xy
bins <- bin2(x = xy,
             ab = matrix(lim,2,2,byrow=TRUE),
             nbin = c(length(xgrid),length(ygrid)))
foo <- bins$nc[,ncol(bins$nc):1,drop=F] # flip horizontally
empUD <- as.vector(foo)/sum(foo)

trueUD <- values(rsfRaster)/sum(values(rsfRaster))

ind <- which(empUD>0)
plotlim <- range(trueUD[ind], empUD[ind])
plot(trueUD[ind], empUD[ind], log="xy")
abline(0,1)
