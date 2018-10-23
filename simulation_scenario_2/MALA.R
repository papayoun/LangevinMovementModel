
# Load packages
library(raster)
library(numDeriv)
library(fields)
library(expm)
library(ggplot2)
library(viridis)
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

# Store covariates in an array 
# (note that each layer is rotated, as required by interp.surface)
covarray <- array(NA, dim = c(length(xgrid),length(ygrid),length(covlist)))
for(i in 1:length(covlist))
    covarray[,,i] <- t(apply(as.matrix(covlist[[i]]),2,rev))

# Define artificial RSF
beta <- c(2, 4) # true parameter values (used in simulations)
rsfRaster <- 0
for(i in 1:length(covlist))
    rsfRaster <- rsfRaster + beta[i]*covlist[[i]]
rsfRaster <- exp(rsfRaster)

###################
## Simulate data ##
###################
set.seed(1)
Tmax <- 100
alldt <- c(0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1)
speed <- 0.5

allsim <- list()
t0 <- Sys.time()
for(iter in 1:length(alldt)) {
    cat("Iteration",iter,"\n")
    dt <- alldt[iter]
    time <- seq(0, Tmax, by = dt)    
    allsim[[iter]] <- simLangMH(beta=beta, speed=speed, time=time, xy0=c(0,0), 
                                xgrid=xgrid, ygrid=ygrid, covarray=covarray)
    print(Sys.time()-t0)
}

# Derive acceptance rates
lapply(allsim, function(sim) sim$acc/(sim$acc+sim$rej))
