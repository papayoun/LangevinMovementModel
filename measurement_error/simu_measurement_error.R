
library(raster)
library(numDeriv)
library(fields)
library(expm)
source("../simCov.R")
source("../simLang.R")
source("../utility.R")
source("kalman.R")

###############################
## Define covariates and RSF ##
###############################
# Simulate covariate fields
set.seed(1)
lim <- c(-50,50,-50,50) # limits of map
res <- 1 # grid resolution
rho <- c(5,10) # correlation coefficients
nbCov <- 2
covlist <- list()
for(i in 1:nbCov)
    covlist[[i]] <- simCov(rho = rho[i], lim = lim, res = res)

# Coordinates of centres of grid cells
xgrid <- seq(lim[1]+res/2,lim[2]-res/2,by=res)
ygrid <- seq(lim[3]+res/2,lim[4]-res/2,by=res)

# # Include squared distance to centre of map as covariate
# xygrid <- expand.grid(xgrid,ygrid)
# mu <- c((lim[1]+lim[2])/2,(lim[3]+lim[4])/2) # centre of map
# dist2 <- (xygrid[,1]-mu[1])^2+(xygrid[,2]-mu[2])^2
# dist2 <- dist2/max(dist2) # scale to [0,1]
# dist2XYZ <- cbind(xygrid,dist2)
# covlist[[nbCov+1]] <- rasterFromXYZ(dist2XYZ)
# nbCov <- nbCov + 1

# Store covariates in an array 
# (note that each layer is rotated, as required by interp.surface)
covarray <- array(NA, dim = c(length(xgrid),length(ygrid),length(covlist)))
for(i in 1:length(covlist))
    covarray[,,i] <- t(apply(as.matrix(covlist[[i]]),2,rev))

# Define artificial RSF
beta <- c(-2, 3) # true parameter values (used in simulations)
rsfRaster <- 0
for(i in 1:length(covlist))
    rsfRaster <- rsfRaster + beta[i]*covlist[[i]]
rsfRaster <- exp(rsfRaster)
plot(rsfRaster)
plot(log(rsfRaster))

###################
## Simulate data ##
###################
ntrack <- 10
dt <- 0.1
Tmax <- 150
speed <- 1.5
data <- NULL
time <- seq(0, Tmax, by = dt)

set.seed(1)
for(zoo in 1:ntrack) {
    track <- simLang(beta=beta, speed=speed, time=time, xy0=runif(2, -10, 10), 
                     xgrid=xgrid, ygrid=ygrid, covarray=covarray)
    data <- rbind(data, cbind(rep(zoo,nrow(track)),track,time))
}
colnames(data) <- c("ID","x","y","time")

ID <- data[,"ID"]
xy <- data[,c("x","y")] + rnorm(2*nrow(data), 0, 0.5)
time <- data[,"time"]

# Evaluate covariate gradients at observed locations
gradarray <- covGrad(xy, xgrid, ygrid, covarray)

# Optimise log-likelihood
beta0 <- c(0,0)
speed0 <- 2
h0 <- 1
par <- c(beta0, log(speed0), log(h0))
fit <- nlminb(start = par, objective = kalman, xy = xy, time = time, ID=ID,
              gradarray = gradarray, control = list(trace=1))
