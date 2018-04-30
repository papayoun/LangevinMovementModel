
# Load packages
library(raster)
library(numDeriv)
library(fields)
library(expm)
# Load other source files
source("simCov.R")
source("RSF.R")
source("utility.R")
source("simLang.R")
source("nllkLang.R")
source("eulerLSE.R")

###############################
## Define covariates and RSF ##
###############################
# Simulate covariate fields
set.seed(1)
lim <- c(-50,50,-50,50) # limits of map
res <- 1 # grid resolution
rho <-   c(5,10) # correlation coefficients
nbCov <- 2
covlist <- list()
for(i in 1:nbCov)
    covlist[[i]] <- simCov(rho = rho[i], lim = lim, res = res)

# Coordinates of centres of grid cells
xgrid <- seq(lim[1]+res/2,lim[2]-res/2,by=res)
ygrid <- seq(lim[3]+res/2,lim[4]-res/2,by=res)

# Include squared distance to centre of map as covariate
xygrid <- expand.grid(xgrid,ygrid)
mu <- c((lim[1]+lim[2])/2,(lim[3]+lim[4])/2) # centre of map
dist2 <- (xygrid[,1]-mu[1])^2+(xygrid[,2]-mu[2])^2
dist2 <- dist2/max(dist2) # scale to [0,1]
dist2XYZ <- cbind(xygrid,dist2)
covlist[[nbCov+1]] <- rasterFromXYZ(dist2XYZ)

# Store covariates in an array 
# (note that each layer is rotated, as required by interp.surface)
covarray <- array(NA, dim = c(length(xgrid),length(ygrid),length(covlist)))
for(i in 1:length(covlist))
    covarray[,,i] <- t(apply(as.matrix(covlist[[i]]),2,rev))

# Define artificial RSF
beta <- c(1, -1, -10) # true parameter values (used in simulations)
rsfRaster <- 0
for(i in 1:length(covlist))
    rsfRaster <- rsfRaster + beta[i]*covlist[[i]]
rsfRaster <- exp(rsfRaster)
plot(rsfRaster)
plot(log(rsfRaster))

###################
## Simulate data ##
###################

ntrack <- 2
dt <- 1
Tmax <- 500
time <- seq(0, Tmax, by = dt)
data <- NULL

set.seed(1)
for(zoo in 1:ntrack) {
    track <- simLang(beta=beta, time=time, xy0=runif(2, -10, 10), 
                     xgrid=xgrid, ygrid=ygrid, covarray=covarray)
    data <- rbind(data, cbind(rep(zoo,nrow(track)),track,time))
}
colnames(data) <- c("ID","x","y","time")

# # save full simulated data set before thinning
# alldata <- data
# # thin
# nobs <- nrow(alldata)
# ind <- sort(sample(1:nrow(alldata),size=nobs,replace=FALSE))
# data <- alldata[ind,]

ID <- data[,"ID"]
xy <- data[,c("x","y")]
time <- data[,"time"]

plot(log(rsfRaster))
lapply(split(data.frame(xy), ID), function(z) {
  points(z,type="o",cex=0.4); 
  points(z[c(1, nrow(z)), ], col = c("blue", "red"), pch = 20)})

# Evaluate covariate gradients at observed locations
gradarray <- covGrad(xy, xgrid, ygrid, covarray)

######################################
## Fit UD with Euler discretization ##
######################################
# Via numerical MLE
# fit <- nlminb(start = beta, objective = nllkLang, xy = xy, time = time, ID=ID,
#               gradarray = gradarray, control = list(trace=1))

# Via least squares
lse <- eulerLSE(ID=ID, time=time, xy=xy, gradarray=gradarray)

# 95% CI
cbind(lse$Bhat-1.96*sqrt(diag(lse$var)),lse$Bhat,lse$Bhat+1.96*sqrt(diag(lse$var)))

# Compute estimated RSF
betaMLE <- lse$Bhat
rsfRasterMLE <- 0
for(i in 1:length(covlist))
    rsfRasterMLE <- rsfRasterMLE + betaMLE[i]*covlist[[i]]
rsfRasterMLE <- exp(rsfRasterMLE)
plot(rsfRasterMLE)
plot(rsfRaster)
# Plot estimated utilisation vs true utilisation for each grid cell 
# (should align with identity line)
plot(values(rsfRaster)/sum(values(rsfRaster)),
     values(rsfRasterMLE)/sum(values(rsfRasterMLE)))
abline(0,1,col=2)
