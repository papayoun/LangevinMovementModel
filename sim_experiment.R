
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

###############################
## Define covariates and RSF ##
###############################
# Simulate covariate fields
set.seed(1)
lim <- c(-100,100,-100,100) # limits of map
res <- 1 # grid resolution
rho <- c(10,100) # correlation coefficient
nbCov <- 2
covlist <- list()
for(i in 1:nbCov)
    covlist[[i]] <- simCov(rho=rho[i],lim=lim,res=res)

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

# Store covariates in an array (note that each layer is rotated, as required by interp.surface)
covarray <- array(NA,dim=c(length(xgrid),length(ygrid),length(covlist)))
for(i in 1:length(covlist))
    covarray[,,i] <- t(apply(as.matrix(covlist[[i]]),2,rev))

# Define artificial RSF
beta <- c(2,3,-7)
rsfRaster <- 0
for(i in 1:length(covlist))
    rsfRaster <- rsfRaster + beta[i]*covlist[[i]]
rsfRaster <- exp(rsfRaster)

###################
## Simulate data ##
###################
dt <- 0.1
Tmax <- 1e4
alltime <- seq(0,Tmax,by=dt)
set.seed(1)
allxy <- simLang(beta=beta, time=alltime, xgrid=xgrid, ygrid=ygrid, covarray=covarray)

# thin to emulate tracking data
nbObs <- 1000
keep <- sort(sample(1:nrow(allxy),size=nbObs,replace=FALSE))
xy <- allxy[keep,]
time <- alltime[keep]

plot(rsfRaster)
points(xy,type="o",cex=0.4)

###########################
## Fit with Euler method ##
###########################
# Evaluate covariate gradients at observed locations
gradarray <- covGrad(xy, xgrid, ygrid, covarray)
# Optimise log-likelihood
fit <- nlminb(start=rep(0,length(beta)), objective=nllkLang, xy=xy, time=time,
              gradarray=gradarray, control=list(trace=1))

###########################
## Fit with Ozaki method ## 
###########################
# Evaluate covariate hessian at observed locations
hessarray <- covHessian(xy, xgrid, ygrid, covarray)
# Optimise log-likelihood
fitOz <- nlminb(start=rep(0,length(beta)), objective=nllkLang, xy=xy, time=time,
                gradarray=gradarray, hessarray=hessarray, method="ozaki", 
                control=list(trace=1))
