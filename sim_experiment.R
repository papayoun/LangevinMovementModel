
# Load packages
library(numDeriv)
library(fields)
# Load other source files
source("resources.R")
source("RSF.R")
source("utility.R")
source("simLang.R")
source("nllkLang.R")

# List of covariate rasters
covlist <- list(cov1rast,cov2rast)
# Extract limits of map, and x/y resolutions
lim <- as.vector(extent(covlist[[1]]))
res <- c(xres(covlist[[1]]),yres(covlist[[1]]))
xgrid <- seq(lim[1]+res[1]/2,lim[2]-res[1]/2,by=res[1])
ygrid <- seq(lim[3]+res[2]/2,lim[4]-res[2]/2,by=res[2])

# Store covariates in an array (note that each layer is rotated, as required by interp.surface)
covarray <- array(NA,dim=c(length(xgrid),length(ygrid),length(covlist)))
for(i in 1:length(covlist))
    covarray[,,i] <- t(apply(as.matrix(covlist[[i]]),2,rev))

# Define artificial RSF
beta <- c(4,-2)
rsfRaster <- exp(beta[1]*covlist[[1]]+beta[2]*covlist[[2]])
plot(rsfRaster)

# Simulate data
set.seed(1)
dt <- 1
xy <- simLang(nbObs=500, beta=beta, h=dt, xgrid=xgrid, ygrid=ygrid, covarray=covarray)
points(xy,type="l")

# Evaluate covariate gradients at observed locations
gradarray <- array(NA,dim=c(nrow(xy),2,2))
for(i in 1:dim(covarray)[3])
    gradarray[,,i] <- apply(xy,1,function(x) grad(interpCov,x,xgrid=xgrid,ygrid=ygrid,covmat=covarray[,,i]))

# Optimise log-likelihood
fit <- nlminb(start=c(0,0),objective=nllkLang,xy=xy,gradarray=gradarray,dt=dt,
              control=list(trace=1))
