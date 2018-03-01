
# Load packages
library(numDeriv)
library(fields)
library(ggplot2)
# Load other source files
source("resources.R")
source("RSF.R")
source("utility.R")
source("simLang.R")
source("nllkLang.R")
mypal <- rev(RColorBrewer::brewer.pal(11,"RdYlBu"))

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
map <- data.frame(x=coordinates(rsfRaster)[,1],y=coordinates(rsfRaster)[,2],z=values(rsfRaster))

# Simulate data
set.seed(1)
nbObs <- 500
time <- cumsum(rgamma(nbObs,10,10))
xy <- simLang(nbObs=nbObs, beta=beta, time=time, xgrid=xgrid, ygrid=ygrid, covarray=covarray)
xydf <- data.frame(x=xy[,1],y=xy[,2])

# Plot simulated data on artificial RSF
ggplot(map,aes(x,y)) + geom_raster(aes(fill=z)) + coord_equal() +
    scale_fill_gradientn(colours=mypal,name="RSF") +
    geom_point(aes(x,y),xydf,size=0.3) + geom_path(aes(x,y),xydf,size=0.3)

# Evaluate covariate gradients at observed locations
gradarray <- covGrad(xy, xgrid, ygrid, covarray)

# Optimise log-likelihood
fit <- nlminb(start=c(0,0),objective=nllkLang,xy=xy,time=time,gradarray=gradarray,
              control=list(trace=1))
