
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

# Simulating Data with tiny dt --------------------------------------------
# Delta <- 10^(-2)
# simTimes <- seq(0, 500, by = Delta)
# xySim <- simLang(beta = beta, time = simTimes,
#                  xgrid=xgrid, ygrid=ygrid, covarray=covarray)
# WholeData <- cbind(xySim, simTimes)
# colnames(WholeData) <- c("x", "y", "t")
# write.table(WholeData, file = "ReferenceSimulatedDataSet.txt", col.names = T, row.names = F)
WholeData <- read.table("ReferenceSimulatedDataSet.txt", header = T)
nbObs <- 500
Sel <- sort(sample(1:nrow(WholeData), size = nbObs, replace = F))
time <- WholeData$t[Sel]
xydf <- WholeData[Sel, c("x", "y")]
xy <- as.matrix(xydf)
rownames(xy) <- NULL
# xy <- simLang(beta=beta, time=time, xgrid=xgrid, ygrid=ygrid, covarray=covarray)



# Plot simulated data on artificial RSF
ggplot(map,aes(x,y)) + geom_raster(aes(fill=z)) + coord_equal() +
    scale_fill_gradientn(colours=mypal,name="RSF") +
    geom_point(aes(x,y),xydf,size=0.3) + geom_path(aes(x,y),xydf,size=0.3)

# Evaluate covariate gradients at observed locations
gradarray <- covGrad(xy, xgrid, ygrid, covarray)
jacobarray <- covHessian(xy, xgrid, ygrid, covarray)
# Optimise log-likelihood
fit <- nlminb(start=c(0,0),objective=nllkLang,xy=xy,time=time,gradarray=gradarray,
              control=list(trace=1))
fitOz <- nlminb(start=c(0,0),objective=nllkLang,xy=xy,time=time,gradarray=gradarray,
                jacobarray = jacobarray, method = "ozaki",
              control=list(trace=1))
