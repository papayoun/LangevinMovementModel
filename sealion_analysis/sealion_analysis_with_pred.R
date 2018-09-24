
library(raster)
library(numDeriv)
library(fields)
library(expm)
library(viridis)
library(ggplot2)
library(gridExtra)
source("utility.R")
source("eulerLSE.R")
####################
## Prepare tracks ##
####################
# load regularised track (from wilsonEtAl_with_pred.R)
tracks <- read.csv("SSLpreddat.csv")
ID <- as.integer(tracks$ID)
ID[ID==14809] <- 1
ID[ID==15136] <- 2
ID[ID==15137] <- 3
time <- as.POSIXct(tracks$time)
time <- as.numeric(time)
time <- (time-min(time))/3600
xy <- matrix(c(tracks$x, tracks$y)/1000, ncol=2) # convert to km

# for plots
xydf <- data.frame(x=xy[,1], y=xy[,2])

#####################
## Load covariates ##
#####################
# download rasters from Github URL
download.file(paste0(baseURL,"habitat/aleut_habitat.grd"), destfile="aleut_habitat.grd")
download.file(paste0(baseURL,"habitat/aleut_habitat.gri"), destfile="aleut_habitat.gri")
hbfull <- brick("aleut_habitat.grd", values=TRUE)
covlist <- list(bathy = hbfull$bathy,
                slope = hbfull$slope,
                d2site = hbfull$d2site,
                d2shelf = hbfull$d2shelf)

# convert to km
for(i in 1:length(covlist)) {
    extent(covlist[[i]]) <- extent(c(xmin(covlist[[i]]), xmax(covlist[[i]]), 
                                     ymin(covlist[[i]]), ymax(covlist[[i]]))/1000)
    projection(covlist[[i]]) <- gsub("units=m", "units=km", projection(covlist[[i]]))
}

ncov <- length(covlist)
# resample covariates to the same grid
for(i in 2:ncov)
    covlist[[i]] <- resample(covlist[[i]],covlist[[1]])

# crop covariates to area of interest
border <- 30
lim <- c(min(xy[,1])-border,max(xy[,1])+border,min(xy[,2])-border,max(xy[,2])+border)
covlist <- lapply(covlist, crop, y=extent(lim))
lim <- as.vector(extent(covlist[[1]]))
res <- res(covlist[[1]])

# Store covariates in an array (note that each layer is rotated, 
# as required by interp.surface)
xgrid <- seq(lim[1]+res[1]/2, lim[2]-res[1]/2, by=res[1])
ygrid <- seq(lim[3]+res[2]/2, lim[4]-res[2]/2, by=res[2])
covarray <- array(NA,dim=c(length(xgrid),length(ygrid),length(covlist)))
for(i in 1:length(covlist))
    covarray[,,i] <- t(apply(as.matrix(covlist[[i]]),2,rev))

#####################
## Plot covariates ##
#####################
covnames <- c("Bathymetry", "Slope", "Distance to SSL sites", "Distance to shelf")
covplot <- list()
for(i in 1:ncov) {
    mycov <- covlist[[i]]
    covmap <- data.frame(coordinates(mycov),val=values(mycov))
    covplot[[i]] <- ggplot(covmap,aes(x,y)) + geom_raster(aes(fill=val)) +
        coord_equal() + scale_fill_viridis(name=bquote(paste("c"[.(i)]))) +
        ggtitle(covnames[i])
}
covplot[[ncov+1]] <- ncov/2 # "ncol" in grid.arrange call
names(covplot)[ncov+1] <- "ncol"
do.call("grid.arrange",covplot)

###################
## Model fitting ##
###################
# Evaluate covariate gradients at observed locations
system.time(gradarray <- covGrad(xy, xgrid, ygrid, covarray))

# Fit model
selectedID <- ID %in% c(1, 2, 3)
selectedCovariate <- (1:4) %in% c(1, 2, 3, 4)
lse <- eulerLSE(ID = ID[selectedID], time = time[selectedID], xy = xy[selectedID, ], 
                gradarray = gradarray[selectedID, , selectedCovariate, drop = F], withspeed = T)

significantCovariate <- which(apply(lse$betaHat95CI, 1, function(x) diff(sign(x)) == 0))
# Compute estimated RSF
betaMLE <- lse$betaHat
if(length(significantCovariate > 0)){
    rsfRasterMLE <- Reduce("+", lapply(significantCovariate, function(i){
        betaMLE[i] * covlist[[i]]
    }))
}
rsfRasterMLE <- exp(rsfRasterMLE)

# plot estimated RSF
covmap <- data.frame(coordinates(rsfRasterMLE),
                     val = values(rsfRasterMLE) / sum(values(rsfRasterMLE)))
ggplot(covmap, aes(x,y)) + geom_raster(aes(fill = val)) +
    coord_equal() + scale_fill_viridis(name = expression(pi)) +
    geom_point(aes(x,y), data = xydf, size = 0.3)
