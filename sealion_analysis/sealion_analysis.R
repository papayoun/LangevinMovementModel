
library(raster)
library(numDeriv)
library(fields)
library(expm)
library(viridis)
library(ggplot2)
library(gridExtra)
source("../utility.R")
source("../eulerLSE.R")

####################
## Prepare tracks ##
####################
# load track from Github URL
baseURL <- "https://raw.githubusercontent.com/kenady/ctmcUD-MEE/master/"
tracksLL <- read.csv(url(paste0(baseURL,"sea_lion_telemetry.csv")))
tracksLL <- subset(tracksLL, !is.na(longitude) & !is.na(latitude))

utmcoord <- SpatialPoints(tracksLL[,c("longitude","latitude")], 
                          proj4string=CRS("+proj=longlat"))
llcoord <- spTransform(utmcoord, CRS("+init=epsg:3338"))
tracks <- data.frame(ID = tracksLL$Deploy_ID,
                     x = attr(llcoord,"coords")[,1]/1000,
                     y = attr(llcoord,"coords")[,2]/1000)

ID <- as.integer(tracks$ID)
ID[ID==14809] <- 1
ID[ID==15136] <- 2
ID[ID==15137] <- 3
xy <- matrix(c(tracks$x,tracks$y),ncol=2)
time <- as.POSIXct(tracksLL$GMT, format="%m/%d/%Y %H:%M", tz="GMT")
time <- as.numeric(time)
time <- (time-min(time))/3600

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

# negative values in Fig 2 of Wilson et al.
values(covlist$slope) <- -values(covlist$slope)
values(covlist$d2site) <- -values(covlist$d2site)
values(covlist$d2shelf) <- -values(covlist$d2shelf)

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
gradarray <- covGrad(xy, xgrid, ygrid, covarray)

# Fit model
lse <- eulerLSE(ID=ID, time=time, xy=xy, gradarray=gradarray)

# 95% CI
cbind(lse$est[1:ncov]-1.96*sqrt(diag(lse$var)),
      lse$est[1:ncov],
      lse$est[1:ncov]+1.96*sqrt(diag(lse$var)))

# Compute estimated RSF
betaMLE <- lse$est[1:ncov]
rsfRasterMLE <- 0
for(i in 1:length(covlist))
    rsfRasterMLE <- rsfRasterMLE + betaMLE[i]*covlist[[i]]
rsfRasterMLE <- exp(rsfRasterMLE)

# plot estimated RSF
covmap <- data.frame(coordinates(rsfRasterMLE),
                     val=values(rsfRasterMLE)/sum(values(rsfRasterMLE)))
ggplot(covmap,aes(x,y)) + geom_raster(aes(fill=log(val))) +
    coord_equal() + scale_fill_viridis(name=expression(log(pi))) +
    geom_point(aes(x,y), data=xydf, size=0.3)
