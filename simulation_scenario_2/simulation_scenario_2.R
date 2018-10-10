
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

############################
## Plot covariates and UD ##
############################
ggopts <- theme(axis.title = element_text(size=20), axis.text = element_text(size=20), 
                legend.title = element_text(size=25), legend.text = element_text(size=20), 
                legend.key.height=unit(3,"line"))

c1map <- as.data.frame(cbind(coordinates(covlist[[1]]), val=values(covlist[[1]])))
c1plot <- ggplot(c1map, aes(x,y)) + geom_raster(aes(fill=val)) +
    coord_equal() + scale_fill_viridis(name=expression(c[1])) + ggopts

c2map <- as.data.frame(cbind(coordinates(covlist[[2]]), val=values(covlist[[2]])))
c2plot <- ggplot(c2map, aes(x,y)) + geom_raster(aes(fill=val)) +
    coord_equal() + scale_fill_viridis(name=expression(c[2])) + ggopts

UDmap <- as.data.frame(cbind(coordinates(rsfRaster),
                             val=values(rsfRaster)/sum(values(rsfRaster))))
UDplot <- ggplot(UDmap, aes(x,y)) + geom_raster(aes(fill=val)) +
    coord_equal() + scale_fill_viridis(name=expression(pi)) + ggopts

pdf("sim2c1.pdf", width = 7, height = 6)
plot(c1plot)
dev.off()
pdf("sim2c2.pdf", width = 7, height = 6)
plot(c2plot)
dev.off()
pdf("sim2UD.pdf", width = 7, height = 6)
plot(UDplot)
dev.off()

###################
## Simulate data ##
###################
Tmax <- 250
dt <- 0.01
ntrack <- 200
speed <- 0.5
time <- seq(0, Tmax, by = dt)
alldata <- NULL

set.seed(1)
zoo <- 1
t0 <- Sys.time()
while(zoo <= ntrack) {
    cat("Simulating track ",zoo,"/",ntrack,"...\n",sep="")
    tryCatch({
        track <- simLang(beta=beta, speed=speed, time=time, xy0=runif(2, -20, 20), 
                         xgrid=xgrid, ygrid=ygrid, covarray=covarray)
        alldata <- rbind(alldata, cbind(rep(zoo,nrow(track)),track,time))
        zoo <- zoo + 1
    }, error = function(e) {
        cat("Simulated track reached limit of map. Started new track.\n")
    })
    cat("Time elapsed:",Sys.time()-t0,"\n")
}
colnames(alldata) <- c("ID","x","y","time")

##############
## Estimate ##
##############
thin <- c(1,2,5,10,25,50,100)
allpar <- matrix(NA, length(thin), 3)
allvar <- matrix(NA, length(thin), 2)
allCI <- list()
for(i in 1:length(thin)) {
    cat("Iteration",i,"/",length(thin),"\n")
    # thin
    ind <- seq(1, nrow(alldata), by=thin[i])
    thindata <- alldata[ind,]
    
    # keep 50000 locations
    data <- NULL
    for(id in unique(thindata[,"ID"]))
        data <- rbind(data, thindata[which(thindata[,"ID"]==id)[1:250],])
    
    ID <- data[,"ID"]
    xy <- data[,c("x","y")]
    time <- data[,"time"]
    
    # derive covariate gradients at observed locations
    gradarray <- covGrad(xy, xgrid, ygrid, covarray)

    # Fit UD with Euler discretization
    lse <- eulerLSE(ID=ID, time=time, xy=xy, gradarray=gradarray)
    allpar[i,] <- c(lse$betaHat, lse$gammaHat)
    allvar[i,] <- diag(lse$betaHatCovariance)
    allCI[[i]] <- lse$betaHat95CI
    
    # # Fit by numerical MLE
    # par0 <- c(0, 0, log(1))
    # mod <- optim(par=par0, fn=nllkLang, xy=xy, time=time, ID=ID, gradarray=gradarray,
    #              control=list(trace=1), hessian=TRUE)
    # 
    # allpar[i,] <- c(mod$par[1:2], exp(mod$par[3]))
    # allvar[i,] <- diag(solve(mod$hessian))
}

# Plot beta estimates
pdf("sim2beta.pdf", width=8, height=4)
par(mfrow=c(1,2), mar=c(5,4,1,1)+0.1)

plot(thin*dt, allpar[,1], log="x", xlab="interval", ylab=expression(beta[1]), ylim=c(0,6))
lCI <- allpar[,1] - 1.96*sqrt(allvar[,1])
uCI <- allpar[,1] + 1.96*sqrt(allvar[,1])
segments(x0=thin*dt, y0=lCI, x1=thin*dt, y1=uCI)
abline(h=0, lty=2)
abline(h=2, lty=2, col=2)

plot(thin*dt, allpar[,2], log="x", xlab="interval", ylab=expression(beta[2]), ylim=c(0,6))
lCI <- allpar[,2] - 1.96*sqrt(allvar[,2])
uCI <- allpar[,2] + 1.96*sqrt(allvar[,2])
segments(x0=thin*dt, y0=lCI, x1=thin*dt, y1=uCI)
abline(h=0, lty=2)
abline(h=4, lty=2, col=2)
dev.off()

# Plot gamma estimates
pdf("sim2gamma.pdf", width=4, height=4)
par(mfrow=c(1,1), mar=c(5,5,1,1)+0.1)
plot(thin*dt, allpar[,3], log="x", xlab="interval", ylab=expression(gamma^2), ylim=c(0.48,0.52))
abline(h=0.5, lty=2, col=2)
dev.off()

##############################
## Irregular time intervals ##
##############################
set.seed(1)
# keep 50000 locations
thin <- 50 # thinning factor
data <- NULL
for(id in unique(alldata[,"ID"])) {
    ind <- sort(sample(1:(thin*250), size=250))
    data <- rbind(data, alldata[which(alldata[,"ID"]==id)[ind],])
}

ID <- data[,"ID"]
xy <- data[,c("x","y")]
time <- data[,"time"]

# derive covariate gradients at observed locations
gradarray <- covGrad(xy, xgrid, ygrid, covarray)

# Fit UD with Euler discretization
lse <- eulerLSE(ID=ID, time=time, xy=xy, gradarray=gradarray)
