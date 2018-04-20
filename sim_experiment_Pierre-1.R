
#' Ce script génère deux covariables aléatoires, et calcule la distance au centre (au carré)
#' comme autre covariable. Tu peux les visualiser avec plot(covlist[[i]]) (pour i=1,2,3).
#' Ensuite, la distribution d'utilisation est définie par une RSF, i.e. exp(beta%*%covariables),
#' où je choisis beta = (2,4,-6). Tu peux la visualiser avec plot(rsfRaster).
#' Je simule 1000 positions avec un schéma d'Euler, sur cette distribution d'utilisation. Finalement,
#' j'utilise la pseudo-vraisemblance basée sur le schéma d'Euler, au même pas de temps que celui de
#' la simulation, pour estimer les paramètres du modèle. Donc, en l'occurrence, j'utilise la "vraie"
#' vraisemblance du processus générateur, ce qui devrait nous permettre de retrouver les paramètres
#' sans problème. Pourtant, c'est pas vraiment le cas.

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
beta <- c(3,4,-.7) # true parameter values (used in simulations)
rsfRaster <- 0
for(i in 1:length(covlist))
    rsfRaster <- rsfRaster + beta[i]*covlist[[i]]
rsfRaster <- exp(rsfRaster)

plot(log(rsfRaster))

###################
## Simulate data ##
###################
ntrack <- 2
dt <- 1
Tmax <- 500
time <- seq(0,Tmax,by=dt)
set.seed(1)
xy <- NULL
ID <- rep(1:ntrack,each=Tmax/dt+1)
for(zoo in 1:ntrack) {
    track <- simLang(beta=beta, time=time, xy0=runif(2,-10,10), 
                     xgrid=xgrid, ygrid=ygrid, covarray=covarray)
    xy <- rbind(xy, track)
}
time <- rep(time,ntrack)

plot(log(rsfRaster))
points(xy,type="o",cex=0.4)

###########################
## Fit with Euler method ##
###########################
# Evaluate covariate gradients at observed locations
gradarray <- covGrad(xy, xgrid, ygrid, covarray)
# Optimise log-likelihood
fit <- nlminb(start = beta, objective = nllkLang, xy = xy, time = time, ID=ID,
              gradarray = gradarray, control = list(trace=1))

# Compute estimated RSF
betaMLE <- fit$par
rsfRasterMLE <- 0
for(i in 1:length(covlist))
    rsfRasterMLE <- rsfRasterMLE + betaMLE[i]*covlist[[i]]
rsfRasterMLE <- exp(rsfRasterMLE)

# Plot estimated utilisation vs true utilisation for each grid cell 
# (should align with identity line)
plot(values(rsfRaster)/sum(values(rsfRaster)),
     values(rsfRasterMLE)/sum(values(rsfRasterMLE)))
abline(0,1,col=2)


#############################################
## Fit with Euler method using lm approach ##
#############################################
d <- dim(gradarray)[3]
n <- dim(gradarray)[1]

i0 <- which(ID[-1]!=ID[-n])
i2 <- c(i0, n) # last obs of each track

##  design matrix construction
X <- 0.5*rbind(gradarray[-i2,1,], gradarray[-i2,2,])
TT <- diag(rep(diff(time)[-i0],2))

TZ <- matrix(apply(xy, 2, function(c) diff(c)[-i0]), ncol=1)

## estimation
XTTX    <- t(X) %*% TT %*% X
XTTXinv <- solve(XTTX)
(Bhat  <- XTTXinv %*% t(X) %*% TZ)

diag(XTTXinv)

# 95% CI
cbind(Bhat-1.96*sqrt(diag(XTTXinv)),Bhat,Bhat+1.96*sqrt(diag(XTTXinv)))
