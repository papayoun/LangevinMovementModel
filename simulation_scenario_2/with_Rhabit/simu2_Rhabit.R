
# Load packages
library(raster)
library(ggplot2)
library(viridis)
library(Rhabit)
library(parallel)
set.seed(1)

#######################
## Define covariates ##
#######################
# Generate two random covariates
lim <- c(-1, 1, -1, 1)*50
resol <- 1
ncov <- 2
covlist <- list()
for(i in 1:ncov) {
    covlist[[i]] <- simSpatialCov(lim = lim, nu = 5, rho = 10, sigma2 = 0.1, 
                                  resol = resol, raster_like = TRUE)
}

# Compute utilisation distribution
beta <- c(4,2)
UD <- getUD(covariates=covlist, beta=beta)
UDrast <- rasterFromXYZ(rasterToGGplot(UD))

###################
## Simulate data ##
###################
Tmax <- 250
dt <- 0.01
time <- seq(0, Tmax, by = dt)
ntrack <- 200
speed <- 0.5

t0 <- Sys.time()

# Check whether data file already exists
if(!file.exists("LMMsimu.RData")) {
    # Time grids
    alltimes <- list()
    for(i in 1:ntrack)
        alltimes[[i]] <- time
    
    # Generate tracks
    alldat <- mclapply(alltimes, function(times) {
        simLangevinMM(beta = beta, gamma2 = speed, times = times, 
                      loc0 = runif(2, -20, 20), cov_list = covlist)
    }, mc.cores = 3)
    
    # Add ID column
    for(zoo in 1:ntrack)
        alldat[[zoo]] <- cbind(ID = rep(zoo, length(time)), alldat[[zoo]])
    
    save(alldat, file="LMMsimu.RData")    
} else {
    load("LMMsimu.RData")
}

t1 <- Sys.time()
print(t1 - t0)

##############
## Estimate ##
##############
thin <- c(1,2,5,10,25,50,100)
allfits <- list()
for(i in 1:length(thin)) {
    cat("Iteration",i,"/",length(thin),"\n")
    
    # Thin and keep 250*200 observations
    thinDatList <- lapply(alldat, function(dat) {
        ind <- thin[i]*(1:250)
        dat[ind,]
    })
    # Stack all tracks into one data frame
    thinDat <- as.matrix(do.call(rbind.data.frame, thinDatList))
    
    ID <- thinDat[,1]
    locs <- thinDat[,2:3]
    times <- thinDat[,4]

    # Derive covariate gradients at observed locations
    gradarray <- covGradAtLocs(locs = locs, cov_list = covlist)

    # Fit Langevin model
    allfits[[i]] <- langevinUD(locs=locs, times=times, ID = ID, grad_array=gradarray)
}

####################
## Plot estimates ##
####################
beta1 <- data.frame(dt = thin*dt,
                    est = unlist(lapply(allfits, function(fit) fit$betaHat[1])),
                    lci = unlist(lapply(allfits, function(fit) fit$betaHat95CI[1,1])),
                    uci = unlist(lapply(allfits, function(fit) fit$betaHat95CI[1,2])))
beta2 <- data.frame(dt = thin*dt,
                    est = unlist(lapply(allfits, function(fit) fit$betaHat[2])),
                    lci = unlist(lapply(allfits, function(fit) fit$betaHat95CI[2,1])),
                    uci = unlist(lapply(allfits, function(fit) fit$betaHat95CI[2,2])))
gamma2 <- data.frame(dt = thin*dt,
                     est = unlist(lapply(allfits, function(fit) fit$gamma2Hat)),
                     lci = unlist(lapply(allfits, function(fit) fit$betaHat95CI[3,1])),
                     uci = unlist(lapply(allfits, function(fit) fit$betaHat95CI[3,2])))

ggplot(beta1, aes(dt, est)) + geom_point() + scale_x_continuous(trans="log", breaks=beta1$dt) +
    geom_segment(aes(x=dt, y=lci, xend=dt, yend=uci)) +
    geom_hline(yintercept = c(0, beta[1]), lty=2, col=c(1,2)) +
    xlab("Time interval") + ylab(expression(beta[1]))

ggplot(beta2, aes(dt, est)) + geom_point() + scale_x_continuous(trans="log", breaks=beta2$dt) +
    geom_segment(aes(x=dt, y=lci, xend=dt, yend=uci)) +
    geom_hline(yintercept = c(0, beta[2]), lty=2, col=c(1,2)) +
    xlab("Time interval") + ylab(expression(beta[2]))

ggplot(gamma2, aes(dt, est)) + geom_point() + scale_x_continuous(trans="log", breaks=gamma2$dt) +
    geom_segment(aes(x=dt, y=lci, xend=dt, yend=uci)) +
    geom_hline(yintercept = speed, lty=2, col=2) + coord_cartesian(ylim=c(0.49,0.51)) +
    xlab("Time interval") + ylab(expression(gamma^2))
