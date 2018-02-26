
##' This script generates two raster covariates
library(gstat)
library(raster)
set.seed(123)

# spatial grid
xy <- expand.grid(seq(-49.5,49.5,1), seq(-49.5,49.5,1))
colnames(xy) <- c("x","y")

# define and predict two resource distributions
field <- gstat(formula=z~1, locations=~x+y, dummy=TRUE, beta=1, 
               model=vgm(psill=100, range=100, model='Exp'), nmax=20)
cov <- predict(field, newdata=xy, nsim=2)

# scale to be in [0,1]
cov[,3] <- (cov[,3]-min(cov[,3]))/(max(cov[,3])-min(cov[,3]))
cov[,4] <- (cov[,4]-min(cov[,4]))/(max(cov[,4])-min(cov[,4]))

# transform to rasters
cov1rast <- rasterFromXYZ(cov[,c("x","y","sim1")])
cov2rast <- rasterFromXYZ(cov[,c("x","y","sim2")])
