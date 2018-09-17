
## The code below is taken from:
## Wilson, K., Hanks, E., & Johnson, D. (2018). 
## Estimating animal utilization densities using continuous‐time Markov chain models.
## Methods in Ecology and Evolution, 9(5), 1232-1240.
## (github.com/kenady/ctmcUD-MEE)

library(crawl) # devtools::install_github("NMML/crawl") for latest version
library(ctmcmove)
library(gdistance)
library(lubridate)
library(raster)
library(dplyr)
library(purrr)
library(readr)

# 1) Load the data
#------------------------------------------------
ssl_data <- readr::read_csv("sea_lion_telemetry.csv") %>% mutate(GMT=mdy_hm(GMT)) # telemetry data
aleut_hab <- brick("habitat/aleut_habitat.grd", values=TRUE) %>% stack()   # habitat covariates
land <- 1.0*(aleut_hab$bathy>=0)

# Project telemetry data
dat_toproj <- filter(ssl_data, !is.na(latitude))
coordinates(dat_toproj) <- ~longitude+latitude
proj4string(dat_toproj) <- CRS("+proj=longlat")
dat_toproj <- spTransform(dat_toproj, CRS("+init=epsg:3338")) %>% 
  as.data.frame() %>% rename(x=longitude, y=latitude) %>% 
  dplyr::select(Deploy_ID, GMT, x, y)
ssl_data = full_join(ssl_data, dat_toproj, by=c("Deploy_ID", "GMT"))

###################################
# Fit CTMC UD to first animal    ##
###################################
# temp_dat <- dplyr::filter(ssl_data, Deploy_ID==14809) %>% arrange(GMT) #pull out one sea lion
temp_dat <- ssl_data %>% arrange(GMT)

#set Argos error & add a constraint parameter for estimating Argos error
# in CTCRW imputation model
temp_dat$argos_class = factor(temp_dat$argos_class, levels=c("3","2","1","0","A","B"))

# Define 'haul' variable for determining when the animal is hauled out. These will
# be removed in the CTMC analysis
temp_dat %>%  mutate(
                  haul = (
                    DryTime==1 & c(0,diff(DryTime==1)==0) |
                      DryTime==0 & c(0,diff(DryTime==1)==-1)
                    )
                  ) -> temp_dat

## Initial state values from CTCRW model in crawl (for imputation of raw telemetry)
initial <- list(a=c(temp_dat$x[1],0, temp_dat$y[1],0),
                P=diag(c(10000^2,5400^2,10000^2,5400^2)))
# Fixed parameter values for CTCRW imputation model
fixPar = c(log(250), log(500), log(1500), rep(NA,5), 0) 
# Lower bounds for parameter estimates
constr = list(lower=c(rep(log(1500),3),rep(-Inf,2)), upper=rep(Inf,5))
# Check that parameters are as expected
crawl::displayPar( mov.model=~1, err.model=list(x=~argos_class-1),data=temp_dat,
            activity=~I(1-DryTime),fixPar=fixPar) #with dry time

# 2) Run crawl model (Fit CTCRW model)
# --------------------------------------------------------------------
set.seed(123)
mi_fit <- crawl::crwMLE(
  mov.model=~1, err.model=list(x=~argos_class-1), activity=~I(1-DryTime),
  data=temp_dat, coord=c("x","y"), Time.name="GMT",
  initial.state=initial, fixPar=fixPar,
  constr=constr, #prior=ln.prior,
  method="L-BFGS-B",
  control=list(maxit=2000, trace=1, REPORT=10),
  initialSANN=list(maxit=250, temp=10, trace=1, REPORT=10))

# 3) Predict locations at times of observations
# --------------------------------------------------------------------
predloc <- crwPredict(mi_fit)
whichObs <- which(!is.na(temp_dat$longitude))
predloc <- predloc[whichObs,]

preddat <- data.frame(ID = predloc$Deploy_ID, 
                      time = temp_dat$GMT[whichObs],
                      x = predloc$mu.x,
                      y = predloc$mu.y)

# write.csv(preddat, file="~/git/LangevinMovementModel/sealion_analysis/SSLpreddat.csv", row.names=FALSE)
