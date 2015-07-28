#############################################################################
##################################### ReTS ##################################
#############################################################################

## File: main.R

## INFO: Object based implementation of
## an eTS by Angelov & Filev (2004) in
## R, hence ReTS.
## TYPE: on-line, evolving, unsupervised


## DISTRIBUTED UNDER THE GNU AFFERO GENERAL PUBLIC LICENSE
## FOR MORE INFORMATION SEE license.txt or http://fsf.org

## Repository: http://triangle1.net/repo/ReTS.git/.git
## Author: Max van Rooijen (max@triangle1.net)
## Version: 1.0

# Copyright (C) 2015  Max van Rooijen
#   
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as
# published by the Free Software Foundation, either version 3 of the
# License, or (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Affero General Public License for more details (included
# as license.txt).


########### Misc requires ############
require(ggplot2)
require(deSolve)
require(ggthemes)
require(frbs)
require(forecast)
require(astsa)
require(RJDBC)
require(plyr)
require(dplyr)
require(ggplot2)
require(lubridate)
require(zoo) # for interpolation
require(lmtest)
require(Hmisc)
require(stargazer)
require(gridExtra)

########### Sourced files ############
source("ReTS/miscFunctions.R")
source("ReTS/abstracts.R")
source("ReTS/ets.plus.R")
source("ReTS/ets.basic.R")
source("ReTS/dataLoader.R")
source("ReTS/theming.R")

################ DATA ################
timeSpan <- 4000
predPeriod <- 50

# examples
# tf <- data.radioFreq(1000, predPeriod)
# tf <- data.MackeyGlass(timeSpan, predPeriod, TRUE, FALSE, FALSE) # input; timespan, periods to predict ahead,
                                                                 #        0 demand periods, sine wave, random noise
# tf <- data.squareWave(20,8)
tf <- data.riverflow(815, 2)

############## Implementation ##############
ptm <- proc.time()

#  init
alg <- ReTS.alg.ets$new(dim = (length(tf[1,]) - 1))
alg$init(.3, as.numeric(tf[1, 1:(length(tf[1,]) - 1)]), 5) #last value is omega

pb <- txtProgressBar(1, length(tf[,1]), style = 3) # set progressbar 
for (i in 2:length(tf[,1])) {
  # set kW
  alg$k <- i
  
  # load next sample
  x.dat <- tf[i,1:alg$dim]
  y.dat <- tf[i,(alg$dim+1)]
  
  # update sample
  alg$sample$update(x.dat, i)
  
  # update point potential var for plotting
  pointpotential <- c(pointpotential, alg$sample$pot)
  
  # update cluster potentials
  alg$updateClusters(x.dat)
  
  # Test potential values and proximity
  g <- alg$proximityTest(x.dat)
  if (g[1] != 0 && !is.na(g[2])) {
    # cluster in close proximity of previous cluster
    # so update existing cluster
    alg$modCluster(g[2], alg$sample$x)
  } else if (g[1] != 0){
    # add new cluster and update consequent part
    alg$addCluster(x.dat, alg$sample$pot)
  }
  
  alg$updateFiringDegrees(x.dat)
  alg$normalizeFiringLvls()
  alg$calcNuMatrix(x.dat)  
  alg$updatewRLS(y.dat)
  alg$calcTheta()

  # update progress
  setTxtProgressBar(pb, i)
}

# close pb
close(pb)
proc.time() - ptm
############## end implementation ##############
