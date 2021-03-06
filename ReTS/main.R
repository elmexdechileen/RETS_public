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

## Repository: https://triangle1.net/repo/ReTS.git/
## Author: Max van Rooijen (max@triangle1.net)
## Version: 1.1

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
require(plyr)
require(dplyr)
require(zoo) # for interpolation

########### Sourced files ############
source("ReTS/miscFunctions.R")
source("ReTS/abstracts.R")
source("ReTS/ets.plus.R")
source("ReTS/ets.basic.R")
source("ReTS/dataLoader.R")

################ DATA ################
timeSpan <- 4000
predPeriod <- 50

# examples
# tf <- data.radioFreq(1000, predPeriod)
 tf <- data.MackeyGlass(timeSpan, predPeriod, TRUE, FALSE, FALSE) # input; timespan, periods to predict ahead,
                                                                 #        0 demand periods, sine wave, random noise
# tf <- data.squareWave(20,8)
tf <- data.riverflow(815, 2)

############## Implementation ##############
ptm <- proc.time()

#  init
prediction <- 0
alg <- ReTS.alg.etsplus$new(dim = (length(tf[1,]) - 1))
alg$init(1,
         as.numeric(tf[1, 1:(length(tf[1,]) - 1)]),
         5,
         0.5)

pb <- txtProgressBar(1, length(tf[,1]), style = 3) # set progressbar 

for (i in 2:length(tf[,1])) {
  # set kW
  alg$k <- i
  
  # load next sample
  x.dat <- tf[i,1:alg$dim]
  y.dat <- tf[i,(alg$dim+1)]

  # update sample
  alg$sample$update(x.dat, i)
  
  # update cluster potentials
  alg$updateClusters(x.dat)
  
  # Test potential values and proximity
  g <- alg$proximityTest(x.dat)
  if (g[1] == 1 && !is.na(g[2])) {
    # cluster in close proximity of previous cluster
    # so update existing cluster
    alg$modCluster(g[2], alg$sample$x, alg$rho)
  } else if (g[1] == 1){
    # add new cluster and update consequent part
    alg$addCluster(x.dat, alg$sample$pot)
  } else if (g[1] == 0 && !is.na(g[2])) {
    # if inclose proximity only assign point as support
    # for closest cluster. So mod cluster by assigning
    # existing cluster center to update (bad code though)
    alg$modCluster(g[2], alg$clusters[[g[2]]]$z, alg$rho)
  }
  
  
  # update consequent
  alg$updateFiringDegrees(x.dat)
  alg$normalizeFiringLvls()
  alg$calcNuMatrix(x.dat)  
  alg$updatewRLS(y.dat)
  alg$calcTheta()

  # store prediction
  prediction <- rbind(prediction,
                      c(i, sum(t(alg$nu)*alg$theta)))
  
  # update progress
  setTxtProgressBar(pb, i)
}

# close pb and print processing time
close(pb)
proc.time() - ptm
############## end implementation ##############
