#############################################################################
##################################### ReTS ##################################
#############################################################################

## File: ets.plus.R

## Implements the version with variable cluster radii per dimension
## based on the very nice work by Angelov and Zhou (2012);
## Angelov, P. (2012). Evolving fuzzy systems.
##Computational Complexity: Theory, Techniques, and Applications, 1053-1065.


## DISTRIBUTED UNDER THE GNU AFFERO GENERAL PUBLIC LICENSE
## FOR MORE INFORMATION SEE license.txt or http://fsf.org

## Repository: https://triangle1.net/repo/ReTS.git
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

ReTS.cluster.etsplus <- setRefClass("ReTS.cluster.etsplus",
                                    fields = list(n = "numeric",
                                                  S = "numeric",         # support
                                                  age = "numeric",
                                                  r = "numeric",         # radius vector
                                                  prevx = "numeric",
                                                  pi = "matrix",
                                                  cv = "matrix",
                                                  prevLambda = "numeric",
                                                  sigma = "numeric",     # local scatter
                                                  M.mean = "numeric",     # temporary mean
                                                  M.2 = "numeric"
                                    ),
                                    contains = "ReTS.cluster",
                                    methods = list(
                                      update = function(zk, k) {
                                        # update potential
                                        pot <<- ((k-1)*pot)/(k-2+pot+pot*sum(ReTS.2ddist(z, zk)^2))
                                        # set previous datapoint
                                        prevx <<- zk
                                      },
                                      calcScatter = function(q) {
                                        # A numerically stable algorithm for the sample 
                                        # variance is given below. It also computes the mean. 
                                        # This algorithm is due to Knuth,[5] who cites Welford,[6]
                                        # and it has been thoroughly analyzed.
                                        # cumulative distance per dimension
                                        
                                        if (mode(q) != "numeric"){
                                          stop("Please use numeric values, not" + mode(q))
                                        }
                                        
                                        if(.self$S > 2){
                                          t.delta <- q - M.mean
                                          M.mean <<- t.delta / .self$S
                                          # make sure this is positive
                                          M.2 <<- abs(t.delta * (q - M.mean))
                                          
                                          sigma <<- sqrt(M.2 / (.self$S - 1))
                                        }
                                      },
                                      updateRadii = function(t.rho) {
                                        if (length(.self$sigma) > 0 ) {
                                          r <<- t.rho * .self$r + (1 - t.rho) * .self$sigma
                                        }
                                      },
                                      calcFiringDegree = function(q) {
                                        t.vdist <- ReTS.2ddist(.self$z, q)
                                        mu <<- exp(((-4/.self$r) * ((t.vdist)^2)))
                                        tau <<- prod(.self$mu)
                                      },
                                      updatewRLS = function(y) {
                                        t.ext <- rbind(1, matrix(prevx))# calculate extended data vector
                                        pi <<- .self$pi +  (.self$cv %*% t.ext %*% .self$lambda) * as.numeric(y - (t(t.ext) %*% .self$pi))
                                        cv <<- .self$cv - (.self$lambda*.self$cv %*% t.ext %*% t(t.ext) %*% .self$cv) /
                                          as.numeric(1+.self$lambda * t(t.ext) %*% .self$cv %*% t.ext)
                                      }
                                    )
)

ReTS.alg.etsplus <- setRefClass("ReTS.alg.etsplus",
                                fields = list(
                                  globalRadius = "numeric",
                                  nu = "matrix",
                                  dim = "numeric",
                                  Omega = "numeric",
                                  rho = "numeric"
                                ),
                                contains = 'ReTS.alg',
                                methods = list(
                                  mindist = function() {
                                    ls <- NULL
                                    for (v in 1:length(clusters)) {
                                      ls <- rbind(ls, ReTS.2ddist(clusters[[v]]$z, sample$x))
                                    }
                                    
                                    return(which.min(ls))
                                  },
                                  proximityTest = function(q) {
                                    # Compares the potential of the sample with all clusters
                                    # in the array 'clusters' and optionally calculates the  
                                    # distance to the nearest cluster.
                                    #
                                    # Arguments: current input vector zk known as q
                                    # Output:  Vector with two items;
                                    #          item1;  1 if Potential of sample is higher than
                                    #                  potential of cluster
                                    #          item2;  if 1 is satisfied and point is within range
                                    #                  of a cluster, it returns the index of the 
                                    #                  cluster. NA otherwise.
                                    if (length(q) <= 0) {
                                      stop("Please provide input variable.")
                                    } else if (!is.numeric(q)) {
                                      q <- as.numeric(q)
                                    }
                                    t.pot <- NULL
                                    t.coor <- NULL
                                    
                                    for (v in 1:length(clusters)) {
                                      t.pot <- rbind(t.pot, clusters[[v]]$pot)
                                      t.coor <- rbind(t.coor, clusters[[v]]$z)
                                    }
                                    t.dist <- ReTS.vectorDist(t.coor, q)                          
                                    
                                    if (sample$pot > max(t.pot)) {
                                      t.proximity <- max((sample$pot/max(t.pot)) - (min(t.dist)/clusters[[1]]$r))
                                      if (t.proximity > 1){
                                        # if in close proximity and higher potential
                                        # signal to change cluster center
                                        return(c(1, which.min(t.dist)))
                                      } else {
                                        # if higher potential but not in proximity
                                        # signal to create new cluster
                                        return(c(1, NA)) 
                                      } 
                                    } else {
                                      # If the potential is not higher assign
                                      # data point to cluster
                                      return(c(0, which.min(t.dist)))
                                    }
                                  },
                                  modCluster = function(id, q, t.rho) {
                                    # Desc: This actually moves the cluster center
                                    # to the new location and updates support (S)
                                    # value.
                                    if (!is.numeric(q)){
                                      q <- as.numeric(q)
                                    }
                                    # update location 
                                    clusters[[id]]$z <<- q
                                    
                                    # update support
                                    clusters[[id]]$S <<- clusters[[id]]$S + 1
                                    
                                    # update scatter
                                    clusters[[id]]$calcScatter(q)
                                    
                                    # update radius
                                    clusters[[id]]$updateRadii(t.rho)
                                    
                                  },
                                  addCluster = function(q, pp) {
                                    if (!is.numeric(q)){
                                      q <- as.numeric(q)
                                    }
                                    # calculate regression vector by taking the avg.
                                    # angelov 2k4
                                    
                                    # get weighted pi
                                    t.wpi <- 0
                                    t.totalScatter <- 0
                                    for (v in 1:length(clusters)) {
                                      t.wpi <- cbind(t.wpi, (clusters[[v]]$lambda * clusters[[v]]$pi) )
                                      
                                      if (length(clusters[[v]]$sigma) > 0){
                                        t.totalScatter <- t.totalScatter + clusters[[v]]$sigma 
                                      }
                                    }
                            
                                    t.clust <- ReTS.cluster.etsplus$new(z = q,
                                                                        r = rep( 1, .self$dim),
                                                                        pot = pp, prevx = q,
                                                                        cv = diag(dim+1)*Omega,
                                                                        tau = 0, lambda = 0, prevLambda = 0,
                                                                        mu = rep(0, dim),
                                                                        S = 1,
                                                                        #sigma = rep(1, 8),
                                                                        sigma = (t.totalScatter / length(t.totalScatter)),
                                                                        M.2 = rep(0, .self$dim),
                                                                        M.mean = rep(0, .self$dim),
                                                                        pi = as.matrix(rowSums(t.wpi)))
                                    clusters <<- c(clusters, t.clust)
                                  },
                                  normalizeFiringLvls = function() {
                                    # from angelov 2004 p. 485
                                    
                                    ls <- NULL
                                    nf <- NULL
                                    
                                    for (v in 1:length(clusters)) {
                                      ls <- c(ls, clusters[[v]]$tau)
                                    }
                                    t.sum <- sum(ls)
                                    
                                    for (v in 1:length(clusters)) {
                                      clusters[[v]]$lambda <<- ls[v]/t.sum
                                    }                                
                                  },
                                  calcNuMatrix = function(q) {
                                    t.matrix <- NULL
                                    q <- as.vector(q)
                                    q <- as.numeric(c(1, q)) #note... cbind or c?
                                    for (v in 1:length(clusters)) {
                                      t.matrix <- cbind(t.matrix, clusters[[v]]$lambda * q) 
                                    }
                                    nu <<- (as.matrix(t.matrix))
                                  }, 
                                  updatewRLS = function(y) {
                                    # updates all clusters
                                    for (v in 1:length(clusters)) {
                                      clusters[[v]]$updatewRLS(y)
                                    }
                                  },
                                  calcScatter = function(y) {
                                    # updates ALL clusters
                                    # should not be used as I believe it 
                                    # only should change when a point is 
                                    # assigned to a cluster.
                                    for (v in 1:length(clusters)) {
                                      clusters[[v]]$calcScatter(y)
                                    }
                                  },
                                  calcTheta = function() {
                                    t.tt <- NULL
                                    # fills theta var
                                    for (v in 1:length(clusters)) {
                                      t.tt <- rbind(t.tt, t(clusters[[v]]$pi))
                                    }
                                    theta <<- t.tt
                                  },
                                  init = function(t.globalRadius, samp, t.omega, t.rho){
                                    # function takes first sample and initializes 
                                    # systema
                                    rho <<- t.rho
                                    k <<- 2
                                    Omega <<- t.omega
                                    globalRadius <<- t.globalRadius
                                    t.clt = ReTS.cluster.etsplus$new(z = as.numeric(samp), r = rep(globalRadius, .self$dim),
                                                                     pot = 1,
                                                                     prevx = samp,  cv = diag(dim+1)*Omega,
                                                                     pi = matrix(0, nrow=dim+1, ncol=1),
                                                                     tau = 0, mu = rep(0, dim), lambda = 0,
                                                                     M.2 = rep(0, .self$dim),
                                                                     M.mean = rep(0, .self$dim), S = 1,
                                                                     sigma = rep(1, .self$dim))
                                    clusters <<- list(t.clt)
                                    
                                    t.smpl <- ReTS.sample$new(x = as.numeric(samp),
                                                              sigma = 0, upsilon = 0,
                                                              beta = c(0,0), nu = 0,
                                                              pot = 0)
                                    sample <<- t.smpl                                   
                                  })
)
