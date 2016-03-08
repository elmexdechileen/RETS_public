#############################################################################
##################################### ReTS ##################################
#############################################################################

## File: abstracts.R
## Contains abstracts for clusters and samples


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

ReTS.impl <- setRefClass("ReTS.impl",
                         fields = list(dataset = "data.frame"),
                         methods = list(
                           sequence = function(){
                             return(NULL)
                           })
)


## BEGIN Abstract classes ##
ReTS.sample <- setRefClass("ReTS.sample",
                           fields = list(x = "numeric", # note this is z_(k-1)
                                         sigma = "numeric",
                                         upsilon = "numeric",
                                         beta = "numeric",
                                         nu = "numeric",
                                         pot = "numeric"),
                           methods = list(
                             update = function(zk, k) {
                               if(!is.numeric(zk)){
                                 zk <- as.numeric(zk)
                               }
                               upsilon <<- sum(zk^2)
                               sigma <<- sigma + upsilon
                               beta <<- beta + zk
                               nu <<- sum(x * beta)
                               x <<- zk
                               pot <<- (k-1)/((k-1)*(upsilon+1)+sigma-2*nu)                           
                               
                             })
)


ReTS.alg <- setRefClass("ReTS.alg",
                        fields = list(k = "numeric",
                                      theta = "matrix",
                                      clusters = "list", # or ReTS.cluster?
                                      sample = "ReTS.sample"),
                        methods = list(
                          init = function(zk) {
                            
                          },
                          updateClusters = function(q){
                            # test and convert if not numeric
                            if(!is.numeric(q)){
                              q <- as.numeric(q)
                            }
                            
                            for (v in 1:length(clusters)) {
                              clusters[[v]]$update(q, .self$k)
                            }
                          },
                          updateFiringDegrees = function(q) {
                            # test and convert if not numeric
                            if(!is.numeric(q)){
                              q <- as.numeric(q)
                            }
                            
                            for (v in 1:length(clusters)) {
                              clusters[[v]]$calcFiringDegree(q)
                            }
                          }))

ReTS.cluster <- setRefClass("ReTS.cluster",
                            fields = list(z = "numeric",
                                          r = "numeric",
                                          pot = "numeric",
                                          mu = "numeric", # firing degree, 
                                          tau = "numeric",
                                          lambda = "numeric"),
                            methods = list()
)

## END Abstract classes ##