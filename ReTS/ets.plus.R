#############################################################################
##################################### ReTS ##################################
#############################################################################

## File: ets.plus.R

## eTS variant with evolving clusters


## DISTRIBUTED UNDER THE GNU AFFERO GENERAL PUBLIC LICENSE
## FOR MORE INFORMATION SEE license.txt or http://fsf.org

## Repository: http://triangle1.net/repo/ReTS.git
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
                                                  x.mean = "numeric",
                                                  x.sigma = "numeric"),
                                    contains = "ReTS.cluster",
                                    methods = list(
                                      update= function(zk, k) {
                                        ## update cluster intermediate values
                                        # store prev average and update mean
                                        prevmean <- x.mean
                                        x.mean <<- ((n - 1)/n)*x.mean + (zk/n)
                                        # update variance (sigma)
                                        x.sigma <<- (x.sigma + prevmean^2 - x.mean^2 +
                                                       (zk^2 - x.mean-prevmean^2)/(n+1))
                                        
                                        # update radius. NOTE: IMPLIES RHO == .5 (Eq. 17-18)
                                        r <<- (r+x.sigma)/2
                                        
                                        ## update potential 
                                        pot <<- ((k-1) * pot)/((k-1)+(k-2)*((1/pot - 1) +
                                                                              + sum(ReTS.cosdist(z, zk)) ))
                                      }
                                    )
)

ReTS.alg.etsplus <- setRefClass("ReTS.alg.etsplus",
                                fields = list(),
                                contains = 'ReTS.alg',
                                methods = list()
)
