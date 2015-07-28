#############################################################################
##################################### ReTS ##################################
#############################################################################

## File: miscFunctions.R

## Misc. functions (mainly distance functions)

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

ReTS.cosdist <- function(a, b) {
  # function to calculate the cosine distance
  # between two points or vectors
  # a higher number indicates higher correlation
  return(1-(sum(a*b)/sqrt(sum(a^2)*sum(b^2))))
}

ReTS.dist <- function(x1, x2) {
  sqrt(sum((x1 - x2) ^ 2))
} 

ReTS.vectorDist <- function(vect, x) {
  t.rt <- NULL
  if (length(vect[1,]) != length(x)) {
    stop("Please use input arguments of same length!")
  }
  for (v in 1:length(vect[,1])) {
    t.rt <- rbind(t.rt, sqrt(sum((vect[v,] - x) ^ 2)))
  }
  return(t.rt)
} 


ReTS.2ddist <- function(x, y) {
  t.rt <- NULL
  if (length(x) != length(y)) {
    stop("Please use input arguments of same length!")
  }
  for (v in 1:length(x)) {
    t.rt <- c(t.rt, abs(x[v]-y[v]))
  }
  return(t.rt)
} 
