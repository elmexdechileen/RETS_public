#############################################################################
##################################### ReTS ##################################
#############################################################################

## File: dataLoader.R

## Dataloader written to load standardized datasets


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


data.radioFreq <- function(timeSpan, predPeriod) {
  y <- read.csv("datasets/monthly-mean-discharge-in-cubic.csv")
  names(y) <- c("index", "y")
  y <- as.integer(y$y)
  
  # normalize
  y <- (y - min(y, na.rm = TRUE)) /( max(y, na.rm = TRUE)-min(y, na.rm = TRUE))
  ori.y <- y[1:timeSpan] # store original y
  
  x18 <- head(y, -(17 + predPeriod))
  x12 <- head(y, -(11 + predPeriod))
  x6 <- head(y, -(5 + predPeriod))
  x1 <- head(y, -(1 + predPeriod))
  
  x12 <- tail(x12, -(6))
  x6 <- tail(x6, -(12))
  x1 <- tail(x1, -(16))
  y <- tail(y, -(17 + predPeriod))
  
  tf2 <- cbind(x18, x12, x6, x1, y)
  tf2 <- as.data.frame(na.omit(tf2))
  return(as.data.frame(na.omit(tf2)))
}

data.riverflow <- function(timeSpan, predPeriod) {
#   @book{hipel1994time,
#         title={Time series modelling of water resources and environmental systems},
#         author={Hipel, Keith W and McLeod, A Ian},
#         year={1994},
#         publisher={Elsevier}
#   }
#   Monthly riverflow in cms, Oostanaula River at resaca, GA, 1893 – 1960
  
  
  y <- read.csv("datasets/monthly-riverflow-in-cms-oostana.csv")
  names(y) <- c("index", "y")
  y <- as.integer(y$y)
  
  # normalize
  y <- (y - min(y, na.rm = TRUE)) /( max(y, na.rm = TRUE)-min(y, na.rm = TRUE))
  ori.y <- y[1:timeSpan] # store original y
  
  x18 <- head(y, -(17 + predPeriod))
  x12 <- head(y, -(11 + predPeriod))
  x6 <- head(y, -(5 + predPeriod))
  x1 <- head(y, -(1 + predPeriod))
  
  x12 <- tail(x12, -(6))
  x6 <- tail(x6, -(12))
  x1 <- tail(x1, -(16))
  y <- tail(y, -(17 + predPeriod))
  
  tf2 <- cbind(x18, x12, x6, x1, y)
  tf2 <- as.data.frame(na.omit(tf2))
  return(as.data.frame(na.omit(tf2)))
}


data.MackeyGlass <- function(timeSpan, predPeriod, zero, sine, jitter) {
  mackeyGlass <- function(t, y, parms, tau) {
    tlag <- t - tau
    
    if (tlag <= 0)
      ylag <- 0.5
    else
      ylag <- lagvalue(tlag)
    
    dy <- 0.2 * ylag * 1/(1+ylag^10) - 0.1 * y
    list(c(dy))
  }
  
  tmax = timeSpan * 1.5
  yinit = 1.2
  times <- seq(from = 0, to = tmax, by = 1)
  y <- dede(y = yinit, times = times, func = mackeyGlass, parms = NULL, tau = 17)
  
  y <- y[seq(from=1, to=nrow(y), by=1), 2]
  
  if(sine) {
    t.t <- seq(from = 1, to = length(y), by = 1)
    t.sine <- sin(t.t/500)*.2
    y <- y + t.sine
  }
  
  # normalize
  y <- (y - min(y, na.rm = TRUE)) /
    ( max(y, na.rm = TRUE)-min(y, na.rm = TRUE))
  
  if(zero) {
    y <- pmax(0, (y - mean(y)))
  }
  
  if(jitter) {
    y <- y + runif(length(y), min = -.1, max = .1)
  }
  
  
  ori.y <- y[1:timeSpan] # store original y
  
  x18 <- head(y, -(17 + predPeriod))
  x12 <- head(y, -(11 + predPeriod))
  x6 <- head(y, -(5 + predPeriod))
  x1 <- head(y, -(1 + predPeriod))
  
  x12 <- tail(x12, -(6))
  x6 <- tail(x6, -(12))
  x1 <- tail(x1, -(16))
  y <- tail(y, -(17 + predPeriod))
  
  tf2 <- cbind(x18, x12, x6, x1, y)
  tf2 <- as.data.frame(tf2[1:timeSpan,])
}

data.oversample <- function(tf, n) {
  index <- c(1:length(tf[,1])) * n
  tf$index <- index
  t.index <- data.frame("index" = c(1:(length(tf[,1])*10)))
  t.tf <- merge(t.index, tf, by = c("index"), all.x = TRUE)
  t.tf$PMI <- na.approx(t.tf$PMI, na.rm="FALSE")
  t.tf$EURIBOR <- na.approx(t.tf$EURIBOR, na.rm="FALSE")
  t.tf$PHLX <- na.approx(t.tf$PHLX, na.rm="FALSE")
  t.tf$ASML <- na.approx(t.tf$ASML, na.rm="FALSE")
  t.tf$pred <- na.approx(t.tf$pred, na.rm="FALSE")
  t.tf$laggedActual <- na.approx(t.tf$laggedActual, na.rm="FALSE")
  t.tf$actual <- na.approx(t.tf$actual, na.rm="FALSE")
  t.tf <- t.tf[10:length(t.tf[,1]),]
  #t.tf$index <- c(1:length(t.tf[,1]))
  return(t.tf)
}

data.squareWave <- function(n, predPeriod) {
  y <- rep(c(rep(0,50), rep(1,50)), n)
  x18 <- head(y, -(17 + predPeriod))
  x12 <- head(y, -(11 + predPeriod))
  x6 <- head(y, -(5 + predPeriod))
  x1 <- head(y, -(1 + predPeriod))
  
  x12 <- tail(x12, -(6))
  x6 <- tail(x6, -(12))
  x1 <- tail(x1, -(16))
  y <- tail(y, -(17 + predPeriod))
  
  tf2 <- cbind(x18, x12, x6, x1, y)
  tf2 <- as.data.frame(na.omit(tf2))
  return(as.data.frame(na.omit(tf2)))
}
