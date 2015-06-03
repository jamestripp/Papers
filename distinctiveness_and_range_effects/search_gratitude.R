# File:           search_gratitude.R
# Version:        1
# Last changed:   Wednesday 3rd June 2015
# Purpose:        Grid search of sdbs, sgems and sdbs + range 
# Author:         Dr James Tripp
# Copyright:      (C) James Tripp 2015

rm(list = ls())

#import packages
library(foreach)
library(doParallel)
library(boot)
options("scipen"=100, "digits"=4)
#setup parallel backend to use 8 processors
cl<-makeCluster(8)
registerDoParallel(cl)

#read data
data <- read.csv(file = 'data_gratitude.csv')
names(data)[1] <- 'ppt'
ppt.max <- length(unique(data$ppt)) #run analysis up to this participant

data.1 <- data.frame(ppt=data$ppt, distribution= data$distribution)
data.1 <- cbind(data.1, data[,12:(ncol(data)-6)])

library(reshape)
data.m <- melt(data.1, id=c("ppt", "distribution"))
data.m$stimuli <- NaN
data.m$question <- NaN
data.u <- subset(data.m, distribution=='Unimodal')
data.b <- subset(data.m, distribution=='bimodal')

stim.u <- c(9,23,27,30,33,36,39,42,45,49,63)
stim.b <- c(9,12,15,19,23,36,49,53,57,60,63)

for (stim.idx in 1:length(stim.u)){
  term <- paste('_', as.character(stim.idx), sep='')  
  data.u$stimuli[grep(pattern = term, x = data.u$variable)] <- stim.u[stim.idx]
}

for (stim.idx in 1:length(stim.b)){
  term <- paste('_', as.character(stim.idx), sep='')  
  data.b$stimuli[grep(pattern = term, x = data.b$variable)] <- stim.b[stim.idx]
}

data <- rbind(data.u, data.b)
data$question <- substr(data$variable, 1, 2)
names(data)[4] <- 'response'
data$question <- as.factor(data$question)

write.table('', file="wood-data-appended.csv", sep=",")

#generic functions
rep.row = function(x, n)
{
  matrix(rep(x, each = n), nrow = n)
}

rep.col = function(x, n)
{
  matrix(rep(x, each = n), ncol = n, byrow = TRUE)
}

scale.to.data = function (prediction, data)
{
  (prediction * (max(data) - min(data))) + min(data)
}

scale.transform <- function(value, min, max){
  value * (max - min) + min
}

#generic sdbs functions
sdbs.calculate.distance <- function(stimuli)
{
  n = length(stimuli)
  abs(rep.row(stimuli, n) - t(rep.row(stimuli, n)))
}

sdbs.calculate.similarity <- function(distance, confusability)
{
  similarity = exp(distance * ((-1) * confusability))
  similarity
}

sdbs.calculate.discriminability <- function(similarity)
{
  discrim <- similarity / rep.col(rowSums(similarity,1), nrow(similarity))
  discrim
}

sdbs.calculate.threshold <- function(discriminability, slope, threshold)
{
  1/ (1 + exp (-slope * (discriminability - threshold)))
}

sdbs.calculate.overall <- function(thresholded_probability)
{ 
  1 - apply(1 - thresholded_probability, 1, prod)
}

sdbs.predictions.simple <- function(stim, c, s, t)
{
  sdbs.calculate.overall(sdbs.calculate.threshold(sdbs.calculate.discriminability(sdbs.calculate.similarity(sdbs.calculate.distance(stim), c)), s, t))
}

sdbs.predictions.dbs <- function(recall.probability)
{
  if (sum(recall.probability) == 0)
  {
    rep(0, length(recall.probability))
  } else {
    sum_lower <- c (0, cumsum (recall.probability))
    sum_lower <- sum_lower[1: (length (sum_lower) - 1)]
    sum_higher <- c (rev (cumsum (rev (recall.probability))), 0)
    sum_higher <- sum_higher[2:length (sum_higher)] 
    sdbs <- sum_lower/ (sum_lower + sum_higher)# +.000000001) #added constant to prevent NaN from 0/0
    sdbs[is.nan(sdbs)] <- 0 #turn NaN values into zero. This should prevent error
    sdbs
  }
}

sdbs.prediction.sdbs <- function(stimuli, c, s, t)
{
  pred <- sdbs.predictions.dbs(sdbs.predictions.simple(stimuli, c, s, t))
  pred
}

get.likelihood <- function(prediction, data, std)
{
  #get intervals on cdf for ordinal responses
  low <- data - .5
  high <- data + .49
  
  lhood.upper = pnorm(high, prediction, sd = std)
  lhood.lower = pnorm(low, prediction, sd = std)
  
  lhood <- lhood.upper - lhood.lower
  lhood <- lhood + 0.000000000000000000001 #avoid the inf values from log transforming 0
  lnL <- -2 * log(lhood)
  rec <- data.frame(data, prediction, lhood.upper, lhood.lower, lhood, lnL)
  lnL
}

sdbs.wrapper.data <- function(df, c, s, std)
{
  #model likelihood
  c <- (1/(1+exp(-c)))*100 #transformation of t to keep between 0 and 1
  s <- (1/(1+exp(-c)))*100
  pred <- sdbs.prediction.sdbs(df$stimuli, c, s, t=0.5)
  scale.pred <- scale.to.data(pred, df$response)
  print.rmsd(scale.pred, df$response)
  get.likelihood(scale.pred, df$response, std)  
}

print.rmsd <- function(pred, data){
  diff <- as.character(abs(pred-data))
  diff <- data.frame(response= data, p=pred, d = diff, model=model)
  s <- sum(sqrt(0.5 * sum(abs(pred-data)^2)))
  rmsd <<- s #global variable so we can add this to summary
}

#SDbS_cs
sdbs.wrapper.ppt <- function(par, ppt.data)
{
  fit <- by(ppt.data, ppt.data$question, sdbs.wrapper.data, c = par[1], s = par[2], std = par[3])
  fit <- unlist(fit)
  return(sum(fit))
}

sdbs.wrapper.optim <- function(par, ppt.data)
{
  optim(par = par, fn = sdbs.wrapper.ppt, ppt.data = ppt.data)
}

#create sdbs_cs grid
sdbs.p.min <- c(0.01, 0.01, 0.01)
sdbs.p.max <- c(100, 100, 100)
sdbs.p.n <- c(10, 10, 20)
sdbs.grid <- expand.grid(c.start = seq(from = sdbs.p.min[1], to = sdbs.p.max[1], by = sdbs.p.n[1]), s.start = seq(from = sdbs.p.min[2], to = sdbs.p.max[2], by = sdbs.p.n[2]), sd.start = seq(from = sdbs.p.min[3], to = sdbs.p.max[3], by = sdbs.p.n[3]), participant = seq(from=1, to=ppt.max, by=1))
sdbs.grid <- cbind(sdbs.grid, data.frame(convergence=0, NeglnL=0, c=0, s=0, sd=0))

#multicore sdbs_cs fit
model <<- 'sdbs'
print('multicore sdbs fit')
print(nrow(sdbs.grid))
strt<-Sys.time()
sdbs.search <- foreach(grid.idx = 1:nrow(sdbs.grid), .combine=rbind) %dopar% {
  result <- sdbs.wrapper.optim(par = c(sdbs.grid$c.start[grid.idx], sdbs.grid$s.start[grid.idx], sdbs.grid$sd.start[grid.idx]), subset(data, ppt == sdbs.grid$participant[grid.idx]))
  grid <- sdbs.grid[grid.idx,]
  grid[,c("c", "s", "sd", "convergence", "NeglnL", "rmsd")] <- c(result$par, result$convergence, result$value, rmsd)
  print(grid)
  grid
}

sdbs.time <- Sys.time()-strt

#sdbs_range
sdbs_range.prediction.sdbs <- function(stimuli, c, s, t, w)
{
  pred <- (w * ((stimuli - min(stimuli))/(max(stimuli) - min(stimuli)))) + ((1 - w) * sdbs.predictions.dbs(sdbs.predictions.simple(stimuli, c, s, t)))
  pred
}

sdbs_range.wrapper.data <- function(df, c, s, w, std)
{
  #model likelihood
  c <- (1/(1+exp(-c)))*100 #transformation of t to keep between 0 and 1
  s <- (1/(1+exp(-c)))*100
  w <- 1/(1+exp(-w)) #transformation of w to keep between 0 and 1
  pred <- sdbs_range.prediction.sdbs(df$stimuli, c, s, t=0.5, w)
  scale.pred <- scale.to.data(pred, df$response)
  print.rmsd(scale.pred, df$response)
  get.likelihood(scale.pred, df$response, std)  
}

sdbs_range.wrapper.ppt <- function(par, ppt.data)
{
  fit <- by(ppt.data, ppt.data$question, sdbs_range.wrapper.data, c = par[1], s = par[2], w = par[3], std = par[4])
  fit <- unlist(fit)
  #print(mean(fit))
  #print(fit)
  return(sum(fit))
}

sdbs_range.wrapper.optim <- function(par, ppt.data)
{
  optim(par = par, fn = sdbs_range.wrapper.ppt, ppt.data = ppt.data)#, method = "L-BFGS-B", lower = c(0.01, 0.01, 0.01, 0.01), upper = c(1000, 1000, 1, 1000))
}

#create sdbs_range grid
sdbs_range.p.min <- c(0.01, 0.01, 0.01)
sdbs_range.p.max <- c(100, 100, 100)
sdbs_range.p.n <- c(10, 10, 20)
sdbs_range.grid <- expand.grid(c.start = seq(from = sdbs_range.p.min[1], to = sdbs_range.p.max[1], by = sdbs_range.p.n[1]), 
                               s.start = seq(from = sdbs_range.p.min[2], to = sdbs_range.p.max[2], by = sdbs_range.p.n[2]), 
                               w.start = seq(from = sdbs_range.p.min[3], to = sdbs_range.p.max[3], by = sdbs_range.p.n[3]), 
                               sd.start = seq(from = sdbs_range.p.min[3], to = sdbs_range.p.max[3], by = sdbs_range.p.n[3]), 
                               participant = seq(from=1, to=ppt.max, by=1))

sdbs_range.grid <- cbind(sdbs_range.grid, data.frame(convergence=0, NeglnL=0, c=0, s=0, w=0, sd=0))

#multicore sdbs_range fit
model <<- 'sdbs_range'
print('multicore sdbs_range fit')
print(nrow(sdbs_range.grid))
strt<-Sys.time()
sdbs_range.search <- foreach(grid.idx = 1:nrow(sdbs_range.grid), .combine=rbind) %dopar% {
  #sdbs_range.search <- foreach(grid.idx = 1:10, .combine=rbind) %dopar% {
  result <- sdbs_range.wrapper.optim(par = c(sdbs_range.grid$c.start[grid.idx], sdbs_range.grid$s.start[grid.idx], sdbs_range.grid$w.start[grid.idx], sdbs_range.grid$sd.start[grid.idx]), subset(data, ppt == sdbs_range.grid$participant[grid.idx]))
  grid <- sdbs_range.grid[grid.idx,]
  grid[,c("c", "s", "w", "sd", "convergence", "NeglnL", "rmsd")] <- c(result$par, result$convergence, result$value, rmsd)
  print(grid)
  grid
}

sdbs_range.time <- Sys.time()-strt

#sgems
# raise distance to power gamma
# calculate recall probablility of item
# multiply raised distance by item recall probability
# prediction is the relative rank of the cumulative sum, just like the normal sdbs model

sgems.prediction.sdbs <- function(stimuli, c, s, t, y)
{
  dis <- sdbs.calculate.distance(stimuli)^y
  prob <- rep.row(sdbs.predictions.simple(stimuli, c, s, t), length(stimuli))
  prod <- dis * prob
  
  prod[prod==Inf] <- 0

  lowerAndHigher <- rowSums(prod)
  lower <- prod
  lower[upper.tri(lower)] <- 0
  lower <- rowSums(lower)
  pred <- lower/lowerAndHigher
}

sgems.wrapper.data <- function(df, c, s, y, std)
{
  #model likelihood
  c <- (1/(1+exp(-c)))*100 #transformation of t to keep between 0 and 1
  s <- (1/(1+exp(-c)))*100
  pred <- sgems.prediction.sdbs(df$stimuli, c, s, t=0.5, y)
  #print(pred)
  scale.pred <- scale.to.data(pred, df$response)
  print.rmsd(scale.pred, df$response)
  get.likelihood(scale.pred, df$response, std)  
}

sgems.wrapper.ppt <- function(par, ppt.data)
{
  fit <- by(ppt.data, ppt.data$question, sgems.wrapper.data, c = par[1], s = par[2], y = par[3], std = par[4])
  fit <- unlist(fit)
  #print(mean(fit))
  #print(fit)
  return(sum(fit))
}

sgems.wrapper.optim <- function(par, ppt.data)
{
  optim(par = par, fn = sgems.wrapper.ppt, ppt.data = ppt.data)#, method = "L-BFGS-B", lower = c(0.01, 0.01, 0.01, 0.01), upper = c(1000, 1000, 1, 1000))
}

sgems.p.min <- c(0.01, 0.01, 0.01)
sgems.p.max <- c(50, 100, 100)
sgems.p.n <- c(10, 10, 20)
sgems.grid <- expand.grid(c.start = seq(from = sgems.p.min[1], to = sgems.p.max[1], by = sgems.p.n[1]), 
                               s.start = seq(from = sgems.p.min[2], to = sgems.p.max[2], by = sgems.p.n[2]), 
                               y.start = seq(from = sgems.p.min[3], to = sgems.p.max[3], by = sgems.p.n[3]), 
                               sd.start = seq(from = sgems.p.min[3], to = sgems.p.max[3], by = sgems.p.n[3]), 
                               participant = seq(from=1, to=ppt.max, by=1))

sgems.grid <- cbind(sgems.grid, data.frame(convergence=0, NeglnL=0, c=0, s=0, y=0, sd=0))

#multicore sdbs_range fit
model <<- 'sgems'
print('multicore sgems fit')
print(nrow(sgems.grid))
strt<-Sys.time()
sgems.search <- foreach(grid.idx = 1:(nrow(sgems.grid)), .combine=rbind) %dopar% {
  result <- sgems.wrapper.optim(par = c(sgems.grid$c.start[grid.idx], sgems.grid$s.start[grid.idx], sgems.grid$y.start[grid.idx], sgems.grid$sd.start[grid.idx]), subset(data, ppt == sgems.grid$participant[grid.idx]))
  grid <- sgems.grid[grid.idx,]
  grid[,c("c", "s", "y", "sd", "convergence", "NeglnL", "rmsd")] <- c(result$par, result$convergence, result$value, rmsd)
  print(grid)
  grid
}

#transform c, s and w values back to the value used to create best model predictions
sdbs.search$c <- (1/(1+exp(-sdbs.search$c)))*100
sdbs.search$s <- (1/(1+exp(-sdbs.search$s)))*100
sdbs_range.search$c <- (1/(1+exp(-sdbs_range.search$c)))*100
sdbs_range.search$s <- (1/(1+exp(-sdbs_range.search$s)))*100
sdbs_range.search$w <- (1/(1+exp(-sdbs_range.search$w)))
sgems.search$c <- (1/(1+exp(-gems.search$c)))*100
sgems.search$s <- (1/(1+exp(-gems.search$s)))*100
