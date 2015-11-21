# multivariate point BF test file 
rm(list=ls())
library(ggplot2)
library(spatialkernel)

# number of grid points along each dimension
n.steps <- 20
max.x <- 1
max.y <- 10

#bivariate uniform samples
n <- 10000
x <- runif(n, min=0, max=max.x)
y <- runif(n, min=0, max=max.y)
data <- matrix(c(x,y),nrow=n)

ggplot(as.data.frame(data), aes(V1, V2)) +
  geom_bin2d()

h <- 0.2 # bandwidth parameter of how much the data should be smoothed

# set up matrix of points to test
xpts <- seq(from = 0, to = 1, length = n.steps)
ypts <- seq(from = 0, to = 10, length = n.steps)
gpts <- as.matrix(expand.grid(xpts = xpts, ypts = ypts))

# matrix of boundaries of the distribution
boundaries <- matrix(c(c(0,0),c(0,10),c(1,10),c(1,0)), nrow=4, byrow=T)

k.est <- lambdahat(data, h, gpts=gpts, poly=boundaries, edge=T)

d <- data.frame(x = k.est$gpts[,1], y = k.est$gpts[,2], z = (k.est$lambda/sum(k.est$lambda)))

ggplot(d, aes(x, y)) +
  geom_tile(aes(fill = z), colour = 'white') +
  scale_fill_gradient(low = "white", high = "steelblue")

A <- (max.x * max.y)/(n.steps-1)^2

key.index <- which(gpts[,1]==1 & gpts[,2]==0)
k.est$gpts[key.index,]
post.dens <- k.est$lambda[key.index]/(A*sum(k.est$lambda))
prior.dens <- 1/10
BF <- post.dens/prior.dens
