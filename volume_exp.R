library(gmat)
library(clusterGeneration)
library(randcorr)

N <- 10000; p <- 50

sample <- rmh(N = N, p = p, h = 1000, eps = 0.5)
sample <- gmat::vectorize(sample)
plot_elliptope(x = sample)
pairs(x = sample, lwd = 1 , pch  =20, cex  = 0.3, asp = 1, labels = c("x", "y", "z"))
nrow(sub_sample(sample, c(-0.5,0.5))) / N


sample <- rpolar(N = N, p = p)
sample <- gmat::vectorize(sample)
plot_elliptope(x = sample)
pairs(x = sample, lwd = 1 , pch  =20, cex  = 0.3, asp = 1, labels = c("x", "y", "z"))
nrow(sub_sample(sample, c(-0.5,0.5))) / N


sample <- ronion(N = N, p = p)
sample <- gmat::vectorize(sample)
plot_elliptope(x = sample)
pairs(x = sample, lwd = 1 , pch  =20, cex  = 0.3, asp = 1, labels = c("x", "y", "z"))
nrow(sub_sample(sample, c(-0.5,0.5))) / N


sample <- gmat::diagdom(N = N, p = p )
sample <- array(apply(sample ,MARGIN =  3, cov2cor), dim = dim(sample))
sample <- gmat::vectorize(sample)
plot_elliptope(x = sample)
pairs(x = sample, lwd = 1 , pch  =20, cex  = 0.3, asp = 1, labels = c("x", "y", "z"))
nrow(sub_sample(sample, c(-0.5,0.5))) / N


sample <- gmat::port(N = N, p = p, rfun = rnorm)
sample <- gmat::vectorize(sample)
plot_elliptope(x = sample, sort.fun = NULL)
pairs(x = sample, lwd = 1 , pch  =20, cex  = 0.3, asp = 1, labels = c("x", "y", "z"))
nrow(sub_sample(sample, c(-0.5,0.5))) / N

