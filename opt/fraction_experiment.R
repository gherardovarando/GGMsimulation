library("gmat")
source("plot_utils.R")

dir.create("res_fraction", showWarnings = FALSE)



######## fraction of diagonally dominant matrices generated with port
N <- 1000
M <- 100
ps <- seq(from = 20, to = 100, by = 20)
ds <- c(0.0025, 0.005, 0.025, 0.05, 0.25, 0.5)

res <- matrix(nrow = length(ps), ncol = length(ds), 0, 
              dimnames = list(ps, ds))

for (i in 1:length(ps)){
  for (j in 1:length(ds)){
    p <- ps[i]
    d <- ds[j]
    tmp <- replicate(M, {
      sample <- port(N = N, d = d, p = p, zapzeros = TRUE)
      sum(apply(sample, MARGIN = 3, is.diagdom)) / N  
    })
     res[i, j] <- mean(tmp)
  }
}


write.csv(res, "res_fraction/fraction_port_matrix.csv")

######################  constant dim * dens 

M <- 100
N <- 1000
ps <- seq(from = 20, to = 100, by = 20)
res <- rep(0, length(ps))
for (i in 1:length(ps)){
  p <- ps[i]
  d <- 2 / p
  
  tmp <- replicate(M, {
    sample <- port(N = N, d = d, p = p)
    sum(apply(sample, MARGIN = 3, is.diagdom)) / N  
  })
  res[i] <- mean(tmp)
}
names(res) <- ps
write.csv(res, "res_fraction/fraction_port_const.csv")


res <- rep(0, length(ps))
for (i in 1:length(ps)){
  p <- ps[i]
  d <- 2 / p
  tmp <- replicate(M, {
    sample <- diagdom(N = N, d = d, p = p)
    sample <- array(apply(sample, 3, cov2cor), dim = dim(sample))
    sum(apply(sample, MARGIN = 3, is.diagdom)) / N  
  })
  res[i] <- mean(tmp)
}
names(res) <- ps
write.csv(res, "res_fraction/fraction_diagdom_const.csv")

############ full matrices with different methods

N <- 100000
ps <- 3:7

res <- rep(0, length(ps))
for (i in 1:length(ps)){
    p <- ps[i]
    sample <- ronion(N = N, p = p)
    res[i] <- sum(apply(sample, MARGIN = 3, is.diagdom)) / N
}
names(res) <- ps
write.csv(res, "res_fraction/fraction_onion_full.csv")


res <- rep(0, length(ps))
for (i in 1:length(ps)){
  p <- ps[i]
  sample <- port(N = N, p = p)
  res[i] <- sum(apply(sample, MARGIN = 3, is.diagdom)) / N
}
names(res) <- ps
write.csv(res, "res_fraction/fraction_port_full.csv")



res <- rep(0, length(ps))
for (i in 1:length(ps)){
  p <- ps[i]
  sample <- chol_mh(N = N, p = p, h = 1000, eps = 0.5)
  res[i] <- sum(apply(sample, MARGIN = 3, is.diagdom)) / N
}
names(res) <- ps
write.csv(res, "res_fraction/fraction_cholmh_full.csv")


ps <- 3:20
res <- rep(0, length(ps))
for (i in 1:length(ps)){
  p <- ps[i]
  sample <- diagdom(N = N, p = p, rfun = rnorm)
  sample <- array(apply(sample, 3, cov2cor), dim = dim(sample))
  res[i] <- sum(apply(sample, MARGIN = 3, is.diagdom)) / N
}
names(res) <- ps
write.csv(res, "res_fraction/fraction_cordiagdom_full.csv")

