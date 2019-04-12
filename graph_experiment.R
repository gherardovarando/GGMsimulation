library(gmat)
library(igraph)
library(gRbase)

p <- 50
d <- 0.1

N <- 1000


ug <- rgraph(p,d)
plot(ug)

Adj <- as_adjacency_matrix(ug, sparse = FALSE)

system.time(sample1 <- port_chol(N, ug = ug, h = 10000, eps = 1))

reduced1 <- t(apply(sample1, 3, function(mat){
  return(mat[Adj == 1])
}))

plot(density( reduced1[,1]  ), col = 1, xlim = c(-1,1), ylim = c(0,5))
for (i in 2:(dim(reduced1)[2])){
  lines(density( reduced1[,i] ), col = 1)
}


system.time(sample2 <- port(N, ug = ug))

reduced2 <- t(apply(sample2, 3, function(mat){
  return(mat[Adj == 1])
}))


plot(density( reduced2[,1]  ), col = 1, xlim = c(-1,1), ylim = c(0,5))
for (i in 2:(dim(reduced2)[2])){
  lines(density( reduced2[,i] ), col = 1)
}

