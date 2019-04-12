library(gmat)
source("plot_utils.R")
p <- 5

Adj <- matrix(nrow = p, ncol = p, 0)
for (i in 2:p){
  Adj[i - 1, i] <- 1
}
ix <- c(1)
i <- 2
while (i < (p + 1)){
  ix[i] <- i + 1
  ix[i + 1] <- i
  i <- i + 2
}
ix <- ix[1:p]
ix[p] <- min(p, ix[p])

#ix <- sample(1:p)
#ix <- 1:p
inv <- Matrix::invPerm(ix)
Adj <- Adj[ix, ix]
dag <- igraph::graph_from_adjacency_matrix(Adj)
ug <- igraph::graph_from_adjacency_matrix(Adj,mode = "undirected")
plot(dag)

N <- 1000


########## port 

sample0 <- port(N, ug = ug)

plot(density(sample0[inv[1], inv[2],]), xlim = c(-1,1), ylim = c (0,2), col = 1)
for (i in 3:p){
  lines(density(sample0[inv[i-1], inv[i],]), col = i - 1)
}
legend("topleft", col = 2:p - 1, lty = 1, legend = 2:p ,cex = 0.5)

s0 <- t(apply(sample0, 3, function(mat){
  c(mat[inv[1],inv[2]], mat[inv[2],inv[3]], mat[inv[3],inv[4]])
}))

pairs(s0)
plot_elliptope(s0, pallette = NULL, sort.fun = NULL, col = "black")


########### port_chol
sample1 <- port_chol(N, ug = ug, h = 1000, eps= 1)

plot(density(sample1[inv[1], inv[2],]), xlim = c(-1,1), ylim = c (0,2), col = 1)
for (i in 3:p){
  lines(density(sample1[inv[i-1], inv[i],]), col = i - 1)
}
legend("topleft", col = 2:p - 1, lty = 1, legend = 2:p ,cex = 0.5)

s1 <- t(apply(sample1, 3, function(mat){
  c(mat[inv[1],inv[2]], mat[inv[2],inv[3]], mat[inv[3],inv[4]])
}))

ramp <- colorRamp(c('red', "blue"))
pairs(s1, asp= 1)
plot_elliptope(s1, pallette = NULL, sort.fun = NULL, col = "black")

#### unif

sample2 <- chol_mh(N, dag = dag, h = 10000, eps = 1)

plot(density(sample2[inv[1], inv[2],]), xlim = c(-1,1), ylim = c (0,2), col = 1)
for (i in 3:p){
  lines(density(sample2[inv[i-1], inv[i],]), col = i - 1)
}
legend("topleft", col = 2:p - 1, lty = 1, legend = 2:p ,cex = 0.5)

s2 <- t(apply(sample2, 3, function(mat){
  c(mat[inv[1],inv[2]], mat[inv[2],inv[3]], mat[inv[3],inv[4]])
}))
pairs(s2, asp = 1)
plot_elliptope(s2, pallette = NULL, sort.fun = NULL, col = "black")


#### diagdom

sample3 <- diagdom(N, ug = ug, rfun = rnorm)
sample3 <- array(apply(sample3 ,MARGIN =  3, cov2cor), dim = dim(sample3))


plot(density(sample3[inv[1], inv[2],]), xlim = c(-1,1), ylim = c(0,1.3), 
     main = "", xlab = "")
for (i in 3:p){
  lines(density(sample3[inv[i-1], inv[i],]), col = i)
}


s3 <- t(apply(sample3, 3, function(mat){
  c(mat[inv[1],inv[2]], mat[inv[2],inv[3]], mat[inv[3],inv[4]])
}))

pairs(s3, asp = 1)
plot_elliptope(s3,col = "black", sort.fun = NULL, pallette = NULL)

