library("gmat")

dir.create("res_sizeoffdiag", showWarnings = FALSE)

N <- 10
M <- 10
ps <- seq(from = 100, to = 1000, by = 100)
ds <- c(0.0025, 0.005, 0.025, 0.05, 0.25, 0.5)


res1 <- array(dim = c(length(ps),  length(ds), 6 ), data = 0, 
              dimnames = list(ps, ds, names(summary(runif(10)))))

res2 <- array(dim = c(length(ps),  length(ds), 6 ), data = 0, 
              dimnames = list(ps, ds, names(summary(runif(10)))))

for (i in 1:length(ps)){
  for (j in 1:length(ds)){
    p <- ps[i]
    d <- ds[j]
    temp1 <- temp2 <- array(dim = c(N , M))
    for (k in 1:M){
      ug <- rgraph(p, d)
      sample1 <- port(N = N, d = d, p = p, zapzeros = TRUE, rfun = rnorm)
      sample2 <- diagdom(N = N, d = d, p = p, rfun = rnorm)
      temp1[, k] <- apply(sample1, MARGIN = 3, FUN = function(M){
        max(abs(M[upper.tri(M) & M!=0]))
      })
      temp2[, k] <- apply(sample2, MARGIN = 3, FUN = function(M){
        max(abs(M[upper.tri(M) & M!=0]))
      })
    }
      
    res1[i, j, ] <- summary(c(temp1))
    res2[i, j, ] <- summary(c(temp2))
    
  }
}


write.csv(res1, "res_sizeoffdiag/sizeoff_port_matrix.csv")
write.csv(res2, "res_sizeoffdiag/sizeoff_diagdom_matrix.csv")





Ks <- c(2, 3, 10) 

res1 <- array(dim = c(length(ps),  length(Ks), 6 ), data = 0, 
              dimnames = list(ps, Ks, names(summary(runif(10)))))

res2 <- array(dim = c(length(ps),  length(Ks), 6 ), data = 0, 
              dimnames = list(ps, Ks, names(summary(runif(10)))))

for (i in 1:length(ps)){
  for (j in 1:length(Ks)){
    p <- ps[i]
    d <-  Ks[j] / ( p - 1)
    temp1 <- temp2 <- array(dim = c(N , M))
    for (k in 1:M){
      ug <- rgraph(p, d)
      sample1 <- port(N = N, d = d, p = p, zapzeros = TRUE, rfun = rnorm)
      sample2 <- diagdom(N = N, d = d, p = p, rfun = rnorm)
      temp1[, k] <- apply(sample1, MARGIN = 3, FUN = function(M){
        max(abs(M[upper.tri(M) & M!=0]))
      })
      temp2[, k] <- apply(sample2, MARGIN = 3, FUN = function(M){
        max(abs(M[upper.tri(M) & M!=0]))
      })
    }
    
    res1[i, j, ] <- summary(c(temp1))
    res2[i, j, ] <- summary(c(temp2))
    
  }
}


write.csv(res1, "res_sizeoffdiag/sizeoff_port_K.csv")
write.csv(res2, "res_sizeoffdiag/sizeoff_diagdom_K.csv")


