library(gmat)
library(ggplot2)
source("plot_utils.R")
p <- 50

Adj <- matrix(nrow = p, ncol = p, 0)
for (i in 2:p){
  Adj[i - 1, i] <- 1
}

dag <- igraph::graph_from_adjacency_matrix(Adj)
ug <- igraph::graph_from_adjacency_matrix(Adj,mode = "undirected")
plot(dag)

N <- 10000


########## port 

sample <- port(N, ug = ug)
reduced <- as.data.frame(t(apply(sample, 3, function(M){
  return(M[Adj == 1])
})))
ggplot(stack(reduced)) + geom_density(aes(x = values, group = ind)) + 
  xlab("") + ylab("") +   ggtitle("partial orthogonalization") + 
  ggsave("plot_chain_port.pdf")

########### port_chol
sample <- port_chol(N, ug = ug, h = 1000, eps= 0.5)

reduced <- as.data.frame(t(apply(sample, 3, function(M){
  return(M[Adj == 1])
})))
ggplot(stack(reduced)) + geom_density(aes(x = values, group = ind)) + 
  xlab("") + ylab("") +   ggtitle("uniform") + 
  ggsave("plot_chain_port_chol.pdf")


#### diagdom

sample <- diagdom(N, ug = ug, rfun = rnorm)
sample <- array(apply(sample ,MARGIN =  3, cov2cor), dim = dim(sample))


reduced <- as.data.frame(t(apply(sample, 3, function(M){
  return(M[Adj == 1])
})))
ggplot(stack(reduced)) + geom_density(aes(x = values, group = ind)) + 
  xlab("") + ylab("") +   ggtitle("diagonal dominance") + 
  ggsave("plot_chain_diagdom.pdf")
