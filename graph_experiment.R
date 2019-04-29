library(gmat)
library(igraph)
library(gRbase)
library(ggplot2)
library(gRbase)
p <- 50
d <- 0.05

N <- 10000


ug <- rgraph(p,d)
plot(ug)



## port_chol
Adj <- as_adjacency_matrix(ug, sparse = FALSE)

system.time(sample <- port_chol(N, ug = ug, h = 10000, eps = 0.5))

reduced <- as.data.frame(t(apply(sample, 3, function(mat){
  return(mat[upper.tri(Adj) & Adj == 1])
})))



ggplot(stack(reduced)) + geom_density(aes(x = values, group = ind)) + 
  xlab("") + ylab("") +   ggtitle("uniform + partial orthogonalization") + 
 ggsave("plot_graph_exp_port_chol.pdf")



## port
system.time(sample <- port(N, ug = ug))

reduced <- as.data.frame(t(apply(sample, 3, function(mat){
  return(mat[upper.tri(Adj) & Adj == 1])
})))

ggplot(stack(reduced)) + geom_density(aes(x = values, group = ind)) + 
  xlab("") + ylab("") +   ggtitle("partial orthogonalization") + 
 ggsave("plot_graph_exp_port.pdf")


###diagdom
system.time(sample <- diagdom(N, ug = ug))
sample <- array(dim = dim(sample), data = apply(sample, 3, cov2cor))


reduced <- as.data.frame(t(apply(sample, 3, function(mat){
  return(mat[upper.tri(Adj) & Adj == 1])
})))


ggplot(stack(reduced)) + geom_density(aes(x = values, group = ind)) + 
  xlab("") + ylab("") +   ggtitle("diagonal dominance") + 
 ggsave("plot_graph_exp_diagdom.pdf")


################ chordal graphs


colnames(Adj) <- rownames(Adj) <- 1:p


Adjch <- triangulate(Adj)
colnames(Adjch) <- rownames(Adjch) <- NULL

ugch <- graph_from_adjacency_matrix(Adjch, mode = "undirected")
plot(ugch)


##port
system.time(sample <- port(N, ug = ugch))
reduced <- as.data.frame(t(apply(sample, 3, function(mat){
  return(mat[upper.tri(Adjch) & Adjch == 1])
})))
ggplot(stack(reduced)) + geom_density(aes(x = values, group = ind)) + 
  xlab("") + ylab("") +ggtitle("partial orthogonalization") +
  ggsave("plot_graph_exp_ch_port.pdf")


##port_chol
system.time(sample <- port_chol(N, ug = ugch))
reduced <- as.data.frame(t(apply(sample, 3, function(mat){
  return(mat[upper.tri(Adjch) & Adjch == 1])
})))
ggplot(stack(reduced)) + geom_density(aes(x = values, group = ind)) + 
  ggtitle("") + 
  xlab("") + ylab("") +  ggtitle("uniform") +
  ggsave("plot_graph_exp_ch_port_chol.pdf")

###diagdom
system.time(sample <- diagdom(N, ug = ugch))
sample <- array(dim = dim(sample), data = apply(sample, 3, cov2cor))
reduced <- as.data.frame(t(apply(sample, 3, function(mat){
  return(mat[upper.tri(Adjch) & Adjch == 1])
})))
ggplot(stack(reduced)) + geom_density(aes(x = values, group = ind)) + 
  xlab("") + ylab("") + ggtitle("diagonal dominance") +
  ggsave("plot_graph_exp_ch_diagdom.pdf")

