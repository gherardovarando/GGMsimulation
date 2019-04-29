library(gmat)
library(ggplot2)
library(igraph)
source("plot_utils.R")


Adj <- matrix(nrow = 3, ncol = 3, data = 0)
Adj[1,2] <- 1
Adj[2,3] <- 1
dag <- graph_from_adjacency_matrix(Adj, "directed")
plot(dag)

ug <- as.undirected(dag)
plot(ug)

N <- 5000

sample1 <- port(N, ug = ug)
s1 <- as.data.frame(vectorize(sample1))
ggplot(s1) + 
  geom_point(aes(x = V1, y = V3), size = 0.1) + 
  coord_fixed() +
  xlab("") + 
  ylab("") + 
  ggtitle("partial orthogonalization") +
  ggsave("plot_3var_port.pdf",width = 3.6, height = 3.4)

nrow(sub_sample(s1, bounds = c(-0.5,0.5))) / N



system.time(sample2 <- diagdom(N, ug = ug,rfun = rnorm))
sample2 <- array(dim = dim(sample2), data = apply(sample2, 3, cov2cor))
s2 <- as.data.frame(vectorize(sample2))

ggplot(s2) + 
  geom_point(aes(x = V1, y = V3), size = 0.1) + 
  coord_fixed() +
  xlab("") + 
  ylab("") + 
  ggtitle("diagonal dominance") + 
  ggsave("plot_3var_diagdom.pdf",width = 3.6, height = 3.4)

nrow(sub_sample(s2, bounds = c(-0.5,0.5))) / N


sample3 <- chol_mh(N, dag = dag, h = 1000, eps = 0.5)

s3 <- as.data.frame(vectorize(sample3))
ggplot(s3) + 
  geom_point(aes(x = V1, y = V3), size = 0.1) + 
  coord_fixed() +
  xlab("") + 
  ylab("") + 
  ggtitle("uniform") +
  ggsave("plot_3var_chol_mh.pdf",width = 3.6, height = 3.4)
nrow(sub_sample(s3, bounds = c(-0.5,0.5))) / N

