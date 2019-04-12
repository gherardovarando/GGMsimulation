library(gmat)
library(ggplot2)
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
  xlab("(1,2)") + 
  ylab("(3,2)") + 
  ggsave("plot_3var_port.pdf")

nrow(sub_sample(s1, bounds = c(-0.5,0.5))) / N



system.time(sample2 <- port_chol(N, ug = ug, h = 1000, eps = 0.5))
s2 <- as.data.frame(vectorize(sample2))

ggplot(s2) + 
  geom_point(aes(x = V1, y = V3), size = 0.1) + 
  coord_fixed() +
  xlab("(1,2)") + 
  ylab("(3,2)") + 
  ggsave("plot_3var_port_chol.pdf")

nrow(sub_sample(s2, bounds = c(-0.5,0.5))) / N


sample3 <- chol_mh(N, dag = dag, h = 1000, eps = 0.5)

s3 <- as.data.frame(vectorize(sample3))
ggplot(s3) + 
  geom_point(aes(x = V1, y = V3), size = 0.1) + 
  coord_fixed() +
  xlab("(1,2)") + 
  ylab("(3,2)") + 
  ggsave("plot_3var_chol_mh.pdf")
nrow(sub_sample(s3, bounds = c(-0.5,0.5))) / N

