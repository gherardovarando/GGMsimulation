library("ggplot2")
library("microbenchmark")
source("plot_utils.R")


## Plot time comparison results
p <- seq(from = 100, to = 1000, by = 100)
d <- c(0.0025, 0.005, 0.025, 0.05, 0.25, 0.5)
plot_time(
  p = p, d = d, method = "port", fname = "time_port.pdf",
  dir_name = "res_fast_t", plot_title = "Partial orthogonalization method"
)
#plot_time(
#  p = p, d = d, method = "diagdom", fname = "time_diagdom.pdf",
#  dir_name = "res_t", plot_title = "Diagonal dominance method"
#)
