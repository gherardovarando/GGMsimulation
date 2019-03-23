library(gmat)
library(microbenchmark)

p <- seq(from = 100, to = 1000, by = 100)
d <- c(0.0025, 0.005, 0.025, 0.05, 0.25, 0.5)
N <- 100


exp_name <- outer(p, d, paste, sep = "_")
exp_fname <- matrix(paste0(exp_name, ".rds"), ncol = length(d))

wd <- getwd()

dir_name <- paste0("res_t")
dir.create(paste0(wd, "/", dir_name), showWarnings = FALSE)

for (i in 1:length(p)) {
  for (j in 1:length(d)) {
    saveRDS(microbenchmark(port(N = 100, p = p[i], d = d[j]), times = N),
            file =
              paste0(dir_name, "/t_port_", exp_fname[i, j])
    )
  }
}


