library("ggplot2")
library("RColorBrewer")
library("Matrix")
library("reshape2")

plot_map_reduce <- function(p, d, N, map = function(x) {
                              return(x)
                            },
                            reduce, method, fname, show_sd = FALSE,
                            dir_name = "res", plot_title = "", plot_ylab = "", ...) {
  d_len <- length(d)
  dir_len <- length(dir_name)
  exp_name <- outer(p, d, paste, sep = "_")
  exp_fname <- matrix(paste0(exp_name, ".rds"), ncol = d_len)

  res <- matrix(
    nrow = length(p), ncol = d_len,
    dimnames = list(p = p, d = d)
  )
  res_sd <- matrix(
    nrow = length(p), ncol = d_len,
    dimnames = list(p = p, d = d)
  )


  for (i in 1:length(p)) {
    sample <- array(dim = c(p[i], p[i], N * dir_len))
    for (j in 1:d_len) {
      for (k in 1:dir_len) {
        sample[, , ((k - 1) * N + 1):(k * N)] <-
          readRDS(file = paste0(dir_name[k], "/", method, "_", exp_fname[i, j]))
      }
      mapd_mat <- apply(sample, MARGIN = 3, map, ...)
      res[i, j] <- reduce(mapd_mat)
      res_sd[i, j] <- sd(mapd_mat)
    }
  }

  wd <- getwd()
  dir.create(paste0(wd, "/plot_", dir_len), showWarnings = FALSE)

  palette <- colorRampPalette(colors = c("black", "red"))
  colors <- palette(d_len)

  df <- melt(res)
  df$d <- as.factor(df$d)
  df$sd <- melt(res_sd)$value

  pl <- ggplot(df, aes(x = p, y = value, group = d, color = d)) +
    geom_line() +
    geom_point() +
    theme(text = element_text(size = 20), legend.position = "bottom") +
    scale_color_manual(values = colors) +
    xlab("Number of nodes") +
    ylab(plot_ylab) +
    ggtitle(plot_title)

  if (show_sd == TRUE) {
    pl <- pl + geom_ribbon(aes(ymin = value - sd, ymax = value + sd))
  }

  ggsave(filename = fname, plot = pl, device = "pdf", path = paste0("plot_", dir_len, "/"))
}

plot_map_reduce_cmp <- function(p, d, N, map = function(x) {
                                  return(x)
                                },
                                reduce, fname, show_sd = FALSE,
                                dir_name = "res", plot_title = "", plot_ylab = "", ...) {
  dir_len <- length(dir_name)
  methods <- c("diagdom", "port")
  exp_name <- outer(methods, p, paste, sep = "_")
  exp_fname <- matrix(paste0(exp_name, "_", d, ".rds"),
    ncol = 2,
    dimnames = list(p = p, method = methods), byrow = TRUE
  )


  res <- matrix(
    nrow = length(p), ncol = 2,
    dimnames = list(p = p, method = methods)
  )
  res_sd <- matrix(
    nrow = length(p), ncol = 2,
    dimnames = list(p = p, method = methods)
  )

  for (i in 1:length(p)) {
    sample <- array(dim = c(p[i], p[i], N * dir_len))
    for (m in methods) {
      for (k in 1:dir_len) {
        sample[, , ((k - 1) * N + 1):(k * N)] <-
          readRDS(file = paste0(dir_name[k], "/", exp_fname[i, m]))
      }
      mapd_mat <- apply(sample, MARGIN = 3, map, ...)
      res[i, m] <- reduce(mapd_mat)
      res_sd[i, m] <- sd(mapd_mat)
    }
  }

  wd <- getwd()
  dir.create(paste0(wd, "/plot_", dir_len), showWarnings = FALSE)

  df <- melt(res)
  df$sd <- melt(res_sd)$value

  pl <- ggplot(df, aes(x = p, y = value, group = method)) +
    xlab("Number of nodes") +
    ylab(plot_ylab) +
    ggtitle(plot_title)


  if (show_sd == TRUE) {
    pl <- pl +
      geom_ribbon(aes(ymin = value - sd, ymax = value + sd, fill = method),
        alpha = .3
      ) +
      scale_fill_manual(labels = c("DD", "PO"), values = c("green4", "blue"))
  }

  pl <- pl +
    geom_line(aes(color = method)) +
    geom_point(aes(color = method)) +
    scale_color_manual(labels = c("DD", "PO"), values = c("green4", "blue")) +
    theme(text = element_text(size = 20), legend.position = "bottom")

  ggsave(filename = fname, plot = pl, device = "pdf", path = paste0("plot_", dir_len, "/"))
}

ext_abs_offdiag <- function(mat, fext = max, inverse = FALSE) {
  if (inverse == TRUE) {
    mat <- Matrix::solve(drop0(mat))
  }
  return(fext(0, abs((mat / diag(mat))[upper.tri(mat) & mat != 0])))
}

plot_eigen_freqpol <- function(p, d, N, method, fname, dir_name = "res", bin_fun = nclass.Sturges, plot_title = "", ...) {
  d_len <- length(d)
  r <- length(dir_name)
  exp_fname <- paste0(p, "_", d, ".rds")

  eigen_vals <- array(dim = c(d_len, p * r * N), dimnames = list(d = d))
  sample <- array(dim = c(p, p, r * N))
  for (i in 1:d_len) {
    for (j in 1:r) {
      sample[, , ((j - 1) * N + 1):(j * N)] <-
        readRDS(file = paste0(dir_name[j], "/", method, "_", exp_fname[i]))
    }
    eigen_vals[i, ] <- apply(sample, MARGIN = 3, function(m) {
      return(eigen(m)$values)
    })
  }
  wd <- getwd()
  dir.create(paste0(wd, "/plot_", r), showWarnings = FALSE)

  palette <- colorRampPalette(colors = c("black", "red"))
  colors <- palette(d_len)

  df <- melt(eigen_vals)
  df$d <- as.factor(df$d)

  pl <- ggplot(df, aes(value, group = d, color = d)) +
    geom_freqpoly(bins = bin_fun(df$value)) +
    xlab("Eigenvalue") +
    ylab("Absolute frequency") +
    ggtitle(plot_title) +
    scale_color_manual(values = colors)

  ggsave(filename = fname, plot = pl, device = "pdf", path = paste0("plot_", r, "/"))
}

plot_time <- function(p, d, method, fname, dir_name = "res", plot_title = "", ...) {
  d_len <- length(d)
  exp_name <- outer(p, d, paste, sep = "_")
  exp_fname <- matrix(paste0(exp_name, ".rds"), ncol = d_len)

  res <- matrix(
    nrow = length(p), ncol = d_len,
    dimnames = list(p = p, d = d)
  )

  for (i in 1:length(p)) {
    for (j in 1:length(d)) {
      time <- readRDS(file = paste0(dir_name, "/t_", method, "_", exp_fname[i, j]))
      res[i, j] <- summary(time, unit = "s")[1, 5]
    }
  }

  wd <- getwd()
  dir.create(paste0(wd, "/plot_", dir_name), showWarnings = FALSE)

  palette <- colorRampPalette(colors = c("black", "red"))
  colors <- palette(d_len)

  df <- melt(res)
  df$d <- as.factor(df$d)

  pl <- ggplot(df, aes(x = p, y = value, group = d, color = d)) +
    geom_line() +
    geom_point() +
    theme(text = element_text(size = 20), legend.position = "bottom") +
    scale_color_manual(values = colors) +
    xlab("Number of nodes") +
    ylab("Execution time in seconds") +
    ggtitle(plot_title)

  ggsave(
    filename = fname, plot = pl, device = "pdf",
    path = paste0("plot_", dir_name, "/"),
    width = 7, height = 5
  )
}


is.diagdom <- function(M){
  all(2*diag(M) >= rowSums(abs(M)) )
}


ronion <- function(N, p) {
  sample <- array(dim = c(p, p, N))
  for (i in 1:N){
    sample[, , i] <- clusterGeneration::genPositiveDefMat(dim = p, rangeVar = c(1,1),
                                                          covMethod  ="onion")$Sigma
  }
  return(sample)
}
rvine <- function(N, p) {
  sample <- array(dim = c(p, p, N))
  for (i in 1:N){
    sample[, , i] <- clusterGeneration::genPositiveDefMat(dim = p, rangeVar = c(1,1),
                                                          covMethod  ="c-vine")$Sigma
  }
  return(sample)
}

rpolar <- function(N, p) {
  sample <- array(dim = c(p, p, N))
  for (i in 1:N){
    sample[, , i] <- randcorr::randcorr(p = p)
  }
  return(sample)
}