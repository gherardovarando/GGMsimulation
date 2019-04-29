library("reshape2")
library("RColorBrewer")
library("ggplot2")

plot_kramer <- function(algorithm, dir_name, fname, plot_title, subtitle,
                        stat = "ppv") {

	methods <- c("adalasso", "lasso", "pls", "shrink", "ridge")
	ylab <- c("ppv" = "True discovery rate", "tpr" = "True positive rate")
	N <- seq(25, 200, 25)
	power <- array(dim = c(length(methods), length(N)),
							 dimnames = list(method = methods, N = N))

	for (method in methods) {
		power.raw <- readRDS(paste0(dir_name, "/", stat, ".", method, ".rds"))
		power[method, ] <- apply(X = power.raw, MARGIN = 2, mean)
	}

	palette <- colorRampPalette(colors = c("black", "red"))
	colpal <- brewer.pal(name = "Set1", n = length(N))

	df <- melt(power)
	df$method <- as.factor(df$method)

	pl <- ggplot(df, aes(x = N, y = value, group = method, color = method)) +
		geom_line() +
		geom_point() +
		theme(text = element_text(size = 20), legend.position = "bottom", 
		      legend.title = element_blank()) +
		scale_color_manual(values = colpal) +
		xlab("Sample size") +
		ylab(ylab[stat]) +
		ylim(0, 1) +
		ggtitle(plot_title, subtitle =  subtitle)

	dir.create(path = "plot_kramer/", showWarnings = FALSE)
	ggsave(filename = paste0("plot_kramer/", stat, "_", algorithm, ".pdf"))
}

plot_kramer(algorithm = "diagdom_025", dir_name = "res_kramer_0.25_domdiag",
	plot_title = "diagonal dominance", 
	subtitle = "d = 0.25",
	stat = "ppv")
plot_kramer(algorithm = "diagdom_025", dir_name = "res_kramer_0.25_domdiag",
	plot_title = "diagonal dominance", 
	subtitle = "d = 0.25",
	stat = "tpr")

plot_kramer(algorithm = "port_chol_025", dir_name = "res_kramer_0.25_port_chol",
	plot_title = "uniform + partial orthogonalization", 
	subtitle = "d = 0.25",
	stat = "ppv")
plot_kramer(algorithm = "port_chol_025", dir_name = "res_kramer_0.25_port_chol",
	plot_title = "uniform  + partial orthogonalization",
	subtitle = "d = 0.25",
	stat = "tpr")

plot_kramer(algorithm = "diagdom_005", dir_name = "res_kramer_0.05_domdiag",
            plot_title = "diagonal dominance",
            subtitle = "d = 0.05",
            stat = "ppv")
plot_kramer(algorithm = "diagdom_005", dir_name = "res_kramer_0.05_domdiag",
            plot_title = "diagonal dominance",
            subtitle = "d = 0.05",
            stat = "tpr")

plot_kramer(algorithm = "port_chol_005", dir_name = "res_kramer_0.05_port_chol",
            plot_title = "uniform + partial orthogonalization",
            subtitle = "d = 0.05",
            stat = "ppv")
plot_kramer(algorithm = "port_chol_005", dir_name = "res_kramer_0.05_port_chol",
            plot_title = "uniform + partial orthogonalization",
            subtitle = "d = 0.05",
            stat = "tpr")

