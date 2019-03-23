# Simulation of covariance and concentration graph matrices

This repository contains the files for replicating the experiments described in
the paper

> Córdoba I., Varando G., Bielza C., Larrañaga P. On Gaussian graphical model simulations.  

They are mainly concerned with the method of partial orthogonalization (Córdoba
et al. 2018), implemented in `gmat::port()`, as well as the traditional diagonal
dominance method, implemented in many software packages, and also in
`gmat::diagdom()`.

The experiments in the following paper

> N. Krämer, J. Schäfer, and A.-L. Boulesteix. Regularized estimation of
> large-scale gene association networks using graphical Gaussian models. 
> BMC Bioinformatics, 10(1):384, 2009

have also been used in Córdoba et al. (2018) to validate both approaches, and
the code for its replication is also available in this repository.

## Contents

- `sim_experiment.R`: script that executes both methods for different matrix
  dimensions and sample sizes, saving the generated samples.
- `time_experiment.R`: script that executes both methods for different matrix
  dimensions and sample sizes, measuring and saving their execution time.
- `kramer_experiment.R`: script that replicates the experiments in Krämer and
  Schäfer (2009) whose results are also included in Córdoba et al. (2018).
- `performance.pcor.R`: same as [parcor::performance.pcor](https://github.com/cran/parcor/blob/master/R/performance.pcor.R), but calling `GeneNet::network.test.edges()` instead of `GeneNet::ggm.test.edges()`, which does not exist in the newest version of `GeneNet`. This file can be safely ignored as it will be removed when/if `parcor` is fixed.
- `plot_utils.R`: utility functions for plotting.
- `plot.R`: script that generates the plots describing the results of both the
  simulation and time experiments.
- `plot_kramer.R`: script that generates the plots corresponding to the Kramer
  experiment.
- `opt`: folder containing scripts for running additional experiments. __Work in
  progress__

## Instructions for simulation and time experiments

- R packages required: `doParallel`, `foreach`, `gmat`, `ggplot2`, `Matrix` and
  `reshape2`.
- Run the following commands from a terminal (or source the files on an open R session)
	```bash
	Rscript sim_experiment.R
	Rscript time_experiment.R
	Rscript plot.R
	```
Both the simulation and time experiment are computationally intensive.
