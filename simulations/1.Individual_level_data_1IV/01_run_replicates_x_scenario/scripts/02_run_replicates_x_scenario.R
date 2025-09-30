


library(tidyr)
library(dplyr)
library(ggplot2)
library(cowplot)
library(pheatmap)
library(dichromat)

## Initialize
rm(list = ls())
time0 <- Sys.time()
set.seed(09222025)



################################################################################
#                     1. Simulation study using one IV TODO
################################################################################
# ______________________________________________________________________________
#   Data simulation and estimation of genetic and causal effects with one IV, 
#   across varying ΔβXY, N and N1/N2, q, and ΔβGX. TODO
# ______________________________________________________________________________

## Define dirs
input_dir <- paste(getwd(), "simulations", "one_IV_simulations", "inputs", sep = "/")
out_dir <- paste(getwd(), "simulations", "one_IV_simulations", "outputs", sep = "/")
plot_dir <- paste(getwd(), "simulations", "one_IV_simulations", "plots", sep = "/")

dir.create(input_dir, showWarnings = F)
dir.create(out_dir, showWarnings = F)
dir.create(plot_dir, showWarnings = F)

## Color dict TODO


## Helper functions:
source("simulations/helper_functions/simulate_strata.R")
source("simulations/helper_functions/simulate_genotype.R")
source("simulations/helper_functions/HWE_test.R")
source("simulations/helper_functions/simulate_confounder.R")
source("simulations/helper_functions/simulate_error.R")
source("simulations/helper_functions/simulate_exposure.R")
source("simulations/helper_functions/test_genetic_association.R")
source("simulations/helper_functions/simulate_outcome.R")

