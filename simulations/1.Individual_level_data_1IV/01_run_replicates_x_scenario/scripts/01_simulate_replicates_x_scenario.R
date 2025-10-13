#!/usr/bin/env Rscript

library(tidyr, quietly = T, warn.conflicts = F)
library(dplyr, quietly = T, warn.conflicts = F)
library(ggplot2, quietly = T, warn.conflicts = F)
library(cowplot, quietly = T, warn.conflicts = F)
library(pheatmap, quietly = T, warn.conflicts = F)
library(dichromat, quietly = T, warn.conflicts = F)
library(sessioninfo, quietly = T, warn.conflicts = F)

rm(list = ls())
# Set seed at the beginning of scenario run (outside replication for-loop)
# to ensure within-run replicate variability of a run, but consistent results across reruns
set.seed(09222025)
options(scipen = 0)

# ------------------------------------------------------------------------------
#                   1. Run simulation replicates per scenario 
# ------------------------------------------------------------------------------
#  This script takes a scenario parameters and simulates it n_rep times: 
#  it simulates genotype (G) for 1 IV, an unknown confounder (U), exposure (X), 
#  outcome (Y), and error terms (eX and eY) for N individuals across the two
#  strata of K. It tests for HWE, G-X and G-Y associations within each stratum
#  and estimates per-stratum causal effects.
# ______________________________________________________________________________

## Define dirs
setwd("/localhome/daianna/Stratified_MR_simulation")
input_dir <- paste(getwd(), "simulations", "1.Individual_level_data_1IV", "01_run_replicates_x_scenario", "inputs", sep = "/")
out_dir <- paste(getwd(), "simulations", "1.Individual_level_data_1IV", "01_run_replicates_x_scenario", "outputs", sep = "/")
plot_dir <- paste(getwd(), "simulations", "1.Individual_level_data_1IV", "01_run_replicates_x_scenario", "plots", sep = "/")

dir.create(input_dir, showWarnings = F)
dir.create(out_dir, showWarnings = F)
dir.create(plot_dir, showWarnings = F)

## Input dirs
input_dir00 <- paste(getwd(), "simulations", "1.Individual_level_data_1IV", "00_define_scenarios", "outputs", sep = "/")


## Load helper functions for simulation:
source("simulations/helper_functions/simulate_strata.R")
source("simulations/helper_functions/simulate_genotype.R")
source("simulations/helper_functions/HWE_test.R")
source("simulations/helper_functions/simulate_confounder.R")
source("simulations/helper_functions/simulate_error.R")
source("simulations/helper_functions/simulate_exposure.R")
source("simulations/helper_functions/test_genetic_association.R")
source("simulations/helper_functions/simulate_outcome.R")


##########   Plotting functions (for internal use in script only)   ##########

# 1. Plot phenotype ~ predictor with estimated and true effects in both strata
plot_pheno_vs_predictor <- function(indiv_data, pheno, predictor, out){
  
  labels = c("X" = expression(X),
             "pred_X" = expression(hat(X) == hat(beta[GX])*G),
             "Y" = expression(Y),
             "pred_Y" = expression(hat(Y) == hat(beta[GY])*G),
             "G" = expression(G), 
             "U" = expression(U), 
             "eX" = expression(epsilon[X]), 
             "eY" = expression(epsilon[Y]))
  
  
  p = ggplot(indiv_data, aes(x = get(predictor), y = get(pheno))) + 
    facet_grid(cols = vars(K), 
               labeller = labeller(K = c("1" = "Stratum 1", "2" = "Stratum 2"))) +
    geom_point(color='gray', alpha = 0.3, size = 1) +
    theme_classic() +
    labs(x = labels[predictor], y = labels[pheno]) +
    theme(axis.title = element_text(size = 8),
          axis.text = element_text(size = 7),
          strip.text.x = element_text(size = 8), 
          axis.line = element_line(linewidth = 0.4),
          axis.ticks = element_line(linewidth = 0.3),
          strip.text = element_text(size = 8.1, face = 2),
          strip.background = element_rect(fill="white", color = "gray", linewidth = 0.7),
          panel.spacing = unit(0.4, "lines"))
  
  if(predictor == "G"){
    p <- p + scale_x_continuous(breaks = c(0,1,2), labels = c(0,1,2)) 
    
    if(pheno == "X" | pheno == "pred_X"){
      betas = data.frame("hat_beta" = c(out$hat_BGX1, out$hat_BGX2), "K" = c(1,2))
      betas$beta = c(out$BGX1, out$BGX2)
      
      ## Add β̂ɢxₖ
      p <- p + geom_abline(data = betas, aes(slope = hat_beta, intercept = 0), linewidth = 0.65, alpha = 0.5, color = "orangered2") +
        geom_text(data = betas, aes(label = paste("hat(beta)[GX] ==", signif(hat_beta, digits = 2)), 
                                    x = 1.75, y = min(indiv_data[, pheno]) + ((max(indiv_data[, pheno]) - min(indiv_data[, pheno]))*0.1)), color = "orangered2", size = 3, parse = T) +
        ## Add true βɢxₖ
        geom_abline(data = betas, aes(slope = beta, intercept = 0), linewidth = 0.65, alpha = 0.5, color = "purple3") +
        geom_text(data = betas, aes(label = paste("beta[GX] ==", signif(beta, digits = 2)), 
                                    x = 1.75, y = min(indiv_data[, pheno]) + ((max(indiv_data[, pheno]) - min(indiv_data[, pheno]))*0.2)), color = "purple3", size = 3, parse = T) 
    }
    
    else if(pheno == "Y" | pheno == "pred_Y"){
      betas <- data.frame("hat_beta" = c(out$hat_BGY1, out$hat_BGY2), "K" = c(1,2))
      betas$beta <- c(out$BGY1, out$BGY2)
      
      ## Add estimated effects β̂ɢʏₖ
      p <- p + geom_abline(data = betas, aes(slope = hat_beta, intercept = 0), linewidth = 0.65, alpha = 0.5, color = "orangered2") +
        geom_text(data = betas, aes(label = paste("hat(beta)[GY] ==", signif(hat_beta, digits = 2)), 
                                    x = 1.75, y = min(indiv_data[, pheno]) + ((max(indiv_data[, pheno]) - min(indiv_data[, pheno]))*0.10)), color = "orangered2", size = 3, parse = T) +
        ## Add true βɢʏₖ = βɢxₖ x βxʏₖ 
        geom_abline(data = betas, aes(slope = beta, intercept = 0), linewidth = 0.65, alpha = 0.5, color = "purple3") +
        geom_text(data = betas, aes(label = paste("beta[GY] ==", signif(beta, digits = 2)), 
                                    x = 1.75, y = min(indiv_data[, pheno]) + ((max(indiv_data[, pheno]) - min(indiv_data[, pheno]))*0.2)), color = "purple3", size = 3, parse = T) 
    }
  }
  
  else if(predictor == "U"){
    
    if(pheno == "X"){
      slope = out$BUX
      label = paste("beta[UX] == ", out$BUX)
    }
    else if(pheno == "pred_X"){
      slope = 0
      label = paste("beta[U][hat(X)] == 0")
    }
    else if(pheno == "Y"){
      slope = out$BUY
      label = paste("beta[UY] == ", out$BUY)
    }
    else if(pheno == "pred_Y"){
      slope = 0
      label = paste("beta[U][hat(Y)] == 0")
    }
    
    p <- p + geom_abline(slope = slope, intercept = 0, linewidth = 0.65, alpha = 0.7, color = "gray30") +
      geom_text(x = 0.75, y = min(indiv_data[, pheno]) + ((max(indiv_data[, pheno]) - min(indiv_data[, pheno]))/2.5), label = label, color = 'gray30', size = 3, parse = T)
  }
  
  else if(predictor == "eX"){
    
    if(pheno == "X"){
      slope = 1
      label = paste("beta[epsilon[X]] == 1")
    }
    else if(pheno == "pred_X"){
      slope = 0
      label = paste("beta[epsilon[X]] == 0")
    }
    else if(pheno == "Y"){
      slope = 0
      label = paste("beta[epsilon[X]] == 0")
    }
    else if(pheno == "pred_Y"){
      slope = 0
      label = paste("beta[epsilon[X]] == 0")
    }
    
    p <- p + geom_abline(slope = slope, intercept = 0, linewidth = 0.65, alpha = 0.7, color = "gray30") +
      geom_text(x = 2, y = min(indiv_data[, pheno]) + ((max(indiv_data[, pheno]) - min(indiv_data[, pheno]))/2.5), label = label, color = 'gray30', size = 3, parse = T)
  }
  
  else if(predictor == "eY"){
    
    if(pheno == "X"){
      slope = 0
      label = paste("beta[epsilon[Y]] == 0")
    }
    else if(pheno == "pred_X"){
      slope = 0
      label = paste("beta[epsilon[Y]] == 0")
    }
    else if(pheno == "Y"){
      slope = 1
      label = paste("beta[epsilon[Y]] == 1")
    }
    else if(pheno == "pred_Y"){
      slope = 0
      label = paste("beta[epsilon[Y]] == 0")
    }
    
    p <- p + geom_abline(slope = slope, intercept = 0, linewidth = 0.65, alpha = 0.7, color = "gray30") +
      geom_text(x = 2, y = min(indiv_data[, pheno]) + ((max(indiv_data[, pheno]) - min(indiv_data[, pheno]))/2.5), label = label, color = 'gray30', size = 3, parse = T)
  }
  
  else if(predictor == "X"){
    
    betas = data.frame("hat_beta" = c(out$hat_BXY1, out$hat_BXY2), "K" = c(1,2))
    betas$beta = c(out$BXY1, out$BXY2)
    
    ## Add β̂xʏₖ
    p <- p + geom_abline(data = betas, aes(slope = hat_beta, intercept = 0), linewidth = 0.65, alpha = 0.5, color = "orangered2") +
      geom_text(data = betas, aes(label = paste("hat(beta)[hat(X)][Y] ==", signif(hat_beta, digits = 2)), 
                                  x = max(indiv_data[, predictor]) - ((max(indiv_data[, predictor]) - min(indiv_data[, predictor]))/8),
                                  y = min(indiv_data[, pheno]) + ((max(indiv_data[, pheno]) - min(indiv_data[, pheno]))*0.10)), color = "orangered2", size = 3, parse = T) +
      ## Add true βxʏₖ
      geom_abline(data = betas, aes(slope = beta, intercept = 0), linewidth = 0.65, alpha = 0.5, color = "purple3") +
      geom_text(data = betas, aes(label = paste("beta[hat(X)][Y] ==", signif(beta, digits = 2)), 
                                  x = max(indiv_data[, predictor]) - ((max(indiv_data[, predictor]) - min(indiv_data[, predictor]))/8),
                                  y = min(indiv_data[, pheno]) + ((max(indiv_data[, pheno]) - min(indiv_data[, pheno]))*0.2)), color = "purple3", size = 3, parse = T) 
  }
  
  else if(predictor == "pred_X"){
    
    betas = data.frame("hat_beta" = c(out$hat_BXY1, out$hat_BXY2), "K" = c(1,2))
    betas$beta = c(out$BXY1, out$BXY2)
    
    ## Add β̂xʏₖ
    p <- p + geom_abline(data = betas, aes(slope = hat_beta, intercept = 0), linewidth = 0.65, alpha = 0.5, color = "orangered2") +
      geom_text(data = betas, aes(label = paste("hat(beta)[hat(X)][Y] ==", signif(hat_beta, digits = 2)), 
                                  x = max(indiv_data[, predictor]) - ((max(indiv_data[, predictor]) - min(indiv_data[, predictor]))/8),
                                  y = min(indiv_data[, pheno]) + ((max(indiv_data[, pheno]) - min(indiv_data[, pheno]))*0.10)), color = "orangered2", size = 3, parse = T) +
      ## Add true βxʏₖ
      geom_abline(data = betas, aes(slope = beta, intercept = 0), linewidth = 0.65, alpha = 0.5, color = "purple3") +
      geom_text(data = betas, aes(label = paste("beta[hat(X)][Y] ==", signif(beta, digits = 2)), 
                                  x = max(indiv_data[, predictor]) - ((max(indiv_data[, predictor]) - min(indiv_data[, predictor]))/8), 
                                  y = min(indiv_data[, pheno]) + ((max(indiv_data[, pheno]) - min(indiv_data[, pheno]))*0.2)), color = "purple3", size = 3, parse = T) 
  }
  
  return(p)
}


# 2. Call plot_pheno_vs_predictor() for multiple predictors
plot_pheno_vs_predictors <- function(indiv_data, out, plot_dir_sc_replicate){
  
  phenotypes <- c("X", "pred_X", "Y", "pred_Y")
  predictors <- list("X" = c("G", "U", "eX"), 
                     "pred_X" = c("G", "U", "eX"),
                     "Y" = c("X", "pred_X", "G", "U", "eX", "eY"),
                     "pred_Y" = c("X", "pred_X", "G", "U", "eX", "eY"))
  plots <- list()
  i = 1
  
  for(pheno in phenotypes){
    for(predictor in predictors[[pheno]]){
      plots[[i]] <- plot_pheno_vs_predictor(indiv_data, pheno, predictor, out)
      i = i + 1
    }
  }
  
  plot_grid(plotlist = plots, ncol = 3, align = "vh")
  ggsave(filename = paste0(plot_dir_sc_replicate, "/X_Y_on_predictors.pdf"), width = 13, height = 14)
}


# 3. Plot corr between X/Y and their predictors
cor_pheno_predictors <- function(indiv_data, plot_dir_sc_replicate){
  
  variables <- c("X", "pred_X", "Y", "pred_Y", "G", "U", "eX", "eY")
  labs <- c("X", expression(hat(X) == hat(beta[GX])*G), "Y", expression(hat(Y) == hat(beta[GY])*G), 
            "G", "U", expression(epsilon[X]), expression(epsilon[Y]))
  
  h <- list()
  for(k in 1:2){
    
    data = subset(indiv_data, K == k)
    
    h[[k]] = pheatmap(cor(data[, variables]), 
                      color = colorRampPalette(c('#FFF0F5', '#EE0000'))(50),
                      display_numbers = T, 
                      cluster_rows = F, 
                      cluster_cols = F, 
                      main = paste("Stratum", k), 
                      border_color = "gray20", 
                      number_color = "black", 
                      fontsize = 10,
                      labels_row = labs,
                      labels_col = labs, 
                      angle_col = 90)$gtable
  }
  
  plot_grid(plotlist = h, ncol = 2, align = "h")
  ggsave(file = paste0(plot_dir_sc_replicate, "/Corr_predictors.pdf"), width = 9, height = 4)
}



##########    Main simulation function    ##########

## Parameter definition:
#' @param main_scenario string for main scenario name (one of 00, 01, 10, 11).
#' @param main_scenario_val value for varying main parameter (diff_BGX and/or diff_BXY)
#' @param sub_sce_sub_sce_varying_par string for sub scenario 1 name, according to varying parameter(s).
#' @param sub_sce_varying_par_value string for sub scenario 2 name, according to the value(s) the varying parameter(s) take. 
#' @param N int for total sample size.
#' @param r double for ratio of stratum sample sizes.
#' @param q1 double for MAF of IV in stratum 1.
#' @param q2 double for MAF of IV in stratum 2.
#' @param BGX1 double for true effect of IV on X in stratum 1.
#' @param diff_BGX double for true difference in the effect of IV on X between stratum 2 and 1.
#' @param BXY1 double for true casual effect in stratum 1.
#' @param diff_BXY double for true difference in the causal effect between stratum 2 and 1.
#' @param BUX int for effect of confounder on X (assumed to be the same across strata).
#' @param BUY int for effect of confounder on Y (assumed to be the same across strata).
#' @param replicate int for replicate number.

simulation_indiv_data_1IV <- function(main_scenario, main_scenario_val, sub_sce_varying_par, sub_sce_varying_par_value, N, r, q1, q2, BGX1, diff_BGX, BXY1, diff_BXY, BUX, BUY, replicate){
  
  ## Output subdirs x main_scenario x replicate
  out_dir_sc <- paste(out_dir, main_scenario, main_scenario_val, sub_sce_varying_par, sub_sce_varying_par_value, sep = "/")
  dir.create(out_dir_sc, recursive = T, showWarnings = F)

  ## df to save results (and inputs)
  out <- data.frame("main_scenario" = main_scenario, "main_scenario_val" = main_scenario_val, "sub_sce_varying_par" = sub_sce_varying_par, "sub_sce_varying_par_value" = sub_sce_varying_par_value, 
                    "replicate" = replicate, "N" = N, "r" = r, "q1" = q1, "q2" = q2, 
                    "BGX1" = BGX1, "diff_BGX" = diff_BGX, "BGX2" = BGX1 + diff_BGX, 
                    "BXY1" = BXY1, "diff_BXY" = diff_BXY, "BXY2" = BXY1 + diff_BXY, "BUX" = BUX, "BUY" = BUY)
  
  #  Step 1: Simulate strata: K | N, r
  # ___________________________________________________________________________
  K = simulate_strata(N, r)
  
  ## Confirm N1/N2 = r
  if(round(table(K)["1"] / table(K)["2"], digits = 2) != r){
    message("Created strata don't satisfy r ratio")
    stop()
  }
  ## Stratum sample sizes
  N1 = table(K)["1"] 
  N2 = table(K)["2"]
  
  out <- cbind(out, "N1" = N1, "N2" = N2)
  
  #  Step 2: Simulate genotype for each stratum: Gₖ ~ Binom(2, qₖ) | Nₖ, qₖ
  # ___________________________________________________________________________
  G1 = simulate_genotype(N1, q1)
  G2 = simulate_genotype(N2, q2)
  G = c(G1, G2)
  
  out <- cbind(out, 
               "n_aa_1" = unname(table(G1)["2"]), 
               "n_aA_1" = unname(table(G1)["1"]), 
               "n_AA_1" = unname(table(G1)["0"]),
               "n_aa_2" = unname(table(G2)["2"]), 
               "n_aA_2" = unname(table(G2)["1"]), 
               "n_AA_2" = unname(table(G2)["0"]))
  
  # -----------------------------------------------------------------------
  #  2.1: Chi-squared test for HWE deviations in each stratum and globally
  # -----------------------------------------------------------------------
  
  for(k in c(1, 2, "")){
    
    Gk = get(paste0("G", k))
    n_aa = table(Gk)["2"]
    n_aA = table(Gk)["1"]
    n_AA = table(Gk)["0"]
    
    HWE_res = HWE_test(n_aa, n_aA, n_AA) 
    q_ob <- (n_aA + (2*n_aa))/ (2*length(Gk))
    
      if(k != ""){
    #   print(paste0("Observed MAF in stratum ", k, ": ", q_ob))
    #   print(paste0("P-val for HWE test in stratum ", k, ": ", HWE_res$P.value))
    #   if(HWE_res$P.value <= 0.05){
    #     message(paste0("Warning: SNP deviates from HWE in stratum ", k, "."))
    #   }
      
      out[,paste0("q", k, "_ob")] = q_ob
      out[,paste0("HWE_CHISQ_", k)] = HWE_res$CHISQ
      out[,paste0("HWE_P_", k)] = HWE_res$P.value
      
    } else{
      # print(paste0("Observed MAF in whole population: ", q_ob))
      # print(paste0("P-val for HWE test in whole population: ", HWE_res$P.value))
      # if(HWE_res$P.value <= 0.05){
      #   message("Warning: SNP deviates from HWE in whole population.")
      # }
      
      out[,"q_global_ob"] = q_ob
      out[,"HWE_CHISQ_global"] = HWE_res$CHISQ
      out[,"HWE_P_global"] = HWE_res$P.value
    }
  }
  # -----------------------------------------------------------------------
  
  indiv_data <- data.frame("K" = K, "G" = G)
  
  
  # Step 3: Simulate unknown confounder U ~ Unif(0,1) | N  
  # ___________________________________________________________________
  U = simulate_confounder(N)
  indiv_data <- cbind(indiv_data, "U" = U)
  
  # Step 4: Simulate error terms for X and Y:  εx,εʏ ~ N(0,1) | N  
  # ___________________________________________________________________
  eX = simulate_error(N)
  eY = simulate_error(N)
  
  indiv_data <- cbind(indiv_data, "eX" = eX, "eY" = eY)
  
  # Step 5: Generate exposure Xₖ = βɢxₖ(G) + U + εx
  # ___________________________________________________________________
  ## Define βɢx₂
  BGX2 = BGX1 + diff_BGX
  
  ## Add βɢx to use according to k
  indiv_data$BGX <- c(rep(BGX1, N1), rep(BGX2, N2))
  
  X = simulate_exposure(indiv_data, cols = c("K", "G", "U", "eX", "BGX"), BUX)
  indiv_data <- cbind(indiv_data, "X" = X)
  
  
  # Step 6: Test G effect on X: get β̂ɢxₖ
  # ___________________________________________________________________
  fit_X = test_genetic_association(indiv_data, cols = c("K", "G", "X"))
  
  hat_BGX1 = fit_X$Estimate[1]
  hat_BGX2 = fit_X$Estimate[2]
  se_hat_BGX1 = fit_X$`Std. Error`[1]
  se_hat_BGX2 = fit_X$`Std. Error`[2]
  BGX_t_stat_1 = fit_X$`t value`[1]
  BGX_t_stat_2 = fit_X$`t value`[2]
  BGX_P_1 = fit_X$`Pr(>|t|)`[1]
  BGX_P_2 = fit_X$`Pr(>|t|)`[2]
  
  out <- cbind(out, "hat_BGX1" = hat_BGX1, "hat_BGX2" = hat_BGX2, 
               "se_hat_BGX1" = se_hat_BGX1, "se_hat_BGX2" = se_hat_BGX2, 
               "BGX_t_stat_1" = BGX_t_stat_1, "BGX_t_stat_2" = BGX_t_stat_2, 
               "BGX_P_1" = BGX_P_1, "BGX_P_2" = BGX_P_2)
  
  # ----------------------------------------------------------------------
  #  6.1: Confirm SNP is strong IV in both strata
  # ----------------------------------------------------------------------
  if(any(fit_X$Fstat < 10)){
    message("Warning: weak instrument.")
  }
  out <- cbind(out, "BGX_F_stat_1" = fit_X$Fstat[1],  
                    "BGX_F_stat_2" = fit_X$Fstat[2])
  # ----------------------------------------------------------------------
  
  ## Add pred X: X̂ = β̂ɢxₖ(G)
  indiv_data$hat_BGX <- c(rep(hat_BGX1, N1), rep(hat_BGX2, N2))
  indiv_data$pred_X <- indiv_data$G * indiv_data$hat_BGX
  
  # *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
  #                !!!  Sanity checks  !!!                |
  # *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
  # Sanity check 1:
  # -------------------------------------------------------
  # True vs estimated effect in each stratum |βɢxₖ- β̂ɢxₖ|
  diff_BGX1 = abs(BGX1 - hat_BGX1)
  diff_BGX2 = abs(BGX2 - hat_BGX2)
  
  out <- cbind(out, "diff_BGX1" = diff_BGX1, "diff_BGX2" = diff_BGX2)
  
  # Sanity check 2:
  # -------------------------------------------------------
  # Detected genetic effect difference Δβɢx = βɢx₂ - βɢx₁?
  # Δβ̂ɢx = β̂ɢx₂ - β̂ɢx₁
  hat_diff_BGX = hat_BGX2 - hat_BGX1
  # Difference between true vs estimated delta |Δβɢx - Δβ̂ɢx|
  diff_diff_BGX = abs(diff_BGX - hat_diff_BGX)
  
  out <- cbind(out, "hat_diff_BGX" = hat_diff_BGX, "diff_diff_BGX" = diff_diff_BGX)
  
  # Sanity check 3:
  # -------------------------------------------------------
  # Significant G x K interaction on X?
  
  # Z-stat = (β̂ɢx₂ - β̂ɢx₁) / √se(β̂ɢx₂)² + se(β̂ɢx₁)² 
  Z_diff_BGX = (hat_BGX2 - hat_BGX1) / sqrt((se_hat_BGX2^2) + (se_hat_BGX1^2))
  # P val: Z-stat ~ N(0,1) | Ho
  p_Z_diff_BGX = 2*pnorm(abs(Z_diff_BGX), 0, 1, lower.tail = F)
  
  out <- cbind(out, "Z_diff_BGX" = Z_diff_BGX, 
                    "p_Z_diff_BGX" = p_Z_diff_BGX)
  
  # *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
  
  # Step 7: Generate outcome Yₖ = βxʏₖ(βɢxₖ)(G) + U + εʏ
  # ___________________________________________________________________
  ## Define βxʏ₂
  BXY2 = BXY1 + diff_BXY
  
  ## Add βxʏ to use according to k
  indiv_data$BXY <- c(rep(BXY1, N1), rep(BXY2, N2))
  
  Y = simulate_outcome(indiv_data, cols = c("K", "G", "U", "eY", "BGX", "BXY"), BUX, BUY)
  indiv_data <- cbind(indiv_data, "Y" = Y)
  
  # Step 8: Test G association with Y: get β̂ɢʏₖ
  # ___________________________________________________________________
  fit_Y = test_genetic_association(indiv_data, cols = c("K", "G", "Y"))
  
  hat_BGY1 = fit_Y$Estimate[1]
  hat_BGY2 = fit_Y$Estimate[2]
  se_hat_BGY1 = fit_Y$`Std. Error`[1]
  se_hat_BGY2 = fit_Y$`Std. Error`[2]
  BGY_t_stat_1 = fit_Y$`t value`[1]
  BGY_t_stat_2 = fit_Y$`t value`[2]
  BGY_P_1 = fit_Y$`Pr(>|t|)`[1]
  BGY_P_2 = fit_Y$`Pr(>|t|)`[2]
  
  out <- cbind(out, "hat_BGY1" = hat_BGY1, "hat_BGY2" = hat_BGY2, 
               "se_hat_BGY1" = se_hat_BGY1, "se_hat_BGY2" = se_hat_BGY2, 
               "BGY_t_stat_1" = BGY_t_stat_1, "BGY_t_stat_2" = BGY_t_stat_2, 
               "BGY_P_1" = BGY_P_1, "BGY_P_2" = BGY_P_2)
  
  ## Add Y pred by G: Ŷ =  β̂ɢʏₖ(G) 
  indiv_data$hat_BGY <- c(rep(hat_BGY1, N1), rep(hat_BGY2, N2))
  indiv_data$pred_Y <- indiv_data$G * indiv_data$hat_BGY
  
  # *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
  #                !!!  Sanity checks  !!!                |
  # *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
  # Sanity check 1:
  # -------------------------------------------------------
  # Estimated vs true effects in each stratum |βɢʏₖ- β̂ɢʏₖ|
  BGY1 = BGX1 * BXY1
  BGY2 = BGX2 * BXY2
  diff_BGY1 = abs(BGY1 - hat_BGY1)
  diff_BGY2 = abs(BGY2 - hat_BGY2)
  
  out <- cbind(out, "BGY1" = BGY1, "BGY2" = BGY2, "diff_BGY1" = diff_BGY1, "diff_BGY2" = diff_BGY2)
  
  # Sanity check 2:
  # -------------------------------------------------------
  # True genetic effect difference Δβɢʏ = βɢʏ₂ - βɢʏ₁
  diff_BGY = BGY2 - BGY1
  # Estimated genetic effect difference Δβ̂ɢʏ = β̂ɢʏ₂ - β̂ɢʏ₁
  hat_diff_BGY = hat_BGY2 - hat_BGY1
  # Difference between true vs estimated delta |Δβɢʏ - Δβ̂ɢʏ|
  diff_diff_BGY = abs(diff_BGY - hat_diff_BGY)
  
  out <- cbind(out, "diff_BGY" = diff_BGY, 
               "hat_diff_BGY" = hat_diff_BGY, 
               "diff_diff_BGY" = diff_diff_BGY)
  
  # Sanity check 3:
  # -------------------------------------------------------
  # Significant G x K interaction on Y?
  
  # Z-stat = (β̂ɢʏ₂ - β̂ɢʏ₁) / √se(β̂ɢʏ₂)² + se(β̂ɢʏ₁)² 
  Z_diff_BGY = (hat_BGY2 - hat_BGY1) / sqrt((se_hat_BGY2^2) + (se_hat_BGY1^2))
  # P val: Z-stat ~ N(0,1) | Ho
  p_Z_diff_BGY = 2*pnorm(abs(Z_diff_BGY), 0, 1, lower.tail = F)
  
  out <- cbind(out, "Z_diff_BGY" = Z_diff_BGY, 
                    "p_Z_diff_BGY" = p_Z_diff_BGY)
  
  # *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*

  
  # Step 9: Calculate causal effects: get β̂xʏₖ = β̂ɢʏₖ / β̂ɢxₖ
  # ___________________________________________________________________
  
  hat_BXY1 <- hat_BGY1 / hat_BGX1
  hat_BXY2 <- hat_BGY2 / hat_BGX2
  
  out <- cbind(out, "hat_BXY1" = hat_BXY1, "hat_BXY2" = hat_BXY2)
  
  # *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
  #                   Metrics of interest                 |
  # *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
  # Estimated vs true effects in each stratum |βxʏₖ- β̂xʏₖ|
  # -------------------------------------------------------
  diff_BXY1 = abs(BXY1 - hat_BXY1)
  diff_BXY2 = abs(BXY2 - hat_BXY2)
  
  out <- cbind(out, "diff_BXY1" = diff_BXY1, "diff_BXY2" = diff_BXY2)
  
  # Difference between true vs estimated delta |Δβxʏ - Δβ̂xʏ|
  # ----------------------------------------------------------
  # True causal effect difference: Δβxʏ = βxʏ₂ - βxʏ₁
  # Estimated causal effect difference: Δβ̂xʏ = β̂xʏ₂ - β̂xʏ₁
  hat_diff_BXY = hat_BXY2 - hat_BXY1
  # Difference 
  diff_diff_BXY = abs(diff_BXY - hat_diff_BXY)
  
  out <- cbind(out, "hat_diff_BXY" = hat_diff_BXY, "diff_diff_BXY" = diff_diff_BXY)
  
  # Significant X x K interaction on Y?
  # -------------------------------------------------------
  ## First order from delta expansion for ratio se: se(β̂xʏₖ) = se(β̂ɢʏₖ) / |β̂ɢxₖ|
  se_hat_BXY1 <- se_hat_BGY1 / abs(hat_BGX1)
  se_hat_BXY2 <- se_hat_BGY2 / abs(hat_BGX2)
  
  ## Second order from delta expansion for ratio se:
  # se(β̂xʏₖ) =  √ (se(β̂ɢʏₖ)² / β̂ɢxₖ²) + (β̂ɢʏₖ²*se(β̂ɢxₖ)² / β̂ɢxₖ⁴)
  # se_hat_BXY1_2nd <- sqrt((se_hat_BGY1^2 / hat_BGX1^2) + ((hat_BGY1^2)*(se_hat_BGX1^2) / (hat_BGX1^4)))
  # se_hat_BXY2_2nd <- sqrt((se_hat_BGY2^2 / hat_BGX2^2) + ((hat_BGY2^2)*(se_hat_BGX2^2) / (hat_BGX2^4)))
  
  # Z-stat = (β̂xʏ₂ - β̂xʏ₁) / √se(β̂xʏ₂)² + se(β̂xʏ₁)²
  Z_diff_BXY = (hat_BXY2 - hat_BXY1) / sqrt((se_hat_BXY2^2) + (se_hat_BXY1^2))
  # P val: Z-stat ~ N(0,1) | Ho
  p_Z_diff_BXY = 2*pnorm(abs(Z_diff_BXY), 0, 1, lower.tail = F)

  ## Z-stat and P val based on 2nd order se:
  # Z_diff_BXY_2nd = (hat_BXY2 - hat_BXY1) / sqrt((se_hat_BXY2_2nd^2) + (se_hat_BXY1_2nd^2))
  # p_Z_diff_BXY_2nd = 2*pnorm(abs(Z_diff_BXY_2nd), 0, 1, lower.tail = F)
  
  out <- cbind(out, "se_hat_BXY1" = se_hat_BXY1, 
                    "se_hat_BXY2" = se_hat_BXY2, 
                    "Z_diff_BXY" = Z_diff_BXY,
                    "p_Z_diff_BXY" = p_Z_diff_BXY)
  
  # out <- cbind(out, "se_hat_BXY1_2nd" = se_hat_BXY1_2nd, 
  #                   "se_hat_BXY2_2nd" = se_hat_BXY2_2nd, 
  #                   "Z_diff_BXY_2nd" = Z_diff_BXY_2nd,
  #                   "p_Z_diff_BXY_2nd" = p_Z_diff_BXY_2nd)
  
  if(replicate == 1){
    
    plot_dir_sc_replicate <- paste(plot_dir, main_scenario, main_scenario_val, sub_sce_varying_par, sub_sce_varying_par_value, replicate, sep = "/")
    dir.create(plot_dir_sc_replicate, recursive = T, showWarnings = F)
    
    # Plot: X ~ G, U, εx
    #       X̂ ~ G, U, εx
    #       Y ~ X, X̂, G, U, εx, εʏ
    #       Ŷ ~ X, X̂, G, U, εx, εʏ
    # plot_pheno_vs_predictors(indiv_data, out, plot_dir_sc_replicate)
    # Plot correlation between exposure/outcome and predictors
    cor_pheno_predictors(indiv_data, plot_dir_sc_replicate)
  }
  
  # *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
  
  ## Save indiv data generated
  save(indiv_data, file = paste0(out_dir_sc, "/indiv_data_rep", replicate, ".Rdata"))
  
  return(out)
}



# ---------------   Main script   ---------------  
## (uncomment according to main scenario simulations run)
## All 00 scenarios 
# scenarios <- get(load(paste0(input_dir00, "/scenario.00.Rdata")))

## All 01 scenarios
scenarios <- get(load(paste0(input_dir00, "/scenario.01.Rdata")))


## 100 replicates x scenario
n_rep = 100

## Run for task ID (= scenario row num)
id <- commandArgs(trailingOnly = TRUE)
sim_args <- scenarios[id, ]

## Scenario out 
out_dir_sc <- paste0(out_dir, "/", sim_args[1], "/", sim_args[2], "/", sim_args[3],  "/", sim_args[4])

## Skip if 100 replicates were run already for scenario
# if("scenario_res_across_100replicates.Rdata" %in% list.files(out_dir_sc) & length(list.files(out_dir_sc)) >100){
#   stop("Outputs already generated for scenario.", call. = T)
# } else{
#   
#   scenario_res = vector()
#   for(i in 1:n_rep){
#     scenario_res = rbind(scenario_res, do.call(simulation_indiv_data_1IV, c(sim_args, i)))
#   }
#   
#   ## Save scenario results across replicates
#   save(scenario_res, file = paste0(out_dir_sc, "/scenario_res_across_", n_rep, "replicates.Rdata"))
# }
  

scenario_res = vector()
for(i in 1:n_rep){
  scenario_res = rbind(scenario_res, do.call(simulation_indiv_data_1IV, c(sim_args, i)))
}

## Save scenario results across replicates
save(scenario_res, file = paste0(out_dir_sc, "/scenario_res_across_", n_rep, "replicates.Rdata"))



## Example test
# N = sim_args["N"] %>% as.numeric()
# r = sim_args["r"] %>% as.numeric()
# q1 = sim_args["q1"] %>% as.numeric()
# q2 = sim_args["q2"] %>% as.numeric()
# BGX1 = sim_args["BGX1"] %>% as.numeric()
# diff_BGX = sim_args["diff_BGX"] %>% as.numeric()
# BXY1 = sim_args["BXY1"] %>% as.numeric()
# diff_BXY = sim_args["diff_BXY"] %>% as.numeric()
# BUX = sim_args["BUX"] %>% as.numeric()
# BUY = sim_args["BUY"] %>% as.numeric()
# #
# main_scenario = sim_args["main_scenario"] %>% as.character()
# main_scenario_val = sim_args["main_scenario_val"] %>% as.character()
# sub_sce_varying_par = sim_args["sub_sce_varying_par"] %>% as.character()
# sub_sce_varying_par_value = sim_args["sub_sce_varying_par_value"] %>% as.character()
# 
# replicate = 1







## Reproducibility info
# session_info()

# ─ Session info ─────────────────────────────────────────────────────────────────────────────────────
# setting  value
# version  R version 4.5.1 (2025-06-13)
# os       macOS Monterey 12.7.6
# system   x86_64, darwin20
# ui       RStudio
# language (EN)
# collate  en_US.UTF-8
# ctype    en_US.UTF-8
# tz       Europe/London
# date     2025-10-01
# rstudio  2025.05.1+513 Mariposa Orchid (desktop)
# pandoc   NA
# quarto   1.6.42 @ /private/var/folders/j1/klb6dgpn3dxfdrxfmbb4f0jc0000gp/T/AppTranslocation/7D11FDEE-A995-4F90-9AC5-8C647803C12A/d/RStudio.app/Contents/Resources/app/quarto/bin/quarto
# 
# ─ Packages ─────────────────────────────────────────────────────────────────────────────────────────
# package      * version date (UTC) lib source
# cli            3.6.5   2025-04-23 [1] CRAN (R 4.5.0)
# cowplot      * 1.2.0   2025-07-07 [1] CRAN (R 4.5.1)
# dichromat    * 2.0-0.1 2022-05-02 [1] CRAN (R 4.5.0)
# dplyr        * 1.1.4   2023-11-17 [1] CRAN (R 4.5.0)
# farver         2.1.2   2024-05-13 [1] CRAN (R 4.5.0)
# fs             1.6.6   2025-04-12 [1] CRAN (R 4.5.0)
# generics       0.1.4   2025-05-09 [1] CRAN (R 4.5.0)
# ggplot2      * 4.0.0   2025-09-11 [1] CRAN (R 4.5.1)
# gitcreds       0.1.2   2022-09-08 [1] CRAN (R 4.5.0)
# glue           1.8.0   2024-09-30 [1] CRAN (R 4.5.0)
# gtable         0.3.6   2024-10-25 [1] CRAN (R 4.5.0)
# labeling       0.4.3   2023-08-29 [1] CRAN (R 4.5.0)
# lifecycle      1.0.4   2023-11-07 [1] CRAN (R 4.5.0)
# magrittr       2.0.4   2025-09-12 [1] CRAN (R 4.5.1)
# pheatmap     * 1.0.13  2025-06-05 [1] CRAN (R 4.5.0)
# pillar         1.11.1  2025-09-17 [1] CRAN (R 4.5.1)
# pkgconfig      2.0.3   2019-09-22 [1] CRAN (R 4.5.0)
# purrr          1.1.0   2025-07-10 [1] CRAN (R 4.5.1)
# R6             2.6.1   2025-02-15 [1] CRAN (R 4.5.0)
# RColorBrewer   1.1-3   2022-04-03 [1] CRAN (R 4.5.0)
# rlang          1.1.6   2025-04-11 [1] CRAN (R 4.5.0)
# rstudioapi     0.17.1  2024-10-22 [1] CRAN (R 4.5.0)
# S7             0.2.0   2024-11-07 [1] CRAN (R 4.5.0)
# scales         1.4.0   2025-04-24 [1] CRAN (R 4.5.0)
# sessioninfo  * 1.2.3   2025-02-05 [1] CRAN (R 4.5.0)
# tibble         3.3.0   2025-06-08 [1] CRAN (R 4.5.0)
# tidyr        * 1.3.1   2024-01-24 [1] CRAN (R 4.5.0)
# tidyselect     1.2.1   2024-03-11 [1] CRAN (R 4.5.0)
# usethis      * 3.2.1   2025-09-06 [1] CRAN (R 4.5.1)
# vctrs          0.6.5   2023-12-01 [1] CRAN (R 4.5.0)
# withr          3.0.2   2024-10-28 [1] CRAN (R 4.5.0)
# 
# [1] /Library/Frameworks/R.framework/Versions/4.5-x86_64/Resources/library
# * ── Packages attached to the search path.
# 
# ────────────────────────────────────────────────────────────────────────────────────────────────────

