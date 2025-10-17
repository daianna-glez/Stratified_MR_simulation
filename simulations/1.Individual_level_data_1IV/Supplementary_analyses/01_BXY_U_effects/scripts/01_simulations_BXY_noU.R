
library(tidyr)
library(dplyr)
library(sessioninfo)

rm(list = ls())
set.seed(09222025)

# ------------------------------------------------------------------------------
#                            Supplementary Analyses 
# ------------------------------------------------------------------------------
# ______________________________________________________________________________
#            1.0  Simulations with varying βxʏ₁ and βux effect on Y.
# ______________________________________________________________________________
#  This script repeats stratified MR simulations (scenario 0,0) across 
#  increasing values of βxʏ₁ for zero and non-zero βux on Y.
# ______________________________________________________________________________

input_dir <- paste(getwd(), "simulations", "1.Individual_level_data_1IV", "Supplementary_analyses", "01_BXY_U_effects", "inputs", sep = "/")
out_dir <- paste(getwd(), "simulations", "1.Individual_level_data_1IV", "Supplementary_analyses", "01_BXY_U_effects", "outputs", sep = "/")
plot_dir <- paste(getwd(), "simulations", "1.Individual_level_data_1IV", "Supplementary_analyses", "01_BXY_U_effects", "plots", sep = "/")

dir.create(input_dir, showWarnings = F)
dir.create(out_dir, showWarnings = F)
dir.create(plot_dir, showWarnings = F)

## Input dirs
input_dir00 <- paste(getwd(), "simulations", "1.Individual_level_data_1IV", "00_define_scenarios", "outputs", sep = "/")
input_dir01 <- paste(getwd(), "simulations", "1.Individual_level_data_1IV", "01_run_replicates_x_scenario", "outputs", sep = "/")

## Load helper functions for simulation:
source("simulations/helper_functions/simulate_strata.R")
source("simulations/helper_functions/simulate_genotype.R")
source("simulations/helper_functions/HWE_test.R")
source("simulations/helper_functions/simulate_confounder.R")
source("simulations/helper_functions/simulate_error.R")
source("simulations/helper_functions/simulate_exposure.R")
source("simulations/helper_functions/simulate_outcome.R")
source("simulations/helper_functions/test_genetic_association.R")

## Subset to βxʏ₁ scenarios in 0,0 
scenarios00_BXY1 <- get(load(paste0(input_dir00, "/scenario.00.Rdata"))) %>% dplyr::filter(sub_sce_varying_par %in% "BXY1")

## Add βux on Y values = {0.0  0.4  0.8  1.2  1.6  2.0}
BUXs = seq(from = 0, to = 2, by = 0.4)

scenarios = vector()
for(BUX in BUXs){
  t = scenarios00_BXY1
  t = cbind(t, "BUXonY" = BUX)
  scenarios = rbind(scenarios, t)
}

scenarios$sub_sce_varying_par = "BXY1_BUXonY"
scenarios$sub_sce_varying_par_value = paste0(scenarios$sub_sce_varying_par_value, ",BUXonY=", scenarios$BUXonY)


simulation_BUX_on_Y <- function(sim_args, replicate){
  
  main_scenario = sim_args[1]
  main_scenario_val = sim_args[2]
  sub_sce_varying_par = sim_args[3]
  sub_sce_varying_par_value = sim_args[4]
  N = sim_args[5] %>% as.numeric()
  r = sim_args[6] %>% as.numeric()
  q1 = sim_args[7] %>% as.numeric()
  q2 = sim_args[8] %>% as.numeric()
  BGX1 = sim_args[9] %>% as.numeric()
  diff_BGX = sim_args[10] %>% as.numeric()
  BXY1 = sim_args[11] %>% as.numeric()
  diff_BXY = sim_args[12] %>% as.numeric()
  BUX = sim_args[13] %>% as.numeric()
  BUY = sim_args[14] %>% as.numeric()
  BUX_on_Y = sim_args[15] %>% as.numeric()
  
  ## Output subdirs 
  out_dir_sc <- paste(out_dir, main_scenario, main_scenario_val, sub_sce_varying_par_value, sep = "/")
  dir.create(out_dir_sc, recursive = T, showWarnings = F)
  
  out <- data.frame("main_scenario" = main_scenario, "main_scenario_val" = main_scenario_val, "sub_sce_varying_par" = sub_sce_varying_par, "sub_sce_varying_par_value" = sub_sce_varying_par_value, 
                    "replicate" = replicate, "N" = N, "r" = r, "q1" = q1, "q2" = q2, 
                    "BGX1" = BGX1, "diff_BGX" = diff_BGX, "BGX2" = BGX1 + diff_BGX, 
                    "BXY1" = BXY1, "diff_BXY" = diff_BXY, "BXY2" = BXY1 + diff_BXY, 
                    "BUX" = BUX, "BUY" = BUY, "BUXonY" = BUX_on_Y)
  
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
  
  # Step 5: Generate exposure Xₖ = βɢxₖ(G) + βuxU + εx
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
  
  # Step 7: Generate outcome Yₖ = βxʏₖ(βɢxₖG + βux_on_ʏU) + βuʏU + εʏ
  # ___________________________________________________________________
  ## Define βxʏ₂
  BXY2 = BXY1 + diff_BXY
  
  ## Add βxʏ to use according to k
  indiv_data$BXY <- c(rep(BXY1, N1), rep(BXY2, N2))
  
  Y = simulate_outcome(indiv_data, cols = c("K", "G", "U", "eY", "BGX", "BXY"), BUX_on_Y, BUY)
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
  
  ## Add Y mean and var
  out <- cbind(out, "Y_mean" = mean(Y), "Y_var" = var(Y))
  
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
  se_hat_BXY1_2nd <- sqrt((se_hat_BGY1^2 / hat_BGX1^2) + ((hat_BGY1^2)*(se_hat_BGX1^2) / (hat_BGX1^4)))
  se_hat_BXY2_2nd <- sqrt((se_hat_BGY2^2 / hat_BGX2^2) + ((hat_BGY2^2)*(se_hat_BGX2^2) / (hat_BGX2^4)))
  
  # Z-stat = (β̂xʏ₂ - β̂xʏ₁) / √se(β̂xʏ₂)² + se(β̂xʏ₁)²
  Z_diff_BXY = (hat_BXY2 - hat_BXY1) / sqrt((se_hat_BXY2^2) + (se_hat_BXY1^2))
  # P val: Z-stat ~ N(0,1) | Ho
  p_Z_diff_BXY = 2*pnorm(abs(Z_diff_BXY), 0, 1, lower.tail = F)
  
  ## Z-stat and P val based on 2nd order se:
  Z_diff_BXY_2nd = (hat_BXY2 - hat_BXY1) / sqrt((se_hat_BXY2_2nd^2) + (se_hat_BXY1_2nd^2))
  p_Z_diff_BXY_2nd = 2*pnorm(abs(Z_diff_BXY_2nd), 0, 1, lower.tail = F)
  
  out <- cbind(out, "se_hat_BXY1" = se_hat_BXY1, 
                    "se_hat_BXY2" = se_hat_BXY2, 
                    "Z_diff_BXY" = Z_diff_BXY,
                    "p_Z_diff_BXY" = p_Z_diff_BXY)
  
  out <- cbind(out, "se_hat_BXY1_2nd" = se_hat_BXY1_2nd,
                    "se_hat_BXY2_2nd" = se_hat_BXY2_2nd,
                    "Z_diff_BXY_2nd" = Z_diff_BXY_2nd,
                    "p_Z_diff_BXY_2nd" = p_Z_diff_BXY_2nd)
  
  # *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
  
  return(out)
}



# ---------------   Main script   ---------------  
all_res <- vector()
all_metrics <- vector()
nrep = 100

for(j in 1:nrow(scenarios)){
  
  sim_args = scenarios[j, ]
  
  scenario_res = vector()
  for(i in 1:nrep){
    res = simulation_BUX_on_Y(sim_args, i)
    scenario_res = rbind(scenario_res, res)
  }
  
  ## Save scenario results across replicates
  out_dir_sc <- paste(out_dir, sim_args["main_scenario"], sim_args["main_scenario_val"], sim_args["sub_sce_varying_par_value"], sep = "/")
  save(scenario_res, file = paste0(out_dir_sc, "/scenario_res_across_", nrep, "replicates.Rdata"))
  all_res <- rbind(all_res, scenario_res)
  
  ## Summarize results across replicates
  scenario_metrics <- with(scenario_res, 
                           c("mean_hat_BGX1" = mean(hat_BGX1), "median_hat_BGX1" = median(hat_BGX1), "var_hat_BGX1" = var(hat_BGX1), 
                             "mean_hat_BGX2" = mean(hat_BGX2), "median_hat_BGX2" = median(hat_BGX2), "var_hat_BGX2" = var(hat_BGX2), 
                             
                             "FNR_BGX1" = mean(BGX_P_1 >= 0.05),"FNR_BGX2" = mean(BGX_P_2 >= 0.05),
                             "TPR_BGX1" = mean(BGX_P_1 < 0.05), "TPR_BGX2" = mean(BGX_P_2 < 0.05),
                             
                             "mean_diff_BGX1" = mean(diff_BGX1), "mean_diff_BGX2" = mean(diff_BGX2),
                             "mean_hat_diff_BGX" = mean(hat_diff_BGX), 
                             "mean_diff_diff_BGX" = mean(diff_diff_BGX), 
                             
                             "mean_hat_BGY1" = mean(hat_BGY1), "median_hat_BGY1" = median(hat_BGY1), "var_hat_BGY1" = var(hat_BGY1), 
                             "mean_hat_BGY2" = mean(hat_BGY2), "median_hat_BGY2" = median(hat_BGY2), "var_hat_BGY2" = var(hat_BGY2), 
                             
                             "mean_diff_BGY1" = mean(diff_BGY1), "mean_diff_BGY2" = mean(diff_BGY2),
                             "mean_hat_diff_BGY" = mean(hat_diff_BGY), 
                             "mean_diff_diff_BGY" = mean(diff_diff_BGY), 
                             
                             "mean_hat_BXY1" = mean(hat_BXY1), "median_hat_BXY1" = median(hat_BXY1), "var_hat_BXY1" = var(hat_BXY1), 
                             "mean_hat_BXY2" = mean(hat_BXY2), "median_hat_BXY2" = median(hat_BXY2), "var_hat_BXY2" = var(hat_BXY2), 
                             
                             "mean_diff_BXY1" = mean(diff_BXY1), "mean_diff_BXY2" = mean(diff_BXY2),
                             "mean_hat_diff_BXY" = mean(hat_diff_BXY), 
                             "mean_diff_diff_BXY" = mean(diff_diff_BXY),
                             
                             ## +/- rate for GxK interaction on X
                             "PR_GxK_on_X" = mean(p_Z_diff_BGX < 0.05),
                             "NR_GxK_on_X" = mean(p_Z_diff_BGX >= 0.05),
                             
                             ## +/- rate for G effect on Y in k=1,2
                             "PR_BGY1" = mean(BGY_P_1 < 0.05),
                             "NR_BGY1" = mean(BGY_P_1 >= 0.05),
                             "PR_BGY2" = mean(BGY_P_2 < 0.05),
                             "NR_BGY2" = mean(BGY_P_2 >= 0.05),
                             
                             ## +/- rate for GxK interaction on Y
                             "PR_GxK_on_Y" = mean(p_Z_diff_BGY < 0.05),
                             "NR_GxK_on_Y" = mean(p_Z_diff_BGY >= 0.05),
                             
                             ## +/- rate for XxK interaction on Y (based on 1st term se's)
                             "PR_XxK_on_Y" = mean(p_Z_diff_BXY < 0.05),
                             "NR_XxK_on_Y" = mean(p_Z_diff_BXY >= 0.05),
                             ## +/- rate for XxK interaction on Y (based on 2nd term se's)
                             "PR_XxK_on_Y_2nd" = mean(p_Z_diff_BXY_2nd < 0.05),
                             "NR_XxK_on_Y_2nd" = mean(p_Z_diff_BXY_2nd >= 0.05),
                             
                             "mean_mean_Y" = mean(Y_mean),
                             "mean_var_Y" = mean(Y_var)
                           ))
  
  save(scenario_metrics, file = paste(c(out_dir_sc, "/scenario_summary_metrics.Rdata"), collapse = "/")) 
  all_metrics <- rbind(all_metrics, cbind("main_scenario" = sim_args["main_scenario"], 
                                          "main_scenario_val" = sim_args["main_scenario_val"], 
                                          "sub_sce_varying_par" = sim_args["sub_sce_varying_par"], 
                                          "sub_sce_varying_par_value" = sim_args["sub_sce_varying_par_value"], 
                                          "BXY1" = sim_args["BXY1"],
                                          "BUXonY" = sim_args["BUXonY"],
                                          t(scenario_metrics)))
}


## Plots comparing metrics across subscenarios
plot_metric_across_subscenarios <- function(all_res, all_metrics, metric_y, metric_to_show, plot_boxplots = F){
  
  all_metrics[, metric_to_show] = signif(all_metrics[, metric_to_show], digits = 2)
  
  all_res$x <- all_res$BXY1
  all_res$x  <- factor(all_res$x, levels = unique(all_res$x))
  all_metrics$x <- all_metrics$BXY1
  all_metrics$x <- factor(all_metrics$x, levels = unique(all_metrics$x))
  xvar = "BXY1"
  all_res$x2 <- all_res$BUXonY
  all_res$x2  <- factor(all_res$x2, levels = unique(all_res$x2))
  col_var = "BUXonY"
    
  var_labs = list("BXY1" = expression(beta[XY[1]]),
                  "BUXonY" = expression(beta[UX]*" on "*Y))
  
  metric_to_show_lab = c("TPR_BGX1" = expression("TPR"), 
                         "TPR_BGX2" = expression("TPR"), 
                         "mean_diff_BGX1" = expression(bar(Delta)*beta[GX[1]]*":"), 
                         "mean_diff_BGX2" = expression(bar(Delta)*beta[GX[2]]*":"), 
                         "mean_hat_diff_BGX" = expression(bar(Delta)*hat(beta)[GX]*":"), 
                         "mean_diff_diff_BGX" = expression(bar(Delta)[Delta*beta[GX]]*":"), 
                         "NR_GxK_on_X" = expression("NR" ~ "for" ~ "GxK" ~ "on" ~ "X"),
                         "PR_GxK_on_X" = expression("PR" ~ "for" ~ "GxK" ~ "on" ~ "X"),
                         "PR_BGY1" = expression("PR"), 
                         "PR_BGY2" = expression("PR"), 
                         "NR_BGY1" = expression("NR"), 
                         "NR_BGY2" = expression("NR"), 
                         "mean_diff_BGY1" = expression(bar(Delta)*beta[GY[1]]*":"),
                         "mean_diff_BGY2" = expression(bar(Delta)*beta[GY[2]]*":"),
                         "mean_hat_diff_BGY" = expression(bar(Delta)*hat(beta)[GY]*":"), 
                         "mean_diff_diff_BGY" = expression(bar(Delta)[Delta*beta[GY]]*":"), 
                         "NR_GxK_on_Y" = expression("NR" ~ "for" ~ "GxK" ~ "on" ~ "Y"),
                         "PR_GxK_on_Y" = expression("PR" ~ "for" ~ "GxK" ~ "on" ~ "Y"),
                         "PR_BXY1" = expression("PR"), 
                         "PR_BXY2" = expression("PR"), 
                         "NR_BXY1" = expression("NR"), 
                         "NR_BXY2" = expression("NR"), 
                         "mean_hat_BXY1" = expression(bar(hat(beta))[XY[1]]),
                         "mean_hat_BXY2" = expression(bar(hat(beta))[XY[2]]), 
                         "mean_diff_BXY1" = expression(bar(Delta)*beta[XY[1]]*":"),
                         "mean_diff_BXY2" = expression(bar(Delta)*beta[XY[2]]*":"),
                         "mean_hat_diff_BXY" = expression(bar(Delta)*hat(beta)[XY]*":"), 
                         "mean_diff_diff_BXY" = expression(bar(Delta)[Delta*beta[XY]]*":"),
                         "NR_XxK_on_Y" = expression("NR" ~ "for" ~ "XxK" ~ "on" ~ "Y"),
                         "PR_XxK_on_Y" = expression("PR" ~ "for" ~ "XxK" ~ "on" ~ "Y")
  ) 
  
  yvar_labels = list("hat_BGX1" = expression(hat(beta[GX[1]])), 
                     "hat_BGX2" = expression(hat(beta[GX[2]])), 
                     "se_hat_BGX1" = expression(se(hat(beta[GX[1]]))), 
                     "se_hat_BGX2" = expression(se(hat(beta[GX[2]]))), 
                     "BGX_t_stat_1" = expression(t[hat(beta[GX[1]])]), 
                     "BGX_t_stat_2" = expression(t[hat(beta[GX[2]])]), 
                     "BGX_P_1" = expression(p[hat(beta[GX[1]])]), 
                     "BGX_P_2" = expression(p[hat(beta[GX[2]])]), 
                     "diff_BGX1" = expression(Delta*beta[GX[1]] == abs(beta[GX[1]] - hat(beta[GX[1]]))),
                     "diff_BGX2" = expression(Delta*beta[GX[2]] == abs(beta[GX[2]] - hat(beta[GX[2]]))),
                     "hat_diff_BGX" = expression(Delta*hat(beta[GX]) == hat(beta[GX[2]]) - hat(beta[GX[1]])), 
                     "diff_diff_BGX" = expression(Delta[Delta*beta[GX]] == abs(Delta*beta[GX] - Delta*hat(beta[GX]))),
                     "Z_diff_BGX" = expression(Z[Delta*beta[GX]]),
                     "p_Z_diff_BGX" = expression(p[Z[Delta*beta[GX]]]),
                     "hat_BGY1" = expression(hat(beta[GY[1]])),                
                     "hat_BGY2" = expression(hat(beta[GY[2]])),                      
                     "se_hat_BGY1" = expression(se(hat(beta[GY[1]]))),             
                     "se_hat_BGY2" = expression(se(hat(beta[GY[2]]))),              
                     "BGY_t_stat_1" = expression(t[hat(beta[GY[1]])]),              
                     "BGY_t_stat_2" = expression(t[hat(beta[GY[2]])]),             
                     "BGY_P_1" = expression(p[hat(beta[GY[1]])]),                   
                     "BGY_P_2" = expression(p[hat(beta[GY[2]])]),         
                     "diff_BGY1" = expression(Delta*beta[GY[1]] == abs(beta[GY[1]] - hat(beta[GY[1]]))),                    
                     "diff_BGY2" = expression(Delta*beta[GY[2]] == abs(beta[GY[2]] - hat(beta[GY[2]]))),   
                     "hat_diff_BGY" = expression(Delta*hat(beta[GY]) == hat(beta[GY[2]]) - hat(beta[GY[1]])),             
                     "diff_diff_BGY" = expression(Delta[Delta*beta[GY]] == abs(Delta*beta[GY] - Delta*hat(beta[GY]))),
                     "Z_diff_BGY" = expression(Z[Delta*beta[GY]]),       
                     "p_Z_diff_BGY" = expression(p[Z[Delta*beta[GY]]]),      
                     "hat_BXY1" = expression(hat(beta[XY[1]])),           
                     "hat_BXY2" = expression(hat(beta[XY[2]])),  
                     "se_hat_BXY1" = expression(se(hat(beta[XY[1]]))),             
                     "se_hat_BXY2" = expression(se(hat(beta[XY[2]]))),  
                     "se_hat_BXY1_2nd" = expression(se(hat(beta[XY[1]]))[2*"nd"]),             
                     "se_hat_BXY2_2nd" = expression(se(hat(beta[XY[2]]))[2*"nd"]),  
                     "diff_BXY1" = expression(Delta*beta[XY[1]] == abs(beta[XY[1]] - hat(beta[XY[1]]))),                 
                     "diff_BXY2" = expression(Delta*beta[XY[2]] == abs(beta[XY[2]] - hat(beta[XY[2]]))),  
                     "hat_diff_BXY" = expression(Delta*hat(beta[XY]) == hat(beta[XY[2]]) - hat(beta[XY[1]])),                      
                     "diff_diff_BXY" = expression(Delta[Delta*beta[XY]] == abs(Delta*beta[XY] - Delta*hat(beta[XY]))),
                     "Z_diff_BXY" = expression(Z[Delta*beta[XY]]),       
                     "p_Z_diff_BXY" = expression(p[Z[Delta*beta[XY]]]),
                     "Z_diff_BXY_2nd" = expression(Z[Delta*beta[XY]*2*"nd"]),       
                     "p_Z_diff_BXY_2nd" = expression(p[Z[Delta*beta[XY]]*2*"nd"])
  )
  
  p = ggplot(all_res, aes(x = x, y = get(metric_y))) + 
    geom_point(color='gray', alpha = 0.3, size = 1) +
    theme_classic() +
    facet_grid(.~x2, scales = "fixed", axes = "all_y") +
    labs(x = var_labs[[xvar]], 
         y = yvar_labels[[metric_y]], subtitle = var_labs[[col_var]]) +
    geom_text(x = 0, y = max(all_res[metric_y])+(max(all_res[metric_y]) - min(all_res[metric_y]))/13.5,
              label = as.character(metric_to_show_lab[metric_to_show]), size = 2.1, parse = T) +
    geom_text(data = all_metrics, aes(x = x, y = max(all_res[metric_y])+(max(all_res[metric_y]) - min(all_res[metric_y]))/25, 
                                       label = get(metric_to_show)), size = 1.7, angle = 90, hjust = 1.3) +
    coord_cartesian(clip = "off", ylim = c(min(all_res[metric_y]), max(all_res[metric_y]))) +
    theme(plot.subtitle = element_text(size = 8, hjust = 0.5),
          strip.background = element_blank(),
          panel.spacing.x = unit(0.25, "cm"),
          strip.text.x = element_text(size = 7, margin = margin(0,0,0.25,0, "cm")),
          axis.title.y = element_text(size = 9),
          axis.title.x = element_text(size = 8),
          axis.text = element_text(size = 7),
          axis.line = element_line(linewidth = 0.4),
          axis.ticks = element_line(linewidth = 0.3), 
          plot.margin = unit(c(1.3, 1, 1, 1), "lines"))
  
  ## To show yline
  if(metric_y %in% c("BGX_P_1", "BGX_P_2", "p_Z_diff_BGX", 
                     "BGY_P_1", "BGY_P_2", "p_Z_diff_BGY", 
                     "p_Z_diff_BXY", "p_Z_diff_BXY_2nd")){yintercept = 0.05}
  else if(metric_y == "hat_BGX1"){yintercept = unique(all_res$BGX1)}
  else if(metric_y == "hat_BGX2"){yintercept = unique(all_res$BGX2)}
  else if(metric_y == "hat_diff_BGX"){yintercept = unique(all_res$diff_BGX)}
  else if(metric_y == "hat_BGY1"){yintercept = unique(all_res$BGY1)}
  else if(metric_y == "hat_BGY2"){yintercept = unique(all_res$BGY2)}
  else if(metric_y == "hat_diff_BGY"){yintercept = unique(all_res$diff_BGY)}
  else if(metric_y == "hat_BXY1"){yintercept = unique(all_res$BXY1)}
  else if(metric_y == "hat_BXY2"){yintercept = unique(all_res$BXY2)}
  else if(metric_y == "hat_diff_BXY"){yintercept = unique(all_res$diff_BXY)}
  else{yintercept = NA}
  
  if(length(yintercept) == 1){
    if(!is.na(yintercept)){
      p <- p + geom_hline(yintercept = yintercept, linewidth = 0.5, alpha = 0.8, color = "gray40")
    }
  }
  
  if(plot_boxplots == T){
    p <- p + geom_boxplot(outliers = F, colour = "gray30", fill = NA, width = 0.3, 
                          box.linewidth = 0.2, whisker.linewidth = 0.2, staple.linewidth = 0.2, median.linewidth = 0.4, inherit.aes = T) 
  }
  
  return(p)
}
  

## Plot estimated BGX 
p1 <- plot_metric_across_subscenarios(all_res, all_metrics, "hat_BGX1", "TPR_BGX1", T)
p2 <- plot_metric_across_subscenarios(all_res, all_metrics, "hat_BGX2", "TPR_BGX2", T)
p3 <- plot_metric_across_subscenarios(all_res, all_metrics, "se_hat_BGX1", "TPR_BGX1", T)
p4 <- plot_metric_across_subscenarios(all_res, all_metrics, "se_hat_BGX2", "TPR_BGX2", T)
p5 <- plot_metric_across_subscenarios(all_res, all_metrics, "BGX_t_stat_1", "TPR_BGX1", T)
p6 <- plot_metric_across_subscenarios(all_res, all_metrics, "BGX_t_stat_2", "TPR_BGX2", T)
p7 <- plot_metric_across_subscenarios(all_res, all_metrics, "BGX_P_1", "TPR_BGX1", T)
p8 <- plot_metric_across_subscenarios(all_res, all_metrics, "BGX_P_2", "TPR_BGX2", T)

plot_grid(plotlist = list(p1, p2, p3, p4, p5, p6, p7, p8), ncol = 2, align = "vh")
ggsave(filename = paste0(plot_dir, "/scenario.00/BGXk_estimated_across_replicates_x_subsce.pdf"), width = 33, height = 12)

## Plot true - hat BGXk
p9 <- plot_metric_across_subscenarios(all_res, all_metrics, "diff_BGX1", "mean_diff_BGX1", T)
p10 <- plot_metric_across_subscenarios(all_res, all_metrics, "diff_BGX2", "mean_diff_BGX2", T)

plot_grid(plotlist = list(p9, p10), ncol = 1, align = "vh")
ggsave(filename = paste0(plot_dir, "/scenario.00/diff_BGXk_across_replicates_x_subsce.pdf"), width = 17, height = 6)

## Plot estimated between-strata BGX difference and true vs estimated difference
p11 <- plot_metric_across_subscenarios(all_res, all_metrics, "hat_diff_BGX", "mean_hat_diff_BGX", T)
p12 <- plot_metric_across_subscenarios(all_res, all_metrics, "diff_diff_BGX", "mean_diff_diff_BGX", T)
p13 <- plot_metric_across_subscenarios(all_res, all_metrics, "Z_diff_BGX", "PR_GxK_on_X", T)
p14 <- plot_metric_across_subscenarios(all_res, all_metrics, "p_Z_diff_BGX", "PR_GxK_on_X")

plot_grid(plotlist = list(p11, p12, p13, p14), ncol = 2, align = "vh")
ggsave(filename = paste0(plot_dir, "/scenario.00/between_strata_BGX_diff_across_replicates_x_subsce.pdf"), width = 40, height = 6)

## Plot estimated BGYk across replicates x subscenario 
p15 <- plot_metric_across_subscenarios(all_res, all_metrics, "hat_BGY1", "PR_BGY1", T)
p16 <- plot_metric_across_subscenarios(all_res, all_metrics, "se_hat_BGY1", "PR_BGY1", T)
p17 <- plot_metric_across_subscenarios(all_res, all_metrics, "BGY_t_stat_1", "PR_BGY1", T)
p18 <- plot_metric_across_subscenarios(all_res, all_metrics, "BGY_P_1", "PR_BGY1", T)

p19 <- plot_metric_across_subscenarios(all_res, all_metrics, "hat_BGY2", "PR_BGY2", T)
p20 <- plot_metric_across_subscenarios(all_res, all_metrics, "se_hat_BGY2", "PR_BGY2", T)
p21 <- plot_metric_across_subscenarios(all_res, all_metrics, "BGY_t_stat_2", "PR_BGY2", T)
p22 <- plot_metric_across_subscenarios(all_res, all_metrics, "BGY_P_2", "PR_BGY2", T)

plot_grid(plotlist = list(p15, p16, p17, p18, p19, p20, p21, p22), ncol = 2, align = "vh")
ggsave(filename = paste0(plot_dir, "/scenario.00/BGYk_estimated_across_replicates_x_subsce.pdf"), width = 40, height = 12)

## Plot true - hat BGYk
p23 <- plot_metric_across_subscenarios(all_res, all_metrics, "diff_BGY1", "mean_diff_BGY1", T)
p24 <- plot_metric_across_subscenarios(all_res, all_metrics, "diff_BGY2", "mean_diff_BGY2", T)

plot_grid(plotlist = list(p23, p24), ncol = 2, align = "vh")
ggsave(filename = paste0(plot_dir, "/scenario.00/diff_BGYk_across_replicates_x_subsce.pdf"), width = 40, height = 3)

## Plot estimated between-strata BGY difference and true vs estimated difference
p25 <- plot_metric_across_subscenarios(all_res, all_metrics, "hat_diff_BGY", "mean_hat_diff_BGY", T)
p26 <- plot_metric_across_subscenarios(all_res, all_metrics, "diff_diff_BGY", "mean_diff_diff_BGY", T)
p27 <- plot_metric_across_subscenarios(all_res, all_metrics, "Z_diff_BGY", "NR_GxK_on_Y", T)
p28 <- plot_metric_across_subscenarios(all_res, all_metrics, "p_Z_diff_BGY", "PR_GxK_on_Y")

plot_grid(plotlist = list(p25, p26, p27, p28), ncol = 2, align = "vh")
ggsave(filename = paste0(plot_dir, "/scenario.00/between_strata_BGY_diff_across_replicates_x_subsce.pdf"), width = 40, height = 6)

## Plot estimated BXYk 
p29 <- plot_metric_across_subscenarios(all_res, all_metrics, "hat_BXY1", "mean_hat_BXY1", T)
p30 <- plot_metric_across_subscenarios(all_res, all_metrics, "hat_BXY2", "mean_hat_BXY2", T)
## Plot true - estimated BXYk 
p31 <- plot_metric_across_subscenarios(all_res, all_metrics, "diff_BXY1", "mean_diff_BXY1", T)
p32 <- plot_metric_across_subscenarios(all_res, all_metrics, "diff_BXY2", "mean_diff_BXY2", T)
## Plot estimated BXYk se's
p33 <- plot_metric_across_subscenarios(all_res, all_metrics, "se_hat_BXY1", "mean_hat_BXY1", T)
p34 <- plot_metric_across_subscenarios(all_res, all_metrics, "se_hat_BXY2", "mean_hat_BXY2", T)
p35 <- plot_metric_across_subscenarios(all_res, all_metrics, "se_hat_BXY1_2nd", "mean_hat_BXY1", T)
p36 <- plot_metric_across_subscenarios(all_res, all_metrics, "se_hat_BXY2_2nd", "mean_hat_BXY2", T)

plot_grid(plotlist = list(p29, p30, p31, p32, p33, p34, p35, p36), ncol = 2, align = "vh")
ggsave(filename = paste0(plot_dir, "/scenario.00/BXYk_estimated_across_replicates_x_subsce.pdf"), width = 40, height = 12)


## Plot estimated between-strata BXY diff and true vs estimated diff
p37 <- plot_metric_across_subscenarios(all_res, all_metrics, "hat_diff_BXY", "mean_hat_diff_BXY", T)
p38 <- plot_metric_across_subscenarios(all_res, all_metrics, "diff_diff_BXY", "mean_diff_diff_BXY", T)
## Plot Z score and P val for diff
p39 <- plot_metric_across_subscenarios(all_res, all_metrics, "Z_diff_BXY", "PR_XxK_on_Y", T)
p40 <- plot_metric_across_subscenarios(all_res, all_metrics, "p_Z_diff_BXY", "PR_XxK_on_Y")
p41 <- plot_metric_across_subscenarios(all_res, all_metrics, "Z_diff_BXY_2nd", "PR_XxK_on_Y_2nd", T)
p42 <- plot_metric_across_subscenarios(all_res, all_metrics, "p_Z_diff_BXY_2nd", "PR_XxK_on_Y_2nd")

plot_grid(plotlist = list(p37, p38, p39, p40, p41, p42), ncol = 2, align = "vh")
ggsave(filename = paste0(plot_dir, "/scenario.00/between_strata_BXY_diff_across_replicates_x_subsce.pdf"), width = 40, height = 9)







## Reproducibility info
session_info()

# ─ Session info ───────────────────────────────────────────────────────────────────────────────────────
# setting  value
# version  R version 4.5.1 (2025-06-13)
# os       macOS Monterey 12.7.6
# system   x86_64, darwin20
# ui       RStudio
# language (EN)
# collate  en_US.UTF-8
# ctype    en_US.UTF-8
# tz       Europe/London
# date     2025-10-07
# rstudio  2025.05.1+513 Mariposa Orchid (desktop)
# pandoc   NA
# quarto   1.6.42 @ /private/var/folders/j1/klb6dgpn3dxfdrxfmbb4f0jc0000gp/T/AppTranslocation/7D11FDEE-A995-4F90-9AC5-8C647803C12A/d/RStudio.app/Contents/Resources/app/quarto/bin/quarto
# 
# ─ Packages ───────────────────────────────────────────────────────────────────────────────────────────
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
# stringi        1.8.7   2025-03-27 [1] CRAN (R 4.5.0)
# stringr        1.5.2   2025-09-08 [1] CRAN (R 4.5.1)
# tibble         3.3.0   2025-06-08 [1] CRAN (R 4.5.0)
# tidyr        * 1.3.1   2024-01-24 [1] CRAN (R 4.5.0)
# tidyselect     1.2.1   2024-03-11 [1] CRAN (R 4.5.0)
# usethis      * 3.2.1   2025-09-06 [1] CRAN (R 4.5.1)
# utf8           1.2.6   2025-06-08 [1] CRAN (R 4.5.0)
# vctrs          0.6.5   2023-12-01 [1] CRAN (R 4.5.0)
# withr          3.0.2   2024-10-28 [1] CRAN (R 4.5.0)
# 
# [1] /Library/Frameworks/R.framework/Versions/4.5-x86_64/Resources/library
# * ── Packages attached to the search path.
# 
# ──────────────────────────────────────────────────────────────────────────────────────────────────────





