
library(tidyr)
library(dplyr)
library(sessioninfo)

## Initialize
rm(list = ls())
set.seed(09222025)

################################################################################
#          1. Simulation study based on individual level data for 1 IV
################################################################################

# ------------------------------------------------------------------------------
#                     0.0 Definition of simulation scenarios
# ------------------------------------------------------------------------------
#  Define main scenarios to simulate according to Δβɢx and Δβxʏ, and their 
#  sub-scenarios according to varying parameter(s) and the values they take.
# ______________________________________________________________________________

## Define dirs
scripts_dir <- paste(getwd(), "simulations", "1.Individual_level_data_1IV", "00_define_scenarios", "scripts", sep = "/")
input_dir <- paste(getwd(), "simulations", "1.Individual_level_data_1IV", "00_define_scenarios", "inputs", sep = "/")
out_dir <- paste(getwd(), "simulations", "1.Individual_level_data_1IV", "00_define_scenarios", "outputs", sep = "/")
plot_dir <- paste(getwd(), "simulations", "1.Individual_level_data_1IV", "00_define_scenarios", "plots", sep = "/")

dir.create(scripts_dir, recursive = T, showWarnings = F)
dir.create(input_dir, recursive = T, showWarnings = F)
dir.create(out_dir, recursive = T, showWarnings = F)
dir.create(plot_dir, recursive = T, showWarnings = F)


## Save all scenarios
scenarios <- vector()

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
#  Base model: N = 10k, r = 1, q₁ = q₂ = 0.25, βɢx₁ = 0.5, βxʏ₁ = 0.7
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
N_base = 10000
r_base = 1
q1_base = q2_base = 0.25
BGX1_base = 0.5
BXY1_base = 0.7
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 


# ==============================================================================================
#  Scenario 00: No genetic effect difference Δβɢx = 0 and no causal effect difference Δβxʏ = 0 
# ==============================================================================================
main_scenario =  "scenario.00"
diff_BGX = 0
diff_BXY = 0

## Sub-scenarios for one varying parameter at a time:
# ______________________________________________________________________________________________
#  1. Varying N: N = {10k, 20k, ..., 100k} r = 1, q₁ = q₂ = 0.25, βɢx₁ = 0.5, βxʏ₁ = 0.7
# ______________________________________________________________________________________________
sub_sce_varying_par = "N"
Ns <- seq(from = 10, to = 100, by = 10)*1000 
# N = {10,000  20,000  30,000  40,000  50,000  60,000  70,000  80,000  90,000  100,000}

# ______________________________________________________________________________________________
#  2. Varying r: N = 10k, r = {0.20, ..., 1, ..., 5}, q₁ = q₂ = 0.25, βɢx₁ = 0.5, βxʏ₁ = 0.7
# ______________________________________________________________________________________________
sub_sce_varying_par = "r"
rs <- c(rev(1/(1:5)), 2:5)
# r = {0.20  0.25  0.33  0.50  1.0  2.0  3.0  4.0  5.0}

# ______________________________________________________________________________________________
#  3. Varying q₁ and q₂: N = 10k, r = 1, q₁ = q₂ = {0.05, ..., 0.45}, βɢx₁ = 0.5, βxʏ₁ = 0.7
# ______________________________________________________________________________________________
sub_sce_varying_par = "q1_and_q2"
q1s <- q2s <- seq(from = 0.05, to = 0.5, by = 0.1)
# q₁ = q₂ = {0.05  0.15  0.25  0.35  0.45}

# ______________________________________________________________________________________________
#  4. Varying βɢx₁: N = 10k, r = 1, q₁ = q2 = 0.25, βɢx₁ = {0.01, 0.09,..., 0.73}, βxʏ₁ = 0.7
# ______________________________________________________________________________________________
sub_sce_varying_par = "BGX1"
BGX1s <- seq(from = 0.01, to = 0.8, by = 0.08)
# βɢx₁ = {0.01  0.09  0.17  0.25  0.33  0.41  0.49  0.57  0.65  0.73}

# ______________________________________________________________________________________________
#  5. Varying βxʏ₁: N = 10k, r = 1, q₁ = q2 = 0.25, βɢx₁ = 0.5, βxʏ₁ = {0.1, ..., 1.45}
# ______________________________________________________________________________________________
sub_sce_varying_par = "BXY1"
BXY1s <- seq(from = 0.1, to = 1.5, by = 0.15)
# βxʏ₁ = {0.10  0.25  0.40  0.55  0.70  0.85  1.00  1.15  1.30  1.45}


#  Sub-scenarios for two varying parameters at a time: 
# ______________________________________________________________________________________________
#  6. Varying N and r: N = {10k, ..., 100k} r = {0.20, ..., 1, ..., 5}
#                      q₁ = q₂ = 0.25, βɢx₁ = 0.5, βxʏ₁ = 0.7
# ______________________________________________________________________________________________
sub_sce_varying_par = "N_and_r"
Ns <- seq(from = 10, to = 100, by = 10)*1000 
# N = {10,000  20,000  30,000  40,000  50,000  60,000  70,000  80,000  90,000  100,000}

# _







c = 1 
for (N in Ns){
  scenarios <- rbind(scenarios, 
                     c("scenario" = scenario, "sub_scenario1" = sub_scenario1, "sub_scenario2" = c, 
                       "N" = N, "r" = r_base, "q1" = q1_base, "q2" = q2_base, 
                       "BGX1" = BGX1_base, "diff_BGX" = diff_BGX, "BXY1" = BXY1_base, "diff_BXY" = diff_BXY)) %>% as.data.frame()
  c = c + 1
}

# -------------------------------------------------------------------------------------------------
#  2. Varying r: N = 10k, r = {0.20, ..., 1, ..., 5}, q₁ = q₂ = 0.25, βɢx₁ = 0.5, βxʏ₁ = 0.7
# -------------------------------------------------------------------------------------------------
rs <- c(rev(1/(1:5)), 2:5)
sub_scenario1 = "2"

c = 1 
for (r in rs){
  scenarios <- rbind(scenarios, c(scenario, sub_scenario1, c, 
                                  N_base, "r" = r, q1_base, q2_base, BGX1_base, diff_BGX, BXY1_base, diff_BXY, 
                                  sigma_sq_BGX_base, sigma_sq_BXY_base)) 
  c = c + 1
}

# -------------------------------------------------------------------------------------------------
#  3. Varying q₁ and q₂: N = 10k, r, q₁ = q₂ = {0.05, ..., 0.45}, βɢx₁ = 0.5, βxʏ₁ = 0.7, σβɢx² = σβxʏ² = 0.01
# -------------------------------------------------------------------------------------------------
q1s <- q2s <- seq(from = 0.05, to = 0.5, by = 0.1)
sub_scenario1 = "3"

c = 1 
for (q1 in q1s){
  for(q2 in q2s){
    scenarios <- rbind(scenarios, c(scenario, sub_scenario1, c, 
                                    N_base, r_base, q1, q2, BGX1_base, diff_BGX, BXY1_base, diff_BXY, 
                                    sigma_sq_BGX_base, sigma_sq_BXY_base)) 
    c = c + 1
  }
}
# -------------------------------------------------------------------------------------------------
#  4. Varying βɢx₁: N = 10k, r, q₁ = q2 = 0.25, βɢx₁ = {0.01, 0.09,..., 0.73}, βxʏ₁ = 0.7, σβɢx² = σβxʏ² = 0.01
# -------------------------------------------------------------------------------------------------
BGX1s <- seq(from = 0.01, to = 0.8, by = 0.08)
sub_scenario1 = "4"

c = 1 
for (BGX1 in BGX1s){
  scenarios <- rbind(scenarios, c(scenario, sub_scenario1, c, 
                                  N_base, r_base, q1_base, q2_base, BGX1, diff_BGX, BXY1_base, diff_BXY, 
                                  sigma_sq_BGX_base, sigma_sq_BXY_base)) 
  c = c + 1
  
}

# -------------------------------------------------------------------------------------------------
#  5. Varying βxʏ₁: N = 10k, r, q₁ = q2 = 0.25, βɢx₁ = 0.5, βxʏ₁ = {0.1, ..., 1.45}, σβɢx² = σβxʏ² = 0.01
# -------------------------------------------------------------------------------------------------
BXY1s <- seq(from = 0.1, to = 1.5, by = 0.15)
sub_scenario1 = "5"

c = 1 
for (BXY1 in BXY1s){
  scenarios <- rbind(scenarios, c(scenario, sub_scenario1, c, 
                                  N_base, r_base, q1_base, q2_base, BGX1_base, diff_BGX, BXY1, diff_BXY, 
                                  sigma_sq_BGX_base, sigma_sq_BXY_base)) 
  c = c + 1
}







## Run simulations: 100 replicates x scenario 
scenarios[, 4:11] <- apply(scenarios[, 4:11], 2, as.numeric)

scenario_args = scenarios[1,]


run_simulation_per_scenario <- function(scenario_args, n_replicates = 100){
  
  scenario_res = vector()
  
  for(i in 1:n_replicates){
    plotting_option = if_else(i <= 10, T, F)
    res = do.call(simulation_one_IV, c(scenario_args, replicate = i, plotting = plotting_option))
    
    scenario_res <- rbind(scenario_res, res)
  }
  
  save(scenario_res, file = paste(out_dir, sub_sub_scenario_args["scenario"], 
                                  sub_sub_scenario_args["sub_scenario1"], 
                                  sub_sub_scenario_args["sub_scenario2"], 
                                  "results_across_replicates.Rdata", sep = "/"))
  
  
  ## Compute sub,sub scenario summary metrics across replicates
  
  apply(scenario_res, c("n_aa_1", "n_aA_1", "n_AA_1", 
                        "n_aa_2", "n_aA_2", "n_AA_2", 
                        "q1_ob", "q2_ob", "HWE_CHISQ_1", "HWE_CHISQ_2"))
  
  ## FPR in HWE
  mean(scenario_res$HWE_P_1 < 0.05)
  mean(scenario_res$HWE_P_2 < 0.05)
  mean(scenario_res$HWE_P_global < 0.05)
  
  ## Mean BGXk across replicates
  mean(scenario_res$drawn_BGX1)
  mean(scenario_res$drawn_BGX2)
  
  ## Mean of estimated BGXk across replicates
  mean(scenario_res$hat_BGX1)
  mean(scenario_res$hat_BGX2)
  
  ## Mean F-stat of SNP 
  mean(scenario_res$BGX_F_stat_1)
  mean(scenario_res$BGX_F_stat_2)
  
  
  
  ## FNR for estimated effects in each stratum
  mean(scenario_res$BGX_P_1 >= 0.05)
  mean(scenario_res$BGX_P_2 >= 0.05)
  
  mean(scenario_A1.0_res$diff_BGX1)
  mean(scenario_A1.0_res$diff_BGX2)
  
  ## FPR for estimated effect differences 
  mean(scenario_A1.0_res$p_Z_diff_BGX < 0.05)
  
  
  
  ## Add scenario pars
  scenarios <- rbind(scenarios, 
                     c("scenario" = "scenario_A1.0",
                       "N" = N, "r" = r, "q1" = q1, "q2" = q2, "BGX1" = BGX1, 
                       "diff_BGX" = diff_BGX, "BXY1" = NA, "diff_BXY" = NA))
  
  
  
  
  plots_x_scenario <- function(sim_data, metrics, )
    
    labs = c("diff_BGX1" = "|βGX1 - β̂GX1|",
             "diff_BGX2" = "|βGX2 - β̂GX2|",
             "diff_diff_BGX" = "|ΔβGX - Δβ̂GX|",
             "Z_diff_BGX" = "Z-stat for GxK", 
             "p_Z_diff_BGX" = "Pval for GxK")
  
  ggplot(data = scenario_A1.0, 
         mapping = aes(x = NA, 
                       y = diff_BGX1, 
                       color = NA)) +
    geom_violin(alpha = 0, size = 0.4, color='black', width = 0.7)+
    geom_jitter(width = 0.1, alpha = 0.7, size = 2) +
    geom_boxplot(alpha = 0, size = 0.4, width=0.1, color='black') +
    theme_classic() +
    # scale_color_manual(values = colors) +
    labs(y = labs["diff_BGX1"]) +
    theme(axis.title = element_text(size = (12)),
          axis.text = element_text(size = (11)))
  
}






# ==============================================================================================
#  Scenario 01: No genetic effect difference Δβɢx = 0 and causal effect difference Δβxʏ ≠ 0 
# ==============================================================================================


# ==============================================================================================
#  Scenario 10: No genetic effect difference Δβɢx ≠ 0 and causal effect difference Δβxʏ = 0 
# ==============================================================================================


# ==============================================================================================
#  Scenario 11: No genetic effect difference Δβɢx ≠ 0 and causal effect difference Δβxʏ ≠ 0 
# ==============================================================================================







