
library(tidyr)
library(dplyr)
library(sessioninfo)

## Initialize
rm(list = ls())
set.seed(09222025)

# ------------------------------------------------------------------------------
#                  2.0 Summarize and plot results per scenario 
# ------------------------------------------------------------------------------
#  This script summarizes per-replicate results for each subscenario for a 
#  varying parameter(s) scenario and plots ... 
# ______________________________________________________________________________

## Define dirs
input_dir <- paste(getwd(), "simulations", "1.Individual_level_data_1IV", "02_summarize_and_plot_results", "inputs", sep = "/")
out_dir <- paste(getwd(), "simulations", "1.Individual_level_data_1IV", "02_summarize_and_plot_results", "outputs", sep = "/")
plot_dir <- paste(getwd(), "simulations", "1.Individual_level_data_1IV", "02_summarize_and_plot_results", "plots", sep = "/")

dir.create(input_dir, showWarnings = F)
dir.create(out_dir, showWarnings = F)
dir.create(plot_dir, showWarnings = F)



scenario_args <- c(main = "scenario.00", sub_sce_varying_par = "N")

## Scenario input dir
input_dir_sce <- paste(getwd(), "simulations", "1.Individual_level_data_1IV", "01_run_replicates_x_scenario", "outputs", 
                     scenario_args["main"], 
                     scenario_args["sub_sce_varying_par"], sep = "/")

## Scenario output dir
output_dir_sce <- paste(out_dir, scenario_args["main"], scenario_args["sub_sce_varying_par"], sep = "/")
dir.create(output_dir_sce, recursive = T, showWarnings = F)


## Iterate over all subscenarios 
subscenarios <- scenarios00 %>% 
  dplyr::filter(scenario == scenario_args["main"], sub_sce_varying_par == scenario_args["sub_sce_varying_par"]) %>% 
  pull(sub_sce_varying_par_value)

## Compute summary metrics x subscenario
for(s in subscenarios){
  
  ## Results across 100 replicates
  res = get(load(paste0(input_dir_sce, "/", s, "/scenario_res_across_100replicates.Rdata")))
  
  scenario_metrics <- with(res, 
                           c("mean_q1_ob" = mean(q1_ob), "mean_q2_ob" = mean(q2_ob), "mean_q_global_ob" = mean(q_global_ob), 
                             "FPR_HWE_1" = mean(HWE_P_1 < 0.05), "FPR_HWE_2" = mean(HWE_P_2 < 0.05), "FPR_HWE_global" = mean(HWE_P_global < 0.05),
                             
                             "mean_hat_BGX1" = mean(hat_BGX1), "median_hat_BGX1" = median(hat_BGX1), "var_hat_BGX1" = var(hat_BGX1), 
                             "mean_hat_BGX2" = mean(hat_BGX2), "median_hat_BGX2" = median(hat_BGX2), "var_hat_BGX2" = var(hat_BGX2), 
                             
                             "FNR_BGX1" = mean(BGX_P_1 >= 0.05),"FNR_BGX2" = mean(BGX_P_2 >= 0.05),
                             "TPR_BGX1" = mean(BGX_P_1 < 0.05), "TPR_BGX2" = mean(BGX_P_2 < 0.05),
                             
                             "mean_IV_F_stat_1" = mean(BGX_F_stat_1),"mean_IV_F_stat_2" = mean(BGX_F_stat_2),
                             
                             "mean_diff_BGX1" = mean(diff_BGX1), "mean_diff_BGX2" = mean(diff_BGX2),
                             "mean_hat_diff_BGX" = mean(hat_diff_BGX), 
                             "mean_diff_diff_BGX" = mean(diff_diff_BGX), 
                             
                             "mean_hat_BGY1" = mean(hat_BGY1), "median_hat_BGY1" = median(hat_BGY1), "var_hat_BGY1" = var(hat_BGY1), 
                             "mean_hat_BGY2" = mean(hat_BGY2), "median_hat_BGY2" = median(hat_BGY2), "var_hat_BGY2" = var(hat_BGY2), 
                             
                             "mean_diff_BGY1" = mean(diff_BGY1),"mean_diff_BGY2" = mean(diff_BGY2),
                             "mean_hat_diff_BGY" = mean(hat_diff_BGY), 
                             "mean_diff_diff_BGY" = mean(diff_diff_BGY), 
                             
                             "mean_hat_BXY1" = mean(hat_BXY1), "median_hat_BXY1" = median(hat_BXY1), "var_hat_BXY1" = var(hat_BXY1), 
                             "mean_hat_BXY2" = mean(hat_BXY2), "median_hat_BXY2" = median(hat_BXY2), "var_hat_BXY2" = var(hat_BXY2), 
                             
                             "mean_diff_BXY1" = mean(diff_BXY1), "mean_diff_BXY2" = mean(diff_BXY2),
                             "mean_hat_diff_BXY" = mean(hat_diff_BXY), 
                             "mean_diff_diff_BXY" = mean(diff_diff_BXY)
                             
                             ## Compute se and p(hat BXY1), se and p(hat BXY2), Zdiff and P for hat diff BXY
                           ))
  
  ## False/true positive/negative rate for G x K interaction on X
  col_name_GxK_on_X_pos <- if_else(unique(res$diff_BGX) == 0, "FPR_GxK_on_X", "TPR_GxK_on_X")
  col_name_GxK_on_X_neg <- if_else(unique(res$diff_BGX) == 0, "TNR_GxK_on_X", "FNR_GxK_on_X")
  
  scenario_metrics <- scenario_metrics %>% t() %>% as.data.frame() %>% 
    mutate(!!col_name_GxK_on_X_pos := mean(res$p_Z_diff_BGX < 0.05), 
           !!col_name_GxK_on_X_neg := mean(res$p_Z_diff_BGX >= 0.05))
  
  
  ## False/true positive/negative rate for G effect on Y in stratum 1
  col_name_BGY1_pos <- if_else(unique(res$BGY1) == 0, "FPR_BGY1", "TPR_BGY1")
  col_name_BGY1_neg <- if_else(unique(res$BGY1) == 0, "TNR_BGY1", "FNR_BGY1")
  
  scenario_metrics <- scenario_metrics  %>% 
    mutate(!!col_name_BGY1_pos := mean(res$BGY_P_1 < 0.05), 
           !!col_name_BGY1_neg := mean(res$BGY_P_1 >= 0.05))
  
  ## False/true positive/negative rate for G effect on Y in stratum 2
  col_name_BGY2_pos <- if_else(unique(res$BGY2) == 0, "FPR_BGY2", "TPR_BGY2")
  col_name_BGY2_neg <- if_else(unique(res$BGY2) == 0, "TNR_BGY2", "FNR_BGY2")
  
  scenario_metrics <- scenario_metrics  %>% 
    mutate(!!col_name_BGY2_pos := mean(res$BGY_P_2 < 0.05), 
           !!col_name_BGY2_neg := mean(res$BGY_P_2 >= 0.05))
  
  ## False/true positive/negative rate for G x K interaction on Y
  col_name_GxK_on_Y_pos <- if_else(unique(res$diff_BGY) == 0, "FPR_GxK_on_Y", "TPR_GxK_on_Y")
  col_name_GxK_on_Y_neg <- if_else(unique(res$diff_BGY) == 0, "TNR_GxK_on_Y", "FNR_GxK_on_Y")
  
  scenario_metrics <- scenario_metrics %>% 
    mutate(!!col_name_GxK_on_Y_pos := mean(res$p_Z_diff_BGY < 0.05), 
           !!col_name_GxK_on_Y_neg := mean(res$p_Z_diff_BGY >= 0.05))
  
  assign(paste0("scenario_metrics_", s), scenario_metrics)
  save(list = paste0("scenario_metrics_", s), file = paste0(input_dir_sce, "/", s, "/scenario_summary_metrics.Rdata"))
}

# # TODO 
# if(BXY1 == 0){
#   scenario_metrics <- append(scenario_metrics, "FPR_BXY1" = mean(BXY_P_1 < 0.05),
#                              "TNR_BXY1" = mean(BXY_P_1 >= 0.05))
# } else{
#   scenario_metrics <- append(scenario_metrics, "TPR_BXY1" = mean(BXY_P_1 < 0.05),
#                              "FNR_BXY1" = mean(BXY_P_1 >= 0.05))
# }
# 
# if(BXY2 == 0){
#   scenario_metrics <- append(scenario_metrics, "FPR_BXY2" = mean(BXY_P_2 < 0.05),
#                              "TNR_BXY2" = mean(BXY_P_2 >= 0.05))
# } else{
#   scenario_metrics <- append(scenario_metrics, "TPR_BXY2" = mean(BXY_P_2 < 0.05),
#                              "FNR_BXY2" = mean(BXY_P_2 >= 0.05))
# }
# 
# 
# if(diff_BXY == 0){
#   scenario_metrics <- append(scenario_metrics, "FPR_XxK_on_Y" = mean(p_Z_diff_BXY < 0.05), 
#                              "TNR_XxK_on_Y" = mean(p_Z_diff_BXY >= 0.05))
# } else{
#   scenario_metrics <- append(scenario_metrics, "TPR_XxK_on_Y" = mean(p_Z_diff_BXY < 0.05), 
#                              "FNR_XxK_on_Y" = mean(p_Z_diff_BXY >= 0.05))
# }



## Plot results x subscenario across replicates
all_res <- vector()
all_metrics <- vector()

for(s in subscenarios){
  
  res = get(load(paste0(input_dir_sce, "/", s, "/scenario_res_across_100replicates.Rdata")))
  all_res <- rbind(all_res, res)
  
  scenario_metrics = cbind(scenario = scenario_args["main"], 
                           sub_sce_varying_par = scenario_args["sub_sce_varying_par"],
                           subscenario = s, 
                           get(paste0("scenario_metrics_", s))) 
  
  all_metrics <- rbind(all_metrics, scenario_metrics)
  
}

save(all_res, file = paste0(output_dir_sce, "/all_subscenarios_res_across_100replicates.Rdata"))

rownames(all_metrics) <- NULL
save(all_metrics, file = paste0(output_dir_sce, "/scenario_metrics.Rdata"))


## Compare metrics in all replicates across all subscenarios
plot_metric_across_subscenarios <- function(all_res, all_metrics, metric_y, metrics_to_show){
  
  all_res$x <- stringr::str_split_i(all_res$sub_sce_varying_par_value, "=", 2) 
  all_res$x <- factor(all_res$x, levels = unique(all_res$x))
  
  all_metrics$x <- stringr::str_split_i(all_metrics$subscenario, "=", 2) 
  
  x_labs = list("N" = "Increasing N", 
                "r" = "Ratio of stratum sample sizes N1/N2")
  
  metrics_labels = list("n_aa_1" = "Number of homozygotes for minor allele in stratum 1",
                     "n_aA_1" = "Number of heterozygotes in stratum 1",
                     "n_AA_1" = "Number of homozygotes for major allele in stratum 1",
                     "n_aa_2" = "Number of homozygotes for minor allele in stratum 2",
                     "n_aA_2" = "Number of heterozygotes in stratum 2",
                     "n_AA_2" = "Number of homozygotes for major allele in stratum 2",
                     "HWE_P_1" = "P-value for HWE in stratum 1",
                     "HWE_P_2" = "P-value for HWE in stratum 2", 
                     "HWE_P_global" = "P-value for HWE in both strata", 
                     "HWE_CHISQ_1" = "Chi-square statistic for HWE in stratum 1",
                     "HWE_CHISQ_2" = "Chi-square statistic for HWE in stratum 2",
                     "HWE_CHISQ_global" = "Chi-square statistic for HWE in both strata",
                     "q1_ob" = "Observed MAF in stratum 1",
                     "q2_ob" = "Observed MAF in stratum 2",
                     "q_global_ob" = "Observed MAF in both strata",
                     "hat_BGX1" = expression(hat(beta[GX[1]])), 
                     "hat_BGX2" = expression(hat(beta[GX[2]])), 
                     "se_hat_BGX1" = expression(se(hat(beta[GX[1]]))), 
                     "se_hat_BGX2" = expression(se(hat(beta[GX[2]]))), 
                     "BGX_t_stat_1" = expression(t[hat(beta[GX[1]])]), 
                     "BGX_t_stat_2" = expression(t[hat(beta[GX[2]])]), 
                     "BGX_P_1" = expression(p[hat(beta[GX[1]])]), 
                     "BGX_P_2" = expression(p[hat(beta[GX[2]])]), 
                     "BGX_F_stat_1" = expression(F[hat(beta[GX[1]])]), 
                     "BGX_F_stat_2" = expression(F[hat(beta[GX[2]])]), 
                     "diff_BGX1" = expression(beta[GX[1]] - hat(beta[GX[1]])),
                     "diff_BGX2" = expression(beta[GX[2]] - hat(beta[GX[2]])),
                     "hat_diff_BGX" = expression(hat[delta*beta[GX]]), 
                     "diff_diff_BGX" = expression(delta*beta[GX] - hat[delta*beta[GX]]),
                     "Z_diff_BGX" = expression(Z[delta*beta[GX]]),
                     "p_Z_diff_BGX" = expression(p[Z[delta*beta[GX]]]),
                     "hat_BGY1" = expression(hat(beta[GY[1]])),                
                     "hat_BGY2" = expression(hat(beta[GY[2]])),                      
                     "se_hat_BGY1" = expression(se(hat(beta[GY[1]]))),             
                     "se_hat_BGY2" = expression(se(hat(beta[GY[2]]))),              
                     "BGY_t_stat_1" = expression(t[hat(beta[GY[1]])]),              
                     "BGY_t_stat_2" = expression(t[hat(beta[GY[2]])]),             
                     "BGY_P_1" = expression(p[hat(beta[GY[1]])]),                   
                     "BGY_P_2" = expression(p[hat(beta[GY[2]])]),         
                     "diff_BGY1" = expression(beta[GY[1]] - hat(beta[GY[1]])),                         
                     "diff_BGY2" = expression(beta[GY[2]] - hat(beta[GY[2]])), 
                     "diff_BGY" = expression(beta[GY[2]] - beta[GY[1]]), 
                     "hat_diff_BGY" = expression(hat[delta*beta[GY]]),              
                     "diff_diff_BGY" = expression(delta*beta[GY] - hat[delta*beta[GY]]), 
                     "Z_diff_BGY" = expression(Z[delta*beta[GY]]),              
                     "p_Z_diff_BGY" = expression(p[Z[delta*beta[GY]]]),            
                     "hat_BXY1" = expression(hat(beta[XY[1]])),             
                     "hat_BXY2" = expression(hat(beta[XY[2]])),                   
                     "diff_BXY1" = expression(beta[XY[1]] - hat(beta[XY[1]])),                     
                     "diff_BXY2" = expression(beta[XY[2]] - hat(beta[XY[2]])),      
                     "hat_diff_BXY" = expression(hat[delta*beta[XY]]),                
                     "diff_diff_BXY" = expression(delta*beta[XY] - hat[delta*beta[XY]]) 
                     )
    
  p = ggplot(all_res, aes(x = x, y = get(metric_y))) + 
    geom_point(color='gray', alpha = 0.3, size = 1) +
    theme_classic() +
    labs(x = x_lab, y = metrics_labels[[metric_y]]) +
    theme(axis.title = element_text(size = 8),
          axis.text = element_text(size = 7),
          axis.line = element_line(linewidth = 0.4),
          axis.ticks = element_line(linewidth = 0.3))
  
  if(metric_y %in% c("HWE_P_1", "HWE_P_2", "HWE_P_global", 
                     "BGX_P_1", "BGX_P_2", "p_Z_diff_BGX", 
                     "BGY_P_1", "BGY_P_2", "p_Z_diff_BGY")){
    
    p <- p + geom_hline(yintercept = 0.05, linewidth = 0.65, alpha = 0.8, color = "orangered1")
    p <- p + geom_text(data = all_metrics, aes(x = ))
  }
}






