
library(tidyr)
library(dplyr)
library(sessioninfo)

rm(list = ls())

# ------------------------------------------------------------------------------
#                  2. Summarize and plot results per scenario 
# ------------------------------------------------------------------------------
# ______________________________________________________________________________
#    1.0  Compute and aggregate metrics across replicates per subscenario 
# ______________________________________________________________________________
#  This script summarizes results across replicates per subscenario and
#  aggregates results across subscenarios of a varying parameter scenario. 
# ______________________________________________________________________________

input_dir <- paste(getwd(), "simulations", "1.Individual_level_data_1IV", "02_summarize_and_plot_results", "inputs", sep = "/")
out_dir <- paste(getwd(), "simulations", "1.Individual_level_data_1IV", "02_summarize_and_plot_results", "outputs", sep = "/")
plot_dir <- paste(getwd(), "simulations", "1.Individual_level_data_1IV", "02_summarize_and_plot_results", "plots", sep = "/")

dir.create(input_dir, showWarnings = F)
dir.create(out_dir, showWarnings = F)
dir.create(plot_dir, showWarnings = F)

## Input dirs
input_dir00 <- paste(getwd(), "simulations", "1.Individual_level_data_1IV", "00_define_scenarios", "outputs", sep = "/")
input_dir01 <- paste(getwd(), "simulations", "1.Individual_level_data_1IV", "01_run_replicates_x_scenario", "outputs", sep = "/")


## Scenarios
# all_scenarios <- get(load(paste0(input_dir00, "/scenario.00.Rdata")))
all_scenarios <- get(load(paste0(input_dir00, "/scenario.01.Rdata")))[1:1128,]
scenarios <- all_scenarios[, c("main_scenario", "main_scenario_val", "sub_sce_varying_par")] %>% unique()
sub_scenarios <- all_scenarios[, c("main_scenario", "main_scenario_val", "sub_sce_varying_par","sub_sce_varying_par_value")]


## Confirm res file x subscenario first
for(j in 1:nrow(sub_scenarios)){

  sub_sce = sub_scenarios[j,]
  scenario <- paste(sub_sce, collapse = "/")
  res <- get(load(paste(c(input_dir01, scenario, "scenario_res_across_100replicates.Rdata"), collapse = "/")))
  if(ncol(res) != 77){
    print(ncol(res))
  }
  if(!"scenario_res_across_100replicates.Rdata" %in% list.files(paste(c(input_dir01, scenario), collapse = "/"))){
    print(j)
  }
}

## For agg results across subscenarios
for(j in 1:nrow(scenarios)){
  assign(paste(c("all_res", scenarios[j,]), collapse = "_"), vector())
  assign(paste(c("all_metrics", scenarios[j,]), collapse = "_"), vector())
}


for(j in 1:nrow(sub_scenarios)){
  
  sub_sce = sub_scenarios[j,]
  scenario <- paste(sub_sce, collapse = "/")
  
  ## Verify results were generated for 100 replicates 
  if((list.files(paste(c(input_dir01, scenario), collapse = "/")) %>% grep("indiv_data", .) %>% length()) != 100){
    print(paste0("Not all replicates run for ", paste(sub_sce, collapse = ", "), ". Stopping."))
    stop()
  }
  if(!"scenario_res_across_100replicates.Rdata" %in% list.files(paste(c(input_dir01, scenario), collapse = "/"))){
    print(paste0("Results across 100 replicates not present for ", paste(sub_sce, collapse = ", "), ". Stopping."))
    stop()
  }
  
  ## Extract results across replicates 
  res <- get(load(paste(c(input_dir01, scenario, "scenario_res_across_100replicates.Rdata"), collapse = "/")))

  if(ncol(res) == 75){
   res = cbind(res[,1:11], "BGX2" = res$BGX1 + res$diff_BGX, res[,12:13], "BXY2" = res$BXY1 + res$diff_BXY, res[,14:75])
  }
  
  ## Add se for estimated causal effects including 2nd term from Delta expansion
  res$se_hat_BXY1_2nd <- sqrt((res$se_hat_BGY1^2 / res$hat_BGX1^2) + ((res$hat_BGY1^2)*(res$se_hat_BGX1^2) / (res$hat_BGX1^4)))
  res$se_hat_BXY2_2nd <- sqrt((res$se_hat_BGY2^2 / res$hat_BGX2^2) + ((res$hat_BGY2^2)*(res$se_hat_BGX2^2) / (res$hat_BGX2^4)))
  
  ## Z score and Pval for causal effect diff based on those se's
  res$Z_diff_BXY_2nd <- (res$hat_BXY2 - res$hat_BXY1) / sqrt((res$se_hat_BXY2_2nd^2) + (res$se_hat_BXY1_2nd^2))
  res$p_Z_diff_BXY_2nd = 2*pnorm(abs(res$Z_diff_BXY_2nd), 0, 1, lower.tail = F)
  
  ## Append with other sub scenarios results
  all_res <- get(paste(c("all_res", sub_sce[-4]), collapse = "_"))
  all_res <- rbind(all_res, res)
  assign(paste(c("all_res", sub_sce[-4]), collapse = "_"), all_res)
    
  ## Compute summary metrics x subscenario
  scenario_metrics <- with(res, 
                           c("mean_q1_ob" = mean(q1_ob), "mean_q2_ob" = mean(q2_ob), "mean_q_global_ob" = mean(q_global_ob), 
                             "FPR_HWE_1" = mean(HWE_P_1 < 0.05), "FPR_HWE_2" = mean(HWE_P_2 < 0.05), "FPR_HWE_global" = mean(HWE_P_global < 0.05),
                             
                             "mean_hat_BGX1" = mean(hat_BGX1), "median_hat_BGX1" = median(hat_BGX1), "var_hat_BGX1" = var(hat_BGX1), 
                             "mean_hat_BGX2" = mean(hat_BGX2), "median_hat_BGX2" = median(hat_BGX2), "var_hat_BGX2" = var(hat_BGX2), 
                             
                             "FNR_BGX1" = mean(BGX_P_1 >= 0.05),"FNR_BGX2" = mean(BGX_P_2 >= 0.05),
                             "TPR_BGX1" = mean(BGX_P_1 < 0.05), "TPR_BGX2" = mean(BGX_P_2 < 0.05),
                             
                             "mean_IV_F_stat_1" = mean(BGX_F_stat_1), "mean_IV_F_stat_2" = mean(BGX_F_stat_2),
                             
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
                             "PR_GxK_on_X" = mean(res$p_Z_diff_BGX < 0.05),
                             "NR_GxK_on_X" = mean(res$p_Z_diff_BGX >= 0.05),
                             
                             ## +/- rate for G effect on Y in k=1,2
                             "PR_BGY1" = mean(BGY_P_1 < 0.05),
                             "NR_BGY1" = mean(BGY_P_1 >= 0.05),
                             "PR_BGY2" = mean(BGY_P_2 < 0.05),
                             "NR_BGY2" = mean(BGY_P_2 >= 0.05),
                             
                             ## +/- rate for GxK interaction on Y
                             "PR_GxK_on_Y" = mean(res$p_Z_diff_BGY < 0.05),
                             "NR_GxK_on_Y" = mean(res$p_Z_diff_BGY >= 0.05),
                             
                             ## +/- rate for XxK interaction on Y (based on 1st term se's)
                             "PR_XxK_on_Y" = mean(res$p_Z_diff_BXY < 0.05),
                             "NR_XxK_on_Y" = mean(res$p_Z_diff_BXY >= 0.05),
                             ## +/- rate for XxK interaction on Y (based on 2nd term se's)
                             "PR_XxK_on_Y_2nd" = mean(res$p_Z_diff_BXY_2nd < 0.05),
                             "NR_XxK_on_Y_2nd" = mean(res$p_Z_diff_BXY_2nd >= 0.05)
                           ))
  
  save(scenario_metrics, file = paste(c(input_dir01, scenario, "scenario_summary_metrics.Rdata"), collapse = "/")) 
  
  all_metrics <- get(paste(c("all_metrics", sub_sce[-4]), collapse = "_"))
  to_bind <- data.frame("main_scenario" = sub_sce$main_scenario, 
                "main_scenario_val" = sub_sce$main_scenario_val, 
                "sub_sce_varying_par" = sub_sce$sub_sce_varying_par, 
                "sub_sce_varying_par_value" = sub_sce$sub_sce_varying_par_value, 
                t(scenario_metrics))
  all_metrics <- rbind(all_metrics, to_bind)
  assign(paste(c("all_metrics", sub_sce[-4]), collapse = "_"), all_metrics)
}


## Confirm 100 rep x num sub scenarios for each var par scenario agg
for(j in 1:nrow(scenarios)){
  num_sub <- all_scenarios %>% 
    dplyr::filter(main_scenario == scenarios[j,"main_scenario"], 
                  main_scenario_val == scenarios[j,"main_scenario_val"], 
                  sub_sce_varying_par == scenarios[j, "sub_sce_varying_par"]) %>% nrow()
  
  if(nrow(get(paste(c("all_res", scenarios[j,]), collapse = "_"))) != num_sub*100){print(j)}
  if(nrow(get(paste(c("all_metrics", scenarios[j,]), collapse = "_"))) != num_sub){print(j)}
}

## Save agg
for(j in 1:nrow(scenarios)){
  
  sce = scenarios[j,]
  path <- paste(c(out_dir, sce), collapse = "/")
  dir.create(path, recursive = T, showWarnings = F)
  
  save(list = paste(c("all_res", sce), collapse = "_"), file = paste0(path, "/agg_results_subscenarios_replicates.Rdata"))
  save(list = paste(c("all_metrics", sce), collapse = "_"), file = paste0(path, "/agg_metrics_subscenarios.Rdata"))
}







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
