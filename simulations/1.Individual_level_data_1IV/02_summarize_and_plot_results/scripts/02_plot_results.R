
library(tidyr, quietly = T)
library(dplyr, quietly = T)
library(sessioninfo, quietly = T)

rm(list = ls())

# ______________________________________________________________________________
#    2.0  Plot results across subscenarios of a varying parameter(s) scenario  
# ______________________________________________________________________________
#  This script summarizes results across replicates per subscenario and
#  aggregates results across subscenarios of a varying parameter scenario. 
# ______________________________________________________________________________

input_dir <- paste(getwd(), "simulations", "1.Individual_level_data_1IV", "02_summarize_and_plot_results", "inputs", sep = "/")
out_dir <- paste(getwd(), "simulations", "1.Individual_level_data_1IV", "02_summarize_and_plot_results", "outputs", sep = "/")
plot_dir <- paste(getwd(), "simulations", "1.Individual_level_data_1IV", "02_summarize_and_plot_results", "plots", sep = "/")

## Input dirs
input_dir00 <- paste(getwd(), "simulations", "1.Individual_level_data_1IV", "00_define_scenarios", "outputs", sep = "/")
input_dir01 <- paste(getwd(), "simulations", "1.Individual_level_data_1IV", "01_run_replicates_x_scenario", "outputs", sep = "/")

## Scenarios
all_scenarios <- get(load(paste0(input_dir00, "/scenario.00.Rdata")))
# all_scenarios <- get(load(paste0(input_dir00, "/scenario.01.Rdata")))[1:500,]
scenarios <- all_scenarios[, c("main_scenario", "main_scenario_val", "sub_sce_varying_par")] %>% unique()


## Scenarios
# all_scenarios <- get(load(paste0(input_dir00, "/scenario.00.Rdata")))
# all_scenarios <- get(load(paste0(input_dir00, "/scenario.01.Rdata")))[1:500,]
# scenarios <- all_scenarios[, c("main_scenario", "main_scenario_val", "sub_sce_varying_par")] %>% unique()
# 

## Compare metrics in all replicates across all subscenarios
plot_metric_across_subscenarios <- function(scenario, metric_y, metric_to_show, plot_boxplots = F){
  
  all_res <- get(load(paste(out_dir, paste(scenario, collapse = "/"), "agg_results_subscenarios_replicates.Rdata", sep = "/")))
  all_metrics <- get(load(paste(out_dir, paste(scenario, collapse = "/"), "agg_metrics_subscenarios.Rdata", sep = "/")))
    
  if(metric_y %in% c("BGX_P_1", "BGX_P_2", "BGY_P_1", "BGY_P_2")){
    all_res[, metric_y] <- format( all_res[, metric_y], scientific = T) %>% as.numeric()
  }
  
  if(scenario["sub_sce_varying_par"] %in% c("N", "r", "BGX1", "BXY1")){
    
    all_res$x <- stringr::str_split_i(all_res$sub_sce_varying_par_value, "=", 2) 
    all_res$x <- factor(all_res$x, levels = unique(all_res$x))
    all_metrics$x <- stringr::str_split_i(all_metrics$sub_sce_varying_par_value, "=", 2) 
    
  }
  
  if(scenario["sub_sce_varying_par"] == "q1_q2"){
    all_res$x <- gsub("^ ", "", gsub("q.=", " ", all_res$sub_sce_varying_par_value))
    all_res$x <- factor(all_res$x, levels = unique(all_res$x))
    
    all_metrics$x <- gsub("^ ", "", gsub("q.=", " ",  all_metrics$sub_sce_varying_par_value))
    all_metrics$x <- factor(all_metrics$x, levels = unique(all_metrics$x))
  }
  if(scenario["sub_sce_varying_par"] == "BUX_BUY"){
    all_res$x <- gsub("BUY=", " ", gsub("BUX=", " ", all_res$sub_sce_varying_par_value))
    all_res$x <- factor(all_res$x, levels = unique(all_res$x))
    
    all_metrics$x <- gsub("BUY=", " ", gsub("BUX=", " ", all_metrics$sub_sce_varying_par_value))
    all_metrics$x <- factor(all_metrics$x, levels = unique(all_metrics$x))
  }
  if(grepl("_", scenario["sub_sce_varying_par"]) & scenario["sub_sce_varying_par"] %in% c("q1_q2", "BUX_BUY")){
    
    if(grepl("q1_q2", scenario["sub_sce_varying_par"]) ){
      
    } else{
      cols_to_add_names <- stringr::str_split(scenario["sub_sce_varying_par"], "_") %>% unlist
      cols_to_add <- do.call(rbind, stringr::str_split(all_res[,"sub_sce_varying_par_value"], ",")) %>% as.data.frame()
      colnames(cols_to_add) <- cols_to_add_names
      cols_to_add[,1] <- gsub(paste0(cols_to_add_names[1], "="), "", cols_to_add[,1]) %>% as.numeric()
      cols_to_add[,2] <- gsub(paste0(cols_to_add_names[2], "="), "", cols_to_add[,2]) %>%  as.numeric()
      
      all_res$var1 <- cols_to_add[,1]
      all_res$var2 <- cols_to_add[,2]
    }
    
  }
  
  var_labs = list("N" = "Increasing N", 
                "r" = "Ratio of stratum sample sizes N1/N2", 
                "BGX1" = expression(beta[GX[1]]),
                "BXY1" = expression(beta[XY[1]]), 
                "BUX_BUY" = expression(beta[UX]*" = "*beta[UY]), 
                "q1_q2" = expression(q[1]*", "*q[2]),
                "diff_BGX" = expression(Delta*beta[GX]),
                "diff_BXY" = expression(Delta*beta[XY]))
  
  metric_to_show_lab = c("mean_q1_ob" = expression(bar(q)[ob1]*":"), 
                         "mean_q2_ob" = expression(bar(q)[ob2]*":"),
                         "mean_q_global_ob" = expression(bar(q)[ob]*":"),
                         "FPR_HWE_1" = "FPR:", "FPR_HWE_2" = "FPR:", "FPR_HWE_global" = "FPR:", 
                         "TPR_BGX1" = "TPR:", "TPR_BGX2" = "TPR:", 
                         "mean_IV_F_stat_1" = expression(bar(F)[hat(beta)[GX[1]]]*":"), "mean_IV_F_stat_2" = expression(bar(F)[hat(beta[GX[2]])]*":"), 
                         "mean_diff_BGX1" = expression(bar(Delta)*beta[GX[1]]*":"), 
                         "mean_diff_BGX2" = expression(bar(Delta)*beta[GX[2]]*":"), 
                         "mean_hat_diff_BGX" = expression(bar(Delta)*hat(beta)[GX]*":"), 
                         "mean_diff_diff_BGX" = expression(bar(Delta)[Delta*beta[GX]]*":"), 
                         "NR_GxK_on_X" = "NR for GxK on X:",
                         "PR_GxK_on_X" = "PR for GxK on X:", 
                         "PR_BGY1" = "PR:", 
                         "PR_BGY2" = "PR:", 
                         "NR_BGY1" = "NR:", 
                         "NR_BGY2" = "NR:", 
                         "mean_diff_BGY1" = expression(bar(Delta)*beta[GY[1]]*":"),
                         "mean_diff_BGY2" = expression(bar(Delta)*beta[GY[2]]*":"),
                         "mean_hat_diff_BGY" = expression(bar(Delta)*hat(beta)[GY]*":"), 
                         "mean_diff_diff_BGY" = expression(bar(Delta)[Delta*beta[GY]]*":"), 
                         "NR_GxK_on_Y" = "NR for GxK on Y:", 
                         "PR_GxK_on_Y" = "PR for GxK on Y:",
                         "PR_BXY1" = "PR:", 
                         "PR_BXY2" = "PR:", 
                         "NR_BXY1" = "NR:", 
                         "NR_BXY2" = "NR:", 
                         "mean_hat_BXY1" = expression(bar(hat(beta))[XY[1]]),
                         "mean_hat_BXY2" = expression(bar(hat(beta))[XY[2]]), 
                         "mean_diff_BXY1" = expression(bar(Delta)*beta[XY[1]]*":"),
                         "mean_diff_BXY2" = expression(bar(Delta)*beta[XY[2]]*":"),
                         "mean_hat_diff_BXY" = expression(bar(Delta)*hat(beta)[XY]*":"), 
                         "mean_diff_diff_BXY" = expression(bar(Delta)[Delta*beta[XY]]*":"),
                         "NR_XxK_on_Y" = "NR for XxK on Y:", 
                         "PR_XxK_on_Y" = "PR for XxK on Y:"
  ) 
  
  yvar_labels = list("n_aa_1" = "aa count in stratum 1",
                        "n_aA_1" = "aA count in stratum 1",
                        "n_AA_1" = "AA count in stratum 1",
                        "n_aa_2" = "aa count in stratum 2",
                        "n_aA_2" = "aA count in stratum 2",
                        "n_AA_2" = "AA count in stratum 2",
                        "HWE_P_1" = "P-value for HWE in stratum 1",
                        "HWE_P_2" = "P-value for HWE in stratum 2", 
                        "HWE_P_global" = "P-value for HWE in both strata", 
                        "HWE_CHISQ_1" = expression(chi**2*" statistic for HWE in stratum 1"),
                        "HWE_CHISQ_2" = expression(chi**2*" statistic for HWE in stratum 2"),
                        "HWE_CHISQ_global" = expression(chi**2*" statistic for HWE in both strata"),
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
                        "diff_BXY1" = expression(Delta*beta[XY[1]] == abs(beta[XY[1]] - hat(beta[XY[1]]))),                 
                        "diff_BXY2" = expression(Delta*beta[XY[2]] == abs(beta[XY[2]] - hat(beta[XY[2]]))),  
                        "hat_diff_BXY" = expression(Delta*hat(beta[XY]) == hat(beta[XY[2]]) - hat(beta[XY[1]])),                      
                        "diff_diff_BXY" = expression(Delta[Delta*beta[XY]] == abs(Delta*beta[XY] - Delta*hat(beta[XY]))),
                        "Z_diff_BXY" = expression(Z[Delta*beta[XY]]),       
                        "p_Z_diff_BXY" = expression(p[Z[Delta*beta[XY]]])
  )
  
  # if(metric_y %in% c("BGX_P_1", "BGX_P_2", "BGY_P_1", "BGY_P_2")){
  #   all_res[, metric_y] <- signif( all_res[, metric_y], digits = 2)
  # }
  
  ## Plot for 1 varying par
  if(scenario["sub_sce_varying_par"] %in% c("N", "r", "BGX1", "BXY1", "q1_q2", "BUX_BUY")){
    
    p = ggplot(all_res, aes(x = x, y = get(metric_y))) + 
      geom_point(color='gray', alpha = 0.3, size = 1) +
      theme_classic() +
      labs(x = var_labs[[unlist(scenario["sub_sce_varying_par"])]], y = yvar_labels[[metric_y]], subtitle = metric_to_show_lab[metric_to_show]) +
      geom_text(data = all_metrics, aes(x = x, y = max(all_res[metric_y])+(max(all_res[metric_y]) - min(all_res[metric_y]))/15, 
                                        label = signif(get(metric_to_show), digits = 3)), size = 2.4) +
      coord_cartesian(clip = "off", ylim = c(min(all_res[metric_y]), max(all_res[metric_y]))) +
      theme(plot.subtitle = element_text(size = 8), 
            axis.title.y = element_text(size = 9),
            axis.title.x = element_text(size = 8),
            axis.text = element_text(size = 7),
            axis.line = element_line(linewidth = 0.4),
            axis.ticks = element_line(linewidth = 0.3),
            plot.margin = unit(c(1.5, 1, 1, 1), "lines"))
  }

  
  yintercept = case_when(metric_y %in% c("HWE_P_1", "HWE_P_2", "HWE_P_global", 
                                         "BGX_P_1", "BGX_P_2", "p_Z_diff_BGX", 
                                         "BGY_P_1", "BGY_P_2", "p_Z_diff_BGY", "p_Z_diff_BXY") ~ 0.05, 
                         metric_y == "q1_ob" ~ unique(all_res$q1),
                         metric_y == "q2_ob" ~ unique(all_res$q2),
                         metric_y == "hat_BGX1" ~ unique(all_res$BGX1), 
                         metric_y == "hat_BGX2" ~ unique(all_res$BGX2), 
                         metric_y == "hat_diff_BGX" ~ unique(all_res$diff_BGX), 
                         metric_y == "hat_BGY1" ~ unique(all_res$BGY1), 
                         metric_y == "hat_BGY2" ~ unique(all_res$BGY2), 
                         metric_y == "hat_diff_BGY" ~ unique(all_res$diff_BGY), 
                         metric_y == "hat_BXY1" ~ unique(all_res$BXY1), 
                         metric_y == "hat_BXY2" ~ unique(all_res$BXY2), 
                         metric_y == "hat_diff_BXY" ~ unique(all_res$diff_BXY), .default = NA
  ) %>% unique()
  
  if(length(yintercept) == 1){
    if(!is.na(yintercept)){
      p <- p + geom_hline(yintercept = yintercept, linewidth = 0.5, alpha = 0.8, color = "gray40")
    }
  }
  
  if(plot_boxplots == T){
    p <- p + geom_boxplot(outliers = F, colour = "gray30", fill = NA, width = 0.3, 
                          box.linewidth = 0.2, whisker.linewidth = 0.2, staple.linewidth = 0.2, median.linewidth = 0.4) 
  }
  
  if(unique(all_res$sub_sce_varying_par) == "q1_q2"){
    p <- p + theme(axis.text = element_text(size = 6, angle = 45, hjust = 1))
  }
  
  return(p)
}


## Plot metrics across replicates, across subscenarios
for(scenario in scenarios){
  
  scenario <- scenarios[6, ]
  plot_dir_sce <- paste(plot_dir, paste(scenario, collapse = "/"), sep = "/")
  width <- c("N" = 17, "r" = 17, "q1_q2" = 22, "BGX1" = 17, "BXY1" = 17, "BUX_BUY" = 17)
  
  ## Plot HWE p-vals and chi-square stats 
  p1 <- plot_metric_across_subscenarios(scenario, "HWE_P_1", "FPR_HWE_1")
  p2 <- plot_metric_across_subscenarios(scenario, "HWE_P_2", "FPR_HWE_2")
  p3 <- plot_metric_across_subscenarios(scenario, "HWE_P_global", "FPR_HWE_global")
  p4 <- plot_metric_across_subscenarios(scenario, "HWE_CHISQ_1", "FPR_HWE_1", T)
  p5 <- plot_metric_across_subscenarios(scenario, "HWE_CHISQ_2", "FPR_HWE_2", T)
  p6 <- plot_metric_across_subscenarios(scenario, "HWE_CHISQ_global", "FPR_HWE_global", T)
  
  plot_grid(plotlist = list(p1, p2, p3, p4, p5, p6), ncol = 3, align = "vh")
  ggsave(filename = paste0(plot_dir_sce, "/HWE_metrics_across_replicates_x_subsce.pdf"), width = width[scenario$sub_sce_varying_par], height = 6)
  
  ## Plot ob MAF
  p7 <- plot_metric_across_subscenarios(scenario, "q1_ob", "mean_q1_ob")
  p8 <- plot_metric_across_subscenarios(scenario, "q2_ob", "mean_q2_ob")
  p9 <- plot_metric_across_subscenarios(scenario, "q_global_ob", "mean_q_global_ob")
  
  plot_grid(plotlist = list(p7, p8, p9), ncol = 3, align = "vh")
  ggsave(filename = paste0(plot_dir_sce, "/ob_MAF_across_replicates_x_subsce.pdf"), width = width[scenario$sub_sce_varying_par], height = 3)
  
  ## Plot num of aa, aA, AA individuals 
  p10 <- plot_metric_across_subscenarios(scenario, "n_aa_1", "mean_q1_ob")
  p11 <- plot_metric_across_subscenarios(scenario, "n_aA_1", "mean_q1_ob")
  p12 <- plot_metric_across_subscenarios(scenario, "n_AA_1", "mean_q1_ob")
  p13 <- plot_metric_across_subscenarios(scenario, "n_aa_2", "mean_q2_ob")
  p14 <- plot_metric_across_subscenarios(scenario, "n_aA_2", "mean_q2_ob")
  p15 <- plot_metric_across_subscenarios(scenario, "n_AA_2", "mean_q2_ob")
  
  plot_grid(plotlist = list(p10, p11, p12, p13, p14, p15), ncol = 3, align = "vh")
  ggsave(filename = paste0(plot_dir_sce, "/num_genotype_across_replicates_x_subsce.pdf"), width = width[scenario$sub_sce_varying_par], height = 6)
  
  ## Plot estimated BGX --always show true BGXk and PR (always T)
  p16 <- plot_metric_across_subscenarios(scenario, "hat_BGX1", "TPR_BGX1", T)
  p17 <- plot_metric_across_subscenarios(scenario, "hat_BGX2", "TPR_BGX2", T)
  p18 <- plot_metric_across_subscenarios(scenario, "se_hat_BGX1", "TPR_BGX1", T)
  p19 <- plot_metric_across_subscenarios(scenario, "se_hat_BGX2", "TPR_BGX2", T)
  p20 <- plot_metric_across_subscenarios(scenario, "BGX_t_stat_1", "TPR_BGX1", T)
  p21 <- plot_metric_across_subscenarios(scenario, "BGX_t_stat_2", "TPR_BGX2", T)
  p22 <- plot_metric_across_subscenarios(scenario, "BGX_P_1", "TPR_BGX1", T)
  p23 <- plot_metric_across_subscenarios(scenario, "BGX_P_2", "TPR_BGX2", T)
  
  plot_grid(plotlist = list(p16, p17, p18, p19, p20, p21, p22, p23), ncol = 2, align = "vh")
  ggsave(filename = paste0(plot_dir_sce, "/BGXk_estimated_across_replicates_x_subsce.pdf"), width = width[scenario$sub_sce_varying_par], height = 12)
  
  ## Plot IV strength 
  p24 <- plot_metric_across_subscenarios(scenario, "BGX_F_stat_1", "mean_IV_F_stat_1", T)
  p25 <- plot_metric_across_subscenarios(scenario, "BGX_F_stat_2", "mean_IV_F_stat_2", T)
  
  plot_grid(plotlist = list(p24, p25), ncol = 2, align = "vh")
  ggsave(filename = paste0(plot_dir_sce, "/IV_strength_across_replicates_x_subsce.pdf"), width = width[scenario$sub_sce_varying_par], height = 3)
  
  ## Plot true - hat BGXk
  p26 <- plot_metric_across_subscenarios(scenario, "diff_BGX1", "mean_diff_BGX1", T)
  p27 <- plot_metric_across_subscenarios(scenario, "diff_BGX2", "mean_diff_BGX2", T)
  
  plot_grid(plotlist = list(p26, p27), ncol = 2, align = "vh")
  ggsave(filename = paste0(plot_dir_sce, "/diff_BGXk_across_replicates_x_subsce.pdf"), width = width[scenario$sub_sce_varying_par], height = 3)
  
  ## Plot estimated between-strata BGX difference and true vs estimated difference
  p28 <- plot_metric_across_subscenarios(scenario, "hat_diff_BGX", "mean_hat_diff_BGX", T)
  p29 <- plot_metric_across_subscenarios(scenario, "diff_diff_BGX", "mean_diff_diff_BGX", T)
  p30 <- plot_metric_across_subscenarios(scenario, "Z_diff_BGX", "PR_GxK_on_X", T)
  p31 <- plot_metric_across_subscenarios(scenario, "p_Z_diff_BGX", "PR_GxK_on_X")
  
  plot_grid(plotlist = list(p28, p29, p30, p31), ncol = 2, align = "vh")
  ggsave(filename = paste0(plot_dir_sce, "/between_strata_BGX_diff_across_replicates_x_subsce.pdf"), width = width[scenario$sub_sce_varying_par], height = 6)
  
  ## Plot estimated BGYk across replicates x subscenario 
  p32 <- plot_metric_across_subscenarios(scenario, "hat_BGY1", "PR_BGY1", T)
  p33 <- plot_metric_across_subscenarios(scenario, "se_hat_BGY1", "PR_BGY1", T)
  p34 <- plot_metric_across_subscenarios(scenario, "BGY_t_stat_1", "PR_BGY1", T)
  p35 <- plot_metric_across_subscenarios(scenario, "BGY_P_1", "PR_BGY1", T)
  
  p36 <- plot_metric_across_subscenarios(scenario, "hat_BGY2", "PR_BGY2", T)
  p37 <- plot_metric_across_subscenarios(scenario, "se_hat_BGY2", "PR_BGY2", T)
  p38 <- plot_metric_across_subscenarios(scenario, "BGY_t_stat_2", "PR_BGY2", T)
  p39 <- plot_metric_across_subscenarios(scenario, "BGY_P_2", "PR_BGY2", T)
  
  plot_grid(plotlist = list(p32, p33, p34, p35, p36, p37, p38, p39), ncol = 2, align = "vh")
  ggsave(filename = paste0(plot_dir_sce, "/BGYk_estimated_across_replicates_x_subsce.pdf"), width = width[scenario$sub_sce_varying_par], height = 12)
  
  ## Plot true - hat BGYk
  p40 <- plot_metric_across_subscenarios(scenario, "diff_BGY1", "mean_diff_BGY1", T)
  p41 <- plot_metric_across_subscenarios(scenario, "diff_BGY2", "mean_diff_BGY2", T)
  
  plot_grid(plotlist = list(p40, p41), ncol = 2, align = "vh")
  ggsave(filename = paste0(plot_dir_sce, "/diff_BGYk_across_replicates_x_subsce.pdf"), width = width[scenario$sub_sce_varying_par], height = 3)
  
  
  ## Plot estimated between-strata BGY difference and true vs estimated difference
  p42 <- plot_metric_across_subscenarios(scenario, "hat_diff_BGY", "mean_hat_diff_BGY", T)
  p43 <- plot_metric_across_subscenarios(scenario, "diff_diff_BGY", "mean_diff_diff_BGY", T)
  p44 <- plot_metric_across_subscenarios(scenario, "Z_diff_BGY", "NR_GxK_on_Y", T)
  p45 <- plot_metric_across_subscenarios(scenario, "p_Z_diff_BGY", "PR_GxK_on_Y")
  
  plot_grid(plotlist = list(p42, p43, p44, p45), ncol = 2, align = "vh")
  ggsave(filename = paste0(plot_dir_sce, "/between_strata_BGY_diff_across_replicates_x_subsce.pdf"), width = width[scenario$sub_sce_varying_par], height = 6)
  
  ## Plot estimated BXYk 
  p46 <- plot_metric_across_subscenarios(scenario, "hat_BXY1", "mean_hat_BXY1", T)
  p47 <- plot_metric_across_subscenarios(scenario, "hat_BXY2", "mean_hat_BXY2", T)
  ## Plot true - estimated BXYk 
  p48 <- plot_metric_across_subscenarios(scenario, "diff_BXY1", "mean_diff_BXY1", T)
  p49 <- plot_metric_across_subscenarios(scenario, "diff_BXY2", "mean_diff_BXY2", T)
  ## Plot estimated between-strata BXY diff and true vs estimated diff
  p50 <- plot_metric_across_subscenarios(scenario, "hat_diff_BXY", "mean_hat_diff_BXY", T)
  p51 <- plot_metric_across_subscenarios(scenario, "diff_diff_BXY", "mean_diff_diff_BXY", T)
  
  p52 <- plot_metric_across_subscenarios(scenario, "Z_diff_BXY", "PR_XxK_on_Y", T)
  p53 <- plot_metric_across_subscenarios(scenario, "p_Z_diff_BXY", "PR_XxK_on_Y")
  
  plot_grid(plotlist = list(p46, p47, p48, p49, p50, p51, p52, p53), ncol = 2, align = "vh")
  ggsave(filename = paste0(plot_dir_sce, "/BXYk_estimated_across_replicates_x_subsce.pdf"), width = width[scenario$sub_sce_varying_par], height = 12)
  
}












