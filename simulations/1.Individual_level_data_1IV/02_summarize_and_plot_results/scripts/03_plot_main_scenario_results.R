
library(tidyr)
library(dplyr)
library(ggplot2)
library(cowplot)
library(sessioninfo)

rm(list = ls())

# ______________________________________________________________________________
#   3.0 Plot results across varying parameter(s) scenarios for a main scenario  
# ______________________________________________________________________________
#  This script compares results across the varying parameter(s) scenarios of a
#  main Δβɢx and Δβxʏ scenario. 
# ______________________________________________________________________________

input_dir <- paste(getwd(), "simulations", "1.Individual_level_data_1IV", "02_summarize_and_plot_results", "inputs", sep = "/")
out_dir <- paste(getwd(), "simulations", "1.Individual_level_data_1IV", "02_summarize_and_plot_results", "outputs", sep = "/")
plot_dir <- paste(getwd(), "simulations", "1.Individual_level_data_1IV", "02_summarize_and_plot_results", "plots", sep = "/")

## Input dirs
input_dir00 <- paste(getwd(), "simulations", "1.Individual_level_data_1IV", "00_define_scenarios", "outputs", sep = "/")

## Scenarios
all_scenarios00 <- get(load(paste0(input_dir00, "/scenario.00.Rdata")))
scenarios00 <- all_scenarios00[, c("main_scenario", "main_scenario_val", "sub_sce_varying_par")] %>% unique()

# all_scenarios01 <- get(load(paste0(input_dir00, "/scenario.01.Rdata")))
# scenarios01 <- all_scenarios01[, c("main_scenario", "main_scenario_val", "sub_sce_varying_par")] %>% unique()



## Plot PR/NR for effect estimates and estimated effect differences between strata
plot_PN_rate <- function(data, rate){
  
  ## Par label
  var_labs = list("N" = "paste('Increasing N')", 
                  "r" = "paste('Ratio of stratum sample sizes N1/N2')", 
                  "BGX1" = "beta[GX[1]]",
                  "BXY1" = "beta[XY[1]]", 
                  "BUX_BUY" = "paste(beta[UX]*' = '*beta[UY])", 
                  "q1_q2" = "paste(q[1]*', '*q[2])",
                  "diff_BGX" = "Delta*beta[GX]",
                  "diff_BXY" = "Delta*beta[XY]", 
                  "N_r" = "paste(N*', '*r)", 
                  "N_q1_q2" = "paste(N*', '*q[1]*', '*q[2])")
  
  ## Par color
  var_colors = list("N" = "deeppink", 
                  "r" = "lightsalmon", 
                  "BGX1" = "#5CACEE",
                  "BXY1" = "#7CCD7C", 
                  "BUX_BUY" = "yellow3", 
                  "q1_q2" = "indianred2",
                  "diff_BGX" = "skyblue",
                  "diff_BXY" = "green2", 
                  "N_r" = "#EE9A49", 
                  "N_q1_q2" = "darkorange2")
  
  colors = list("N" = c("10000" = "#EED2EE", "40000" = "#EEAEEE", "70000" = "#CD69C9", "100000" = "#8B4789"), 
                "r" = c("0.2" = "#FFDAB9", "0.5" = "#FFA07A", "1" = "#FF8247", "2" = "#FF6347", "5" = "#CD4F39"), 
                "BGX1" = c("#B2DFEE", "#87CEEB", "#7EC0EE", "#9FB6CD", "#4A708B"),
                "BXY1" = c("#B4EEB4", "#9BCD9B", "#7CCD7C", "#8FBC8F", "#698B69"), 
                "BUX_BUY" = c("#EEE685", "#EEDC82", "#CDC673", "#CDCD00","#8B8B00"), 
                "q1_q2" = c("#FFE4E1", "#EED5D2", "#EEB4B4",
                            "#EED5D2", "#EEB4B4", "#EE799F", 
                            "#EEB4B4", "#EE799F", "#CD6090"))
  
  
  
  ## If rate is T/F
  PR_BGY1_pre <- if_else(unique(data$BGY1) == 0, "F", "T")
  PR_BGY2_pre <- if_else(unique(data$BGY2) == 0, "F", "T")
  NR_BGY1_pre <- if_else(unique(data$BGY1) == 0, "T", "F")
  NR_BGY2_pre <- if_else(unique(data$BGY2) == 0, "T", "F")
  PR_GxK_on_X_pre <- if_else(unique(data$diff_BGX) == 0, "F", "T")
  PR_GxK_on_Y_pre <- if_else(unique(data$diff_BGY) == 0, "F", "T")
  PR_XxK_on_Y_pre <- PR_XxK_on_Y_2nd_pre <- if_else(unique(data$diff_BXY) == 0, "F", "T")
  NR_GxK_on_X_pre <- if_else(unique(data$diff_BGX) == 0, "T", "F")
  NR_GxK_on_Y_pre <- if_else(unique(data$diff_BGY) == 0, "T", "F")
  NR_XxK_on_Y_pre <- NR_XxK_on_Y_2nd_pre <- if_else(unique(data$diff_BXY) == 0, "T", "F")
  
  tf <- c("TPR_BGX1" = '', 
          "TPR_BGX2" = '', 
          "FNR_BGX1" = '', 
          "FNR_BGX2" = '', 
          "PR_BGY1" = PR_BGY1_pre,  
          "PR_BGY2" = PR_BGY2_pre, 
          "NR_BGY1" = NR_BGY1_pre, 
          "NR_BGY2" = NR_BGY2_pre, 
          "PR_GxK_on_X" = PR_GxK_on_X_pre, 
          "PR_GxK_on_Y" = PR_GxK_on_Y_pre, 
          "PR_XxK_on_Y" = PR_XxK_on_Y_pre, 
          "PR_XxK_on_Y_2nd" = PR_XxK_on_Y_2nd_pre,
          "NR_GxK_on_X" = NR_GxK_on_X_pre, 
          "NR_GxK_on_Y" = NR_GxK_on_Y_pre, 
          "NR_XxK_on_Y" = NR_XxK_on_Y_pre, 
          "NR_XxK_on_Y_2nd" = NR_XxK_on_Y_2nd_pre)
  
  ## Rate label
  rate_label <- list("TPR_BGX1" = expression("TPR" ~ "for" ~ beta[GX[1]]), 
                  "TPR_BGX2" = expression("TPR" ~ "for" ~ beta[GX[2]]), 
                  "FNR_BGX1" = expression("FNR" ~ "for" ~ beta[GX[1]]), 
                  "FNR_BGX2" = expression("FNR" ~ "for" ~ beta[GX[2]]), 
                  "PR_BGY1" = paste0("paste(", tf[rate], "*'PR for ', beta[GY[1]])"),
                  "PR_BGY2" = paste0("paste(", tf[rate], "*'PR for ', beta[GY[2]])"),
                  "NR_BGY1" = paste0("paste(", tf[rate], "*'NR for ', beta[GY[1]])"),
                  "NR_BGY2" = paste0("paste(", tf[rate], "*'NR for ', beta[GY[2]])"),
                  "PR_GxK_on_X" = paste0("paste(", tf[rate], "*'PR for ', Delta*beta[GX])"),
                  "PR_GxK_on_Y" = paste0("paste(", tf[rate], "*'PR for ', Delta*beta[GY])"),
                  "PR_XxK_on_Y" = paste0("paste(", tf[rate], "*'PR for ', Delta*beta[XY])"),
                  "PR_XxK_on_Y_2nd" = paste0("paste(", tf[rate], "*'PR for ', Delta*beta[XY]~'(corrected se)')"),
                  "NR_GxK_on_X" = paste0("paste(", tf[rate], "*'NR for ', Delta*beta[GX])"),
                  "NR_GxK_on_Y" = paste0("paste(", tf[rate], "*'NR for ', Delta*beta[GY])"),
                  "NR_XxK_on_Y" = paste0("paste(", tf[rate], "*'NR for ', Delta*beta[XY])"),
                  "NR_XxK_on_Y_2nd" = paste0("paste(", tf[rate], "*'NR for ', Delta*beta[XY]~'(corrected se)')"))
   
  ## Plot scenario settings as title
  title = list("scenario.00" = "paste(Delta*beta[GX]==0~', '~Delta*beta[XY]==0)", 
               "scenario.01" = "paste(Delta*beta[GX]==0~', '~abs(Delta*beta[XY])>0)")
  
  ## If only 1 par varying
  if(unique(stringr::str_count(data$sub_sce_varying_par_value, "=")) == 1 | unique(data$sub_sce_varying_par) %in% c("q1_q2", "BUX_BUY")){
    
    if(unique(data$sub_sce_varying_par) == "q1_q2"){
      data$x = gsub("q1=", "", data$sub_sce_varying_par_value) %>% gsub("q2=", " ", .)
      data$x <- factor(data$x, levels = unique(data$x))
    } 
    else if(unique(data$sub_sce_varying_par) == "BUX_BUY"){
      data$x = gsub("BUX=", "", data$sub_sce_varying_par_value) %>% gsub("BUY=", " ", .)
      data$x <- factor(data$x, levels = unique(data$x))
    }
    else{
      data$x = gsub(".*=", "", data$sub_sce_varying_par_value)
      data$x <- factor(data$x, levels = unique(data$x))
    }
    
    ## Plot rate across values
    p = ggplot(data, aes(x=x, y=get(rate), group=1)) +
      geom_point(size = 1.2, color = var_colors[[unique(data$sub_sce_varying_par)]]) +
      geom_line(colour = var_colors[[unique(data$sub_sce_varying_par)]]) +
      theme_bw() +
      scale_y_continuous(limits = c(0,1)) +
      labs(x = parse(text = var_labs[unique(data$sub_sce_varying_par)]),
           y = parse(text = rate_label[[rate]]),
           title = parse(text = title[[unique(data$main_scenario)]])) +
      theme(plot.title = element_text(size = 9),
            axis.title = element_text(size = 9),
            axis.text = element_text(size = 8), 
            panel.grid = element_line(linewidth = 0.2))
  }
  
  ## If >1 par varying 
  else{
    
    if(unique(stringr::str_count(data$sub_sce_varying_par_value, "=")) == 2){
  
      var1 <- strsplit(unique(data$sub_sce_varying_par), "_")[[1]][1]
      var2 <- strsplit(unique(data$sub_sce_varying_par), "_")[[1]][2]
      data$x <- gsub(paste0(var1, "="), "", data$sub_sce_varying_par_value) %>% gsub(",.*", "", .)
      data$x <- factor(data$x, levels = unique(data$x))
      data$x2 <- gsub(".*=", "", data$sub_sce_varying_par_value)
      data$x2 <- factor(data$x2, levels = unique(data$x2))
    
    } 
    else if(unique(stringr::str_count(data$sub_sce_varying_par_value, "=")) == 3){
      
      if(grepl("^q1_q2", unique(data$sub_sce_varying_par))){
        data$x <- paste0(data$q1, ", ", data$q2)
        data$x <- factor(data$x, levels = unique(data$x))
        var2 <- gsub("q1_q2_", "", unique(data$sub_sce_varying_par))
        data$x2 <- gsub(".*,", "", data$sub_sce_varying_par_value) %>% gsub(paste0(var2, "="), "", .)
        data$x2 <- factor(data$x2, levels = unique(data$x2))
        
      } else{
        var1 = strsplit(unique(data$sub_sce_varying_par), "_q1_q2")[[1]]
        var2 = "q1_q2"
        data$x <- gsub(",q.*", "", data$sub_sce_varying_par_value) %>% gsub(paste0(var1, "="), "", .)
        data$x <- factor(data$x, levels = unique(data$x))
        data$x2 <- paste0(data$q1, ", ", data$q2)
        data$x2 <- factor(data$x2, levels = unique(data$x2))
      }
    }
    
    ## Plot rate across values
    p = ggplot(data, aes(x=x, y=get(rate), color = x2, group = x2)) +
      geom_point(size = 1.2) +
      geom_line(alpha = 0.5) +
      theme_bw() +
      scale_color_manual(values = colors[[var2]]) + 
      scale_y_continuous(limits = c(0,1)) +
      labs(x = parse(text = var_labs[var1]),
           y = parse(text = rate_label[[rate]]),
           title = parse(text = title[[unique(data$main_scenario)]]), 
           color = parse(text = var_labs[var2])) +
      theme(plot.title = element_text(size = 9),
            legend.title = element_text(size = 8),
            legend.text = element_text(size = 7),
            axis.title = element_text(size = 9),
            axis.text = element_text(size = 8), 
            panel.grid = element_line(linewidth = 0.2), 
            legend.key.height = unit(0.4, "cm"),
            legend.key.width  = unit(0.3, "cm"))
    
  }
  
  if(grepl("^q1_q2", unique(data$sub_sce_varying_par)) | grepl("BUX_BUY", unique(data$sub_sce_varying_par))){
    p = p + theme(axis.text.x = element_text(size = 7.5, angle = 45, hjust = 1))
  }
  return(p)
}


## Rates to plot 
rates = c("TPR_BGX1", "FNR_BGX1", "TPR_BGX2", "FNR_BGX2", 
          "PR_BGY1", "NR_BGY1", "PR_BGY2", "NR_BGY2", 
          "PR_GxK_on_X", "PR_GxK_on_Y", "PR_XxK_on_Y", "PR_XxK_on_Y_2nd", 
          "NR_GxK_on_X", "NR_GxK_on_Y", "NR_XxK_on_Y", "NR_XxK_on_Y_2nd")

## Plots x scenario
for(i in 1:nrow(scenarios00)){
  
  ## Collect metrics
  scenario <- scenarios00[i, ]
  plot_dir_sce <- paste(plot_dir, paste(scenario, collapse = "/"), sep = "/")
  width <- c("N" = 21, "r" = 15, "q1_q2" = 22, "BGX1" = 15, "BXY1" = 15, "BUX_BUY" = 15, 
             "N_r" = 20, "N_q1_q2" = 15, "N_BGX1" = 15, "N_BXY1" = 15, "r_q1_q2" = 15,
             "r_BGX1" = 15, "r_BXY1" = 15, "q1_q2_BGX1" = 22, "q1_q2_BXY1" = 22, "BGX1_BXY1" = 15)
  
  metrics <- get(load(paste(out_dir, paste(scenario, collapse = "/"), "agg_metrics_subscenarios.Rdata", sep = "/")))
  all_res <- get(load(paste(out_dir, paste(scenario, collapse = "/"), "agg_results_subscenarios_replicates.Rdata", sep = "/")))
  settings <- all_res[, c("main_scenario", "main_scenario_val", 
                          "sub_sce_varying_par", "sub_sce_varying_par_value", 
                          "N", "r", "q1", "q2", "BGX1", "diff_BGX", "BGX2", 
                          "BXY1", "diff_BXY", "BXY2", 
                          "BGY1", "diff_BGY", "BGY2", 
                          "BUX", "BUY")] %>% unique()
  
  ## Append all scenario settings
  metrics = metrics %>% left_join(settings, by = c("main_scenario", "main_scenario_val", 
                                         "sub_sce_varying_par", "sub_sce_varying_par_value"))
  
  j = 1
  plots <- list()
  for(rate in rates){
    plots[[j]] = plot_PN_rate(metrics, rate)
    j = j + 1
  }
  
  plot_grid(plotlist = plots, ncol = 4, align = "vh")
  ggsave(filename = paste0(plot_dir_sce, "/TF_PN_rates.pdf"), width = width[scenario$sub_sce_varying_par], height = 10)
}







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









