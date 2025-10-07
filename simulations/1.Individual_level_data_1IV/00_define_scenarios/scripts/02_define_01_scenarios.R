
library(tidyr)
library(dplyr)
library(sessioninfo)

rm(list = ls())


# ==============================================================================================
#  2. Scenario 01: No genetic effect difference Δβɢx = 0 and causal effect difference Δβxʏ ≠ 0 
# ==============================================================================================

## Define dirs
scripts_dir <- paste(getwd(), "simulations", "1.Individual_level_data_1IV", "00_define_scenarios", "scripts", sep = "/")
input_dir <- paste(getwd(), "simulations", "1.Individual_level_data_1IV", "00_define_scenarios", "inputs", sep = "/")
out_dir <- paste(getwd(), "simulations", "1.Individual_level_data_1IV", "00_define_scenarios", "outputs", sep = "/")
plot_dir <- paste(getwd(), "simulations", "1.Individual_level_data_1IV", "00_define_scenarios", "plots", sep = "/")


scenarios01 <- vector()
main_scenario =  "scenario.01"
diff_BGX = 0
diff_BXYs = seq(from = -1.7, to = 1.3, by = 0.5) 
# Δβxʏ = {-1.7  -1.2  -0.7  -0.2  0.3  0.8  1.3} -> βxʏ₂ = {-1.0  -0.5  0.0  0.5  1.0  1.5  2.0} (base case)


## Same scenarios as in 00 but multiplied for each Δβxʏ 
load(paste0(out_dir, "/scenario.00.Rdata"), verbose = T)
# Loading objects:
#   scenarios00

scenarios01 <- vector()
for(i in 1:length(diff_BXYs)){
  
  temp <- scenarios00
  temp$diff_BXY <- diff_BXYs[i]
  temp$main_scenario_val <- paste0("diff_BXY=",diff_BXYs[i])
  
  scenarios01 <- rbind(scenarios01, temp)
}

scenarios01$scenario <- "scenario.01"

scenarios01[, 5:14] <- apply(scenarios01[, 5:14] , 2, as.numeric)
save(scenarios01, file = paste0(out_dir, "/scenario.01.Rdata"))







## Reproducibility info
session_info()

# ─ Session info ────────────────────────────────────────────────────────────────────────────────────────────────
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
# ─ Packages ────────────────────────────────────────────────────────────────────────────────────────────────────
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
# ───────────────────────────────────────────────────────────────────────────────────────────────────────────────



