
library(tidyr)
library(dplyr)
library(sessioninfo)

rm(list = ls())


# ==============================================================================================
#  4. Scenario 1,1: Genetic effect difference Δβɢx ≠ 0 and causal effect difference Δβxʏ ≠ 0 
# ==============================================================================================

## Define dirs
scripts_dir <- paste(getwd(), "simulations", "1.Individual_level_data_1IV", "00_define_scenarios", "scripts", sep = "/")
input_dir <- paste(getwd(), "simulations", "1.Individual_level_data_1IV", "00_define_scenarios", "inputs", sep = "/")
out_dir <- paste(getwd(), "simulations", "1.Individual_level_data_1IV", "00_define_scenarios", "outputs", sep = "/")
plot_dir <- paste(getwd(), "simulations", "1.Individual_level_data_1IV", "00_define_scenarios", "plots", sep = "/")


scenarios11 <- vector()
main_scenario =  "scenario.11"
diff_BGXs = seq(from = -0.6, to = 1.4, by = 0.4) 
diff_BXYs = seq(from = -1.7, to = 1.3, by = 0.5) 

# Δβɢx = {-0.6  -0.2  0.2  0.6  1.0  1.4} -> βɢx₂ = {-0.1  0.3  0.7  1.1  1.5  1.9} (base case)
# Δβxʏ = {-1.7  -1.2  -0.7  -0.2  0.3  0.8  1.3} -> βxʏ₂ = {-1.0  -0.5  0.0  0.5  1.0  1.5  2.0} (base case)


## Same scenarios as in 00 for each Δβɢx, Δβxʏ combination
load(paste0(out_dir, "/scenario.00.Rdata"), verbose = T)
# Loading objects:
#   scenarios00

scenarios11 <- vector()
for(j in 1:length(diff_BXYs)){
  
  for(i in 1:length(diff_BGXs)){
    
    temp <- scenarios00
    temp$diff_BXY <- diff_BXYs[j]
    temp$diff_BGX <- diff_BGXs[i]
    temp$main_scenario_val <- paste0("diff_BGX=",diff_BGXs[i], ",diff_BXY=", diff_BXYs[j])
    
    scenarios11 <- rbind(scenarios11, temp)
  }
}

scenarios11$main_scenario <- "scenario.11"

scenarios11[, 5:14] <- apply(scenarios11[, 5:14] , 2, as.numeric)
save(scenarios11, file = paste0(out_dir, "/scenario.11.Rdata"))
write.table(scenarios11, file = paste0(out_dir, "/scenario.11.csv"), col.names = F, row.names = F, quote = F, na = "NA", sep = "\t")







## Reproducibility info
session_info()

# ─ Session info ────────────────────────────────────────────────────────────────────────────────────────────────────
# setting  value
# version  R version 4.5.1 (2025-06-13)
# os       Ubuntu 22.04.5 LTS
# system   x86_64, linux-gnu
# ui       RStudio
# language (EN)
# collate  en_GB.UTF-8
# ctype    en_GB.UTF-8
# tz       Europe/London
# date     2025-10-13
# rstudio  2024.12.1+563 Kousa Dogwood (server)
# pandoc   NA
# quarto   1.6.42 @ /usr/local/bin/quarto
# 
# ─ Packages ────────────────────────────────────────────────────────────────────────────────────────────────────────
# package      * version date (UTC) lib source
# cli            3.6.5   2025-04-23 [1] CRAN (R 4.5.1)
# cowplot      * 1.2.0   2025-07-07 [1] CRAN (R 4.5.1)
# crayon         1.5.3   2024-06-20 [1] CRAN (R 4.5.1)
# dichromat      2.0-0.1 2022-05-02 [1] CRAN (R 4.5.1)
# dplyr        * 1.1.4   2023-11-17 [1] CRAN (R 4.5.1)
# farver         2.1.2   2024-05-13 [1] CRAN (R 4.5.1)
# generics       0.1.4   2025-05-09 [1] CRAN (R 4.5.1)
# ggplot2      * 4.0.0   2025-09-11 [1] CRAN (R 4.5.1)
# glue           1.8.0   2024-09-30 [1] CRAN (R 4.5.1)
# gtable         0.3.6   2024-10-25 [1] CRAN (R 4.5.1)
# labeling       0.4.3   2023-08-29 [1] CRAN (R 4.5.1)
# lifecycle      1.0.4   2023-11-07 [1] CRAN (R 4.5.1)
# magrittr       2.0.4   2025-09-12 [1] CRAN (R 4.5.1)
# pillar         1.11.1  2025-09-17 [1] CRAN (R 4.5.1)
# pkgconfig      2.0.3   2019-09-22 [1] CRAN (R 4.5.1)
# purrr          1.1.0   2025-07-10 [1] CRAN (R 4.5.1)
# R6             2.6.1   2025-02-15 [1] CRAN (R 4.5.1)
# ragg           1.5.0   2025-09-02 [1] CRAN (R 4.5.1)
# RColorBrewer   1.1-3   2022-04-03 [1] CRAN (R 4.5.1)
# rlang          1.1.6   2025-04-11 [1] CRAN (R 4.5.1)
# rstudioapi     0.17.1  2024-10-22 [1] CRAN (R 4.5.1)
# S7             0.2.0   2024-11-07 [1] CRAN (R 4.5.1)
# scales         1.4.0   2025-04-24 [1] CRAN (R 4.5.1)
# sessioninfo  * 1.2.3   2025-02-05 [1] CRAN (R 4.5.1)
# stringi        1.8.7   2025-03-27 [1] CRAN (R 4.5.1)
# stringr        1.5.2   2025-09-08 [1] CRAN (R 4.5.1)
# systemfonts    1.2.3   2025-04-30 [1] CRAN (R 4.5.1)
# textshaping    1.0.3   2025-09-02 [1] CRAN (R 4.5.1)
# tibble         3.3.0   2025-06-08 [1] CRAN (R 4.5.1)
# tidyr        * 1.3.1   2024-01-24 [1] CRAN (R 4.5.1)
# tidyselect     1.2.1   2024-03-11 [1] CRAN (R 4.5.1)
# vctrs          0.6.5   2023-12-01 [1] CRAN (R 4.5.1)
# withr          3.0.2   2024-10-28 [1] CRAN (R 4.5.1)
# 
# [1] /home/daianna/R/x86_64-pc-linux-gnu-library/4.5
# [2] /usr/local/lib/R/site-library
# [3] /usr/lib/R/site-library
# [4] /usr/lib/R/library
# * ── Packages attached to the search path.
# 
# ───────────────────────────────────────────────────────────────────────────────────────────────────────────────────










