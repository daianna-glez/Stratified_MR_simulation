
library(tidyr)
library(dplyr)
library(sessioninfo)

## Initialize
rm(list = ls())


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
scenarios00 <- vector()
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

for (N in Ns){
  scenarios00 <- rbind(scenarios00, c("scenario" = main_scenario, "sub_sce_varying_par" = sub_sce_varying_par, "sub_sce_varying_par_value" = paste0(sub_sce_varying_par,"=", N), 
                                  "N" = N, "r" = r_base, "q1" = q1_base, "q2" = q2_base, 
                                  "BGX1" = BGX1_base, "diff_BGX" = diff_BGX, "BXY1" = BXY1_base, "diff_BXY" = diff_BXY)) %>% as.data.frame()
}

# ______________________________________________________________________________________________
#  2. Varying r: N = 10k, r = {0.20, ..., 1, ..., 5}, q₁ = q₂ = 0.25, βɢx₁ = 0.5, βxʏ₁ = 0.7
# ______________________________________________________________________________________________
sub_sce_varying_par = "r"
rs <- signif(c(rev(1/(1:5)), 2:5), digits = 2)
# r = {0.20  0.25  0.33  0.50  1.0  2.0  3.0  4.0  5.0}

for (r in rs){
  scenarios00 <- rbind(scenarios00, c(main_scenario, sub_sce_varying_par, paste0(sub_sce_varying_par, "=", r), 
                                  N_base, "r" = r, q1_base, q2_base, BGX1_base, diff_BGX, BXY1_base, diff_BXY)) 
}
# ______________________________________________________________________________________________
#  3. Varying q₁ and q₂: N = 10k, r = 1, q₁ = q₂ = {0.05, ..., 0.45}, βɢx₁ = 0.5, βxʏ₁ = 0.7
# ______________________________________________________________________________________________
sub_sce_varying_par = "q1_q2"
q1s <- q2s <- seq(from = 0.05, to = 0.5, by = 0.1)
# q₁ = q₂ = {0.05  0.15  0.25  0.35  0.45}

for (q1 in q1s){
  for (q2 in q2s){
    scenarios00 <- rbind(scenarios00, c(main_scenario, sub_sce_varying_par, paste0("q1=", q1, ",q2=", q2), 
                                    N_base, r_base, "q1" = q1, "q2" = q2, BGX1_base, diff_BGX, BXY1_base, diff_BXY)) 
  }
}

# ______________________________________________________________________________________________
#  4. Varying βɢx₁: N = 10k, r = 1, q₁ = q2 = 0.25, βɢx₁ = {0.01, 0.09,..., 0.73}, βxʏ₁ = 0.7
# ______________________________________________________________________________________________
sub_sce_varying_par = "BGX1"
BGX1s <- seq(from = 0.01, to = 0.8, by = 0.08)
# βɢx₁ = {0.01  0.09  0.17  0.25  0.33  0.41  0.49  0.57  0.65  0.73}

for (BGX1 in BGX1s){
  scenarios00 <- rbind(scenarios00, c(main_scenario, sub_sce_varying_par, paste0(sub_sce_varying_par, "=", BGX1), 
                                  N_base, r_base, q1_base, q2_base, "BGX1" = BGX1, diff_BGX, BXY1_base, diff_BXY)) 
}
# ______________________________________________________________________________________________
#  5. Varying βxʏ₁: N = 10k, r = 1, q₁ = q2 = 0.25, βɢx₁ = 0.5, βxʏ₁ = {0.0, ..., 1.35}
# ______________________________________________________________________________________________
sub_sce_varying_par = "BXY1"
BXY1s <- seq(from = 0.0, to = 1.35, by = 0.15)
# βxʏ₁ = {0.00  0.15  0.30  0.45  0.60  0.75  0.90  1.05  1.20  1.35}

for (BXY1 in BXY1s){
  scenarios00 <- rbind(scenarios00, c(main_scenario, sub_sce_varying_par, paste0(sub_sce_varying_par, "=", BXY1), 
                                  N_base, r_base, q1_base, q2_base, BGX1_base, diff_BGX, "BXY1" = BXY1, diff_BXY)) 
}


#  Sub-scenarios for two varying parameters at a time: 
# ______________________________________________________________________________________________
#  6. Varying N and r: N = {10k, ..., 100k} r = {0.20, ..., 1, ..., 5},
#                      q₁ = q₂ = 0.25, βɢx₁ = 0.5, βxʏ₁ = 0.7
# ______________________________________________________________________________________________
sub_sce_varying_par = "N_r"

for (N in Ns){
  for (r in rs){
    scenarios00 <- rbind(scenarios00, c(main_scenario, sub_sce_varying_par, paste0("N=", N, ",r=", r), 
                                    "N" = N, "r" = r, q1_base, q2_base, BGX1_base, diff_BGX, BXY1_base, diff_BXY)) 
  }
}

# ______________________________________________________________________________________________
#  7. Varying N, q₁ and q₂: N = {10k, ..., 100k}, q₁ = q₂ = {0.05, ..., 0.45}, 
#                           r = 1, βɢx₁ = 0.5, βxʏ₁ = 0.7
# ______________________________________________________________________________________________
sub_sce_varying_par = "N_q1_q2"

for (N in Ns){
  for (q1 in q1s){
    for (q2 in q2s){
      scenarios00 <- rbind(scenarios00, c(main_scenario, sub_sce_varying_par, paste0("N=", N, ",q1=", q1, ",q2=", q2), 
                                      "N" = N, r_base, "q1" = q1, "q2" = q2, BGX1_base, diff_BGX, BXY1_base, diff_BXY)) 
    }
  }
}

# ______________________________________________________________________________________________
#  8. Varying N and βɢx₁: N = {10k, ..., 100k}, βɢx₁ = {0.01, 0.09,..., 0.73}, 
#                         r = 1, q₁ = q₂ = 0.25, βxʏ₁ = 0.7
# ______________________________________________________________________________________________
sub_sce_varying_par = "N_BGX1"

for (N in Ns){
  for (BGX1 in BGX1s){
    scenarios00 <- rbind(scenarios00, c(main_scenario, sub_sce_varying_par, paste0("N=", N, ",BGX1=", BGX1), 
                                    "N" = N, r_base, q1_base, q2_base, "BGX1" = BGX1, diff_BGX, BXY1_base, diff_BXY)) 
  }
}

# ______________________________________________________________________________________________
#  9. Varying N and βxʏ₁: N = {10k, ..., 100k}, βxʏ₁ = {0.1, ..., 1.45}, 
#                         r = 1, q₁ = q₂ = 0.25, βɢx₁ = 0.5
# ______________________________________________________________________________________________
sub_sce_varying_par = "N_BXY1"

for (N in Ns){
  for (BXY1 in BXY1s){
    scenarios00 <- rbind(scenarios00, c(main_scenario, sub_sce_varying_par, paste0("N=", N, ",BXY1=", BXY1), 
                                    "N" = N, r_base, q1_base, q2_base, BGX1_base, diff_BGX, "BXY1" = BXY1, diff_BXY)) 
  }
}

# ______________________________________________________________________________________________
#  10. Varying r, q₁ and q₂: r = {0.20, ..., 1, ..., 5}, q₁ = q₂ = {0.05, ..., 0.45},
#                            N = 10k, βɢx₁ = 0.5, βxʏ₁ = 0.7
# ______________________________________________________________________________________________
sub_sce_varying_par = "r_q1_q2"

for (r in rs){
  for (q1 in q1s){
    for (q2 in q2s){
      scenarios00 <- rbind(scenarios00, c(main_scenario, sub_sce_varying_par, paste0("r=", r, ",q1=", q1, ",q2=", q2), 
                                      N_base, "r" = r, "q1" = q1, "q2" = q2, BGX1_base, diff_BGX, BXY1_base, diff_BXY)) 
    }
  }
}

# ______________________________________________________________________________________________
#  11. Varying r and βɢx₁: r = {0.20, ..., 1, ..., 5}, βɢx₁ = {0.01, 0.09,..., 0.73}, 
#                          N = 10k, q₁ = q₂ = 0.25, βxʏ₁ = 0.7
# ______________________________________________________________________________________________
sub_sce_varying_par = "r_BGX1"

for (r in rs){
  for (BGX1 in BGX1s){
    scenarios00 <- rbind(scenarios00, c(main_scenario, sub_sce_varying_par, paste0("r=", r, ",BGX1=", BGX1), 
                                    N_base, "r" = r, q1_base, q2_base, "BGX1" = BGX1, diff_BGX, BXY1_base, diff_BXY)) 
  }
}

# ______________________________________________________________________________________________
#  12. Varying r and βxʏ₁: r = {0.20, ..., 1, ..., 5}, βxʏ₁ = {0.1, ..., 1.45},  
#                          N = 10k, q₁ = q₂ = 0.25, βɢx₁ = 0.5
# ______________________________________________________________________________________________
sub_sce_varying_par = "r_BXY1"

for (r in rs){
  for (BXY1 in BXY1s){
    scenarios00 <- rbind(scenarios00, c(main_scenario, sub_sce_varying_par, paste0("r=", r, ",BXY1=", BXY1), 
                                    N_base, "r" = r, q1_base, q2_base, BGX1_base, diff_BGX, "BXY1" = BXY1, diff_BXY)) 
  }
}

# ______________________________________________________________________________________________
#  13. Varying q₁, q₂ and βɢx₁: q₁ = q₂ = {0.05, ..., 0.45}, βɢx₁ = {0.01, 0.09,..., 0.73},  
#                               N = 10k, r = 1, βxʏ₁ = 0.7
# ______________________________________________________________________________________________
sub_sce_varying_par = "q1_q2_BGX1"

for (q1 in q1s){
  for (q2 in q2s){
    for(BGX1 in BGX1s){
      scenarios00 <- rbind(scenarios00, c(main_scenario, sub_sce_varying_par, paste0("q1=", q1, ",q2=", q2, ",BGX1=",BGX1), 
                                      N_base, r_base, "q1" = q1, "q2" = q2, "BGX1" = BGX1, diff_BGX, BXY1_base, diff_BXY)) 
    }
  }
}

# ______________________________________________________________________________________________
#  14. Varying q₁, q₂ and βxʏ₁: q₁ = q₂ = {0.05, ..., 0.45}, βxʏ₁ = {0.1, ..., 1.45},  
#                               N = 10k, r = 1, βɢx₁ = 0.5
# ______________________________________________________________________________________________
sub_sce_varying_par = "q1_q2_BXY1"

for (q1 in q1s){
  for (q2 in q2s){
    for(BXY1 in BXY1s){
      scenarios00 <- rbind(scenarios00, c(main_scenario, sub_sce_varying_par, paste0("q1=", q1, ",q2=", q2, ",BXY1=",BXY1), 
                                      N_base, r_base, "q1" = q1, "q2" = q2, BGX1_base, diff_BGX, "BXY1" = BXY1, diff_BXY)) 
    }
  }
}

# ______________________________________________________________________________________________
#  15. Varying βɢx₁ and βxʏ₁: βɢx₁ = {0.01, 0.09,..., 0.73}, βxʏ₁ = {0.1, ..., 1.45},  
#                             N = 10k, r = 1, q₁ = q₂ = 0.25
# ______________________________________________________________________________________________
sub_sce_varying_par = "BGX1_BXY1"

for (BGX1 in BGX1s){
  for(BXY1 in BXY1s){
    scenarios00 <- rbind(scenarios00, c(main_scenario, sub_sce_varying_par, paste0("BGX1=", BGX1, ",BXY1=",BXY1), 
                                    N_base, r_base, q1_base, q2_base, "BGX1" = BGX1, diff_BGX, "BXY1" = BXY1, diff_BXY)) 
  }
}

#  Sub-scenarios for three varying parameters at a time: 
# ______________________________________________________________________________________________
#  16. Varying N, r, q₁ and q₂: N = {10k, ..., 100k}, r = {0.20, ..., 1, ..., 5}, 
#                               q₁ = q₂ = {0.05, ..., 0.45}, βɢx₁ = 0.5, βxʏ₁ = 0.7
# ______________________________________________________________________________________________
sub_sce_varying_par = "N_r_q1_q2"

for (N in Ns){
  for (r in rs){
    for(q1 in q1s){
      for(q2 in q2s){
      scenarios00 <- rbind(scenarios00, c(main_scenario, sub_sce_varying_par, paste0("N=", N, ",r=", r,",q1=", q1, ",q2=", q2), 
                                      "N" = N, "r" = r, "q1" = q1, "q2" = q2, BGX1_base, diff_BGX, BXY1_base, diff_BXY)) 
      }
    }
  }
}

# ______________________________________________________________________________________________
#  17. Varying N, r, and βɢx₁: N = {10k, ..., 100k}, r = {0.20, ..., 1, ..., 5}, 
#                              βɢx₁ = {0.01, 0.09,..., 0.73}, q₁ = q₂ = 0.25, βxʏ₁ = 0.7
# ______________________________________________________________________________________________
sub_sce_varying_par = "N_r_BGX1"

for (N in Ns){
  for (r in rs){
    for(BGX1 in BGX1s){
      scenarios00 <- rbind(scenarios00, c(main_scenario, sub_sce_varying_par, paste0("N=", N, ",r=", r,",BGX1=", BGX1), 
                                        "N" = N, "r" = r, q1_base, q2_base, "BGX1" = BGX1, diff_BGX, BXY1_base, diff_BXY)) 
    }
  }
}

# ______________________________________________________________________________________________
#  18. Varying N, r, and βxʏ₁: N = {10k, ..., 100k}, r = {0.20, ..., 1, ..., 5}, 
#                              βxʏ₁ = {0.1, ..., 1.45}, q₁ = q₂ = 0.25, βɢx₁ = 0.5
# ______________________________________________________________________________________________
sub_sce_varying_par = "N_r_BXY1"

for (N in Ns){
  for (r in rs){
    for(BXY1 in BXY1s){
      scenarios00 <- rbind(scenarios00, c(main_scenario, sub_sce_varying_par, paste0("N=", N, ",r=", r,",BXY1=", BXY1), 
                                      "N" = N, "r" = r, q1_base, q2_base, BGX1_base, diff_BGX, "BXY1" = BXY1, diff_BXY)) 
    }
  }
}

# ______________________________________________________________________________________________
#  19. Varying N, q₁, q₂, and βɢx₁: N = {10k, ..., 100k}, q₁ = q₂ = {0.05, ..., 0.45},
#                                   βɢx₁ = {0.01, 0.09,..., 0.73}, r = 1, βxʏ₁ = 0.7
# ______________________________________________________________________________________________
sub_sce_varying_par = "N_q1_q2_BGX1"

for (N in Ns){
  for (q1 in q1s){
    for(q2 in q2s){
      for(BGX1 in BGX1s){
        scenarios00 <- rbind(scenarios00, c(main_scenario, sub_sce_varying_par, paste0("N=", N, ",q1=", q1,",q2=", q2, ",BGX1=", BGX1), 
                                        "N" = N, r_base, "q1" = q1, "q2" = q2, "BGX1" = BGX1, diff_BGX, BXY1_base, diff_BXY)) 
      }
    }
  }
}

# ______________________________________________________________________________________________
#  20. Varying N, q₁, q₂, and βxʏ₁: N = {10k, ..., 100k}, q₁ = q₂ = {0.05, ..., 0.45},
#                                   βxʏ₁ = {0.1, ..., 1.45}, r = 1, βɢx₁ = 0.5
# ______________________________________________________________________________________________
sub_sce_varying_par = "N_q1_q2_BXY1"

for (N in Ns){
  for (q1 in q1s){
    for(q2 in q2s){
      for(BXY1 in BXY1s){
        scenarios00 <- rbind(scenarios00, c(main_scenario, sub_sce_varying_par, paste0("N=", N, ",q1=", q1,",q2=", q2, ",BXY1=", BXY1), 
                                        "N" = N, r_base, "q1" = q1, "q2" = q2, BGX1_base, diff_BGX, "BXY1" = BXY1, diff_BXY)) 
      }
    }
  }
}

# ______________________________________________________________________________________________
#  21. Varying N, βɢx₁, and βxʏ₁: N = {10k, ..., 100k}, βɢx₁ = {0.01, 0.09,..., 0.73}
#                                 βxʏ₁ = {0.1, ..., 1.45}, r = 1, q₁ = q₂ = 0.25
# ______________________________________________________________________________________________
sub_sce_varying_par = "N_BGX1_BXY1"

for (N in Ns){
  for (BGX1 in BGX1s){
    for(BXY1 in BXY1s){
      scenarios00 <- rbind(scenarios00, c(main_scenario, sub_sce_varying_par, paste0("N=", N, ",BGX1=", BGX1,",BXY1=", BXY1), 
                                        "N" = N, r_base, q1_base, q2_base, "BGX1" = BGX1, diff_BGX, "BXY1" = BXY1, diff_BXY)) 
    }
  }
}

# ______________________________________________________________________________________________
#  22. Varying r, q₁ and q₂, and βɢx₁: r = {0.20, ..., 1, ..., 5}, q₁ = q₂ = {0.05, ..., 0.45},
#                                      βɢx₁ = {0.01, 0.09,..., 0.73}, N = 10k, βxʏ₁ = 0.7
# ______________________________________________________________________________________________
sub_sce_varying_par = "r_q1_q2_BGX1"

for (r in rs){
  for (q1 in q1s){
    for(q2 in q2s){
      for(BGX1 in BGX1s){
        scenarios00 <- rbind(scenarios00, c(main_scenario, sub_sce_varying_par, paste0("r=", r, ",q1=", q1,",q2=", q2, ",BGX1=", BGX1), 
                                          N_base, "r" = r, "q1" = q1, "q2" = q2, "BGX1" = BGX1, diff_BGX, BXY1_base, diff_BXY)) 
      }
    }
  }
}

# ______________________________________________________________________________________________
#  23. Varying r, q₁ and q₂, and βxʏ₁: r = {0.20, ..., 1, ..., 5}, q₁ = q₂ = {0.05, ..., 0.45},
#                                      βxʏ₁ = {0.1, ..., 1.45}, N = 10k, βɢx₁ = 0.5
# ______________________________________________________________________________________________
sub_sce_varying_par = "r_q1_q2_BXY1"

for (r in rs){
  for (q1 in q1s){
    for(q2 in q2s){
      for(BXY1 in BXY1s){
        scenarios00 <- rbind(scenarios00, c(main_scenario, sub_sce_varying_par, paste0("r=", r, ",q1=", q1,",q2=", q2, ",BXY1=", BXY1), 
                                        N_base, "r" = r, "q1" = q1, "q2" = q2, BGX1_base, diff_BGX, "BXY1" = BXY1, diff_BXY)) 
      }
    }
  }
}

# ______________________________________________________________________________________________
#  24. Varying r, βɢx₁, and βxʏ₁: r = {0.20, ..., 1, ..., 5}, βɢx₁ = {0.01, 0.09,..., 0.73},
#                                 βxʏ₁ = {0.1, ..., 1.45}, N = 10k, q₁ = q₂ = 0.25
# ______________________________________________________________________________________________
sub_sce_varying_par = "r_BGX1_BXY1"

for (r in rs){
  for (BGX1 in BGX1s){
    for(BXY1 in BXY1s){
      scenarios00 <- rbind(scenarios00, c(main_scenario, sub_sce_varying_par, paste0("r=", r, ",BGX1=", BGX1,",BXY1=", BXY1), 
                                      N_base, "r" = r, q1_base, q2_base, "BGX1" = BGX1, diff_BGX, "BXY1" = BXY1, diff_BXY)) 
    }
  }
}

# ___________________________________________________________________________________________________
#  25. Varying q₁ and q₂, βɢx₁, and βxʏ₁: q₁ = q₂ = {0.05, ..., 0.45}, βɢx₁ = {0.01, 0.09,..., 0.73},
#                                         βxʏ₁ = {0.1, ..., 1.45}, N = 10k, r = 1
# ___________________________________________________________________________________________________
sub_sce_varying_par = "q1_q2_BGX1_BXY1"

for (q1 in q1s){
  for (q2 in q2s){
    for(BGX1 in BGX1s){
      for(BXY1 in BXY1s){
        scenarios00 <- rbind(scenarios00, c(main_scenario, sub_sce_varying_par, paste0("q1=", q1, ",q2=", q2,",BGX1=", BGX1, ",BXY1=", BXY1), 
                                        N_base, r_base, "q1" = q1, "q2" = q2, "BGX1" = BGX1, diff_BGX, "BXY1" = BXY1, diff_BXY)) 
      }
    }
  }
}

scenarios00[, 4:11] <- apply(scenarios00[, 4:11] , 2, as.numeric)
save(scenarios00, file = paste0(out_dir, "/scenario.00.Rdata"))


# ==============================================================================================
#  Scenario 01: No genetic effect difference Δβɢx = 0 and causal effect difference Δβxʏ ≠ 0 
# ==============================================================================================
scenarios01 <- vector()
main_scenario =  "scenario.01"
diff_BGX = 0
diff_BXYs = seq(from = 0.01, to = 1, by = 0.19)
# βɢx₁ = {0.01  0.20  0.39  0.58  0.77  0.96}



scenarios01 <- rbind(cbind(scenarios00, ""))
  
# ==============================================================================================
#  Scenario 10: No genetic effect difference Δβɢx ≠ 0 and causal effect difference Δβxʏ = 0 
# ==============================================================================================

# ==============================================================================================
#  Scenario 11: No genetic effect difference Δβɢx ≠ 0 and causal effect difference Δβxʏ ≠ 0 
# ==============================================================================================






## Run simulations: 100 replicates x scenario 
scenarios00[, 4:11] <- apply(scenarios00[, 4:11], 2, as.numeric)

scenario_args = scenarios00[1,]


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
  scenarios00 <- rbind(scenarios00, 
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







## Reproducibility info
session_info()

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
