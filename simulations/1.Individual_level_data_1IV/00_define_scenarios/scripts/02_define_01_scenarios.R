
# ==============================================================================================
#  2. Scenario 01: No genetic effect difference Δβɢx = 0 and causal effect difference Δβxʏ ≠ 0 
# ==============================================================================================
scenarios01 <- vector()
main_scenario =  "scenario.01"
diff_BGX = 0
diff_BXYs = seq(from = 0.01, to = 1, by = 0.19)
# βɢx₁ = {0.01  0.20  0.39  0.58  0.77  0.96}



scenarios01 <- rbind(cbind(scenarios00, ""))