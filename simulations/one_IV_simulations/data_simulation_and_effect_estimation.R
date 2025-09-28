
library(tidyr)
library(dplyr)
library(ggplot2)
library(cowplot)
library(pheatmap)
library(dichromat)

## Initialize
rm(list = ls())
time0 <- Sys.time()
set.seed(09222025)


################################################################################
#                     1. Simulation study using one IV
################################################################################
# ______________________________________________________________________________
#   Data simulation and estimation of genetic and causal effects with one IV, 
#   across varying ΔβXY, N and N1/N2, q, and ΔβGX.
# ______________________________________________________________________________

## Define dirs
input_dir <- paste(getwd(), "simulations", "one_IV_simulations", "inputs", sep = "/")
out_dir <- paste(getwd(), "simulations", "one_IV_simulations", "outputs", sep = "/")
plot_dir <- paste(getwd(), "simulations", "one_IV_simulations", "plots", sep = "/")

dir.create(input_dir, showWarnings = F)
dir.create(out_dir, showWarnings = F)
dir.create(plot_dir, showWarnings = F)

## Color dict TODO


## Helper functions:
source("simulations/helper_functions/simulate_strata.R")
source("simulations/helper_functions/simulate_genotype.R")
source("simulations/helper_functions/HWE_test.R")
source("simulations/helper_functions/simulate_confounder.R")
source("simulations/helper_functions/simulate_error.R")
source("simulations/helper_functions/simulate_exposure.R")
source("simulations/helper_functions/test_genetic_association.R")
source("simulations/helper_functions/simulate_outcome.R")

# source("simulations/helper_functions/test_causal_effect.R")


## Definition of used parameters / variables: TODO

# - N: total sample size
# - Nk: stratum k sample size
# - r = N1/N2: sample size ratio
# - qk: minor allele frequency of variant in stratum k

# - ki ∈ {1,2}: stratum individual i belongs to
# - gi ~ Binom(2, q_ki): genotype for individual i
# - ci ~ Unif(0,1): confounder value for individual i
# - εXi (eXi) ~ N(0,1): error in X for individual i
# - εYi (eYi) ~ N(0,1): error in Y for individual i


# - βGXk (BGXk): genetic effect on X in stratum k
# - ΔβGX (diff_BGX): diff in genetic effect on X between strata = βGX2 - βGX1
# - βGYk (BGYk): genetic effect on Y in stratum k
# - ΔβGY (diff_BGY): diff in genetic effect on Y between strata = βGY2 - βGY1
# - βXYk (BXYk): causal effect in stratum k
# - ΔβXY (diff_BXY): diff in causal effect between strata = βXY2 - βXY1

# - xi = βGX_ki(gi) + ci + εXi: exposure for individual i
# - x̂i (pred_xi) = βGX_ki(gi) + ci: predicted exposure for individual i
# - yi = βXY_ki(x̂i) + ci + εYi: outcome for individual i

# - β̂GXk (est_BGXk): estimated βGXk
# - Δβ̂GX (est_diff_BGX) = β̂GX2 - β̂GX1: estimated ΔβGX
# - β̂GYk (est_BGYk): estimated βGYk   
# - Δβ̂GY (est_diff_BGY): β̂GY2 - β̂GY1: estimated ΔβGY
# - β̂XYk (est_BXYk): estimated βXYk
# - Δβ̂XY (est_diff_BXY): β̂XY2 - β̂XY1: estimated ΔβXY





## Simulation with given parameters

## Test 1:
# N = 1000
# r = 1
# q1 = 0.1
# q2 = 0.1
# BGX1 = 0.4
# diff_BGX = 0
# BXY1 = 0.5
# diff_BXY = 0.5
# scenario = "scenario_1"
# replicate = 1
# sigma_sq = 0.05

simulation_one_IV <- function(N, r, q1, q2, BGX1, diff_BGX, sigma_sq_BGX, BXY1, diff_BXY, sigma_sq_BXY, scenario, replicate, plotting){
  
  ## Output subdirs x scenario x replicate
  out_dir_sc <- paste(out_dir, scenario, sep = "/")
  out_dir_sc_replicate <- paste(out_dir_sc, replicate, sep = "/")
  plot_dir_sc <- paste(plot_dir, scenario, sep = "/")
  plot_dir_sc_replicate <- paste(plot_dir_sc, replicate, sep = "/")
  
  dir.create(out_dir_sc_replicate, recursive = T, showWarnings = F)
  dir.create(plot_dir_sc_replicate, recursive = T, showWarnings = F)
  
  ## df to save results (and inputs)
  out <- data.frame("replicate" = replicate, "N" = N, "r" = r, "q1" = q1, "q2" = q2, 
                    "BGX1" = BGX1, "diff_BGX" = diff_BGX, "sigma_sq_BGX" = sigma_sq_BGX,
                    "BXY1" = BXY1, "diff_BXY" = diff_BXY, "sigma_sq_BXY" = sigma_sq_BXY)
  
  #  Step 1: Simulate strata: K | N, r
  # ___________________________________________________________________________
    K = simulate_strata(N, r)
    
    ## Confirm N1/N2 = r
    if(table(K)["1"] / table(K)["2"] != r){
      message("Created strata don't meet r ratio")
      stop()
    }
    ## Stratum sample sizes
    N1 = (r*N)/(r+1)
    N2 = N/(r+1)
  
    out <- cbind(out, "N1" = N1, "N2" = N2)
    
  #  Step 2: Simulate genotype for each stratum: Gk ~ Binom(2, qk) | Nk, qk
  # ___________________________________________________________________________
    G1 = simulate_genotype(N1, q1)
    G2 = simulate_genotype(N2, q2)
    G = c(G1, G2)
  
    out <- cbind(out, 
             "n_aa_1" = unname(table(G1)["2"]), 
             "n_aA_1" = unname(table(G1)["1"]), 
             "n_AA_1" = unname(table(G1)["0"]),
             "n_aa_2" = unname(table(G2)["2"]), 
             "n_aA_2" = unname(table(G2)["1"]), 
             "n_AA_2" = unname(table(G2)["0"]))
    
    # -----------------------------------------------------------------------
    #  2.1: Chi-squared test for HWE deviations in each stratum and globally
    # -----------------------------------------------------------------------

    for(k in c(1, 2, "")){
     
       Gk = get(paste0("G", k))
       n_aa = table(Gk)["2"]
       n_aA = table(Gk)["1"]
       n_AA = table(Gk)["0"]
       
       HWE_res = HWE_test(n_aa, n_aA, n_AA) 
       q_ob <- signif((n_aA + (2*n_aa))/ (2*length(Gk)), digits = 4)
       
       if(k != ""){
         print(paste0("Observed MAF in stratum ", k, ": ", q_ob))
         print(paste0("P-val for HWE test in stratum ", k, ": ", HWE_res$P.value))
         if(HWE_res$P.value <= 0.05){
           message(paste0("Warning: SNP deviates from HWE in stratum ", k, "."))
         }
         
         out[,paste0("q", k, "_ob")] = q_ob
         out[,paste0("HWE_CHISQ_", k)] = HWE_res$CHISQ
         out[,paste0("HWE_P_", k)] = HWE_res$P.value
         
       } else{
         print(paste0("Observed MAF in whole population: ", q_ob))
         print(paste0("P-val for HWE test in whole population: ", HWE_res$P.value))
         if(HWE_res$P.value <= 0.05){
           message("Warning: SNP deviates from HWE in whole population.")
         }
         
         out[,"q_global_ob"] = q_ob
         out[,"HWE_CHISQ_global"] = HWE_res$CHISQ
         out[,"HWE_P_global"] = HWE_res$P.value
       }
    }
    # -----------------------------------------------------------------------
     
    indiv_data <- data.frame("K" = K, "G" = G)
    
    
  # Step 3: Simulate unknown confounder U ~ Unif(0,1) | N  
  # ___________________________________________________________________
    U = simulate_confounder(N)
    indiv_data <- cbind(indiv_data, "U" = U)
    
  # Step 4: Simulate error terms for X and Y: εX,εY  ~ N(0,1) | N  
  # ___________________________________________________________________
    eX = simulate_error(N)
    eY = simulate_error(N)
     
    indiv_data <- cbind(indiv_data, "eX" = eX, "eY" = eY)
     
  # Step 5: Generate exposure Xₖ = βɢxₖ(G) + U + εx
  # ___________________________________________________________________
     
    # -----------------------------------------------------------------
    #  5.1: Draw βɢxₖ from ~ N(βɢxₖ, σβɢx²)
    # -----------------------------------------------------------------
      drawn_BGX1 <- rnorm(1, BGX1, sigma_sq_BGX)
      drawn_BGX2 <- rnorm(1, BGX1 + diff_BGX, sigma_sq_BGX)
      
      out[, "drawn_BGX1"] = drawn_BGX1
      out[, "drawn_BGX2"] = drawn_BGX2
    # -----------------------------------------------------------------
    
    X = simulate_exposure(indiv_data, cols = c("K", "G", "U", "eX"), drawn_BGX1, drawn_BGX2)
    indiv_data <- cbind(indiv_data, "X" = X)
    
  
  # Step 6: Test G effect on X: get β̂ɢxₖ
  # ___________________________________________________________________
    fit_X = test_genetic_association(indiv_data, cols = c("K", "G", "X"), out_dir_sc_replicate)
 
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
  
    # ----------------------------------------------------------------------
    #  6.1: Confirm SNP is strong IV in both strata
    # ----------------------------------------------------------------------
       fit_X$Fstat = fit_X$`t value` ** 2
    
       if(any(fit_X$Fstat < 10)){
         message("Stopping: weak instrument. Simulate different one.")
         stop()
       }
       out <- cbind(out, "BGX_F_stat_1" = fit_X$Fstat[1],  
                         "BGX_F_stat_2" = fit_X$Fstat[2])
    # ----------------------------------------------------------------------
    
    ## Add pred X: X̂ = β̂ɢxₖ(G)
    indiv_data$pred_X <- apply(indiv_data, 1, function(i){if(i["K"] == 1){i["G"]*hat_BGX1} else{i["G"]*hat_BGX2}})
    
  
    # *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
    #                !!!  Sanity checks  !!!                |
    # *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
    # Sanity check 1:
    # -------------------------------------------------------
    # Estimated vs true effects in each stratum |βGXk - β̂GXk|
    diff_BGX1 = abs(drawn_BGX1 - hat_BGX1)
    diff_BGX2 = abs(drawn_BGX2 - hat_BGX2)
    
    out <- cbind(out, "diff_BGX1" = diff_BGX1, "diff_BGX2" = diff_BGX2)
    
    # Sanity check 2:
    # -------------------------------------------------------
    # Detected genetic effect difference ΔβGX = βGX2 - βGX1?
    drawn_diff_BGX = drawn_BGX2 - drawn_BGX1
    # Δβ̂GX = β̂GX2 - β̂GX1
    hat_diff_BGX = hat_BGX2 - hat_BGX1
    # Difference between true vs estimated delta |ΔβGX - Δβ̂GX|
    diff_diff_BGX = abs(drawn_diff_BGX - hat_diff_BGX)
    
    out <- cbind(out, "drawn_diff_BGX" = drawn_diff_BGX, 
                      "hat_diff_BGX" = hat_diff_BGX, 
                      "diff_diff_BGX" = diff_diff_BGX)
    
    # Sanity check 3:
    # -------------------------------------------------------
    # Significant G x K interaction on X?
    
    # Z-stat = (β̂GX2 - β̂GX1) / √se(β̂GX2)² + se(β̂GX1)² 
    Z_diff_BGX = (hat_BGX2 - hat_BGX1) / sqrt(se_hat_BGX2^2 + se_hat_BGX1^2)
    # P val: Z-stat ~ N(0,1) | Ho
    p_Z_diff_BGX = 2*pnorm(abs(Z_diff_BGX), 0, 1, lower.tail = F)
    
    out <- cbind(out, "Z_diff_BGX" = Z_diff_BGX, 
                      "p_Z_diff_BGX" = p_Z_diff_BGX)
    
    # *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*

  # Step 7: Generate outcome Yₖ = βxʏₖ(βɢxₖ)(G) + U + εʏ
  # ___________________________________________________________________
    
    # -----------------------------------------------------------------
    #  7.1: Draw βxʏₖ from ~ N(βxʏₖ, σβxʏ²)
    # -----------------------------------------------------------------
    drawn_BXY1 <- rnorm(1, BXY1, sigma_sq_BXY)
    drawn_BXY2 <- rnorm(1, BXY1 + diff_BXY, sigma_sq_BXY)
    
    out[, "drawn_BXY1"] = drawn_BXY1
    out[, "drawn_BXY2"] = drawn_BXY2
    # -----------------------------------------------------------------
    
    Y = simulate_outcome(indiv_data, cols = c("K", "G", "U", "eY"), 
                         drawn_BGX1, drawn_BGX2, drawn_BXY1, drawn_BXY2)
    indiv_data <- cbind(indiv_data, "Y" = Y)
  
  # Step 8: Test G association with Y: get β̂ɢʏₖ
  # ___________________________________________________________________
    fit_Y = test_genetic_association(indiv_data, cols = c("K", "G", "Y"), out_dir_sc_replicate)
    
    hat_BGY1 = fit_Y$Estimate[1]
    hat_BGY2 = fit_Y$Estimate[2]
    se_hat_BGY1 = fit_Y$`Std. Error`[1]
    se_hat_BGY2 = fit_Y$`Std. Error`[2]
    BGY_t_stat_1 = fit_Y$`t value`[1]
    BGY_t_stat_2 = fit_Y$`t value`[2]
    BGY_P_1 = fit_Y$`Pr(>|t|)`[1]
    BGY_P_2 = fit_Y$`Pr(>|t|)`[2]
    
    out <- cbind(out, "hat_BGY1" = hat_BGY1, "hat_BGY2" = hat_BGY1, 
                 "se_hat_BGY1" = se_hat_BGY1, "se_hat_BGY2" = se_hat_BGY2, 
                 "BGY_t_stat_1" = BGY_t_stat_1, "BGY_t_stat_2" = BGY_t_stat_2, 
                 "BGY_P_1" = BGY_P_1, "BGY_P_2" = BGY_P_2)
    
    ## Add Y pred by G: Ŷ =  β̂ɢʏₖ(G) 
    indiv_data$pred_Y <- apply(indiv_data, 1, function(i){if(i["K"] == 1){i["G"]*hat_BGY1} else{i["G"]*hat_BGY2}})
    

    # *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
    #                !!!  Sanity checks  !!!                |
    # *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
    # Sanity check 1:
    # -------------------------------------------------------
    if(plotting == T){
      # Plot: X ~ G, U, εx
      #       X̂ ~ G, U, εx
      #       Y ~ X, X̂, G, U, εx, εʏ
      #       Ŷ ~ X, X̂, G, U, εx, εʏ
      plot_pheno_vs_predictors(indiv_data, pheno_var = "X", predictors = c("G", "U", "eX"), out, plot_dir_sc_replicate)
      # Plot correlation between exposure/outcome and predictors
      cor_pheno_predictors(indiv_data, variables = c("X", "pred_X", "G", "U", "eX"), plot_dir_sc_replicate)
    }
  return(out)
}



## Plot phenotype ~ predictor with estimated effects in both strata
plot_pheno_vs_predictor <- function(indiv_data, pheno, predictor, out){
  
  labels = c("X" = expression(X),
             "pred_X" = expression(hat(X) == hat(beta[GX])*G),
             "Y" = expression(Y),
             "pred_Y" = expression(hat(Y) == hat(beta[GY])*G),
             "G" = expression(G), 
             "U" = expression(U), 
             "eX" = expression(epsilon[X]), 
             "eY" = expression(epsilon[Y]))
  
  ypos1 = min(indiv_data[, pheno])+0.4
  ypos2 = min(indiv_data[, pheno])+1.1
  
  if(length(grep("pred", pheno))>0){
    ypos1 = min(indiv_data[, pheno])+0.05
    ypos2 = min(indiv_data[, pheno])+0.14
    }
  
  p = ggplot(indiv_data, aes(x = get(predictor), y = get(pheno))) + 
    facet_grid(cols = vars(K), 
               labeller = labeller(K = c("1" = "Stratum 1", "2" = "Stratum 2"))) +
    geom_point(color='gray', alpha = 0.3, size = 1) +
    theme_classic() +
    labs(x = labels[predictor], y = labels[pheno]) +
    theme(axis.title = element_text(size = 8),
          axis.text = element_text(size = 7),
          strip.text.x = element_text(size = 8), 
          axis.line = element_line(linewidth = 0.4),
          axis.ticks = element_line(linewidth = 0.3),
          strip.text = element_text(size = 8.1, face = 2),
          strip.background = element_rect(fill="white", color = "gray", linewidth = 0.7),
          panel.spacing = unit(0.4, "lines"))
  
  if(pheno == "X" & predictor == "G"){
    p <- p + scale_x_continuous(breaks = c(0,1,2), labels = c(0,1,2)) +
      ## Add β̂ɢxₖ
      geom_abline(data = data.frame("hat_beta" = c(out$hat_BGX1, out$hat_BGX2), "K" = c(1,2)), 
                  aes(slope = hat_beta, intercept = 0), linewidth = 0.65, alpha = 0.5, color = "orangered2") +
      geom_text(data = data.frame("hat_beta" = c(out$hat_BGX1, out$hat_BGX2), "K" = c(1,2)), 
                aes(label = paste("hat(beta)[GX] ==", signif(hat_beta, digits = 2)), x = 1.75, y = min(indiv_data[, pheno])+0.4), color = "orangered2", size = 3, parse = T) +
      ## Add true (drawn) βɢxₖ
      geom_abline(data = data.frame("drawn_beta" = c(out$drawn_BGX1, out$drawn_BGX2), "K" = c(1,2)), 
                  aes(slope = drawn_beta, intercept = 0), linewidth = 0.65, alpha = 0.5, color = "purple3") +
      geom_text(data = data.frame("drawn_beta" = c(out$drawn_BGX1, out$drawn_BGX2), "K" = c(1,2)), 
                aes(label = paste("beta[GX] ==", signif(drawn_beta, digits = 2)), x = 1.75, y = min(indiv_data[, pheno])+1.2), color = "purple3", size = 3, parse = T) 
  }
  
  else if(pheno == "Y" & predictor == "G"){
    p <- p + scale_x_continuous(breaks = c(0,1,2), labels = c(0,1,2)) +
      ## Add β̂ɢʏₖ
      geom_abline(data = data.frame("hat_beta" = c(out$hat_BGY1, out$hat_BGY2), "K" = c(1,2)), 
                  aes(slope = hat_beta, intercept = 0), linewidth = 0.65, alpha = 0.5, color = "orangered2") +
      geom_text(data = data.frame("hat_beta" = c(out$hat_BGY1, out$hat_BGY2), "K" = c(1,2)), 
                aes(label = paste("hat(beta)[GY] ==", signif(hat_beta, digits = 2)), x = 1.75, y = min(indiv_data[, pheno])+0.4), color = "orangered2", size = 3, parse = T) +
      ## Add true (drawn) βɢʏₖ = βɢxₖ x βxʏₖ 
      geom_abline(data = data.frame("drawn_beta" = c(out$drawn_BGX1*out$drawn_BXY1, out$drawn_BGX2* out$drawn_BXY2), "K" = c(1,2)), 
                  aes(slope = drawn_beta, intercept = 0), linewidth = 0.65, alpha = 0.5, color = "purple3") +
      geom_text(data = data.frame("drawn_beta" = c(out$drawn_BGX1*out$drawn_BXY1, out$drawn_BGX2* out$drawn_BXY2), "K" = c(1,2)), 
                aes(label = paste("beta[GY] ==", signif(drawn_beta, digits = 2)), x = 1.75, y = min(indiv_data[, pheno])+1.2), color = "purple3", size = 3, parse = T) 
  }
  
  else if(pheno == "pred_X" & predictor == "G"){
    p <- p + scale_x_continuous(breaks = c(0,1,2), labels = c(0,1,2)) +
      ## Add β̂ɢxₖ
      geom_abline(data = data.frame("hat_beta" = c(out$hat_BGX1, out$hat_BGX2), "K" = c(1,2)), 
                  aes(slope = hat_beta, intercept = 0), linewidth = 0.65, alpha = 0.5, color = "orangered2") +
      geom_text(data = data.frame("hat_beta" = c(out$hat_BGX1, out$hat_BGX2), "K" = c(1,2)), 
                aes(label = paste("hat(beta)[GX] ==", signif(hat_beta, digits = 2)), x = 1.75, y = min(indiv_data[, pheno])+0.05), color = "orangered2", size = 3, parse = T) +
      ## Add true (drawn) βɢxₖ
      geom_abline(data = data.frame("drawn_beta" = c(out$drawn_BGX1, out$drawn_BGX2), "K" = c(1,2)), 
                  aes(slope = drawn_beta, intercept = 0), linewidth = 0.65, alpha = 0.5, color = "purple3") +
      geom_text(data = data.frame("drawn_beta" = c(out$drawn_BGX1, out$drawn_BGX2), "K" = c(1,2)), 
                aes(label = paste("beta[GX] ==", signif(drawn_beta, digits = 2)), x = 1.75, y = min(indiv_data[, pheno])+0.14), color = "purple3", size = 3, parse = T) 
  }
  
  else if(pheno == "X" & predictor == "U"){
    p <- p + geom_abline(slope = 1, intercept = 0, linewidth = 0.65, alpha = 0.7, color = "gray30") +
             geom_text(x = 0.75, y = 0.5, label = paste("beta[UX] == 1"), color = 'gray30', size = 3, parse = T)
  }
  
  else if(pheno == "pred_X" & predictor == "U" ){
    p <- p + geom_abline(slope = 0, intercept = 0, linewidth = 0.65, alpha = 0.7, color = "gray30") +
             geom_text(x = 0.75, y = 0.5, label = paste("beta[U][hat(X)] == 0"), color = 'gray30', size = 3, parse = T)
  }
  
  
  else if(pheno == "X" & predictor == "eX"){
    p <- p + geom_abline(slope = 1, intercept = 0, linewidth = 0.65, alpha = 0.7, color = "gray30") +
             geom_text(x = 2, y = 0.5, label = paste("beta[epsilon[X]] == 1"), color = 'gray30', size = 3, parse = T)
  }
  
  else if(pheno == "pred_X" & predictor == "eX"){
    p <- p + geom_abline(slope = 0, intercept = 0, linewidth = 0.65, alpha = 0.7, color = "gray30") +
             geom_text(x = 2, y = 0.5, label = paste("beta[epsilon[X]] == 0"), color = 'gray30', size = 3, parse = T)
  }
  
  else if(pheno == "Y" & predictor == "eX"){
    p <- p + geom_abline(slope = 0, intercept = 0, linewidth = 0.65, alpha = 0.7, color = "gray30") +
      geom_text(x = 2, y = 0.5, label = paste("beta[epsilon[X]] == 0"), color = 'gray30', size = 3, parse = T)
  }
  
  else if(pheno == "X" & predictor == "U"){
  p <- p + geom_abline(slope = 1, intercept = 0, linewidth = 0.65, alpha = 0.7, color = "gray30") +
           geom_text(x = 0.75, y = 0.5, label = paste("beta[UX] == 1"), color = 'gray30', size = 3, parse = T)
  }
  
  else if(pheno == "pred_X" & predictor == "U" ){
    p <- p + geom_abline(slope = 0, intercept = 0, linewidth = 0.65, alpha = 0.7, color = "gray30") +
             geom_text(x = 0.75, y = 0.5, label = paste("beta[U][hat(X)] == 0"), color = 'gray30', size = 3, parse = T)
  }
  
  else if(pheno == "Y" & predictor == "U"){
    p <- p + geom_abline(slope = 1, intercept = 0, linewidth = 0.65, alpha = 0.7, color = "gray30") +
             geom_text(x = 0.75, y = 0.5, label = paste("beta[UY] == 1"), color = 'gray30', size = 3, parse = T)
  }
  
  else if(pheno == "pred_Y" & predictor == "U" ){
    p <- p + geom_abline(slope = 0, intercept = 0, linewidth = 0.65, alpha = 0.7, color = "gray30") +
             geom_text(x = 0.75, y = 0.5, label = paste("beta[U][hat(Y)] == 0"), color = 'gray30', size = 3, parse = T)
  }
  
  else if(pheno == "Y" & predictor == "eY"){
    p <- p + geom_abline(slope = 1, intercept = 0, linewidth = 0.65, alpha = 0.7, color = "gray30") +
      geom_text(x = 2, y = 0.5, label = paste("beta[epsilon[Y]] == 0"), color = 'gray30', size = 3, parse = T)
  }
  
  else if(pheno == "Y" & predictor == "X"){
    ## Add β̂xʏₖ
    p <- p + geom_abline(data = data.frame("hat_beta" = c(out$hat_BXY1, out$hat_BXY2), "K" = c(1,2)), 
                  aes(slope = hat_beta, intercept = 0), linewidth = 0.65, alpha = 0.5, color = "orangered2") +
      geom_text(data = data.frame("hat_beta" = c(out$hat_BXY1, out$hat_BXY2), "K" = c(1,2)), 
                aes(label = paste("hat(beta)[XY] ==", signif(hat_beta, digits = 2)), x = 1.75, y = min(indiv_data[, pheno])+0.05), color = "orangered2", size = 3, parse = T) +
      ## Add true (drawn) βxʏₖ
      geom_abline(data = data.frame("drawn_beta" = c(out$drawn_BXY1, out$drawn_BXY2), "K" = c(1,2)), 
                  aes(slope = drawn_beta, intercept = 0), linewidth = 0.65, alpha = 0.5, color = "purple3") +
      geom_text(data = data.frame("drawn_beta" = c(out$drawn_BXY1, out$drawn_BXY2), "K" = c(1,2)), 
                aes(label = paste("beta[XY] ==", signif(drawn_beta, digits = 2)), x = 1.75, y = min(indiv_data[, pheno])+0.14), color = "purple3", size = 3, parse = T) 
    
    }
  
  
  
  

  return(p)
}


plot_pheno_vs_predictors <- function(indiv_data, predictors, out, plot_dir_sc_replicate){
  
  phenotypes <- c("X", "pred_X", "Y", "pred_Y")
  predictors <- list("X" = c("G", "U", "eX"), 
                     "pred_X" = c("G", "U", "eX"),
                     "Y" = c("X", "pred_X", "G", "U", "eX", "eY"),
                     "pred_Y" = c("X", "pred_X", "G", "U", "eX", "eY"))
  plots <- list()
  i = 1
  
  for(pheno in phenotypes){
    for(predictor in predictors[[pheno]]){
      plots[[i]] <- plot_pheno_vs_predictor(indiv_data, pheno, predictor, out)
      i = i + 1
    }
  }
  
  if(length(plots) == 6){
    plot_grid(plotlist = plots, ncol = 3, align = "vh")
    ggsave(filename = paste0(plot_dir_sc_replicate, "/", pheno_var, "_on_G_U_e", pheno_var, ".pdf"), width = 13, height = 5)
  }
}
  

# Plot corr between X/Y and their predictors
cor_pheno_predictors <- function(indiv_data, variables, plot_dir_sc_replicate){
  
  h <- list()
  for(k in 1:2){
    
    data = subset(indiv_data, K == k)
    
    h[[k]] = pheatmap(cor(data[, variables]), 
             color = colorRampPalette(c('#FFF0F5', '#EE0000'))(50),
             display_numbers = T, 
             cluster_rows = F, 
             cluster_cols = F, 
             main = paste("Stratum", k), 
             border_color = "gray20", 
             number_color = "black", 
             fontsize = 10,
             labels_row = c("X", expression(X == hat(beta[GX])*G), "G", "U", expression(epsilon[X])),
             labels_col = c("X", expression(X == hat(beta[GX])*G), "G", "U", expression(epsilon[X])), 
             angle_col = 90)$gtable
  }
  
  plot_grid(plotlist = h, ncol = 2, align = "h")
  ggsave(file = paste0(plot_dir_sc_replicate, "/Corr_predictors.pdf"), width = 8, height = 3.7)
}






# ______________________________________________________________________________
#        A) Simulations to evaluate estimated genetic effects of one SNP 
#                          on exposure across strata
# ______________________________________________________________________________
# Assess impact of varying N, r, q1, q2, Δq, βGX1, βGX2, and ΔβGX --> 
#                                                         β̂GX1, β̂GX2, and Δβ̂GX 
# ______________________________________________________________________________

scenarios <- vector()

# Base model: 
N = 10000
r = 1
q1 = q2 = 0.25
BGX1 = 0.7
sigma_sq_BGX = sigma_sq_BXY = 0.05
BXY1 = 0.8
diff_BXY = 0

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                Scenario 1: No genetic effect difference ΔβGX = 0
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
diff_BGX = 0

# ------------------------------------------------------------------------------
# 1.0 Base model:  ΔβGX = 0, N = 10,000, r = 1, q1 = q2 = 0.25, βGX1 = 0.7
# ------------------------------------------------------------------------------

scenario_A1.0_res = vector()

for(i in 1:100){
  plotting_option = if_else(i <= 10, T, F)
  res = simulation_one_IV(N, r, q1, q2, BGX1, diff_BGX, BXY1 = BXY1, diff_BXY = diff_BXY, sigma_sq_BGX =  sigma_sq_BGX, 
                          sigma_sq_BXY = sigma_sq_BXY, 
                          scenario = "scenario_A1.0", replicate = i, plotting = plotting_option)
  scenario_A1.0_res <- rbind(scenario_A1.0_res, res)
  
}
save(scenario_A1.0_res, file = paste0(out_dir))
## Mean of drawn BGX1 across replicates
mean(scenario_A1.0_res$drawn_BGX1)
mean(scenario_A1.0_res$drawn_BGX2)

mean(scenario_A1.0_res$hat_BGX1)
mean(scenario_A1.0_res$hat_BGX2)

## Mean F-stat of SNP 
mean(scenario_A1.0_res$BGX_F_stat_1)
mean(scenario_A1.0_res$BGX_F_stat_2)

## FPR in HWE
mean(scenario_A1.0_res$HWE_P_1 < 0.05)
mean(scenario_A1.0_res$HWE_P_2 < 0.05)
mean(scenario_A1.0_res$HWE_P_global < 0.05)

## FNR for estimated effects in each stratum
mean(scenario_A1.0_res$BGX_P_1 >= 0.05)
mean(scenario_A1.0_res$BGX_P_2 >= 0.05)

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

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -





# ------------------------------------------------------------------------------




# ------------------------------------------------------------------------------
# Scenario 1: Genetic effect difference ΔβGX > 0
diff_BGX = seq(from = 0.05, to = 1.0, by = 0.1)
simulation_one_IV(N, r, q1, q2, BGX1, diff_BGX, BXY1 = NULL, diff_BXY = NULL, 
                  scenario = "GWAS_scenario_0", replicate = 1)
# ------------------------------------------------------------------------------



# Scenario 1: Base model and genetic effect difference ΔβGX = 0.1:






# N = 10000, r = 1, q = 0.1, βGX1 = 0.4, 
#             ΔβGX = 0, βXY1 = 0.4, ΔβXY = 0.5
# --------------------------------------------------------------


## Simulate 1 replica x scenario





    