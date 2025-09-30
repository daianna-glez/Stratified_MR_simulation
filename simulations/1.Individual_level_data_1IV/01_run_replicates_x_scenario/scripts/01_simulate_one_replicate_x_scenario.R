
#' @title Individual level data simulation and 1 IV-based MR for a replicate of a scenario. 
#' @description This function runs a simulation for one replicate of a scenario: simulates genotype (G) for 1 IV, an unknown confounder (U), \\
#' exposure (X), outcome (Y), and error terms (eX and eY) for N individuals across the two strata of K. \\
#' It tests for HWE, G-X and G-Y associations within each stratum, and estimates per-stratum causal effects.
#' @param scenario string for main scenario name (one of 00, 01, 10, 11).
#' @param sub_scenario1 string for sub scenario 1 name, according to varying parameter(s).
#' @param sub_scenario2 string for sub scenario 2 name, according to the value(s) the varying parameter(s) take. 
#' @param N int for total sample size.
#' @param r double for ratio of stratum sample sizes.
#' @param q1 double for MAF of IV in stratum 1.
#' @param q2 double for MAF of IV in stratum 2.
#' @param BGX1 double for true effect of IV on X in stratum 1.
#' @param diff_BGX double for true difference in the effect of IV on X between stratum 2 and 1.
#' @param BXY1 double for true casual effect in stratum 1.
#' @param diff_BXY double for true difference in the causal effect between stratum 2 and 1.
#' @param replicate int for replicate number.
#' @return A data.frame with scenario input arguments and output metrics for scenario replicate.  
#' @rdname simulation_indiv_data_1IV
#' @export


simulation_indiv_data_1IV <- function(scenario, sub_scenario1, sub_scenario2, N, r, q1, q2, BGX1, diff_BGX, BXY1, diff_BXY, replicate){
  
  ## Output subdirs x scenario x replicate
  out_dir_sc <- paste(out_dir, scenario, sub_scenario1, sub_scenario2, sep = "/")
  out_dir_sc_replicate <- paste(out_dir_sc, replicate, sep = "/")
  plot_dir_sc <- paste(plot_dir, scenario, sub_scenario1, sub_scenario2, sep = "/")
  plot_dir_sc_replicate <- paste(plot_dir_sc, replicate, sep = "/")
  
  dir.create(out_dir_sc_replicate, recursive = T, showWarnings = F)
  
  ## df to save results (and inputs)
  out <- data.frame("scenario" = scenario, "sub_scenario1" = sub_scenario1, "sub_scenario2" = sub_scenario2, 
                    "replicate" = replicate, "N" = N, "r" = r, "q1" = q1, "q2" = q2, 
                    "BGX1" = BGX1, "diff_BGX" = diff_BGX, "BXY1" = BXY1, "diff_BXY" = diff_BXY)
  
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
  
  #  Step 2: Simulate genotype for each stratum: Gₖ ~ Binom(2, qₖ) | Nₖ, qₖ
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
  
  # Step 4: Simulate error terms for X and Y:  εx,εʏ ~ N(0,1) | N  
  # ___________________________________________________________________
  eX = simulate_error(N)
  eY = simulate_error(N)
  
  indiv_data <- cbind(indiv_data, "eX" = eX, "eY" = eY)
  
  # Step 5: Generate exposure Xₖ = βɢxₖ(G) + U + εx
  # ___________________________________________________________________
  ## Define βɢx₂
  BGX2 = BGX1 + diff_BGX
  
  X = simulate_exposure(indiv_data, cols = c("K", "G", "U", "eX"), BGX1, BGX2)
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
  # True vs estimated effect in each stratum |βɢxₖ- β̂ɢxₖ|
  diff_BGX1 = abs(BGX1 - hat_BGX1)
  diff_BGX2 = abs(BGX2 - hat_BGX2)
  
  out <- cbind(out, "diff_BGX1" = diff_BGX1, "diff_BGX2" = diff_BGX2)
  
  # Sanity check 2:
  # -------------------------------------------------------
  # Detected genetic effect difference Δβɢx = βɢx₂ - βɢx₁?
  # Δβ̂ɢx = β̂ɢx₂ - β̂ɢx₁
  hat_diff_BGX = hat_BGX2 - hat_BGX1
  # Difference between true vs estimated delta |Δβɢx - Δβ̂ɢx|
  diff_diff_BGX = abs(diff_BGX - hat_diff_BGX)
  
  out <- cbind(out, "hat_diff_BGX" = hat_diff_BGX, 
               "diff_diff_BGX" = diff_diff_BGX)
  
  # Sanity check 3:
  # -------------------------------------------------------
  # Significant G x K interaction on X?
  
  # Z-stat = (β̂ɢx₂ - β̂ɢx₁) / √se(β̂ɢx₂)² + se(β̂ɢx₁)² 
  Z_diff_BGX = (hat_BGX2 - hat_BGX1) / sqrt(se_hat_BGX2^2 + se_hat_BGX1^2)
  # P val: Z-stat ~ N(0,1) | Ho
  p_Z_diff_BGX = 2*pnorm(abs(Z_diff_BGX), 0, 1, lower.tail = F)
  
  out <- cbind(out, "Z_diff_BGX" = Z_diff_BGX, 
               "p_Z_diff_BGX" = p_Z_diff_BGX)
  
  # *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
  
  # Step 7: Generate outcome Yₖ = βxʏₖ(βɢxₖ)(G) + U + εʏ
  # ___________________________________________________________________
  ## Define βxʏ₂
  BXY2 = BXY1 + diff_BXY
  
  Y = simulate_outcome(indiv_data, cols = c("K", "G", "U", "eY"), BGX1, BGX2, BXY1, BXY2)
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
  
  out <- cbind(out, "hat_BGY1" = hat_BGY1, "hat_BGY2" = hat_BGY2, 
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
  # Estimated vs true effects in each stratum |βɢʏₖ- β̂ɢʏₖ|
  BGY1 = BGX1 * BXY1
  BGY2 = BGX2 * BXY2
  diff_BGY1 = abs(BGY1 - hat_BGY1)
  diff_BGY2 = abs(BGY2 - hat_BGY2)
  
  out <- cbind(out, "BGY1" = BGY1, "BGY2" = BGY2, 
               "diff_BGY1" = diff_BGY1, "diff_BGY2" = diff_BGY2)
  
  # Sanity check 2:
  # -------------------------------------------------------
  # True genetic effect difference Δβɢʏ = βɢʏ₂ - βɢʏ₁
  diff_BGY = BGY2 - BGY1
  # Estimated genetic effect difference Δβ̂ɢʏ = β̂ɢʏ₂ - β̂ɢʏ₁
  hat_diff_BGY = hat_BGY2 - hat_BGY1
  # Difference between true vs estimated delta |Δβɢʏ - Δβ̂ɢʏ|
  diff_diff_BGY = abs(diff_BGY - hat_diff_BGY)
  
  out <- cbind(out, "diff_BGY" = diff_BGY, 
               "hat_diff_BGY" = hat_diff_BGY, 
               "diff_diff_BGY" = diff_diff_BGY)
  
  # Sanity check 3:
  # -------------------------------------------------------
  # Significant G x K interaction on Y?
  
  # Z-stat = (β̂ɢʏ₂ - β̂ɢʏ₁) / √se(β̂ɢʏ₂)² + se(β̂ɢʏ₁)² 
  Z_diff_BGY = (hat_BGY2 - hat_BGY1) / sqrt(se_hat_BGY2^2 + se_hat_BGY1^2)
  # P val: Z-stat ~ N(0,1) | Ho
  p_Z_diff_BGY = 2*pnorm(abs(Z_diff_BGY), 0, 1, lower.tail = F)
  
  out <- cbind(out, "Z_diff_BGY" = Z_diff_BGY, 
               "p_Z_diff_BGY" = p_Z_diff_BGY)
  
  # *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
  
  
  # Step 9: Calculate causal effects: get β̂xʏₖ
  # ___________________________________________________________________
  
  hat_BXY1 <- hat_BGY1 / hat_BGX1
  hat_BXY2 <- hat_BGY2 / hat_BGX2
  
  out <- cbind(out, "hat_BXY1" = hat_BXY1, "hat_BXY2" = hat_BXY2)
  
  # *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
  #                   Metrics of interest                 |
  # *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
  # Estimated vs true effects in each stratum |βxʏₖ- β̂xʏₖ|
  # -------------------------------------------------------
  diff_BXY1 = abs(BXY1 - hat_BXY1)
  diff_BXY2 = abs(BXY2 - hat_BXY2)
  
  out <- cbind(out, "diff_BXY1" = diff_BXY1, "diff_BXY2" = diff_BXY2)
  
  # # Difference between true vs estimated delta |Δβxʏ - Δβ̂xʏ|
  # ----------------------------------------------------------
  # True causal effect difference: Δβxʏ = βxʏ₂ - βxʏ₁
  # Estimated causal effect difference: Δβ̂xʏ = β̂xʏ₂ - β̂xʏ₁
  hat_diff_BXY = hat_BXY2 - hat_BXY1
  # Difference 
  diff_diff_BXY = abs(diff_BXY - hat_diff_BXY)
  
  out <- cbind(out, "hat_diff_BXY" = hat_diff_BXY, 
               "diff_diff_BXY" = diff_diff_BXY)
  
  # Significant X x K interaction on Y? se?????
  # -------------------------------------------------------
  # Z-stat = (β̂xʏ₂ - β̂xʏ₁) / √se(β̂xʏ₂)² + se(β̂xʏ₁)² 
  # Z_diff_BXY = (hat_BXY2 - hat_BXY1) / sqrt(se_hat_BXY2^2 + se_hat_BXY1^2)
  # # P val: Z-stat ~ N(0,1) | Ho
  # p_Z_diff_BXY = 2*pnorm(abs(Z_diff_BXY), 0, 1, lower.tail = F)
  # 
  # out <- cbind(out, "Z_diff_BXY" = Z_diff_BXY, 
  #              "p_Z_diff_BXY" = p_Z_diff_BXY)
  
  if(replicate <= 3){
    
    dir.create(plot_dir_sc_replicate, recursive = T, showWarnings = F)
    
    # Plot: X ~ G, U, εx
    #       X̂ ~ G, U, εx
    #       Y ~ X, X̂, G, U, εx, εʏ
    #       Ŷ ~ X, X̂, G, U, εx, εʏ
    plot_pheno_vs_predictors(indiv_data, out, plot_dir_sc_replicate)
    # Plot correlation between exposure/outcome and predictors
    cor_pheno_predictors(indiv_data, plot_dir_sc_replicate)
  }
  
  # *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
  
  ## Save indiv data generated
  save(indiv_data, file = paste0(out_dir_sc_replicate, "/indiv_data.Rdata"))
  
  return(out)
}



## Plot phenotype ~ predictor with estimated and true effects in both strata
plot_pheno_vs_predictor <- function(indiv_data, pheno, predictor, out){
  
  labels = c("X" = expression(X),
             "pred_X" = expression(hat(X) == hat(beta[GX])*G),
             "Y" = expression(Y),
             "pred_Y" = expression(hat(Y) == hat(beta[GY])*G),
             "G" = expression(G), 
             "U" = expression(U), 
             "eX" = expression(epsilon[X]), 
             "eY" = expression(epsilon[Y]))
  
  
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
  
  if(predictor == "G"){
    p <- p + scale_x_continuous(breaks = c(0,1,2), labels = c(0,1,2)) 
    
    if(pheno == "X" | pheno == "pred_X"){
      betas = data.frame("hat_beta" = c(out$hat_BGX1, out$hat_BGX2), "K" = c(1,2))
      betas$beta = c(out$BGX1, out$BGX2)
      
      ## Add β̂ɢxₖ
      p <- p + geom_abline(data = betas, aes(slope = hat_beta, intercept = 0), linewidth = 0.65, alpha = 0.5, color = "orangered2") +
        geom_text(data = betas, aes(label = paste("hat(beta)[GX] ==", signif(hat_beta, digits = 2)), 
                                    x = 1.75, y = min(indiv_data[, pheno]) + ((max(indiv_data[, pheno]) - min(indiv_data[, pheno]))*0.1)), color = "orangered2", size = 3, parse = T) +
        ## Add true βɢxₖ
        geom_abline(data = betas, aes(slope = beta, intercept = 0), linewidth = 0.65, alpha = 0.5, color = "purple3") +
        geom_text(data = betas, aes(label = paste("beta[GX] ==", signif(beta, digits = 2)), 
                                    x = 1.75, y = min(indiv_data[, pheno]) + ((max(indiv_data[, pheno]) - min(indiv_data[, pheno]))*0.2)), color = "purple3", size = 3, parse = T) 
    }
    
    else if(pheno == "Y" | pheno == "pred_Y"){
      betas <- data.frame("hat_beta" = c(out$hat_BGY1, out$hat_BGY2), "K" = c(1,2))
      betas$beta <- c(out$BGY1, out$BGY2)
      
      ## Add estimated effects β̂ɢʏₖ
      p <- p + geom_abline(data = betas, aes(slope = hat_beta, intercept = 0), linewidth = 0.65, alpha = 0.5, color = "orangered2") +
        geom_text(data = betas, aes(label = paste("hat(beta)[GY] ==", signif(hat_beta, digits = 2)), 
                                    x = 1.75, y = min(indiv_data[, pheno]) + ((max(indiv_data[, pheno]) - min(indiv_data[, pheno]))*0.10)), color = "orangered2", size = 3, parse = T) +
        ## Add true βɢʏₖ = βɢxₖ x βxʏₖ 
        geom_abline(data = betas, aes(slope = beta, intercept = 0), linewidth = 0.65, alpha = 0.5, color = "purple3") +
        geom_text(data = betas, aes(label = paste("beta[GY] ==", signif(beta, digits = 2)), 
                                    x = 1.75, y = min(indiv_data[, pheno]) + ((max(indiv_data[, pheno]) - min(indiv_data[, pheno]))*0.2)), color = "purple3", size = 3, parse = T) 
    }
  }
  
  else if(predictor == "U"){
    
    if(pheno == "X"){
      slope = 1
      label = paste("beta[UX] == 1")
    }
    else if(pheno == "pred_X"){
      slope = 0
      label = paste("beta[U][hat(X)] == 0")
    }
    else if(pheno == "Y"){
      slope = 1
      label = paste("beta[UY] == 1")
    }
    else if(pheno == "pred_Y"){
      slope = 0
      label = paste("beta[U][hat(Y)] == 0")
    }
    
    p <- p + geom_abline(slope = slope, intercept = 0, linewidth = 0.65, alpha = 0.7, color = "gray30") +
      geom_text(x = 0.75, y = min(indiv_data[, pheno]) + ((max(indiv_data[, pheno]) - min(indiv_data[, pheno]))/2.5), label = label, color = 'gray30', size = 3, parse = T)
  }
  
  else if(predictor == "eX"){
    
    if(pheno == "X"){
      slope = 1
      label = paste("beta[epsilon[X]] == 1")
    }
    else if(pheno == "pred_X"){
      slope = 0
      label = paste("beta[epsilon[X]] == 0")
    }
    else if(pheno == "Y"){
      slope = 0
      label = paste("beta[epsilon[X]] == 0")
    }
    else if(pheno == "pred_Y"){
      slope = 0
      label = paste("beta[epsilon[X]] == 0")
    }
    
    p <- p + geom_abline(slope = slope, intercept = 0, linewidth = 0.65, alpha = 0.7, color = "gray30") +
      geom_text(x = 2, y = min(indiv_data[, pheno]) + ((max(indiv_data[, pheno]) - min(indiv_data[, pheno]))/2.5), label = label, color = 'gray30', size = 3, parse = T)
  }
  
  else if(predictor == "eY"){
    
    if(pheno == "X"){
      slope = 0
      label = paste("beta[epsilon[Y]] == 0")
    }
    else if(pheno == "pred_X"){
      slope = 0
      label = paste("beta[epsilon[Y]] == 0")
    }
    else if(pheno == "Y"){
      slope = 1
      label = paste("beta[epsilon[Y]] == 1")
    }
    else if(pheno == "pred_Y"){
      slope = 0
      label = paste("beta[epsilon[Y]] == 0")
    }
    
    p <- p + geom_abline(slope = slope, intercept = 0, linewidth = 0.65, alpha = 0.7, color = "gray30") +
      geom_text(x = 2, y = min(indiv_data[, pheno]) + ((max(indiv_data[, pheno]) - min(indiv_data[, pheno]))/2.5), label = label, color = 'gray30', size = 3, parse = T)
  }
  
  else if(predictor == "X"){
    
    betas = data.frame("hat_beta" = c(out$hat_BXY1, out$hat_BXY2), "K" = c(1,2))
    betas$beta = c(out$BXY1, out$BXY2)
    
    ## Add β̂xʏₖ
    p <- p + geom_abline(data = betas, aes(slope = hat_beta, intercept = 0), linewidth = 0.65, alpha = 0.5, color = "orangered2") +
      geom_text(data = betas, aes(label = paste("hat(beta)[hat(X)][Y] ==", signif(hat_beta, digits = 2)), 
                                  x = max(indiv_data[, predictor]) - ((max(indiv_data[, predictor]) - min(indiv_data[, predictor]))/8),
                                  y = min(indiv_data[, pheno]) + ((max(indiv_data[, pheno]) - min(indiv_data[, pheno]))*0.10)), color = "orangered2", size = 3, parse = T) +
      ## Add true βxʏₖ
      geom_abline(data = betas, aes(slope = beta, intercept = 0), linewidth = 0.65, alpha = 0.5, color = "purple3") +
      geom_text(data = betas, aes(label = paste("beta[hat(X)][Y] ==", signif(beta, digits = 2)), 
                                  x = max(indiv_data[, predictor]) - ((max(indiv_data[, predictor]) - min(indiv_data[, predictor]))/8),
                                  y = min(indiv_data[, pheno]) + ((max(indiv_data[, pheno]) - min(indiv_data[, pheno]))*0.2)), color = "purple3", size = 3, parse = T) 
  }
  
  else if(predictor == "pred_X"){
    
    betas = data.frame("hat_beta" = c(out$hat_BXY1, out$hat_BXY2), "K" = c(1,2))
    betas$beta = c(out$BXY1, out$BXY2)
    
    ## Add β̂xʏₖ
    p <- p + geom_abline(data = betas, aes(slope = hat_beta, intercept = 0), linewidth = 0.65, alpha = 0.5, color = "orangered2") +
      geom_text(data = betas, aes(label = paste("hat(beta)[hat(X)][Y] ==", signif(hat_beta, digits = 2)), 
                                  x = max(indiv_data[, predictor]) - ((max(indiv_data[, predictor]) - min(indiv_data[, predictor]))/8),
                                  y = min(indiv_data[, pheno]) + ((max(indiv_data[, pheno]) - min(indiv_data[, pheno]))*0.10)), color = "orangered2", size = 3, parse = T) +
      ## Add true βxʏₖ
      geom_abline(data = betas, aes(slope = beta, intercept = 0), linewidth = 0.65, alpha = 0.5, color = "purple3") +
      geom_text(data = betas, aes(label = paste("beta[hat(X)][Y] ==", signif(beta, digits = 2)), 
                                  x = max(indiv_data[, predictor]) - ((max(indiv_data[, predictor]) - min(indiv_data[, predictor]))/8), 
                                  y = min(indiv_data[, pheno]) + ((max(indiv_data[, pheno]) - min(indiv_data[, pheno]))*0.2)), color = "purple3", size = 3, parse = T) 
  }
  
  return(p)
}

plot_pheno_vs_predictors <- function(indiv_data, out, plot_dir_sc_replicate){
  
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
  
  plot_grid(plotlist = plots, ncol = 3, align = "vh")
  ggsave(filename = paste0(plot_dir_sc_replicate, "/X_Y_on_predictors.pdf"), width = 13, height = 14)
}

## Plot corr between X/Y and their predictors
cor_pheno_predictors <- function(indiv_data, plot_dir_sc_replicate){
  
  variables <- c("X", "pred_X", "Y", "pred_Y", "G", "U", "eX", "eY")
  labs <- c("X", expression(hat(X) == hat(beta[GX])*G), "Y", expression(hat(Y) == hat(beta[GY])*G), 
            "G", "U", expression(epsilon[X]), expression(epsilon[Y]))
  
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
                      labels_row = labs,
                      labels_col = labs, 
                      angle_col = 90)$gtable
  }
  
  plot_grid(plotlist = h, ncol = 2, align = "h")
  ggsave(file = paste0(plot_dir_sc_replicate, "/Corr_predictors.pdf"), width = 9, height = 4)
}














