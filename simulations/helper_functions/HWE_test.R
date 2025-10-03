
#' @title Test for Hardy-Weinberg Equilibrium
#' @description This function performs a chi-squared test to assess whether the observed genotype frequencies of a bi-allelic SNP deviate from Hardy–Weinberg equilibrium.
#' @param n_aa observed number of homozygote individuals for minor allele "a"
#' @param n_Aa observed number of heterozygote individuals 
#' @param n_AA observed number of heterozygote individuals for major allele "A"
#' @return A vector with computed $\chi^2$ test statistic and its p-value for 1 degree of freedom.
#' @rdname HWE_test
#' @export 
#' 

HWE_test = function(n_aa, n_aA, n_AA){

  # Step 1: get observed allele frequencies
  N = n_aa + n_aA+ n_AA
  q_obs = (2*n_aa + n_aA)/(2*N)
  p_obs = (2*n_AA + n_aA)/(2*N)
  
  # Step 2: compute expected genotype frequencies
  aa_exp <- q_obs^2
  aA_exp <- 2*q_obs*p_obs
  AA_exp <- p_obs^2
  
  # Step 3: get observed genotype frequencies
  aa_obs <- n_aa / N
  aA_obs <- n_aA / N
  AA_obs <- n_AA / N
  
  # Step 4: calculate chi^2 test statistic
  X_sq = (((aa_obs - aa_exp)^2 /aa_exp) + 
          ((aA_obs - aA_exp)^2 / aA_exp) +
          ((AA_obs - AA_exp)^2 / AA_exp)) * N
  
  ## Step 5: compute p-value from the chi^2 distribution with 
  #          df = # categories - # constraints 
  #             = # genotypes - 
  #               1 df by fixed N (constraining obs genotype counts) - 
  #               1 df by fixed q (determining exp genotype counts) = 3 - 1 - 1
  pval = pchisq(X_sq, df=1, lower.tail = F)
  
  res = data.frame(CHISQ = unname(X_sq), P.value = unname(pval))
  
  return(res)
  
}

# -----------------------------
## Test positive control 
# q = 0.3
# p = 1 - q
# N = 10000
# 
# n_aa = q^2 * N
# n_aA = 2*q*p * N
# n_AA = p^2 * N
# 
# HWE_test(n_aa, n_Aa, n_AA)
# # CHISQ P.value 
# #     0       1 

# -----------------------------
## Test negative control
# n_aa2 = n_aa - 100
# n_aA2 = n_aA + 20
# n_AA2 = n_AA + 80
# 
# HWE_test(n_aa2, n_Aa2, n_AA2)
# #       CHISQ      P.value 
# #    4.29622708 0.03819702 



