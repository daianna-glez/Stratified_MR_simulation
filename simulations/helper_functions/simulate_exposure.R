
#' @title Simulate exposure
#' @description This function simulates the exposure for all N individuals, according to the stratum they belong to.
#' @param data data.frame with columns for stratifying variable, genotype, confounder, error, and βɢxₖ across all N individuals.
#' @param cols vector with column names for stratifying variable (K), genotype (G), confounder (U), error (eX), and βɢxₖ (BGX). 
#' @param BUX confounder effect on X (assumed to be the same in both strata).
#' @return A numeric vector containing the exposure for all N individuals.  
#' @rdname simulate_exposure
#' @export


simulate_exposure = function(data, cols, BUX){
 
  K = data[, cols[1]]
  G = data[, cols[2]]
  U = data[, cols[3]]
  eX = data[, cols[4]]
  BGX = data[, cols[5]]
  
  X = (BGX*G) + (BUX*U) + eX

  return(X)
}
