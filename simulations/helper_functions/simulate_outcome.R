
#' @title Simulate outcome
#' @description This function simulates the outcome for all N individuals, according to their genetically proxied X and the stratum they belong to.
#' @param data data.frame with columns for stratifying variable, genotype, confounder, error, βɢxₖ, and βxʏₖ across all N individuals.
#' @param cols vector with column names for stratifying variable (K), genotype (G), confounder (U), error (eY), βɢxₖ (BGX), and βxʏₖ (BXY). 
#' @param BUX confounder effect on X (assumed to be the same in both strata).
#' @param BUY confounder effect on Y (assumed to be the same in both strata).
#' @return A numeric vector containing the outcome for all N individuals.  
#' @rdname simulate_outcome
#' @export

simulate_outcome = function(data, cols, BUX, BUY){
  
  K = data[, cols[1]]
  G = data[, cols[2]]
  U = data[, cols[3]]
  eY = data[, cols[4]]
  BGX = data[, cols[5]]
  BXY = data[, cols[6]]
  
  Y = BXY*((BGX*G) + (BUX*U)) + (BUY*U) + eY
  return(Y)
}
