
#' @title Simulate exposure
#' @description This function simulates the exposure for all N individuals, according to the stratum they belong to.
#' @param data data.frame with columns for stratifying variable, genotype, confounder, and error terms across all N individuals.
#' @param cols vector with column names for stratifying variable (K), genotype (G), confounder (U), and error (eX). 
#' @param BGX1 true effect of the variant on X in stratum 1.
#' @param BGX2 true effect of the variant on X in stratum 2.
#' @param BUX confounder effect on X (assumed to be the same in both strata).
#' @return A numeric vector containing the exposure for all N individuals.  
#' @rdname simulate_exposure
#' @export


simulate_exposure = function(data, cols, BGX1, BGX2, BUX){
 
  K = data[, cols[1]]
  G = data[, cols[2]]
  U = data[, cols[3]]
  eX = data[, cols[4]]
  
  X = vector()
  for(i in seq_along(K)){
    
    if(K[i] == 1){
      xi = (BGX1*G[i]) + (BUX*U[i]) + eX[i]
    } else{
      xi = (BGX2*G[i]) + (BUX*U[i]) + eX[i]
    }
    X[i] = xi
  }
  
  return(X)
}
