
#' @title Simulate exposure
#' @description This function simulates the exposure for all N individuals, according to the stratum they belong to.
#' @param data data.frame with columns for stratifying variable, genotype, confounder, and error terms across all N individuals.
#' @param cols vector with column names for stratifying variable (K), genotype (G), confounder (U), and errors (eX). 
#' @param drawn_BGX1 drawn effect of the variant on the exposure from ~ N(βɢx₁, σβɢx²) in stratum 1.
#' @param drawn_BGX2 drawn effect of the variant on the exposure from ~ N(βɢx₂, σβɢx²) in stratum 2.
#' @param BUX optional beta for confounder effect on X. Defaults to 1.
#' @return A numeric vector containing the exposure for all N individuals.  
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  hist(simulate_exposure(indiv_data, c("K", "G", "U", "eX"), 0.4, 0.39))
#'  }
#' }
#' @rdname simulate_exposure
#' @export


simulate_exposure = function(data, cols, drawn_BGX1, drawn_BGX2, BUX = 1){
 
  K = data[, cols[1]]
  G = data[, cols[2]]
  U = data[, cols[3]]
  eX = data[, cols[4]]
  
  X = vector()
  for(i in seq_along(K)){
    
    if(K[i] == 1){
      xi = (drawn_BGX1*G[i]) + (BUX*U[i]) + eX[i]
    } else{
      xi = (drawn_BGX2*G[i]) + (BUX*U[i]) + eX[i]
    }
    X[i] = xi
  }
  
  return(X)
}
