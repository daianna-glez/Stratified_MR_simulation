
#' @title Simulate exposure
#' @description This function simulates the exposure for all N individuals, according to the stratum they belong to.
#' @param K vector for stratifying variable across all N individuals.
#' @param G vector with genotype across all N individuals.
#' @param U vector for confounder across all N individuals.
#' @param eX vector with error term for X across all N individuals.
#' @param BGX1 effect of the variant on the exposure in stratum 1.
#' @param BGX2 effect of the variant on the exposure in stratum 2.
#' @param BCX optional beta for confounder effect on X. Defaults to 1.
#' @return A numeric vector containing the exposure for all N individuals.  
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  hist(simulate_exposure(K, G, U, eX, 0.4, 0.39))
#'  }
#' }
#' @rdname simulate_exposure
#' @export


simulate_exposure = function(K, G, U, eX, drawn_BGX1, drawn_BGX2, BUX = 1){
  
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
