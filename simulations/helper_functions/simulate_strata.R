
#' @title Simulate strata
#' @description This function stratifies N individuals into two strata considering their sample size ratio.
#' @param N Total number of individuals. 
#' @param r Ratio of stratum sample sizes (N1/N2). 
#' @return K, a binary vector indicating the stratum each individual was assigned to.  
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  table(simulate_strata(1000, 0.33))
#'  }
#' }
#' @rdname simulate_strata
#' @export


simulate_strata = function(N, r){
  
  stratum_1 <- rep(1, round((r*N)/(r+1)))
  stratum_2 <- rep(2, round(N/(r+1)))
  
  K <- c(stratum_1, stratum_2)
  
  if(length(K) != N){
    message("Stopping: not total N in strata generation.")
  }
  
  return(K)
}


