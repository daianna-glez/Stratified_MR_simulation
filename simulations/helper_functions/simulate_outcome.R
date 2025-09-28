
#' @title Simulate outcome
#' @description This function simulates the outcome for all N individuals, according to their genetically proxied X and the stratum they belong to.
#' @param data data.frame with columns for stratifying variable, genotype, confounder, and error terms across all N individuals.
#' @param cols vector with column names for stratifying variable (K), genotype (G), confounder (U), and errors (eY). 
#' @param drawn_BGX1 drawn effect of the variant on the exposure from ~ N(βɢx₁, σβɢx²) in stratum 1.
#' @param drawn_BGX2 drawn effect of the variant on the exposure from ~ N(βɢx₂, σβɢx²) in stratum 2.
#' @param drawn_BXY1 drawn causal effect of the exposure on the outcome from ~ N(βxʏ₁, σβxʏ²) in stratum 1.
#' @param drawn_BXY2 drawn causal effect of the exposure on the outcome from ~ N(βxʏ₂, σβxʏ²) in stratum 2.
#' @param BUY optional beta for confounder effect on Y. Defaults to 1.
#' @return A numeric vector containing the outcome for all N individuals.  
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  hist(simulate_outcome(indiv_data, c("K", "G", "U", "eY"), 0.4, 0.39))
#'  }
#' }
#' @rdname simulate_outcome
#' @export

simulate_outcome = function(data, cols, drawn_BGX1, drawn_BGX2, drawn_BXY1, drawn_BXY2, BUY = 1){
  
  K = data[, cols[1]]
  G = data[, cols[2]]
  U = data[, cols[3]]
  eY = data[, cols[4]]
  
  Y = vector()
  
  for(i in seq_along(K)){
    
    if(K[i] == 1){
      yi = (drawn_BGX1*drawn_BXY1*G[i]) + (BUY*U[i]) + eY[i]
    } else{
      yi = (drawn_BGX2*drawn_BXY2*G[i]) + (BUY*U[i]) + eY[i]
    }
    Y[i] = yi
  }
  
  return(Y)
}
