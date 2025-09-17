
#' @title Simulate genotype data 
#' @description This function simulates genotypes in Hardy-Weinberg equilibrium for a bi-allelic SNP of MAF q across N individuals.
#' @param N Total number of individuals. 
#' @param q Frequency for the minor allele of the SNP.
#' @param centered Should the minor allele dosages be centered to 0 or not? Defaults to F.
#' @return A numeric vector containing the minor allele dosages for all N individuals.  
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  table(simulate_genotype(100, 0.25, F))
#'  }
#' }
#' @rdname simulate_genotype
#' @export


simulate_genotype = function(N, q, centered = F){
  
  ## Confirm q corresponds to minor allele
  if (!(q <= 0.5)) stop("MAF cannot be >0.5. Consider switching the minor allele.")
  
  G = rbinom(n = N, 
             size = 2, 
             prob = q)
  
  if(centered == T){
    
    ## Observed q and p=1-q
    q_obs = (table(G)["1"] + (2*table(G)["2"])) / (2*N)
    p_obs = (table(G)["1"] + (2*table(G)["0"])) / (2*N)
    
    ## Genotype mean and variance (approx if large N)
    mu_G = 2*q_obs
    sigma_G = 2*q_obs*p_obs
    
    ## Centered G
    G =  G - mu_G
    
    # Scaled G?
    # G =  G/sqrt(sigma_G)
  }
  
  return(G)
  
}
  