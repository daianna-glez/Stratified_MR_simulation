
#' @title Test genetic association with a trait (exposure X or outcome Y) per stratum.
#' @description This function fits a linear model to assess the association of a genetic variant with a continuous trait in each stratum.
#' @param data data.frame with columns for stratifying variable, genotype, and trait across all N individuals.
#' @param cols vector with column names for stratifying variable (K), genotype (G), and trait (P). 
#' @param out_dir_sc_replicate output path to save list with model fit results per stratum. 
#' @return A data.frame containing the estimated genetic effect β̂, se, t-statistic, and p-value, in each stratum. 
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  test_genetic_association(indiv_data, c("K", "G", "X"), out_dir_sc_replicate)
#'  }
#' }
#' @rdname test_genetic_association
#' @export


test_genetic_association = function(data, cols, out_dir_sc_replicate){
  
  ## Define col names
  K = cols[1]
  G = cols[2]
  P = cols[3]
  
  ## Linear model
  formula <- as.formula(paste(P, "~",  G))
  
  ## Fit separately per stratum
  res <- by(data, data[, K], function(k) lm(formula, data = k)) %>% as.list()
  save(res, file = paste0(out_dir_sc_replicate, "/stratum_fit_res_G_on_", P, ".Rdata"))
  
  ## Data frame with estimated genetic effects and associated stats
  df = do.call(rbind, lapply(res, function(k){summary(k)$coefficients["G",]})) %>% as.data.frame()
  df$K = c(1, 2)
 
  return(df)
}
