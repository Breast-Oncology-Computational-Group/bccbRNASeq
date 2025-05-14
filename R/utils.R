#' Get PAM50 genes
#'
#' @return vector with gene symbols of PAM50 genes
#' @export
get_pam50_genes <- function() {
  return(PAM50genes$hugo_symbol)
}


#' Get LM22 genes
#'
#' @return vector with gene symbols of genes in the LM22 signature matrix
#' @export
get_lm22_genes <- function() {
  return(LM22genes)
}

#' Get UNC ER+ proportion
#'
#' @return numeric value of the UNC ER+ proportion cohort
#' @export
get_unc_erpos <- function() {
  return(unc_erpos)
}

