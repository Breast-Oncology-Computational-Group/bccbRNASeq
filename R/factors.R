#' Get factors for PAM50 normalization
#'
#' @param x  Numeric matrix of gene expression. Rows are genes, columns are samples
#'
#' @return A list of two vectors of factors per row:
#' \itemize{
#'  \item{maf} {ma factors}
#'  \item{mif} {mi factors}
#' }
#'
qfactors <- function(x) {

  q = 0.05
  stopifnot("x must be a numeric matrix" = rlang::is_bare_numeric(x) & length(dim(x)) == 2)
  list(
    maf =  apply(x, 1, function(x) {  stats::quantile(x, probs=1-(q/2), na.rm=TRUE)}),
    mif = apply(x, 1, function(x) {  stats::quantile(x, probs=q/2, na.rm=TRUE)})
  )
}

#' Get factors
#'
#' @param x  Numeric matrix of gene expression. Rows are genes, columns are samples
#' @param qs Named vector for quantile probabilities.
#' Rownames in q and names in qs must match
#' @return Named numeric vector with x row quantile factors
#'
sfactors <- function(x, qs) {
  stopifnot("x must be a numeric matrix" = rlang::is_bare_numeric(x) & length(dim(x)) == 2,
            "qs must be a numeric vector" = rlang::is_bare_double(qs) & is.null(dim(qs)),
            "qs must have names" = rlang::is_named(qs),
            "rownames in x and names in qs must match" = all(names(qs) %in% rownames(x)))

  mfs <- lapply(names(qs), function(r) {
    stats::quantile(x[r,], probs=qs[r], na.rm=TRUE, names = FALSE)
  })
  names(mfs) <- names(qs)
  return(mfs)
}


#' Get lower and upper quantile and decile factors per rows
#'
#' @param x Numeric matrix of gene expression. Rows are genes, columns are samples
#'
#' @return A list of four vectors of factors per row:
#' \itemize{
#'  \item{lq} {Lower quantile values}
#'  \item{uq} {Upper quantile values}
#'  \item{lq10} {Lower decile values}
#'  \item{uq10} {Upper decile values}
#' }
lufactors <- function(x) {

  stopifnot("x must be a numeric matrix" = is.numeric(x) & length(dim(x)) == 2)

  list(
    lq = apply(x, 1, stats::quantile, probs=0.25, na.rm=TRUE),
    uq = apply(x, 1, stats::quantile, probs=0.75, na.rm=TRUE),
    lq10 = apply(x, 1, stats::quantile, probs=0.10, na.rm=TRUE),
    uq10 = apply(x, 1, stats::quantile, probs=0.9, na.rm=TRUE)
  )
}
