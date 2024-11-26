#' Get deciles/quantiles expression classification with cohort rebalancing
#'
#' @inheritParams pam50_rescaled
#' @return Character matrix with the following categories: "Upper Decile",
#' "Upper Quartile", "Lower Decile", "Lower Quartile"
#' @export
#'
classification_rescaled <- function(x, er_status, n_iter = 100, seed = 37, cores = getOption("mc.cores", 2L)) {

  stopifnot("x must be a numeric matrix" = rlang::is_bare_numeric(x) & length(dim(x)) == 2,
            "er_status must be a logical vector" = rlang::is_logical(er_status),
            "er_status vector must have names" = rlang::is_named(er_status),
            "names in er_status vector must match columns in x" = all(names(er_status) %in% colnames(x)),
            "n_iter must be an scalar integer" = is_scalar_numinteger(n_iter),
            "seed must be an scalar integer" = is_scalar_numinteger(seed),
            "cores must be an scalar integer" = is_scalar_numinteger(cores))

  factor_medians <- iterative_luf(x, er_status, n_iter, seed, cores)

  class_expr <- apply(x, 2, function(x) {
    dplyr::case_when(
      x > factor_medians$uq10 ~ "Upper Decile",
      x > factor_medians$uq ~ "Upper Quartile",
      x < factor_medians$lq10 ~ "Lower Decile",
      x < factor_medians$lq ~ "Lower Quartile",
      .default = ""
    )
  })
  rownames(class_expr) <- rownames(x)

  return(class_expr)
}
