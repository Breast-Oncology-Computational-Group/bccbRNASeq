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

  factor_medians <- lufactors_median(x, er_status, unc_erpos, n_iter, seed, cores)

  class_expr <- apply(x, 2, function(y) {
    dplyr::case_when(
      y > factor_medians$uq10 ~ "Upper Decile",
      y > factor_medians$uq ~ "Upper Quartile",
      y < factor_medians$lq10 ~ "Lower Decile",
      y < factor_medians$lq ~ "Lower Quartile",
      .default = ""
    )
  }) |>
    matrix(nrow = nrow(x), ncol = ncol(x))
  rownames(class_expr) <- rownames(x)
  colnames(class_expr) <- colnames(x)

  return(class_expr)
}
