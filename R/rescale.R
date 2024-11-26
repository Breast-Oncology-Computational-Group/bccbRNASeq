is_scalar_numinteger <- function(n) {
  return(is.numeric(n) && isTRUE(all.equal(n, as.integer(n))) &&
           rlang::is_scalar_integer(as.integer(n)))
}


#' Get number of ER negative samples to achieve a balanced cohort
#' according to J Clin Oncol. 2020 Dec 10; 38(35): 4184–4193.
#'
#' @param er_status Logical vector indicating ER positivity. NAs are treated as
#' Unknown or Indeterminate status
#'
#' @return Number of ER negative additional samples required. If set has enough negative
#' samples, the function returns 0.
#' @export
#'
#' @examples
#' er_neg_missing(sample(c(TRUE, FALSE, NA), size = 50, replace = TRUE, prob = c(0.9, 0.07, 0.03)))
er_neg_missing <- function(er_status) {

  rlang::check_required(er_status, "er_status")
  stopifnot("er_status must be a logical vector" = rlang::is_logical(er_status))
  unc_er_ratio <- 126/195 # 64.6% ER+

  pos_n <- sum(er_status, na.rm = TRUE)
  num_n <- sum(!is.na(er_status))
  negm_n <- as.integer(max(round(pos_n / unc_er_ratio - num_n), 0))
  ### If Positive missing, warning

  return(negm_n)
}

#' Function to add extra columns in a matrix by sampling from
#' a character vector of column names
#'
#' @param x Numeric matrix
#' @param column_names Column names to be sampled
#' @param n Number of sampled columns
#'
#' @return A numeric matrix of columns in x plus n sampled columns from column_names
expand_matrix <- function(x, column_names, n = length(column_names)) {

  rlang::check_required(column_names, "column_names")
  rlang::check_required(x)

  stopifnot("x must be a numeric matrix" = rlang::is_bare_numeric(x) & length(dim(x)) == 2,
            "column_names must be a character vector" = rlang::is_character(column_names),
            "values in column_names must match columns in x" = all(column_names %in% colnames(x)),
            "n must be an scalar integer" = is_scalar_numinteger(n))

  if(n == 0) {
    return(x)
  }
  added_cols <- sample(column_names, size = n, replace = TRUE)

  extra_cols <- x[, added_cols, drop = FALSE]
  colnames(extra_cols) <-  paste0(colnames(extra_cols), "_DUP")
  return(cbind(x, extra_cols))
}

## Function factory to get factors for multiple iterations in parallel
## ADD TESTS
iterative_factors <- function(fun_factor) {
  function(x, er_status, n_iter, seed, cores) {

    stopifnot("x must be a numeric matrix" = rlang::is_bare_numeric(x) & length(dim(x)) == 2,
              "er_status must be a logical vector" = rlang::is_logical(er_status),
              "er_status vector must have names" = rlang::is_named(er_status),
              "names in er_status vector must match columns in x" = all(names(er_status) %in% colnames(x)),
              "n_iter must be an scalar integer" = is_scalar_numinteger(n_iter),
              "seed must be an scalar integer" = is_scalar_numinteger(seed),
              "cores must be an scalar integer" = is_scalar_numinteger(cores))

    ## To use set.seed with parallel::mclapply
    ## https://stackoverflow.com/questions/30456481/controlling-seeds-with-mclapply
    RNGkind("L'Ecuyer-CMRG")
    set.seed(seed)
    negm <- er_neg_missing(er_status)
    negc <- names(Filter(isFALSE, er_status))
    factors <- parallel::mclapply(FUN = function(i) {
      rmatrix <- expand_matrix(x, negc, negm)
      fun_factor(rmatrix)
    }, X = 1:n_iter, mc.cores = cores)

    return(factors_median(factors))
  }
}

## Function factory to get factors for multiple iterations for testing
iterative_factors_noparallel <- function(fun_factor) {
  function(x, er_status, n_iter, seed) {

    stopifnot("x must be a numeric matrix" = rlang::is_bare_numeric(x) & length(dim(x)) == 2,
              "er_status must be a logical vector" = rlang::is_logical(er_status),
              "er_status vector must have names" = rlang::is_named(er_status),
              "names in er_status vector must match columns in x" = all(names(er_status) %in% colnames(x)),
              "n_iter must be an scalar integer" = is_scalar_numinteger(n_iter),
              "seed must be an scalar integer" = is_scalar_numinteger(seed))


    stopifnot("ER status vector must have names" = rlang::is_named(er_status))
    set.seed(seed)
    negm <- er_neg_missing(er_status)
    negc <- names(Filter(isFALSE, er_status))
    factors <- lapply(1:n_iter, function(i) {
        rmatrix <- expand_matrix(x, negc, negm)
        fun_factor(rmatrix)
    })
    return(factors_median(factors))
  }
}

# Function to get multiple q (ma, mi) factors with cohort
# balancing
iterative_qf <- iterative_factors(qfactors)

# Function to get multiple lower and upper quartile and decile
# factors with cohort balancing
iterative_luf <- iterative_factors(lufactors)


iterative_qf_np <- iterative_factors_noparallel(qfactors)
iterative_luf_np <- iterative_factors_noparallel(lufactors)
