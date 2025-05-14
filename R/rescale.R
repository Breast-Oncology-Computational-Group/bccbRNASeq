is_scalar_numinteger <- function(n) {
  return(is.numeric(n) && isTRUE(all.equal(n, as.integer(n))) &&
           rlang::is_scalar_integer(as.integer(n)))
}


#' Get number of negative samples to achieve a balanced cohort
#' From J Clin Oncol. 2020 Dec 10; 38(35): 4184–4193 for ER balance
#'
#' @param status Logical vector for status. NAs are treated as
#' Unknown or Indeterminate status
#' @param pos_ratio Numeric between 0 and 1 indicating the proportion of
#' positive samples wanted, without taking NAs into account
#' @return Number of extra negative samples required. If negative, extra positive
#' samples are needed
#' @export
#'
#' @examples
#' v  <- sample(c(TRUE, FALSE, NA), size = 50, replace = TRUE,
#'  prob = c(0.9, 0.07, 0.03))
#' negative_samples_needed(v, 0.8)
negative_samples_needed <- function(status, pos_ratio=126/195) {

  rlang::check_required(status, "er_status")
  stopifnot("status must be a logical vector" = rlang::is_logical(status),
            "pos_ratio must be numeric" = rlang::is_scalar_double(pos_ratio),
            "pos_ratio outside [0,1]" = pos_ratio > 0 & pos_ratio < 1)
  status <- status[!is.na(status)]
  pos_n <- sum(status)
  negm_n <- as.integer(round(pos_n /pos_ratio - length(status)))

  if(negm_n < 0) {
    negm_n <- as.integer(round(((pos_ratio)/(1-pos_ratio)) * length(status) - ( pos_n/(1-pos_ratio))))
    negm_n <- negm_n*(-1)
  }
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

# This function takes a list where each element is also list with same number of
# elements, each one of them a vector with the same number of elements
factors_median <- function(factors) {
  fs_names <- names(factors[[1]])
  mfs <- lapply(fs_names, FUN = function(fn) {
    apply(do.call("cbind", lapply(factors, "[[", fn)), 1, median)
  })
  names(mfs) <- fs_names
  return(mfs)
}

## Function factory to get factors for multiple iterations in parallel
## ADD TESTS
rescale_factors_median <- function(fun_factor) {
  function(x, status, pos_ratio, n_iter, seed, cores, ...) {

    stopifnot("x must be a numeric matrix" = rlang::is_bare_numeric(x) & length(dim(x)) == 2,
              "status must be a logical vector" = rlang::is_logical(status),
              "status vector must have names" = rlang::is_named(status),
              "names in status vector must match columns in x" = all(names(status) %in% colnames(x)),
              "n_iter must be an scalar integer" = is_scalar_numinteger(n_iter),
              "seed must be an scalar integer" = is_scalar_numinteger(seed),
              "cores must be an scalar integer" = is_scalar_numinteger(cores))

    ## To use set.seed with parallel::mclapply
    ## https://stackoverflow.com/questions/30456481/controlling-seeds-with-mclapply
    RNGkind("L'Ecuyer-CMRG")
    set.seed(seed)
    status <- status[!is.na(status)]
    negm <- negative_samples_needed(status, pos_ratio)

    if(negm > 0) {
      sample_names <- names(Filter(isFALSE, status))
      message("Adding ",  negm, " negative samples for rescale")
    } else {
      sample_names <- names(Filter(isTRUE, status))
      negm <- negm*1
      message("Adding ",  negm, " positive samples for rescale")
    }

    factors <- parallel::mclapply(FUN = function(i) {
      # Factors should be calculated with the non NA status columns
      rmatrix <- expand_matrix(x[, names(status)], sample_names, negm)
      fun_factor(rmatrix, ...)
    }, X = 1:n_iter, mc.cores = cores)

    message("Finished ", n_iter, " factors iterations")
    return(factors_median(factors))
  }
}

## Function factory to get factors for multiple iterations for testing
rescale_factors_median_noparallel <- function(fun_factor) {
  function(x, status, pos_ratio, n_iter, seed, ...) {

    stopifnot("x must be a numeric matrix" = rlang::is_bare_numeric(x) & length(dim(x)) == 2,
              "er_status must be a logical vector" = rlang::is_logical(status),
              "er_status vector must have names" = rlang::is_named(status),
              "names in er_status vector must match columns in x" = all(names(status) %in% colnames(x)),
              "n_iter must be an scalar integer" = is_scalar_numinteger(n_iter),
              "seed must be an scalar integer" = is_scalar_numinteger(seed))


    stopifnot("ER status vector must have names" = rlang::is_named(status))
    set.seed(seed)
    negm <- negative_samples_needed(status)

    if(negm > 0) {
      sample_names <- names(Filter(isFALSE, status))
    } else {
      sample_names <- names(Filter(isTRUE, status))
      negm <- negm*1
    }

    factors <- lapply(1:n_iter, function(i) {
        rmatrix <- expand_matrix(x[, names(status)], sample_names, negm)
        fun_factor(rmatrix, ...)
    })
    return(factors_median(factors))
  }
}

# Function to get q (ma, mi) factors with cohort
# balancing
qfactors_median <- rescale_factors_median(qfactors)

# Function to get lower and upper quartile and decile
# factors with cohort balancing
lufactors_median <- rescale_factors_median(lufactors)

# Function to get factors with user specified probabi
sfactors_median <- rescale_factors_median(sfactors)

qfactors_median_np <- rescale_factors_median_noparallel(qfactors)
lufactors_median_np <- rescale_factors_median_noparallel(lufactors)
sfactors_median_np <- rescale_factors_median_noparallel(sfactors)

#' Rescale expression matrix with ma/mi factors
#'
#' @param x Numeric matrix with expression values
#' @param status Logical vector for status. NAs are treated as
#' Unknown or Indeterminate status
#' @param pos_ratio Numeric between 0 and 1 indicating the proportion of
#' positive samples wanted, without taking NAs into account
#' @param n_iter Number of iterations for sampling
#' @param seed Random seed
#' @param cores Number of cores for parallel processing
#'
#' @return Numeric rescaled matrix as a balance cohort
#' @export
rescale_expr <- function(x, status, pos_ratio, n_iter = 1000,
                         seed = 37, cores = getOption("mc.cores", 2L)) {

  factor_medians <- qfactors_median(x, status, pos_ratio, n_iter, seed, cores)

  x_res <- apply(x, 2, function(x) {
    X <- (x - factor_medians$mif) / (factor_medians$maf - factor_medians$mif)
    X <- (X - 0.5)*2
  })

  return(as.matrix(x_res))
}

#' Rescale expression matrix with user provided factors
#'
#' @param x Numeric matrix with expression values
#' @param status Logical vector for status. NAs are treated as
#' Unknown or Indeterminate status
#' @param pos_ratio Numeric between 0 and 1 indicating the proportion of
#' positive samples wanted, without taking NAs into account
#' @param qs Named vector for quantile probabilities.
#' Rownames in q and names in qs must match
#' @param n_iter Number of iterations for sampling
#' @param seed Random seed
#' @param cores Number of cores for parallel processing
#'
#' @return Numeric rescaled matrix as a balance cohort
#' @export
rescale_expr_median_center <- function(x, status, pos_ratio, qs, n_iter = 1000,
                                  seed = 37, cores = getOption("mc.cores", 2L)) {

  fms <- sfactors_median(x, status, pos_ratio, n_iter, seed, cores, qs)
  fms <- unlist(fms, use.names = TRUE)
  message(paste0(fms, collapse = ","))
  x_res <- apply(x, 2, function(s) s - fms)

  return(as.matrix(x_res))
}


