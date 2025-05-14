#' Filter matrix to PAM50 genes
#'
#' @param x Numeric matrix with expression values
#'
#' @return Numeric matrix subsetted for PAM50 genes
subset_pam50 <- function(x) {

  id_matches <- sapply(colnames(PAM50genes), function(id) {
    length(intersect(rownames(x), PAM50genes[[id]]))
  })

  if(any(id_matches == 50)) {
    col_match <- names(which(id_matches == 50))
    message("Found PAM50 genes matches with ", col_match)
  } else {
    stop("Rownames must include complete matches for PAM50 genes in emsembl_id, entrez_id or hugo_symbol")
  }

  if(col_match != "hugo_symbol") {
    x <- x[PAM50genes[[col_match]], ]
    rownames(x) <- PAM50genes[["hugo_symbol"]]
  } else {
    x <- x[PAM50genes[["hugo_symbol"]], ]
  }
  return(x)
}


#' Get genefu predictions
#'
#' @param x Numeric matrix with expression values
#'
#' @return Genefu object with PAM50 predictions
genefu_pam50 <- function(x) {

  if(!rlang::is_installed("genefu")) {
    stop("Please install genefu from Bioconductor")
  }

  x_pam50 <- subset_pam50(x)
  annot <- PAM50genes[ ,c("entrez_id", "hugo_symbol"), drop = FALSE]
  rownames(annot) <- PAM50genes[["hugo_symbol"]]
  colnames(annot) <- c("EntrezGene.ID", "probe")

  genefu_pam50_model <- utils::data('pam50', package = "genefu", envir = environment())
  pam50_preds <- genefu::intrinsic.cluster.predict(sbt.model = get(genefu_pam50_model),
                                                   data = t(x_pam50),
                                                   annot = annot,
                                                   do.mapping = TRUE,
                                                   mapping = annot,
                                                   verbose = TRUE)
  return(pam50_preds)
}

#' Get PAM50 correlations
#'
#' @param pam50_preds Genefu object
#'
#' @return Data Frame with PAM50 correlations
#' @export
#'
pam50_mixed <- function(pam50_preds){

  # extract mixed subtypes, probabilities, and correlations
  PAM50_mixed_th <- 0.66
  min_PAM50_mixed <- 0.1
  mixed_PAM50_ratio <- 0.5
  threshold_name <- as.character(100 * PAM50_mixed_th)


  pam50_mixed <- pam50_preds$subtype.proba |>
    as.data.frame() |>
    tibble::rownames_to_column("sample_id") |>
    tidyr::pivot_longer(cols = -"sample_id", names_to = 'subtype', values_to = 'proba') |>
    dplyr::group_by(.data$sample_id) |>
    tidyr::nest() |>
    dplyr::mutate(subt_probs = purrr::map(data, ~.x |>
                                            dplyr::pull(proba, name = subtype) |>
                                            sort(decreasing = TRUE)),
                  PAM50_subtype = purrr::map(.data$subt_probs, ~names(.x[1])),
                  PAM50_subtype_mixed = purrr::map(.data$subt_probs, ~{
                    if (.x[1] > PAM50_mixed_th) {
                      return(names(.x[1]))
                    } else {
                      mixed_subt <- lapply(2:length(.x), function(i) {
                        if(.x[i] > min_PAM50_mixed | .x[i]/.x[2] > mixed_PAM50_ratio) {
                          return(names(.x[i]))
                        }
                        return(NA)
                      })
                      mixed_subt <- unlist(mixed_subt)[!is.na(mixed_subt)]
                      return(paste(c(names(.x[1]), mixed_subt), collapse = '_'))
                    }
                  })) |>
    dplyr::select("sample_id", "PAM50_subtype", "PAM50_subtype_mixed", "subt_probs") |>
    tidyr::unnest_wider("subt_probs") |>
    dplyr::rename_with(.fn = ~paste0('pr_', .), .cols = -c("sample_id", tidyr::starts_with("pam50"))) |>
    tidyr::unnest(c("PAM50_subtype", "PAM50_subtype_mixed")) |>
    dplyr::ungroup()

  # merge with correlation
  pam50_mixed <- pam50_mixed |>
    dplyr::left_join(pam50_preds$cor |>
                       as.data.frame() |>
                       tibble::rownames_to_column('sample_id') |>
                       dplyr::rename_with(.fn = ~paste0('corr_', .), .cols = -c("sample_id")) |>
                       dplyr::mutate(max_corr = apply(pam50_preds$cor, FUN = max, 1)), by = "sample_id") |>
    dplyr::mutate(PAM50_subtype = ifelse(.data$max_corr < 0.1, "NC", .data$PAM50_subtype),
                  PAM50_subtype_mixed = ifelse(.data$max_corr < 0.1, "NC", .data$PAM50_subtype_mixed))

  return(pam50_mixed)
}

#' Get pam50 predictions with ER status rescaling
#' This method requires a mixed ER+/ER- cohort. Rescaling is performed to match
#' the ER+ proportion in the UNC cohort.
#'
#' @param x Numeric matrix with expression values
#' @param er_status Logical vector indicating ER positivity. NAs are treated as
#' Unknown or Indeterminate status
#' @param n_iter Number of iterations for sampling
#' @param seed Random seed
#' @param cores Number of cores for parallel processing
#'
#' @return Data frame with PAM50 predictions for each sample and its corresponding subtype correlations
#' @export
pam50_rescaled <- function(x, er_status, n_iter = 1000, seed = 37,
                           cores = getOption("mc.cores", 2L)) {
  stopifnot("x must be a numeric matrix" = rlang::is_bare_numeric(x) & length(dim(x)) == 2,
            "er_status must be a logical vector" = rlang::is_logical(er_status),
            "er_status vector must have names" = rlang::is_named(er_status),
            "names in er_status vector must match columns in x" = all(names(er_status) %in% colnames(x)),
            "n_iter must be an scalar integer" = is_scalar_numinteger(n_iter),
            "seed must be an scalar integer" = is_scalar_numinteger(seed),
            "cores must be an scalar integer" = is_scalar_numinteger(cores))

  rx <- rescale_expr(x, er_status, unc_erpos, n_iter, seed, cores)
  gpam50 <- genefu_pam50(rx)
  return(pam50_mixed(gpam50))
}



#' Get pam50 predictions using gene centering.
#' This method is intended for single receptor cohorts.
#'
#' @param x Numeric matrix with expression values
#' @param receptor_status Receptor status of the samples in the cohort. Accepted
#' values are: ERpos_HER2neg, HER2pos_ERneg, HER2pos_ERpos, TNBC
#'
#' @returns Data frame with PAM50 predictions for each sample and its corresponding subtype correlations
#' @export
pam50_receptor_gene_centering <- function(x, receptor_status) {
  stopifnot("x must be a numeric matrix" = rlang::is_bare_numeric(x) & length(dim(x)) == 2)
  rlang::arg_match(receptor_status, values = c("ERpos_HER2neg", "HER2pos_ERneg", "HER2pos_ERpos", "TNBC"))

  maf_mif <- qfactors(x)
  qs <- PAM50sigma[,receptor_status]

  x_res <- apply(x, 2, function(x) {
    X <- (x - maf_mif$mif) / (maf_mif$maf - maf_mif$mif)
    X <- (X - 0.5)*2
  })

  x_res_pam50 <- subset_pam50(x_res)

  sfs <- unlist(sfactors(x_res_pam50, qs))
  x_res_pam50 <- x_res_pam50[names(sfs), ]
  x_res_pam50 <- apply(x_res_pam50, 2, function(s) s - sfs)
  gpam50 <- genefu_pam50(x_res_pam50)
  return(pam50_mixed(gpam50))
}



