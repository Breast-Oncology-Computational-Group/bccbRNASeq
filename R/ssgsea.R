#' Get single sample GSEA for a list of samples
#'
#' @param norm_counts Numeric matrix with normalized counts. Rownames are gene symbols, columns are samples
#' @param set_list List of signatures for GSEA
#' @param seed Random seed
#' @param cores Number of cores for parallel processing
#'
#' @return Data frame with GSEA outputs from fgsea::fgsea. The data frame contains the following columns:
##' \itemize{
##'  \item{"signature_e"}{}
##'  \item{"NES"}{Normalized enrichment scores}
##'  \item{"ES"}{Enrichment score}
##'  \item{"pval"}{}
##'  \item{"padj"}{}
##'  \item{"size"}{}
##'  \item{"leadingEdge"}{Genes in the leadingEdge vector}
##'  \item{"sample"}{Sample id}
##' }
ss_gsea_par <- function(norm_counts, set_list, seed = 37, cores = getOption("mc.cores", 2L)) {

  if(!rlang::is_installed("fgsea")) {
    stop("Please install fgsea from Bioconductor")
  }

  stopifnot("norm_counts must be a numeric matrix" = rlang::is_bare_numeric(norm_counts) & length(dim(norm_counts)) == 2,
            "norm_counts must have rownames" = !is.null(rownames(norm_counts)),
            "seed must be an scalar integer" = is_scalar_numinteger(seed),
            "cores must be an scalar integer" = is_scalar_numinteger(cores))

  fgsea_list <- parallel::mclapply(colnames(norm_counts), function (sample) {
      vstat <- norm_counts[, sample]
      names(vstat) <- rownames(norm_counts)
      set.seed(seed)
      results <- tryCatch(
        {
          fgsea::fgsea(pathways = set_list, stats = vstat,
                       eps = 1e-8, minSize  = 10, maxSize  = 1000, nPermSimple = 10000, scoreType="pos")
        }, error = function(cond) {
          message(paste0("Error running fgsea for sample: ", sample))
          message("fgsea error:")
          message(conditionMessage(cond))
          NA
        }
      )
      results <- tibble::tibble(results)
      results |>
        dplyr::mutate(sample_id = sample)

   }, mc.cores = cores)


  fgsea_list <- dplyr::bind_rows(fgsea_list)
  fgsea_list <- fgsea_list |>
    dplyr::mutate(leadingEdge = purrr::map_chr(.data$leadingEdge,
                                               \(x) paste(x, collapse = ","))) |>
    dplyr::rename("signature" = "pathway")
  return(fgsea_list)
}

#' Get single sample GSEA for a list of samples
#'
#' @param norm_counts Numeric matrix with normalized counts. Rownames are gene symbols, columns are samples
#' @param set_list List of signatures for GSEA
#' @param seed Random seed
#'
#' @return Data frame with GSEA outputs from fgsea::fgsea. The data frame contains the following columns:
##' \itemize{
##'  \item{"signature_e"}{}
##'  \item{"NES"}{Normalized enrichment scores}
##'  \item{"ES"}{Enrichment score}
##'  \item{"pval"}{}
##'  \item{"padj"}{}
##'  \item{"size"}{}
##'  \item{"leadingEdge"}{Genes in the leadingEdge vector}
##'  \item{"sample"}{Sample id}
##' }
#' @export
ss_gsea <- function(norm_counts, set_list, seed = 37) {

  if(!rlang::is_installed("fgsea")) {
    stop("Please install fgsea from Bioconductor")
  }

  stopifnot("norm_counts must be a numeric matrix" = rlang::is_bare_numeric(norm_counts) & length(dim(norm_counts)) == 2,
            "norm_counts must have rownames" = !is.null(rownames(norm_counts)),
            "seed must be an scalar integer" = is_scalar_numinteger(seed))

  fgsea_list <- lapply(colnames(norm_counts), function (sample) {
    vstat <- norm_counts[, sample]
    names(vstat) <- rownames(norm_counts)
    set.seed(seed)
    results <- tryCatch(
      {
        fgsea::fgsea(pathways = set_list, stats = vstat,
            eps = 1e-8, minSize  = 10, maxSize  = 1000, nPermSimple = 10000, scoreType="pos")
      }, error = function(cond) {
        message(paste0("Error running fgsea for sample: ", sample))
        message("fgsea error:")
        message(conditionMessage(cond))
        NA
      }
    )
    results <- tibble::tibble(results)
    results |>
      dplyr::mutate(sample_id = sample)

  })

  fgsea_list <- dplyr::bind_rows(fgsea_list)
  fgsea_list <- fgsea_list |>
    dplyr::mutate(leadingEdge = purrr::map_chr(.data$leadingEdge,
                                               \(x) paste(x, collapse = ","))) |>
    dplyr::rename("signature" = "pathway")
  return(fgsea_list)
}


#' Get expression signature values for selected sets
#'
#' @param norm_counts Numeric matrix with normalized counts. Rownames are gene symbols, columns are samples
#' @param signature_set Character vector of signature sets to calculate ssGSEA. Values include: "hallmarks", "breast", "all"
#' @param signature_names Character vector of signature names to calculate ssGSEA.
#' @param seed Random seed
#'
#' @return Data frame with GSEA outputs from fgsea::fgsea. The data frame contains the following columns:
##' \itemize{
##'  \item{"signature"}{}
##'  \item{"NES"}{Normalized enrichment scores}
##'  \item{"ES"}{Enrichment score}
##'  \item{"pval"}{}
##'  \item{"padj"}{}
##'  \item{"size"}{}
##'  \item{"leadingEdge"}{Genes in the leadingEdge vector}
##'  \item{"sample"}{Sample id}
##'  \item{"signature_set"}{Name of the signature set}
##' }
#' @export
#'
bccb_expression_signatures <- function(norm_counts, signature_set = "all", signature_names = NULL,
                                       seed = 37) {

  rlang::check_exclusive(signature_set, signature_names, .require = FALSE)

  stopifnot("norm_counts must be a numeric matrix" = rlang::is_bare_numeric(norm_counts) & length(dim(norm_counts)) == 2,
            "norm_counts must have rownames" = !is.null(rownames(norm_counts)),
            "seed must be an scalar integer" = is_scalar_numinteger(seed))

  rlang::arg_match(signature_set, values = c("hallmark", "c2cgp_breast", "tirosh", "all"))
  if(!is.null(signature_names) && !all(signature_names %in% bccb_signature_names())) {
    stop("`signature_names` must include values in the bccb signature names set")
  }

  if(!is.null(signature_names)) {
    signature_results <- ss_gsea(norm_counts,
                                 unlist(rlang::set_names(bccb_signatures, NULL), recursive = FALSE)[signature_names],
                                 seed = seed)

  } else if(signature_set != "all") {
    message("Getting ssGSEA values for ", length(bccb_signatures[[signature_set]]), " signatures in ", signature_set, " set")
    signature_results <- ss_gsea(norm_counts, bccb_signatures[[signature_set]],
                                 seed = seed)
  } else {
    signature_results <- lapply(names(bccb_signatures), function(ss) {
      results <- ss_gsea(norm_counts, bccb_signatures[[ss]], seed = seed)
      results <- results |>
        dplyr::mutate("signature_set" = ss)
    })
    signature_results <- dplyr::bind_rows(signature_results)
  }
  return(signature_results)
}


#' Get a list of available signature names in the bccb_signatures set
#'
#' @return Character vector with all expression signatures in the bccb_signatures set
#' @export
#'
bccb_signature_names <- function() {
  return(unlist(lapply(bccb_signatures, names), use.names = FALSE))
}
