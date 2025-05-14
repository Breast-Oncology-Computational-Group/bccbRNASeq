#' Get list of expressed genes plus genes of interest
#'
#' @param counts Numeric matrix of gene counts. Rownames are gene symbols, columns are samples
#' @param min_count Minimum number of counts
#' @param min_total_count Minimum number of total counts
#'
#' @return A character vector of expressed genes.
expressed_genes <- function(counts, min_count, min_total_count) {

  if(!rlang::is_installed("edgeR")) {
    stop("Please install edgeR from Bioconductor")
  }

  expressed <- suppressMessages(edgeR::filterByExpr(counts, min.count = min_count,
                    min.total.count = min_total_count, large.n = 10, min.prop = 0.7))

  message("Found ", sum(expressed), " expressed genes")
  return(union(names(expressed[expressed == TRUE]), bccb_genes))
}

#' Get genes from our curated list missing from a counts matrix
#'
#' @param counts Numeric matrix of gene counts. Rownames are gene symbols, columns are samples
#'
#' @return Character vector of gene symbols
missing_genes <- function(counts) {
  stopifnot("counts must be a numeric matrix" = rlang::is_bare_numeric(counts) & length(dim(counts)) == 2)
  return(sort(bccb_genes[!bccb_genes %in% rownames(counts)]))
}


#' Filter counts and normalized counts (ex. tpms)  for expressed genes
#'
#' @inheritParams expressed_genes
#' @param genes_keep Character vector of genes to keep in the matrix regarding their counts
#' @param norm_counts Numeric matrix with normalized counts (ex. tpms). Should have the same dimensions as counts.
#' @return  If norm_counts are not provided, it returns a numeric matrix with expressed genes. If norm_counts are provided, it returns a  List of size two with counts matrix and norm_counts filtered by expressed genes with genes in the BCCB curated list and
#' genes passed in the genes_keep parameter.
#'
#' This function uses the filterByExpr function from the edgeR Bioconductor package,
#' and re-adds genes of interest regarding of their expression values.
#' These include PAM50 genes, genes in commonly used
#' expression signatures and oncogenic genes.
#'
#' @export
subset_expressed_genes <- function(counts, min_count = 100, min_total_count = 150, genes_keep = NULL, norm_counts = NULL) {

  stopifnot("counts must be a numeric matrix" = rlang::is_bare_numeric(counts) & length(dim(counts)) == 2,
            "counts must have rownames" = !is.null(rownames(counts)),
            "min_count must be an scalar integer" = is_scalar_numinteger(min_count),
            "min_total_count must be an scalar integer" = is_scalar_numinteger(min_total_count),
            "genes_keep must be a character vector" = is.null(genes_keep) || rlang::is_character(genes_keep),
            "norm_counts must be a numeric matrix" =  is.null(norm_counts) || (rlang::is_bare_numeric(norm_counts) && length(dim(norm_counts)) == 2),
            "counts and norm_counts should have the same dimensions" =  is.null(norm_counts) || (all(dim(counts) == dim(norm_counts))),
            "counts and norm_counts should have the same rownames" =  is.null(norm_counts) ||
              all(sort(rownames(counts)) == sort(rownames(norm_counts))))

  not_found <- length(missing_genes(counts))
  if(not_found > 1) {
    warning(not_found, "/", length(bccb_genes), " BCCB genes were not found in the matrix")
  }

  not_found <- sum(!genes_keep %in% rownames(counts))
  if(not_found > 1) {
    warning(not_found, "/", length(genes_keep), " genes_keep were not found in the matrix")
  }

  e_genes <- expressed_genes(counts, min_count = min_count, min_total_count = min_total_count)
  e_genes <- intersect(rownames(counts), union(e_genes, genes_keep))

  if(is.null(norm_counts)) {
    return(counts[e_genes, ])
  } else{
    return(list(counts = counts[e_genes, ], norm_counts = norm_counts[e_genes, ]))
  }

}


