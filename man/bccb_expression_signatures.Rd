% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ssgsea.R
\name{bccb_expression_signatures}
\alias{bccb_expression_signatures}
\title{Get expression signature values for selected sets}
\usage{
bccb_expression_signatures(
  norm_counts,
  signature_set = "all",
  signature_names = NULL,
  seed = 37
)
}
\arguments{
\item{norm_counts}{Numeric matrix with normalized counts. Rownames are gene symbols, columns are samples}

\item{signature_set}{Character vector of signature sets to calculate ssGSEA. Values include: "hallmarks", "breast", "all"}

\item{signature_names}{Character vector of signature names to calculate ssGSEA.}

\item{seed}{Random seed}
}
\value{
Data frame with GSEA outputs from fgsea::fgsea. The data frame contains the following columns:
\itemize{
\item{"signature"}{}
\item{"NES"}{Normalized enrichment scores}
\item{"ES"}{Enrichment score}
\item{"pval"}{}
\item{"padj"}{}
\item{"size"}{}
\item{"leadingEdge"}{Genes in the leadingEdge vector}
\item{"sample"}{Sample id}
\item{"signature_set"}{Name of the signature set}
}
}
\description{
Get expression signature values for selected sets
}
