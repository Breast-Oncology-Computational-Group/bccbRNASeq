% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rescale.R
\name{expand_matrix}
\alias{expand_matrix}
\title{Function to add extra columns in a matrix by sampling from
a character vector of column names}
\usage{
expand_matrix(x, column_names, n = length(column_names))
}
\arguments{
\item{x}{Numeric matrix}

\item{column_names}{Column names to be sampled}

\item{n}{Number of sampled columns}
}
\value{
A numeric matrix of columns in x plus n sampled columns from column_names
}
\description{
Function to add extra columns in a matrix by sampling from
a character vector of column names
}
