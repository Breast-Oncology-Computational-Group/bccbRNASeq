library(mockery)

nm_matrix <- function(n, m, min = 0, max = 1) {
  matrix(runif(n*m, min = min, max = max), nrow = n, ncol = m,
         dimnames = list(paste0("r", 1:n), NULL))
}

named_nm_matrix <- function(n, m, min = 0, max = 1) {
  nm <- nm_matrix(n, m, min = min, max = max)
  colnames(nm) <- paste0("sample", 1:m)
  rownames(nm) <- paste0("gene", 1:n)
  return(nm)
}

log_vector <- function(n, ns = NULL, prob = c(0.75, 0.20, 0.5)) {
  lv <- sample(c(TRUE, FALSE, NA), n, replace = TRUE, prob = prob)
  names(lv) <- ns
  return(lv)
}
