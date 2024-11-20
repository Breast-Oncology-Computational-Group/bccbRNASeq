library(mockery)

nm_matrix <- function(n, m) {
  matrix(runif(n*m), nrow = n, ncol = m,
         dimnames = list(paste0("r", 1:n), NULL))
}

named_nm_matrix <- function(n, m) {
  nm <- nm_matrix(n, m)
  colnames(nm) <- paste0("sample", 1:m)
  return(nm)
}

log_vector <- function(n, ns = NULL, prob = c(0.75, 0.20, 0.5)) {
  lv <- sample(c(TRUE, FALSE, NA), n, replace = TRUE, prob = prob)
  names(lv) <- ns
  return(lv)
}
