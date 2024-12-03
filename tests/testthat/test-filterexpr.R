test_that("expressed_genes throws an error if counts is not provided", {
  expect_error(expressed_genes())
})

test_that("missing_genes throws an error if counts is not provided", {
  expect_error(missing_genes())
})

test_that("missing_genes throws an error if counts is not numeric", {
  expect_error(missing_genes(counts = letters[1:5]), "counts must be a numeric matrix")
  expect_error(missing_genes(counts = matrix(rep(TRUE, 12), nrow = 3)), "counts must be a numeric matrix")
})

test_that("expressed_genes returns correct message", {
  nm <- named_nm_matrix(10, 15, min = 1000, max = 1500)
  nm <- rbind(nm, named_nm_matrix(5, 15))
  expect_message(expressed_genes(nm, min_count = 100, min_total_count = 150), "Found 10 expressed genes")
})

test_that("subset_expressed_genes throws an error if counts is not provided", {
  expect_error(subset_expressed_genes())
})

test_that("subset_expressed_genes throws an error if counts is not numeric", {
  expect_error(subset_expressed_genes(counts = letters[1:5]), "counts must be a numeric matrix")
  expect_error(subset_expressed_genes(counts = matrix(rep(TRUE, 12), nrow = 3)), "counts must be a numeric matrix")
})

test_that("subset_expressed_genes throws an error if min_count is not an scalar integer", {
  nm <- nm_matrix(5, 4)
  expect_error(subset_expressed_genes(counts = nm,  min_count = "a"),
               "min_count must be an scalar integer")
  expect_error(subset_expressed_genes(counts = nm, runif(3)),
               "min_count must be an scalar integer")
})

test_that("subset_expressed_genes throws an error if min_total_count is not an scalar integer", {
  nm <- nm_matrix(5, 4)
  expect_error(subset_expressed_genes(counts = nm,  min_total_count = "a"),
               "min_total_count must be an scalar integer")
  expect_error(subset_expressed_genes(counts = nm, min_total_count = runif(3)),
               "min_total_count must be an scalar integer")
})

test_that("subset_expressed_genes throws an error if norm_counts is not numeric", {
  nm <- named_nm_matrix(10, 15, min = 1000, max = 1500)
  expect_error(subset_expressed_genes(counts = nm, norm_counts = letters[1:5]), "norm_counts must be a numeric matrix")
  expect_error(subset_expressed_genes(counts = nm, norm_counts = matrix(rep(TRUE, 12), nrow = 3)), "norm_counts must be a numeric matrix")
})

test_that("subset_expressed_genes throws an error if genes_keep is not character", {
  nm <- named_nm_matrix(10, 15, min = 1000, max = 1500)
  expect_error(subset_expressed_genes(counts = nm, genes_keep = runif(10)), "genes_keep must be a character vector")
  expect_error(subset_expressed_genes(counts = nm, genes_keep = rep(TRUE, 12)), "genes_keep must be a character vector")
})

test_that("subset_expressed_genes throws error if counts and norm_counts have different dimensions", {
  counts <- nm_matrix(5, 4)
  norm_counts <- nm_matrix(3, 8)
  expect_error(subset_expressed_genes(counts = counts, norm_counts = norm_counts), "counts and norm_counts should have the same dimensions")
})

test_that("subset_expressed_genes throws warning when missing BCCB genes found", {
  n <- sample.int(100, size = 1) + 1
  nm <- named_nm_matrix(100, 15)
  egm <- mock(sample(rownames(nm), size = 50))
  mgm <- mock(sample(letters, size = n, replace = TRUE))

  with_mocked_bindings({
    expect_warning(subset_expressed_genes(nm), paste0(n, "/", length(bccb_genes), " BCCB genes were not found in the matrix. Are you using gene symbols?"))
  }, expressed_genes = egm, missing_genes = mgm)
})

test_that("subset_expressed_genes throws warning when missing genes found", {
  n <- sample.int(100, size = 1) + 2
  nm <- named_nm_matrix(100, 15)
  genes <- sample(letters, 5)
  mgm <- mock(c())
  egm <- mock(sample(rownames(nm), size = 50))
  with_mocked_bindings({
    expect_warning(subset_expressed_genes(nm, genes_keep = genes),
                   paste0(length(genes), "/", length(genes), " genes_keep were not found in the matrix"))
  }, expressed_genes = egm, missing_genes = mgm)
})

test_that("subset_expressed_genes throws error if counts and norm_counts have different dimension names", {
  counts <- named_nm_matrix(5, 8)
  norm_counts <- named_nm_matrix(5, 8)
  genes <- paste0("GENE", 1:5)
  rownames(counts) <- genes
  rownames(norm_counts) <- c("GENEA", sample(genes)[1:4])
  expect_error(subset_expressed_genes(counts = counts, norm_counts = norm_counts), "counts and norm_counts should have the same rownames")
})

test_that("subset_expressed_genes filters matrices correctly", {
  n <- sample.int(20, size = 1) +1
  genes <- paste0("GENE", 1:n)
  nn <- sample.int(20, size = 1) +1
  genesn <-  paste0("GENE11", 1:nn)
  egm <- mock(genes, cycle = TRUE)
  mgm <- mock(c())
  counts <- named_nm_matrix(n+nn, 10)
  rownames(counts) <- c(genes, genesn)
  norm_counts <- named_nm_matrix(n+nn, 10)
  rownames(norm_counts) <- c(genes, genesn)
  with_mocked_bindings({
    filtered <- subset_expressed_genes(counts, norm_counts = norm_counts, genes_keep = genesn)
    expect_equal(filtered, list(counts = counts, norm_counts = norm_counts))
  }, expressed_genes = egm, missing_genes = mgm)
})

test_that("subset_expressed_genes filters count correctly", {
  n <- sample.int(20, size = 1) +1
  genes <- paste0("GENE", 1:n)
  nn <- sample.int(20, size = 1) +1
  genesn <-  paste0("GENE11", 1:nn)
  egm <- mock(genes, cycle = TRUE)
  mgm <- mock(c())

  counts <- named_nm_matrix(n+nn, 10)
  rownames(counts) <- c(genes, genesn)
  with_mocked_bindings({
    filtered <- subset_expressed_genes(counts, genes_keep = genesn)
    expect_equal(filtered, counts)
  }, expressed_genes = egm, missing_genes = mgm)
})

# TO DO: test function calls
