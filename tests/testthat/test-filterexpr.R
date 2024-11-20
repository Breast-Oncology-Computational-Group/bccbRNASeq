test_that("expressed_genes throws an error if counts is not provided", {
  expect_error(expressed_genes())
})

test_that("expressed_genes throws an error if counts is not numeric", {
  expect_error(expressed_genes(counts = letters[1:5]), "counts must be a numeric matrix")
  expect_error(expressed_genes(counts = matrix(rep(TRUE, 12), nrow = 3)), "counts must be a numeric matrix")
})

test_that("expressed_genes throws an error if min_count is not an scalar integer", {
  nm <- nm_matrix(5, 4)
  expect_error(expressed_genes(counts = nm,  min_count = "a"),
               "min_count must be an scalar integer")
  expect_error(expressed_genes(counts = nm, runif(3)),
               "min_count must be an scalar integer")
})

test_that("expressed_genes throws an error if min_total_count is not an scalar integer", {
  nm <- nm_matrix(5, 4)
  expect_error(expressed_genes(counts = nm,  min_total_count = "a"),
               "min_total_count must be an scalar integer")
  expect_error(expressed_genes(counts = nm, min_total_count = runif(3)),
               "min_total_count must be an scalar integer")
})

test_that("expressed_genes throws warning when missing genes found", {
  not_in_genes <- sample(bccb_genes, sample.int(20, size = 1))
  in_genes <- bccb_genes[!bccb_genes %in% not_in_genes]
  diff <- nrow(sample_expr) - length(in_genes)
  rownames(sample_expr) <- sample(c(paste0("gene_", 1:diff), in_genes))
  expect_warning(expressed_genes(sample_expr), paste0(length(not_in_genes), "/", length(bccb_genes), " BCCB genes were not found in the matrix"))
})

test_that("missing_genes throws an error if counts is not provided", {
  expect_error(missing_genes())
})

test_that("missing_genes throws an error if counts is not numeric", {
  expect_error(missing_genes(counts = letters[1:5]), "counts must be a numeric matrix")
  expect_error(missing_genes(counts = matrix(rep(TRUE, 12), nrow = 3)), "counts must be a numeric matrix")
})

test_that("missing_genes returns correct number of missing genes", {
  not_in_genes <- sample(bccb_genes, sample.int(20, size = 1))
  in_genes <- bccb_genes[!bccb_genes %in% not_in_genes]
  diff <- nrow(sample_expr) - length(in_genes)
  rownames(sample_expr) <- sample(c(paste0("gene_", 1:diff), in_genes))
  expect_equal(sort(not_in_genes), missing_genes(sample_expr))
})

test_that("filter_matrix throws an error if counts is not provided", {
  expect_error(filter_matrix())
})

test_that("filter_matrix throws an error if counts is not numeric", {
  expect_error(filter_matrix(counts = letters[1:5]), "counts must be a numeric matrix")
  expect_error(filter_matrix(counts = matrix(rep(TRUE, 12), nrow = 3)), "counts must be a numeric matrix")
})

test_that("filter_matrix throws an error if min_count is not an scalar integer", {
  nm <- nm_matrix(5, 4)
  expect_error(filter_matrix(counts = nm,  min_count = "a"),
               "min_count must be an scalar integer")
  expect_error(filter_matrix(counts = nm, runif(3)),
               "min_count must be an scalar integer")
})

test_that("filter_matrix throws an error if min_total_count is not an scalar integer", {
  nm <- nm_matrix(5, 4)
  expect_error(filter_matrix(counts = nm,  min_total_count = "a"),
               "min_total_count must be an scalar integer")
  expect_error(filter_matrix(counts = nm, min_total_count = runif(3)),
               "min_total_count must be an scalar integer")
})

test_that("filter_matrix throws warning when missing genes found", {
  n <- sample.int(20, size = 1)
  genes <- paste0("GENE", 1:n)
  gm <- mock(sample(rownames(sample_expr), size = 50))

  with_mocked_bindings({
    expect_warning(filter_matrix(sample_expr, genes_keep = genes),
                   paste0(length(genes), "/", length(genes), " genes_keep were not found in the matrix"))
  }, expressed_genes = gm)
})

test_that("filter_matrix throws an error if norm_counts is not numeric", {
  expect_error(filter_matrix(counts = sample_expr, norm_counts = letters[1:5]), "norm_counts must be a numeric matrix")
  expect_error(filter_matrix(counts = sample_expr, norm_counts = matrix(rep(TRUE, 12), nrow = 3)), "norm_counts must be a numeric matrix")
})

test_that("filter_matrix throws an error if genes_keep is not character", {
  expect_error(filter_matrix(counts = sample_expr, genes_keep = runif(10)), "genes_keep must be a character vector")
  expect_error(filter_matrix(counts = sample_expr, genes_keep = rep(TRUE, 12)), "genes_keep must be a character vector")
})

test_that("filter_matrix filters matrices correctly warning when missing genes found", {
  n <- sample.int(20, size = 1) +1
  genes <- paste0("GENE", 1:n)
  gm <- mock(sample(rownames(sample_expr), size = 50))

  with_mocked_bindings({
    expect_warning(filter_matrix(sample_expr, genes_keep = genes),
                   paste0(length(genes), "/", length(genes), " genes_keep were not found in the matrix"))
  }, expressed_genes = gm)
})

test_that("filter_matrix throws error if counts and norm_counts have different dimensions", {
  counts <- nm_matrix(5, 4)
  norm_counts <- nm_matrix(3, 8)
  expect_error(filter_matrix(counts = counts, norm_counts = norm_counts), "counts and norm_counts should have the same dimensions")
})

test_that("filter_matrix throws error if counts and norm_counts have different dimension names", {
  counts <- named_nm_matrix(5, 8)
  norm_counts <- named_nm_matrix(5, 8)
  genes <- paste0("GENE", 1:5)
  rownames(counts) <- genes
  rownames(norm_counts) <- c("GENEA", sample(genes)[1:4])
  expect_error(filter_matrix(counts = counts, norm_counts = norm_counts), "counts and norm_counts should have the same rownames")
})

test_that("filter_matrix filters matrices correctly", {
  n <- sample.int(20, size = 1) +1
  genes <- paste0("GENE", 1:n)
  gm <- mock(genes, cycle = TRUE)
  nn <- sample.int(20, size = 1) +1
  genesn <-  paste0("GENE11", 1:nn)

  counts <- named_nm_matrix(n+nn, 10)
  rownames(counts) <- c(genes, genesn)
  norm_counts <- named_nm_matrix(n+nn, 10)
  rownames(norm_counts) <- c(genes, genesn)
  with_mocked_bindings({
    filtered <- filter_matrix(counts, norm_counts = norm_counts, genes_keep = genesn)
    expect_equal(filtered, list(counts = counts, norm_counts = norm_counts))
  }, expressed_genes = gm, )
})

test_that("filter_matrix filters count correctly", {
  n <- sample.int(20, size = 1) +1
  genes <- paste0("GENE", 1:n)
  gm <- mock(genes, cycle = TRUE)
  nn <- sample.int(20, size = 1) +1
  genesn <-  paste0("GENE11", 1:nn)

  counts <- named_nm_matrix(n+nn, 10)
  rownames(counts) <- c(genes, genesn)
  with_mocked_bindings({
    filtered <- filter_matrix(counts, genes_keep = genesn)
    expect_equal(filtered, list(counts = counts, norm_counts = NULL))
  }, expressed_genes = gm, )
})


