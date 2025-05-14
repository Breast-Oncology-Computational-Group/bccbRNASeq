test_that("negative_samples_needed throws an error if status not supplied", {
  expect_error(negative_samples_needed())
})

test_that("negative_samples_needed throws an error if status is not a logical vector", {
  expect_error(negative_samples_needed(LETTERS), "status must be a logical vector")
  expect_error(negative_samples_needed(NULL), "status must be a logical vector")
  expect_error(negative_samples_needed(runif(50)), "status must be a logical vector")
})

test_that("negative_samples_needed throws an error if pos_ratio is not an scalar number", {
  expect_error(negative_samples_needed(log_vector(50), LETTERS), "pos_ratio must be numeric")
  expect_error(negative_samples_needed(log_vector(50), TRUE), "pos_ratio must be numeric")
  expect_error(negative_samples_needed(log_vector(50), runif(50)), "pos_ratio must be numeric")
})

test_that("negative_samples_needed throws an error if pos_ratio is not between 0 and 1", {
  expect_error(negative_samples_needed(log_vector(50), 0), "pos_ratio outside \\[0,1\\]")
  expect_error(negative_samples_needed(log_vector(50), 1), "pos_ratio outside \\[0,1\\]")
  expect_error(negative_samples_needed(log_vector(50), as.numeric(sample.int(1, n =100))), "pos_ratio outside \\[0,1\\]")
  expect_error(negative_samples_needed(log_vector(50), -1*runif(1)), "pos_ratio outside \\[0,1\\]")
})


test_that("negative_samples_needed returns a single positive value",  {
  nm <- negative_samples_needed(sample(c(TRUE, FALSE, NA), size = 50, replace = T, prob = c(0.8, 0.15, 0.05)))
  expect_equal(length(nm), 1)
})

test_that("negative_samples_needed returns correct number for negatively unbalanced cohort",  {
  nm <- negative_samples_needed(sample(c(rep(TRUE, 80), rep(FALSE, 20), rep(NA, 15))))
  expect_equal(length(nm), 1)
  expect_equal(nm, 24)
})

test_that("negative_samples_needed returns cero for balanced cohort",  {
  nm <- negative_samples_needed(sample(c(rep(TRUE, 80), rep(FALSE, 20), rep(NA, 15))), 0.8)
  expect_equal(length(nm), 1)
  expect_equal(nm, 0)
})

test_that("negative_samples_needed returns correct number for negatively unbalanced cohort with
          given ratio",  {
  nm <- negative_samples_needed(sample(c(rep(TRUE, 80), rep(FALSE, 20), rep(NA, 15))), 0.1)
  expect_equal(length(nm), 1)
  expect_equal(nm, 700)
})

test_that("negative_samples_needed returns correct number for positively unbalanced cohort with
          given ratio",  {
  nm <- negative_samples_needed(sample(c(rep(TRUE, 80), rep(FALSE, 20), rep(NA, 15))), 0.9)
  expect_equal(length(nm), 1)
  expect_equal(nm, -100)
})

test_that("expand_matrix throws an error if x is not provided", {
  expect_error(expand_matrix(column_names = letters[1:5]))
})

test_that("expand_matrix throws an error if x is not numeric", {
  expect_error(expand_matrix(x = letters[1:5], column_names = letters[1:5]), "x must be a numeric matrix")
})

test_that("expand_matrix throws an error if column_names is not provided", {
  expect_error(expand_matrix(x = nm_matrix(5, 4)))
})

test_that("expand_matrix throws an error if column_names is not a character vector", {
  expect_error(expand_matrix(x = nm_matrix(5, 4), column_names = runif(10)),
               "column_names must be a character vector")
})

test_that("expand_matrix throws an error if column_names are not in colnames(x)", {
  nm <- named_nm_matrix(5, 4)

  expect_error(expand_matrix(x = nm, column_names = c(sample(colnames(nm), 2), "extra")),
               "values in column_names must match columns in x")
})

test_that("expand_matrix throws an error if n is not an scalar integer", {
  nm <- named_nm_matrix(5, 4)
  cnames <-  sample(colnames(nm), 2)
  expect_error(expand_matrix(x = nm, column_names = cnames, NULL),
               "n must be an scalar integer")
  expect_error(expand_matrix(x = nm, column_names = cnames, "a"),
               "n must be an scalar integer")
  expect_error(expand_matrix(x = nm, column_names = cnames, runif(3)),
               "n must be an scalar integer")
})

test_that("expand_matrix adds correct number of columns", {
  nm <- named_nm_matrix(sample.int(100, 1), sample.int(100, 1))
  n_extra <- sample.int(40, 1)
  extracs <-  sample(colnames(nm), sample.int(ncol(nm), 1))
  em <- expand_matrix(x = nm, column_names = extracs, n_extra)
  expect_equal(ncol(em), ncol(nm) + n_extra)
  expect_equal(sum(grepl("_DUP", colnames(em))), n_extra)
})

test_that("expand_matrix adds one column correctly", {
  nm <- named_nm_matrix(sample.int(100, 1), sample.int(100, 1))
  extracs <-  sample(colnames(nm), sample.int(ncol(nm), 1))
  em <- expand_matrix(x = nm, column_names = extracs, 1L)
  expect_equal(ncol(em), ncol(nm) + 1)
  expect_equal(sum(grepl("_DUP", colnames(em))), 1)
})

test_that("expand_matrix adds zero columns correctly", {
  nm <- named_nm_matrix(sample.int(100, 1), sample.int(100, 1))
  em <- expand_matrix(x = nm, column_names = colnames(nm), 0L)
  expect_equal(ncol(em), ncol(nm))
  expect_equal(sum(grepl("_DUP", colnames(em))), 0)
})

test_that("qfactors_median calls negative_samples_needed", {
  ngm <- mock(sample.int(100, 1))
  x <- named_nm_matrix(5, 100)
  er <- log_vector(100, colnames(x))
  with_mocked_bindings(code = {
    qfactors_median(x, er, runif(1), 10, seed = 37, cores = 2)
    expect_called(ngm, 1)
  }, negative_samples_needed = ngm)
})

test_that("iterative_luf calls negative_samples_needed", {
  ngm <- mock(sample.int(100, 1))
  x <- named_nm_matrix(5, 100)
  er <- log_vector(100, colnames(x))
  with_mocked_bindings(code = {
    qfactors_median(x, er, runif(1), 10, seed = 37, cores = 2)
    expect_called(ngm, 1)
  }, negative_samples_needed = ngm)
})

##
test_that("qfactors_median_np calls expand_matrix correct number of times", {
  x <- named_nm_matrix(5, 100)
  er <- log_vector(100, colnames(x))
  nn <- sample.int(100, 1)
  nc <- sample.int(20, 1)
  ngm <- mock(nn)
  em <- mock(x, cycle = TRUE)
  with_mocked_bindings(code = {
    qfactors_median_np(x, er, runif(1), n_iter =  nc, seed = 37)
    expect_called(em, nc)
  }, negative_samples_needed = ngm, expand_matrix = em)
})

##
test_that("qfactors_median_np calls expand_matrix with correct arguments", {
  x <- named_nm_matrix(5, 100)
  er <- log_vector(100, colnames(x))
  nn <- sample.int(100, 1)
  nc <- sample.int(20, 1)
  ngm <- mock(nn)
  em <- mock(x, cycle = TRUE)
  cargs <- lapply(1:nc, function(i) list(x, names(which(!er)), nn))
  with_mocked_bindings(code = {
    qfactors_median_np(x, er, runif(1), n_iter = nc, seed = 37)
    expect_called(em, nc)
    margs <- mock_args(em)
    expect_equal(margs, cargs)
  }, negative_samples_needed = ngm, expand_matrix = em)
})

