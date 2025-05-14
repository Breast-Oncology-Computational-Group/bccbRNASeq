zeroes <- matrix(rep(0, times = 12), ncol = 4)

test_that("qfactors throws error if x is not a numeric matrix", {
  expect_error(qfactors(NULL), "x must be a numeric matrix")
  expect_error(qfactors(matrix(letters[1:12], ncol = 4)),  "x must be a numeric matrix")
})

test_that("qfactors returns a list with two elements", {
  qf <- qfactors(nm_matrix(sample.int(100, 1), sample.int(100, 1)))
  expect_equal(length(qf), 2)
})

test_that("qfactors returns a list of vectors with correct number of elements", {
  n = sample.int(100, 1)
  qf <- qfactors(nm_matrix(n, sample.int(100, 1)))
  expect_equal(length(qf[[1]]), n)
  expect_equal(length(qf[[2]]), n)
})

test_that("sfactors throws error if x is not a numeric matrix", {
  expect_error(sfactors(NULL), "x must be a numeric matrix")
  expect_error(sfactors(matrix(letters[1:12], ncol = 4)),  "x must be a numeric matrix")
})

test_that("sfactors throws error if qs is not a numeric vector", {
  expect_error(sfactors(nm_matrix(20, 30), LETTERS), "qs must be a numeric vector")
  expect_error(sfactors(nm_matrix(20, 30), nm_matrix(20, 30)),  "qs must be a numeric vector")
  expect_error(sfactors(nm_matrix(20, 30), log_vector(50)),  "qs must be a numeric vector")
})

test_that("sfactors throws error if x and qs don't share rownames/names", {
  qs <- runif(10)
  expect_error(sfactors(nm_matrix(20, 30), qs), "qs must have names")
  names(qs) <- paste0("x", 1:10)
  expect_error(sfactors(nm_matrix(20, 30), qs), "rownames in x and names in qs must match")
})

test_that("sfactors returns a list with same number of elements as qs", {
  nrs <- sample.int(100, 1)
  nm <- nm_matrix(nrs, 20)
  qs <- runif(nrs)
  names(qs) <- rownames(nm)
  sfs <- sfactors(nm, qs)
  expect_equal(length(sfs), nrs)
  expect_equal(names(sfs), rownames(nm))
})

test_that("lufactors throws error if x is not a numeric matrix", {
  expect_error(lufactors(NULL), "x must be a numeric matrix")
  expect_error(lufactors(matrix(letters[1:12], ncol = 4)),  "x must be a numeric matrix")
})

test_that("lufactors returns a list with four elements", {
  luf <- lufactors(nm_matrix(sample.int(100, 1), sample.int(100, 1)))
  expect_equal(length(luf), 4)
})

test_that("lufactors returns a list of vectors with correct number of elements", {
  n = sample.int(100, 1)
  luf <- lufactors(nm_matrix(n, sample.int(100, 1)))
  expect_equal(length(luf[[1]]), n)
  expect_equal(length(luf[[2]]), n)
})

