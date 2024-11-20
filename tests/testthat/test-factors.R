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

