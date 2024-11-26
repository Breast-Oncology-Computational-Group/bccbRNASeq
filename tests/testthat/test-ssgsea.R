test_that("bccb_expression_signatures throws error if both signature_set and signature_names are provided", {
  expr <- named_nm_matrix(25, 25)
  expect_error(bccb_expression_signatures(expr, signature_set = "all",
                                          signature_names = c("HALLMARK")),
               "Exactly one of `signature_set` or `signature_names` must be supplied.")
})

test_that("bccb_expression_signatures throws an error if seed is not an scalar integer", {
  expr <- named_nm_matrix(25, 25)
  expect_error(bccb_expression_signatures(expr, seed = "a"),
               "seed must be an scalar integer")
  expect_error(bccb_expression_signatures(expr, seed = runif(3)),
               "seed must be an scalar integer")
})

test_that("bccb_expression_signatures throws error if signature_names have names not in the signature set", {
  expr <- named_nm_matrix(25, 25)
  signatures <- c('RTK_ACT_UP', month.name[10])
  expect_error(bccb_expression_signatures(expr, signature_names = signatures),
               "`signature_names` must include values in the bccb signature names set")
})
