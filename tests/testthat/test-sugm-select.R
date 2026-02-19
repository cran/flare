make_sugm_est <- function(seed = 1) {
  set.seed(seed)
  x <- matrix(rnorm(96), nrow = 24, ncol = 4)
  sugm(x, nlambda = 3, method = "clime", verbose = FALSE, standardize = FALSE, max.ite = 2000)
}

test_that("sugm.select cv stores selected loss metadata", {
  est <- make_sugm_est(11)

  sel <- sugm.select(est, criterion = "cv", fold = 3, verbose = FALSE)

  expect_s3_class(sel, "select")
  expect_identical(sel$criterion, "cv")
  expect_true(!is.null(sel$loss))
  expect_true(sel$loss %in% c("sugm.likelihood", "sugm.tracel2"))
})

test_that("sugm.select stars falls back to largest lambda when threshold is not reached", {
  est <- make_sugm_est(12)

  sel <- sugm.select(est, criterion = "stars", rep.num = 3, stars.thresh = 1, verbose = FALSE)

  expect_s3_class(sel, "select")
  expect_equal(sel$opt.index, est$nlambda)
})

test_that("sugm.select validates stars hyper-parameters", {
  est <- make_sugm_est(13)

  expect_error(
    sugm.select(est, criterion = "stars", rep.num = 0, verbose = FALSE),
    "rep.num"
  )
  expect_error(
    sugm.select(est, criterion = "stars", rep.num = 2, stars.thresh = 1.1, verbose = FALSE),
    "stars.thresh"
  )
  expect_error(
    sugm.select(est, criterion = "stars", rep.num = 2, stars.subsample.ratio = 1, verbose = FALSE),
    "stars.subsample.ratio"
  )
})

test_that("sugm rejects non-positive or non-finite lambda values", {
  set.seed(14)
  x <- matrix(rnorm(96), nrow = 24, ncol = 4)

  est <- NULL
  utils::capture.output({
    est <- sugm(x, lambda = c(0.2, 0, 0.1), method = "clime", verbose = FALSE)
  })
  expect_null(est)
})
