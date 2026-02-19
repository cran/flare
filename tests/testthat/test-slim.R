test_that("slim handles vector y with missing values", {
  set.seed(1)
  x <- matrix(rnorm(80), nrow = 20, ncol = 4)
  y <- rnorm(20)
  y[c(3, 7)] <- NA_real_

  fit <- NULL
  utils::capture.output({
    fit <- slim(x, y, nlambda = 2, method = "lasso", verbose = FALSE, max.ite = 2000)
  })

  expect_s3_class(fit, "slim")
  expect_equal(nrow(fit$X), 18L)
  expect_equal(nrow(fit$Y), 18L)
  expect_equal(ncol(fit$beta), fit$nlambda)
  expect_false(anyNA(fit$X))
  expect_false(anyNA(fit$Y))
})

test_that("predict.slim validates input shape and returns prediction matrix", {
  set.seed(2)
  x <- matrix(rnorm(120), nrow = 30, ncol = 4)
  y <- rnorm(30)

  fit <- slim(x, y, nlambda = 3, method = "lasso", verbose = FALSE, max.ite = 2000)

  bad_newdata <- matrix(rnorm(30), nrow = 10, ncol = 3)
  bad_pred <- NULL
  utils::capture.output({
    bad_pred <- predict.slim(fit, bad_newdata, lambda.idx = 1, Y.pred.idx = 1)
  })
  expect_null(bad_pred)

  good_pred <- NULL
  utils::capture.output({
    good_pred <- predict.slim(fit, x, lambda.idx = 1:2, Y.pred.idx = 1:2)
  })
  expect_type(good_pred, "list")
  expect_equal(dim(good_pred$Y.pred), c(nrow(x), 2L))
})

test_that("slim rejects non-positive or non-finite lambda values", {
  set.seed(3)
  x <- matrix(rnorm(60), nrow = 15, ncol = 4)
  y <- rnorm(15)

  fit <- NULL
  utils::capture.output({
    fit <- slim(x, y, lambda = c(0.2, 0, 0.1), method = "lasso", verbose = FALSE)
  })
  expect_null(fit)
})
