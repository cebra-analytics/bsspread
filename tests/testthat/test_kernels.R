context("Kernels")

test_that("initializes with kernel generators", {
  expect_silent(kernal_gen <- Kernels())
  expect_is(kernal_gen, "Kernels")
  expect_equal(get("multiplier",
                   envir = environment(kernal_gen$get_beta_function)), 1)
  expect_silent(kernal_gen <- Kernels(multiplier = 0.3))
  expect_equal(get("multiplier",
                   envir = environment(kernal_gen$get_beta_function)), 0.3)
  expect_is(kernal_gen$get_beta_function, "function")
  expect_is(kernal_gen$get_cauchy_function, "function")
  expect_is(kernal_gen$get_exp_function, "function")
  expect_is(kernal_gen$get_gaussian_function, "function")
  expect_is(kernal_gen$get_lognormal_function, "function")
  expect_is(kernal_gen$get_weibull_function, "function")
  expect_is(kernal_gen$get_lookup_function, "function")
})

test_that("generates various kernel functions", {
  kernal_gen <- Kernels(multiplier = 0.3)
  expect_silent(beta_function <-
                  kernal_gen$get_beta_function(alpha = 2, beta = 4, upper = 5))
  expect_is(beta_function, "function")
  expect_equal(beta_function(1:5), 0.3*(2/(2 + 4))*stats::dbeta((1:5)/5, 2, 4))
  expect_silent(cauchy_function <- kernal_gen$get_cauchy_function(scale = 2))
  expect_is(cauchy_function, "function")
  expect_equal(cauchy_function(1:5), 2*0.3*2*stats::dcauchy(1:5, scale = 2))
  expect_silent(exp_function <- kernal_gen$get_exp_function(mean = 2))
  expect_is(exp_function, "function")
  expect_equal(exp_function(1:5), 0.3*2*stats::dexp(1:5, rate = 1/2))
  expect_silent(gaussian_function <- kernal_gen$get_gaussian_function(sd = 2))
  expect_is(gaussian_function, "function")
  expect_equal(gaussian_function(1:5), 2*0.3*2*stats::dnorm(1:5, sd = 2))
  expect_silent(lognormal_function <-
                  kernal_gen$get_lognormal_function(mean = 3, sd = 2))
  expect_is(lognormal_function, "function")
  expect_equal(lognormal_function(1:5),
               0.3*3*stats::dlnorm(1:5, meanlog = log(3^2/sqrt(3^2 + 2^2)),
                                   sdlog = sqrt(log(1 + 2^2/3^2))))
  expect_silent(weibull_function <- kernal_gen$get_weibull_function(shape = 1, scale = 2))
  expect_is(weibull_function, "function")
  expect_equal(weibull_function(1:5), 0.3*2*stats::dweibull(1:5, shape = 1,
                                                            scale = 2))
})

test_that("generates various kernel functions", {
  ## generates table look-up function
  kernal_gen <- Kernels(multiplier = 0.3)
  lookup_table <- data.frame(value = 1:4, pr = 1/(1:4))
  expect_silent(lookup_function <-
                  kernal_gen$get_lookup_function(table = lookup_table))
  expect_is(lookup_function, "function")
  expected <- 0.3*rep(lookup_table$pr, each = 2)
  expected <- c(expected[1], rep(expected, 2), expected[8])
  expected <- rowMeans(array(expected, c(9, 2)))
  expect_equal(lookup_function((1:9)/2), expected)
})
