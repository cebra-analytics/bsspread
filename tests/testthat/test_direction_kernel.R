context("DirectionKernel")

test_that("initializes with kernel generators", {
  expect_silent(kernel_gen <- DirectionKernel())
  expect_is(kernel_gen, "DirectionKernel")
  expect_s3_class(kernel_gen, "Kernels")
  expect_is(kernel_gen$get_beta_function, "function")
  expect_is(kernel_gen$get_lookup_function, "function")
})

test_that("generates beta kernel direction function", {
  expect_silent(kernel_gen <- DirectionKernel(multiplier = 0.3))
  expect_silent(beta_function <- kernel_gen$get_beta_function(
    alpha = 2, beta = 4))
  expect_is(beta_function, "function")
  expect_equal(beta_function((0:36)*10),
               0.3*(2/(2 + 4))*stats::dbeta(((0:36)*10)/360, 2, 4))
  expect_error(beta_function <- kernel_gen$get_beta_function(
    alpha = 2, beta = 4, shift = 370),
    "The Beta direction kernel shift parameter must be <= 360 degrees.")
  expect_silent(beta_function <- kernel_gen$get_beta_function(
    alpha = 2, beta = 4, shift = 120))
  expect_equal(beta_function((0:36)*10),
               0.3*(2/(2 + 4))*stats::dbeta(c(24:35, 0:24)*10/360, 2, 4))
  expect_silent(beta_function <- kernel_gen$get_beta_function(
    alpha = 2, beta = 4))
  expect_silent(kernel_gen <- DirectionKernel(
    direction_type = "trigonometric", orientation = "from"))
  expect_silent(beta_function <- kernel_gen$get_beta_function(
    alpha = 2, beta = 4, shift = 0))
  expect_equal(beta_function((0:36)*10),
               (2/(2 + 4))*stats::dbeta(c(18:35, 0:18)*10/360, 2, 4))
  expect_silent(beta_function <- kernel_gen$get_beta_function(
    alpha = 2, beta = 4, shift = 40))
  expect_equal(beta_function((0:36)*10),
               (2/(2 + 4))*stats::dbeta(c(14:35, 0:14)*10/360, 2, 4))
  expect_silent(kernel_gen <- DirectionKernel(
    direction_type = "navigational", orientation = "to"))
  expect_silent(beta_function <- kernel_gen$get_beta_function(
    alpha = 2, beta = 4, shift = 0))
  expect_equal(beta_function((0:36)*10),
               (2/(2 + 4))*stats::dbeta(c(9:0, 35:9)*10/360, 2, 4))
  expect_silent(beta_function <- kernel_gen$get_beta_function(
    alpha = 2, beta = 4, shift = 40))
  expect_equal(beta_function((0:36)*10),
               (2/(2 + 4))*stats::dbeta(c(5:0, 35:5)*10/360, 2, 4))
  expect_silent(kernel_gen <- DirectionKernel(
    direction_type = "navigational", orientation = "from"))
  expect_silent(beta_function <- kernel_gen$get_beta_function(
    alpha = 2, beta = 4, shift = 0))
  expect_equal(beta_function((0:36)*10),
               (2/(2 + 4))*stats::dbeta(c(27:0, 35:27)*10/360, 2, 4))
  expect_silent(beta_function <- kernel_gen$get_beta_function(
    alpha = 2, beta = 4, shift = 40))
  expect_equal(beta_function((0:36)*10),
               (2/(2 + 4))*stats::dbeta(c(23:0, 35:23)*10/360, 2, 4))
})

test_that("generates table look-up direction function", {
  expect_silent(kernel_gen <- DirectionKernel(
    direction_type = "trigonometric", orientation = "to"))
  lookup_table <- data.frame(direction = (0:12)*30,
                             prob = pmax(0, round(sin((0:12)*pi/6), 4)))[-5,]
  lookup_table$prob[5] <- 0
  expect_silent(
    lookup_function <- kernel_gen$get_lookup_function(table = lookup_table))
  expect_is(lookup_function, "function")
  expected <- c(round(sin((0:3)*pi/6), 4), 0.5, rep(0, 8))
  expect_equal(lookup_function((0:12)*30), expected[1:13])
  expect_silent(kernel_gen <- DirectionKernel(
    direction_type = "trigonometric", orientation = "from"))
  expect_silent(
    lookup_function <- kernel_gen$get_lookup_function(table = lookup_table))
  expect_equal(lookup_function((0:12)*30), expected[c(7:12, 1:7)])
  expect_silent(kernel_gen <- DirectionKernel(
    direction_type = "navigational", orientation = "to"))
  expect_silent(
    lookup_function <- kernel_gen$get_lookup_function(table = lookup_table))
  expect_equal(lookup_function((0:12)*30), expected[c(4:1, 12:4)])
  expect_silent(kernel_gen <- DirectionKernel(
    direction_type = "navigational", orientation = "from"))
  expect_silent(
    lookup_function <- kernel_gen$get_lookup_function(table = lookup_table))
  expect_equal(lookup_function((0:12)*30), expected[c(10:1, 12:10)])
})
