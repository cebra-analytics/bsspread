context("Welford collation")

# Mirrors Results$collate() / finalize() for scalar summaries.
results_welford_sequential <- function(x) {
  n <- 0L
  mean <- 0
  m2 <- 0
  for (xi in x) {
    n <- n + 1L
    prev_mean <- mean
    mean <- prev_mean + (xi - prev_mean) / n
    m2 <- m2 + (xi - prev_mean) * (xi - mean)
  }
  list(
    mean = mean,
    sd = if (n > 1L) sqrt(m2 / (n - 1L)) else NA_real_
  )
}

# Process batches sequentially but update one sample at a time with global r.
# Mirrors parallel replicate merge: replicates may finish in any order/batching,
# but each is still folded into the running summary with n_B = 1.
welford_global_r_batches <- function(batches) {
  n_done <- 0L
  mean <- 0
  m2 <- 0
  for (b in batches) {
    for (xi in b) {
      n_done <- n_done + 1L
      prev_mean <- mean
      mean <- prev_mean + (xi - prev_mean) / n_done
      m2 <- m2 + (xi - prev_mean) * (xi - mean)
    }
  }
  list(
    mean = mean,
    sd = if (n_done > 1L) sqrt(m2 / (n_done - 1L)) else NA_real_
  )
}

expect_close_stats <- function(obs, ref_mean, ref_sd, tol = 1e-12) {
  expect_equal(obs$mean, ref_mean, tolerance = tol)
  expect_equal(obs$sd, ref_sd, tolerance = tol)
}

test_that("online Welford matches reference mean and sd", {
  set.seed(42)
  x <- runif(1000)
  ref_mean <- mean(x)
  ref_sd <- sd(x)
  grp <- sample(LETTERS, 1000, replace = TRUE)
  x_batched <- split(x, grp)

  expect_close_stats(results_welford_sequential(x), ref_mean, ref_sd)
  expect_close_stats(welford_global_r_batches(x_batched), ref_mean, ref_sd)
})

test_that("one-at-a-time merge is invariant to arrival order", {
  set.seed(1)
  x <- rnorm(500)
  ref_mean <- mean(x)
  ref_sd <- sd(x)

  by_letter <- split(x, sample(letters, length(x), replace = TRUE))
  by_index <- split(x, (seq_along(x) - 1L) %% 16L)

  expect_close_stats(welford_global_r_batches(by_letter), ref_mean, ref_sd)
  expect_close_stats(welford_global_r_batches(by_index), ref_mean, ref_sd)
})

test_that("online Welford matches reference for matrix stacks", {
  set.seed(7)
  n_cells <- 50L
  n_reps <- 32L
  mat <- matrix(rnorm(n_cells * n_reps), nrow = n_cells, ncol = n_reps)
  ref_mean <- rowMeans(mat)
  ref_sd <- apply(mat, 1, sd)

  # Fold columns in shuffled order, one replicate at a time (n_B = 1).
  col_order <- sample(seq_len(n_reps))
  merged_mean <- numeric(n_cells)
  merged_m2 <- numeric(n_cells)
  merged_n <- 0L

  for (j in col_order) {
    x <- mat[, j]
    merged_n <- merged_n + 1L
    prev_mean <- merged_mean
    merged_mean <- prev_mean + (x - prev_mean) / merged_n
    merged_m2 <- merged_m2 + (x - prev_mean) * (x - merged_mean)
  }

  merged_sd <- sqrt(merged_m2 / (merged_n - 1L))
  expect_equal(merged_mean, ref_mean, tolerance = 1e-12)
  expect_equal(merged_sd, ref_sd, tolerance = 1e-12)
})
