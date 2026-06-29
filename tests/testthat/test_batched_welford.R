context("batched Welford")

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

# Chan et al. merge of two groups (n, mean, M2).
merge_welford_groups <- function(n_a, mean_a, m2_a, n_b, mean_b, m2_b) {
  n <- n_a + n_b
  delta <- mean_b - mean_a
  list(
    n = n,
    mean = mean_a + delta * n_b / n,
    m2 = m2_a + m2_b + delta^2 * n_a * n_b / n
  )
}

batch_m2 <- function(x) {
  n <- length(x)
  mu <- mean(x)
  list(n = n, mean = mu, m2 = sum((x - mu)^2))
}

# Merge arbitrary batches (e.g. split(x, grp)) via Chan's formula.
welford_from_batches <- function(batches) {
  state <- list(n = 0L, mean = 0, m2 = 0)
  for (b in batches) {
    if (!length(b)) next
    bs <- batch_m2(b)
    if (state$n == 0L) {
      state <- bs
    } else {
      state <- merge_welford_groups(
        state$n, state$mean, state$m2,
        bs$n, bs$mean, bs$m2
      )
    }
  }
  list(
    mean = state$mean,
    sd = if (state$n > 1L) sqrt(state$m2 / (state$n - 1L)) else NA_real_
  )
}

# Process batches sequentially but update one sample at a time with global r.
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

test_that("batched Chan merge matches reference mean and sd", {
  set.seed(42)
  x <- runif(1000)
  ref_mean <- mean(x)
  ref_sd <- sd(x)
  grp <- sample(LETTERS, 1000, replace = TRUE)
  x_batched <- split(x, grp)

  expect_close_stats(welford_from_batches(x_batched), ref_mean, ref_sd)
  expect_close_stats(results_welford_sequential(x), ref_mean, ref_sd)
  expect_close_stats(welford_global_r_batches(x_batched), ref_mean, ref_sd)
})

test_that("batched merge is invariant to batch partition", {
  set.seed(1)
  x <- rnorm(500)
  ref_mean <- mean(x)
  ref_sd <- sd(x)

  by_letter <- split(x, sample(letters, length(x), replace = TRUE))
  by_index <- split(x, (seq_along(x) - 1L) %% 16L)

  expect_close_stats(welford_from_batches(by_letter), ref_mean, ref_sd)
  expect_close_stats(welford_from_batches(by_index), ref_mean, ref_sd)
})

test_that("vectorised per-batch merge matches reference for matrix stacks", {
  set.seed(7)
  n_cells <- 50L
  n_reps <- 32L
  mat <- matrix(rnorm(n_cells * n_reps), nrow = n_cells, ncol = n_reps)
  ref_mean <- rowMeans(mat)
  ref_sd <- apply(mat, 1, sd)

  batch_size <- 8L
  batch_ids <- split(seq_len(n_reps), ((seq_len(n_reps) - 1L) %/% batch_size))

  merged_mean <- numeric(n_cells)
  merged_m2 <- numeric(n_cells)
  merged_n <- 0L

  for (cols in batch_ids) {
    block <- mat[, cols, drop = FALSE]
    n_b <- length(cols)
    mean_b <- rowMeans(block)
    m2_b <- rowSums((block - mean_b)^2)

    if (merged_n == 0L) {
      merged_n <- n_b
      merged_mean <- mean_b
      merged_m2 <- m2_b
    } else {
      n <- merged_n + n_b
      delta <- mean_b - merged_mean
      merged_m2 <- merged_m2 + m2_b + (delta^2) * merged_n * n_b / n
      merged_mean <- merged_mean + delta * n_b / n
      merged_n <- n
    }
  }

  merged_sd <- sqrt(merged_m2 / (merged_n - 1L))
  expect_equal(merged_mean, ref_mean, tolerance = 1e-12)
  expect_equal(merged_sd, ref_sd, tolerance = 1e-12)
})
