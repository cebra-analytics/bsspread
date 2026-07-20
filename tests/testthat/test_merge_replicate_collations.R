context("merge_replicate_collations")

test_that("parallel merge passes Welford merge count, not replicate index", {
  recorded <- list()
  results <- list(
    collate = function(r, tm, n) {
      recorded[[length(recorded) + 1L]] <<- as.integer(r)
    }
  )
  out_r2 <- list(
    r = 2L,
    collations = list(
      list(r = 2L, tm = 0L, n = 1L),
      list(r = 2L, tm = 1L, n = 1L)
    )
  )
  out_r1 <- list(
    r = 1L,
    collations = list(
      list(r = 1L, tm = 0L, n = 2L)
    )
  )

  reps_merged <- merge_replicate_collations(
    out_r2, results, sim_env = NULL, reps_merged = 0L
  )
  expect_equal(reps_merged, 1L)
  expect_equal(unlist(recorded), c(1L, 1L))

  reps_merged <- merge_replicate_collations(
    out_r1, results, sim_env = NULL, reps_merged = reps_merged
  )
  expect_equal(reps_merged, 2L)
  expect_equal(unlist(recorded), c(1L, 1L, 2L))
})
