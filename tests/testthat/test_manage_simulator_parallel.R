context("Simulator parallel replicates")

TEST_INPUTS_DIR <- test_path("test_inputs")

make_spread_simulator <- function(replicates = 4L, parallel_cores = NULL) {
  TEST_DIRECTORY <- TEST_INPUTS_DIR
  template <- terra::rast(file.path(TEST_DIRECTORY, "greater_melb.tif"))
  region <- Region(template)
  initial_n <- rep(0, region$get_locations())
  initial_n[5922] <- 10
  population_model <- UnstructPopulation(region, growth = 1.5)
  initializer <- Initializer(
    initial_n,
    region = region,
    population_model = population_model
  )
  dispersal <- Dispersal(
    region,
    population_model,
    proportion = 1,
    max_distance = 1000
  )
  Simulator(
    region,
    time_steps = 3,
    replicates = replicates,
    parallel_cores = parallel_cores,
    initializer = initializer,
    population_model = population_model,
    dispersal_models = list(dispersal)
  )
}

psock_worker_init <- local({
  test_inputs_dir <- TEST_INPUTS_DIR
  function(sim_env) {
    template <- terra::rast(file.path(test_inputs_dir, "greater_melb.tif"))
    region <- Region(template)
    initial_n <- rep(0, region$get_locations())
    initial_n[5922] <- 10
    population_model <- UnstructPopulation(region, growth = 1.5)
    initializer <- Initializer(
      initial_n,
      region = region,
      population_model = population_model
    )
    dispersal <- Dispersal(
      region,
      population_model,
      proportion = 1,
      max_distance = 1000
    )
    sim_env$region <- region
    sim_env$population_model <- population_model
    sim_env$initializer <- initializer
    sim_env$dispersal_models <- list(dispersal)
    sim_env$impacts <- list()
    sim_env$actions <- list()
    sim_env$user_function <- NULL
    sim_env$continued_incursions <- initializer$continued_incursions()
    invisible(NULL)
  }
})

test_that("parallel PSOCK replicates match serial with per-replicate seeding", {
  skip_if_not(.Platform$OS.type == "unix", "PSOCK cluster requires Unix")
  simulator <- make_spread_simulator(replicates = 4L)
  expect_silent(
    res_serial <- simulator$run(
      random_seed = 100L,
      per_replicate_seed = TRUE
    )
  )
  simulator <- make_spread_simulator(replicates = 4L)
  res_parallel <- suppressMessages(simulator$run(
    parallel_replicates = TRUE,
    replicate_workers = 2L,
    random_seed = 100L,
    per_replicate_seed = TRUE,
    worker_init = psock_worker_init
  ))
  expect_equal(
    res_serial$get_list(),
    res_parallel$get_list()
  )
})

test_that("PSOCK parallel replicates require worker_init", {
  skip_if_not(.Platform$OS.type == "unix", "PSOCK test requires Unix cluster")
  simulator <- make_spread_simulator(replicates = 2L)
  expect_error(
    simulator$run(
      parallel_replicates = TRUE,
      replicate_workers = 2L
    ),
    "worker_init"
  )
})

test_that("replicate_workers defaults to min(parallel_cores, replicates)", {
  skip_if_not(.Platform$OS.type == "unix", "PSOCK cluster requires Unix")
  simulator <- make_spread_simulator(replicates = 6L, parallel_cores = 4L)
  res <- suppressMessages(simulator$run(
    parallel_replicates = TRUE,
    replicate_workers = 4L,
    random_seed = 1L,
    worker_init = psock_worker_init
  ))
  expect_is(res, "Results")
})

test_that("timestep_callback fires once per time step in serial runs", {
  simulator <- make_spread_simulator(replicates = 1L)
  calls <- 0L
  timestep_cb <- function(tm, r, t0, t1, t2, t3, t4, t5,
                          n, gc_time_prev, collations = NULL) {
    calls <<- calls + 1L
    gc_time_prev
  }
  suppressMessages(simulator$run(
    timestep_callback = timestep_cb
  ))
  expect_equal(calls, 3L)
})

test_that("parallel sim_env has no results before workers start", {
  skip_if_not(.Platform$OS.type == "unix", "PSOCK cluster requires Unix")
  simulator <- make_spread_simulator(replicates = 2L)
  merge_cb <- function(phase, sim_env, ...) {
    if (identical(phase, "before_pool") || identical(phase, "pool_ready")) {
      expect_null(sim_env$results)
    }
    invisible(NULL)
  }
  suppressMessages(simulator$run(
    parallel_replicates = TRUE,
    replicate_workers = 2L,
    worker_init = psock_worker_init,
    parallel_merge_callback = merge_cb
  ))
})

test_that("parallel_merge_callback fires for each merged replicate", {
  skip_if_not(.Platform$OS.type == "unix", "PSOCK cluster requires Unix")
  simulator <- make_spread_simulator(replicates = 3L)
  phases <- character()
  merge_cb <- function(phase, sim_env, reps_merged, reps_total,
                       rep_outputs = NULL, out = NULL, ...) {
    phases <<- c(phases, phase)
    invisible(NULL)
  }
  suppressMessages(simulator$run(
    parallel_replicates = TRUE,
    replicate_workers = 2L,
    worker_init = psock_worker_init,
    parallel_merge_callback = merge_cb
  ))
  expect_true("before_pool" %in% phases)
  expect_true("pool_ready" %in% phases)
  expect_true("after_merge" %in% phases)
  expect_lt(match("before_pool", phases), match("pool_ready", phases))
  expect_lt(match("pool_ready", phases),
            match(phases[grepl("^received r=", phases)][1L], phases))
  expect_equal(sum(grepl("^received r=", phases)), 3L)
  expect_equal(sum(grepl("^merged r=", phases)), 3L)
})

test_that("parallel_merge_callback receives pool wall timing on merge", {
  skip_if_not(.Platform$OS.type == "unix", "PSOCK cluster requires Unix")
  simulator <- make_spread_simulator(replicates = 3L)
  merged_timing <- list()
  merge_cb <- function(phase, sim_env, reps_merged, reps_total,
                       wall_s = NA_real_, avg_s_per_rep = NA_real_,
                       rep_outputs = NULL, out = NULL, ...) {
    if (grepl("^merged r=", phase)) {
      merged_timing[[length(merged_timing) + 1L]] <<- list(
        reps_merged = reps_merged,
        wall_s = wall_s,
        avg_s_per_rep = avg_s_per_rep
      )
    }
    invisible(NULL)
  }
  suppressMessages(simulator$run(
    parallel_replicates = TRUE,
    replicate_workers = 2L,
    worker_init = psock_worker_init,
    parallel_merge_callback = merge_cb
  ))
  expect_length(merged_timing, 3L)
  expect_true(all(vapply(merged_timing, function(x) {
    is.finite(x$wall_s) && x$wall_s >= 0
  }, logical(1L))))
  expect_true(all(vapply(merged_timing, function(x) {
    is.finite(x$avg_s_per_rep) && x$avg_s_per_rep > 0
  }, logical(1L))))
  expect_equal(vapply(merged_timing, `[[`, integer(1L), "reps_merged"), 1:3)
})

test_that("parallel runs attach parallel_stats to results", {
  skip_if_not(.Platform$OS.type == "unix", "PSOCK cluster requires Unix")
  simulator <- make_spread_simulator(replicates = 2L)
  res <- suppressMessages(simulator$run(
    parallel_replicates = TRUE,
    replicate_workers = 2L,
    worker_init = psock_worker_init
  ))
  stats <- attr(res, "parallel_stats", exact = TRUE)
  expect_is(stats, "list")
  expect_equal(stats$reps, 2L)
  expect_equal(stats$backend, "PSOCK persistent")
  expect_true(is.numeric(stats$wall_s))
})
