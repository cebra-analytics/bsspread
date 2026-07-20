#' Shared worker state for parallel replicate execution
#'
#' Populated on the parent before \code{makeCluster(PSOCK)}. Socket workers
#' repopulate terra-backed state via \code{worker_init}.
#'
#' @noRd
manage_parallel_worker_state <- new.env(parent = emptyenv())

#' Active worker state (reads \code{.GlobalEnv} copy after clusterExport)
#'
#' @noRd
parallel_worker_state <- function() {
  if (exists("manage_parallel_worker_state", envir = .GlobalEnv,
             inherits = FALSE)) {
    get("manage_parallel_worker_state", envir = .GlobalEnv)
  } else {
    manage_parallel_worker_state
  }
}

#' Scalar-only sim_env for PSOCK export (terra objects cannot be serialised)
#'
#' @noRd
parallel_lean_sim_env <- function(sim_env) {
  lean <- new.env(parent = emptyenv())
  for (name in c(
    "time_steps", "step_duration", "step_units", "collation_steps",
    "replicates", "result_stages"
  )) {
    lean[[name]] <- sim_env[[name]]
  }
  lean
}

#' Prepare worker_state before cluster creation
#'
#' @noRd
parallel_prepare_worker_state <- function(sim_env) {
  manage_parallel_worker_state$timestep_callback <- sim_env$timestep_callback
  manage_parallel_worker_state$dispersal_callback <- sim_env$dispersal_callback
  manage_parallel_worker_state$sim_env <- parallel_lean_sim_env(sim_env)
  invisible(NULL)
}

#' Restore callbacks on worker sim_env (when set on worker_state)
#'
#' Rebuild logging helpers on the worker; never attach serialised parent
#' callbacks (nested closures break after PSOCK export).
#'
#' @noRd
parallel_attach_worker_callbacks <- function(sim_env) {
  worker_state <- parallel_worker_state()
  wants_timestep <- !is.null(worker_state$timestep_callback)
  wants_dispersal <- !is.null(worker_state$dispersal_callback)
  if (wants_timestep || wants_dispersal) {
    if (exists("parallel_psock_worker_ensure_timestep_logging",
               mode = "function", inherits = TRUE)) {
      parallel_psock_worker_ensure_timestep_logging()
    }
  }
  sim_env$timestep_callback <- NULL
  sim_env$dispersal_callback <- NULL
  if (wants_timestep &&
      exists("log_timestep", mode = "function", inherits = TRUE)) {
    sim_env$timestep_callback <- get("log_timestep", inherits = TRUE)
  }
  if (wants_dispersal &&
      exists("log_dispersal_model", mode = "function", inherits = TRUE)) {
    sim_env$dispersal_callback <- get("log_dispersal_model", inherits = TRUE)
  }
  invisible(sim_env)
}

#' Force region and dispersal models to serial inner parallelism
#'
#' @noRd
force_serial_inner_parallel <- function(sim_env) {
  if (is.function(sim_env$region$set_cores)) {
    sim_env$region$set_cores(1)
  }
  if (length(sim_env$dispersal_models)) {
    for (i in seq_along(sim_env$dispersal_models)) {
      if (is.function(sim_env$dispersal_models[[i]]$set_cores)) {
        sim_env$dispersal_models[[i]]$set_cores(1)
      }
    }
  }
  invisible(NULL)
}

#' Elapsed pool wall time and average seconds per merged replicate
#'
#' \code{rep_wall_start} is set after cluster init and first-batch dispatch,
#' immediately before the parent blocks on \code{recvOneResult}.
#'
#' @noRd
parallel_replicate_timing <- function(rep_wall_start, reps_merged) {
  if (is.null(rep_wall_start) || is.na(reps_merged) || reps_merged <= 0L) {
    return(list(wall_s = NA_real_, avg_s_per_rep = NA_real_))
  }
  wall_s <- as.numeric(difftime(Sys.time(), rep_wall_start, units = "secs"))
  list(
    wall_s = wall_s,
    avg_s_per_rep = wall_s / reps_merged
  )
}

#' Invoke \code{parallel_merge_callback} with optional pool timing
#'
#' \code{phase} is one of \code{before_pool}, \code{pool_ready},
#' \code{received r=N}, \code{merged r=N}, or \code{after_merge}.
#'
#' @noRd
parallel_invoke_merge_callback <- function(parallel_merge_callback,
                                           phase,
                                           sim_env,
                                           results,
                                           reps_merged,
                                           reps_total,
                                           rep_wall_start = NULL,
                                           rep_outputs = NULL,
                                           out = NULL) {
  if (!is.function(parallel_merge_callback)) {
    return(invisible(NULL))
  }
  timing <- parallel_replicate_timing(rep_wall_start, reps_merged)
  parallel_merge_callback(
    phase = phase,
    sim_env = sim_env,
    results = results,
    reps_merged = reps_merged,
    reps_total = reps_total,
    wall_s = timing$wall_s,
    avg_s_per_rep = timing$avg_s_per_rep,
    rep_outputs = rep_outputs,
    out = out
  )
  invisible(NULL)
}

#' Merge deferred replicate collations into parent results
#'
#' @noRd
merge_replicate_collations <- function(out,
                                       results,
                                       sim_env,
                                       parallel_merge_callback = NULL,
                                       reps_merged = 0L,
                                       reps_total = NA_integer_,
                                       rep_wall_start = NULL) {
  parallel_invoke_merge_callback(
    parallel_merge_callback = parallel_merge_callback,
    phase = sprintf("received r=%d", out$r),
    sim_env = sim_env,
    results = results,
    reps_merged = reps_merged,
    reps_total = reps_total,
    rep_wall_start = rep_wall_start,
    rep_outputs = list(out),
    out = out
  )
  # collate()'s first argument is the number of replicates merged so far
  # (Welford count), not the replicate index. Completion order varies in parallel.
  merge_count <- reps_merged + 1L
  for (col in out$collations) {
    results$collate(merge_count, col$tm, col$n)
  }
  reps_merged <- reps_merged + 1L
  parallel_invoke_merge_callback(
    parallel_merge_callback = parallel_merge_callback,
    phase = sprintf("merged r=%d", out$r),
    sim_env = sim_env,
    results = results,
    reps_merged = reps_merged,
    reps_total = reps_total,
    rep_wall_start = rep_wall_start,
    rep_outputs = NULL,
    out = out
  )
  reps_merged
}

#' Format a remote worker error for the parent process
#'
#' @noRd
parallel_worker_error_message <- function(x) {
  if (!inherits(x, "try-error")) {
    return(as.character(x))
  }
  cond <- attr(x, "condition", exact = TRUE)
  if (!is.null(cond) && inherits(cond, "condition")) {
    return(conditionMessage(cond))
  }
  ch <- as.character(x)
  if (length(ch) && nzchar(ch[[1L]])) {
    return(paste(ch, collapse = "\n"))
  }
  "unknown worker error (empty try-error)"
}

#' Worker entry point for parallel replicate execution
#'
#' @noRd
manage_parallel_worker <- function(r) {
  worker_state <- parallel_worker_state()
  sim_env <- worker_state$sim_env
  parallel_attach_worker_callbacks(sim_env)
  random_seed <- worker_state$random_seed
  per_replicate_seed <- worker_state$per_replicate_seed
  if (isTRUE(per_replicate_seed) && !is.null(random_seed)) {
    rep_seed <- random_seed + as.integer(r) - 1L
    set.seed(rep_seed)
    message(sprintf(
      "Parallel replicate %d: set.seed(%d) (pid %d)",
      r, rep_seed, Sys.getpid()
    ))
  } else {
    message(sprintf(
      "Parallel replicate %d: (pid %d)",
      r, Sys.getpid()
    ))
  }
  force_serial_inner_parallel(sim_env)
  tryCatch(
    run_one_replicate(r, sim_env, defer_collate = TRUE),
    error = function(e) {
      stop(sprintf(
        "replicate %d on pid %d: %s",
        r, Sys.getpid(), conditionMessage(e)
      ), call. = FALSE)
    }
  )
}

#' Initialise a persistent PSOCK cluster for replicate execution
#'
#' @noRd
parallel_init_cluster <- function(n_workers,
                                  worker_init = NULL,
                                  psock_exports = NULL,
                                  psock_export_envir = .GlobalEnv) {
  if (!is.function(worker_init)) {
    stop(
      paste(
        "PSOCK parallel replicates require worker_init(sim_env);",
        "terra-backed objects cannot be serialised to socket workers."
      ),
      call. = FALSE
    )
  }

  cl <- parallel::makeCluster(n_workers, type = "PSOCK", outfile = "")

  assign(".parallel_worker_init", worker_init, envir = psock_export_envir)
  on.exit(
    rm(list = ".parallel_worker_init", envir = psock_export_envir),
    add = TRUE
  )
  export_vars <- c(
    "manage_parallel_worker_state",
    "parallel_worker_state",
    "manage_parallel_worker",
    "run_one_replicate",
    "force_serial_inner_parallel",
    "parallel_attach_worker_callbacks"
  )
  parallel::clusterExport(cl, export_vars, envir = environment())
  parallel::clusterExport(cl, ".parallel_worker_init",
                          envir = psock_export_envir)

  psock_vars <- unique(psock_exports)
  if (length(psock_vars)) {
    parallel::clusterExport(cl, psock_vars, envir = psock_export_envir)
  }

  parallel::clusterEvalQ(cl, {
    worker_state <- parallel_worker_state()
    if (exists("parallel_psock_sanitize_globals", mode = "function",
               inherits = TRUE)) {
      parallel_psock_sanitize_globals()
    }
    .parallel_worker_init(worker_state$sim_env)
    parallel_attach_worker_callbacks(worker_state$sim_env)
    force_serial_inner_parallel(worker_state$sim_env)
  })

  pids <- tryCatch(
    parallel::clusterCall(cl, function() Sys.getpid()),
    error = function(e) integer()
  )
  if (length(pids)) {
    attr(cl, "worker_pids") <- as.integer(unlist(pids))
  }
  cl
}

#' Stop a parallel cluster; force-kill workers when needed
#'
#' @noRd
parallel_stop_cluster <- function(cl, force = FALSE) {
  if (is.null(cl)) {
    return(invisible(NULL))
  }
  pids <- attr(cl, "worker_pids", exact = TRUE)
  tryCatch(parallel::stopCluster(cl), error = function(e) NULL)
  if (length(pids) && .Platform$OS.type == "unix") {
    pids <- unique(as.integer(pids))
    pids <- pids[is.finite(pids) & pids > 1L]
    if (length(pids)) {
      tryCatch(tools::pskill(pids), error = function(e) NULL)
    }
  }
  invisible(NULL)
}

#' Run replicates on a persistent PSOCK worker pool
#'
#' @noRd
run_parallel_replicates <- function(sim_env,
                                    results = NULL,
                                    results_factory = NULL,
                                    n_workers,
                                    random_seed = NULL,
                                    per_replicate_seed = TRUE,
                                    worker_init = NULL,
                                    psock_exports = NULL,
                                    psock_export_envir = .GlobalEnv,
                                    parallel_merge_callback = NULL) {
  n_workers <- min(as.integer(n_workers), sim_env$replicates)
  replicate_seq <- seq_len(sim_env$replicates)
  reps_total <- length(replicate_seq)

  if (!reps_total) {
    return(list(
      wall_s = 0,
      reps = 0L,
      cores = n_workers,
      backend = "PSOCK persistent"
    ))
  }

  force_serial_inner_parallel(sim_env)

  parallel_invoke_merge_callback(
    parallel_merge_callback = parallel_merge_callback,
    phase = "before_pool",
    sim_env = sim_env,
    results = results,
    reps_merged = 0L,
    reps_total = reps_total,
    rep_wall_start = NULL,
    rep_outputs = NULL,
    out = NULL
  )

  manage_parallel_worker_state$random_seed <- random_seed
  manage_parallel_worker_state$per_replicate_seed <- per_replicate_seed
  parallel_prepare_worker_state(sim_env)

  cl <- parallel_init_cluster(
    n_workers,
    worker_init = worker_init,
    psock_exports = psock_exports,
    psock_export_envir = psock_export_envir
  )

  parallel_invoke_merge_callback(
    parallel_merge_callback = parallel_merge_callback,
    phase = "pool_ready",
    sim_env = sim_env,
    results = results,
    reps_merged = 0L,
    reps_total = reps_total,
    rep_wall_start = NULL,
    rep_outputs = NULL,
    out = NULL
  )

  if (is.null(results)) {
    if (!is.function(results_factory)) {
      parallel_stop_cluster(cl)
      stop(
        "Parallel replicates require results or results_factory.",
        call. = FALSE
      )
    }
    results <- results_factory()
  }

  cleanup_cluster_pool <- function(force = FALSE) {
    parallel_stop_cluster(cl, force = force)
    manage_parallel_worker_state$sim_env <- NULL
    manage_parallel_worker_state$random_seed <- NULL
    manage_parallel_worker_state$per_replicate_seed <- NULL
    manage_parallel_worker_state$timestep_callback <- NULL
    manage_parallel_worker_state$dispersal_callback <- NULL
    invisible(NULL)
  }

  on.exit(cleanup_cluster_pool(force = FALSE), add = TRUE)

  pending_reps <- as.list(replicate_seq)
  active_rep <- rep(NA_integer_, length(cl))
  reps_merged <- 0L

  for (i in seq_along(cl)) {
    if (!length(pending_reps)) {
      break
    }
    r <- pending_reps[[1L]]
    pending_reps <- pending_reps[-1L]
    active_rep[i] <- r
    parallel:::sendCall(cl[[i]], manage_parallel_worker, list(r = r))
  }

  rep_wall_start <- Sys.time()
  tryCatch({
    while (any(!is.na(active_rep))) {
      res <- parallel:::recvOneResult(cl)
      worker_i <- res$node
      if (inherits(res$value, "try-error")) {
        worker_err <- parallel_worker_error_message(res$value)
        stop(sprintf(
          "Parallel replicates: replicate %d FAILED on worker %d:\n%s",
          active_rep[worker_i],
          worker_i,
          worker_err
        ), call. = FALSE)
      }
      reps_merged <- merge_replicate_collations(
        res$value,
        results,
        sim_env,
        parallel_merge_callback = parallel_merge_callback,
        reps_merged = reps_merged,
        reps_total = reps_total,
        rep_wall_start = rep_wall_start
      )
      active_rep[worker_i] <- NA_integer_

      if (length(pending_reps)) {
        r <- pending_reps[[1L]]
        pending_reps <- pending_reps[-1L]
        active_rep[worker_i] <- r
        parallel:::sendCall(cl[[worker_i]], manage_parallel_worker, list(r = r))
      }
    }

    parallel_invoke_merge_callback(
      parallel_merge_callback = parallel_merge_callback,
      phase = "after_merge",
      sim_env = sim_env,
      results = results,
      reps_merged = reps_total,
      reps_total = reps_total,
      rep_wall_start = rep_wall_start,
      rep_outputs = NULL,
      out = NULL
    )

    list(
      wall_s = as.numeric(Sys.time() - rep_wall_start, units = "secs"),
      reps = reps_total,
      cores = n_workers,
      backend = "PSOCK persistent"
    )
  }, interrupt = function(cond) {
    cleanup_cluster_pool(force = TRUE)
    message("Parallel replicates interrupted; workers stopped.", call. = FALSE)
    if (!interactive()) {
      quit(save = "no", status = 130, runLast = FALSE)
    }
    invokeRestart("abort")
  })
}
