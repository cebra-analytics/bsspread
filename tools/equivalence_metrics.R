# Helper utilities for cross-branch behavioural equivalence checks.
#
# Designed for manual workflows where branch checkout/install/run is done
# outside this script (e.g. inside a Docker container with mounted volume).

eq_init <- function(run_label,
                    git_sha = NA_character_,
                    branch = NA_character_,
                    scenario = NA_character_,
                    out_dir = ".",
                    save_occupancy = TRUE) {
  if (!dir.exists(out_dir)) {
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  }
  list(
    run_label = as.character(run_label),
    git_sha = as.character(git_sha),
    branch = as.character(branch),
    scenario = as.character(scenario),
    out_dir = normalizePath(out_dir, winslash = "/", mustWork = FALSE),
    save_occupancy = isTRUE(save_occupancy),
    metrics = list(),
    occupancy = list(),
    prev_occ = NULL
  )
}

eq_total_pop <- function(n) {
  if (is.logical(n)) {
    return(sum(n, na.rm = TRUE))
  }
  sum(as.matrix(n), na.rm = TRUE)
}

eq_occ_indices <- function(n) {
  n_mat <- as.matrix(n)
  if (ncol(n_mat) == 1L) {
    which(n_mat[, 1L] > 0)
  } else {
    which(rowSums(n_mat > 0, na.rm = TRUE) > 0)
  }
}

eq_capture_step <- function(state,
                            n,
                            timestep,
                            replicate = 1L,
                            grow_s = NA_real_,
                            dispersal_s = NA_real_,
                            impacts_s = NA_real_,
                            actions_s = NA_real_,
                            other_s = NA_real_,
                            total_s = NA_real_,
                            rng_check = NA_character_) {
  stopifnot(is.list(state), !is.null(state$metrics))

  occ <- eq_occ_indices(n)
  occ_n <- length(occ)
  new_occ_n <- if (is.null(state$prev_occ)) {
    occ_n
  } else {
    length(setdiff(occ, state$prev_occ))
  }

  rec <- data.frame(
    run_label = state$run_label,
    git_sha = state$git_sha,
    branch = state$branch,
    scenario = state$scenario,
    replicate = as.integer(replicate),
    timestep = as.integer(timestep),
    occupied_n = as.integer(occ_n),
    new_occupied_n = as.integer(new_occ_n),
    total_pop = as.numeric(eq_total_pop(n)),
    total_area_occupied = NA_real_,
    grow_s = as.numeric(grow_s),
    dispersal_s = as.numeric(dispersal_s),
    impacts_s = as.numeric(impacts_s),
    actions_s = as.numeric(actions_s),
    other_s = as.numeric(other_s),
    total_s = as.numeric(total_s),
    rng_check = as.character(rng_check),
    stringsAsFactors = FALSE
  )

  state$metrics[[length(state$metrics) + 1L]] <- rec
  if (isTRUE(state$save_occupancy)) {
    key <- paste0("r", as.integer(replicate), "_t", as.integer(timestep))
    state$occupancy[[key]] <- occ
  }
  state$prev_occ <- occ
  state
}

eq_write_outputs <- function(state) {
  stopifnot(is.list(state), !is.null(state$metrics), !is.null(state$out_dir))
  if (!length(state$metrics)) {
    stop("No metrics captured in state.", call. = FALSE)
  }

  metrics_df <- do.call(rbind, state$metrics)
  metrics_path <- file.path(state$out_dir, paste0(state$run_label, "_metrics.csv"))
  utils::write.csv(metrics_df, metrics_path, row.names = FALSE)

  occ_path <- file.path(state$out_dir, paste0(state$run_label, "_occupancy_sets.rds"))
  if (isTRUE(state$save_occupancy)) {
    saveRDS(state$occupancy, occ_path)
  }

  list(
    metrics_path = metrics_path,
    occupancy_path = if (isTRUE(state$save_occupancy)) occ_path else NA_character_
  )
}
