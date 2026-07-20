#' Simulation state for Simulator runs
#'
#' @noRd
manage_sim_env <- function(region,
                           time_steps,
                           step_duration,
                           step_units,
                           collation_steps,
                           replicates,
                           result_stages,
                           initializer,
                           population_model,
                           dispersal_models,
                           impacts,
                           actions,
                           user_function) {
  sim_env <- new.env(parent = emptyenv())
  sim_env$region <- region
  sim_env$time_steps <- time_steps
  sim_env$step_duration <- step_duration
  sim_env$step_units <- step_units
  sim_env$collation_steps <- collation_steps
  sim_env$replicates <- replicates
  sim_env$result_stages <- result_stages
  sim_env$initializer <- initializer
  sim_env$population_model <- population_model
  sim_env$dispersal_models <- dispersal_models
  sim_env$impacts <- impacts
  sim_env$actions <- actions
  sim_env$user_function <- user_function
  sim_env$continued_incursions <- NULL
  sim_env$results <- NULL
  sim_env$timestep_callback <- NULL
  sim_env$dispersal_callback <- NULL
  sim_env
}

#' Continued-incursions closure for a simulation run
#'
#' @noRd
manage_setup_continued_incursions <- function(sim_env) {
  sim_env$continued_incursions <- sim_env$initializer$continued_incursions()
  invisible(sim_env)
}

#' Create Results for a simulation run
#'
#' @noRd
manage_create_results <- function(sim_env) {
  Results(
    sim_env$region,
    sim_env$population_model,
    impacts = sim_env$impacts,
    actions = sim_env$actions,
    time_steps = sim_env$time_steps,
    step_duration = sim_env$step_duration,
    step_units = sim_env$step_units,
    collation_steps = sim_env$collation_steps,
    replicates = sim_env$replicates,
    combine_stages = sim_env$result_stages
  )
}

#' Initialise continued incursions and Results on sim_env
#'
#' @noRd
manage_setup_results <- function(sim_env) {
  manage_setup_continued_incursions(sim_env)
  sim_env$results <- manage_create_results(sim_env)
  invisible(sim_env)
}

#' Run one simulation replicate and collate into sim_env$results
#'
#' When \code{defer_collate = TRUE}, collation steps are returned as
#' \code{list(r = r, collations = ...)} for parallel worker merge on the parent.
#'
#' Optional \code{sim_env$timestep_callback} and \code{sim_env$dispersal_callback}
#' support per-timestep observability (e.g. platform benchmark logging).
#'
#' @noRd
run_one_replicate <- function(r, sim_env, defer_collate = FALSE) {
  region <- sim_env$region
  time_steps <- sim_env$time_steps
  initializer <- sim_env$initializer
  population_model <- sim_env$population_model
  dispersal_models <- sim_env$dispersal_models
  impacts <- sim_env$impacts
  actions <- sim_env$actions
  user_function <- sim_env$user_function
  continued_incursions <- sim_env$continued_incursions
  results <- if (defer_collate) {
    NULL
  } else {
    sim_env$results
  }
  timestep_callback <- sim_env$timestep_callback
  dispersal_callback <- sim_env$dispersal_callback
  collation_buf <- if (defer_collate) {
    # Mutable buffer avoids <<- (unsafe when PSOCK exports pollute .GlobalEnv).
    buf <- new.env(parent = emptyenv())
    buf$collations <- list()
    buf
  } else {
    NULL
  }
  gc_time_prev <- if (is.function(timestep_callback)) {
    sum(gc.time())
  } else {
    NULL
  }

  result_collate <- function(tm, n) {
    if (defer_collate) {
      i <- length(collation_buf$collations) + 1L
      collation_buf$collations[[i]] <- list(
        r = r,
        tm = tm,
        n = n
      )
    } else {
      results$collate(r, tm, n)
    }
  }

  # Initialize population array
  n <- initializer$initialize()

  # Set diffusion attributes when spatially implicit (single patch)
  if (region$spatially_implicit()) {

    # Diffusion model
    if (any(sapply(dispersal_models,
                   function(dm) inherits(dm, "Diffusion")))) {
      idx <- which(sapply(dispersal_models,
                          function(dm) inherits(dm, "Diffusion")))[1]
      attr(n, "initial_n") <- n
      attr(n, "diffusion_rate") <-
        dispersal_models[[idx]]$get_diffusion_rate()
      attr(n, "diffusion_radius") <- 0
    }

    # Area spread model
    if (any(sapply(dispersal_models,
                   function(dm) inherits(dm, "AreaSpread")))) {
      capacity <- population_model$get_capacity()
      capacity_area <- attr(capacity, "area")
      if (population_model$get_type() == "stage_structured") {
        stages <- population_model$get_capacity_stages()
        attr(n, "spread_area") <-
          sum(n[, stages]) / as.numeric(capacity) * capacity_area
      } else { # unstructured
        attr(n, "spread_area") <- n / as.numeric(capacity) * capacity_area
      }
    }
  }

  # Calculate impacts
  if (length(impacts)) {

    # Calculate each impact
    for (i in seq_along(impacts)) {
      n <- impacts[[i]]$calculate(n, 0)
    }

    # Update recovery delays where applicable
    for (i in seq_along(impacts)) {
      n <- impacts[[i]]$update_recovery_delay(n, 0)
    }

    # Apply any dynamically linked impacts to capacity
    population_model$set_capacity_mult(n)
  }

  # Apply actions
  if (length(actions)) {
    for (i in seq_along(actions)) {
      n <- actions[[i]]$apply(n, 0)
    }
  }

  # Initial results (t = 0)
  result_collate(0L, n)

  # Time steps
  for (tm in seq_len(time_steps)) {

    if (is.function(timestep_callback)) {
      t0 <- Sys.time()
    }

    # Population growth
    n <- population_model$grow(n, tm)

    if (is.function(timestep_callback)) {
      t1 <- Sys.time()
    }

    # Dispersal for each spread vector
    if (length(dispersal_models)) {

      # Pack into list of original, remaining and relocated populations
      n <- dispersal_models[[1]]$pack(n)

      # Perform dispersal for each spread vector
      for (i in seq_along(dispersal_models)) {
        if (is.function(dispersal_callback) || is.function(timestep_callback)) {
          t_dm <- Sys.time()
        }
        n <- dispersal_models[[i]]$disperse(n, tm)
        if (is.function(dispersal_callback)) {
          dispersal_callback(
            i,
            as.numeric(Sys.time() - t_dm, units = "secs"),
            r = if (defer_collate) r else NA_integer_
          )
        }
      }

      # Unpack population array from separated list
      n <- dispersal_models[[1]]$unpack(n)
    }

    if (is.function(timestep_callback)) {
      t2 <- Sys.time()
    }

    # Calculate impacts
    if (length(impacts)) {

      # Calculate each impact
      for (i in seq_along(impacts)) {
        n <- impacts[[i]]$calculate(n, tm)
      }

      # Update recovery delays where applicable
      for (i in seq_along(impacts)) {
        n <- impacts[[i]]$update_recovery_delay(n, tm)
      }

      # Apply any dynamically linked impacts to capacity
      population_model$set_capacity_mult(n)
    }

    if (is.function(timestep_callback)) {
      t3 <- Sys.time()
    }

    # Apply actions
    if (length(actions)) {

      # Clear attributes
      for (i in seq_along(actions)) {
        n <- actions[[i]]$clear_attributes(n)
      }

      # Apply sequentially
      for (i in seq_along(actions)) {
        n <- actions[[i]]$apply(n, tm)
      }
    }

    if (is.function(timestep_callback)) {
      t4 <- Sys.time()
    }

    # User-defined function
    if (is.function(user_function)) {
      n_attr <- attributes(n) # get attributes
      if (length(formals(user_function)) == 3) {
        n <- user_function(n, r, tm)
      } else { # previously just n
        n <- user_function(n)
      }
      if (length(n_attr)) {
        for (i in seq_along(n_attr)) { # restore attributes
          if (!names(n_attr[i]) %in% names(attributes(n))) {
            attr(n, names(n_attr[i])) <- n_attr[[i]]
          }
        }
      }
    }

    # Collate results
    result_collate(tm, n)

    # Continued incursions
    if (is.function(continued_incursions)) {
      n <- continued_incursions(tm, n)
    }

    if (is.function(timestep_callback)) {
      t5 <- Sys.time()
      gc_time_prev <- timestep_callback(
        tm = tm,
        r = r,
        t0 = t0,
        t1 = t1,
        t2 = t2,
        t3 = t3,
        t4 = t4,
        t5 = t5,
        n = n,
        gc_time_prev = gc_time_prev,
        collations = if (defer_collate) collation_buf$collations else NULL
      )
    }

  } # time steps

  if (defer_collate) {
    return(list(r = r, collations = collation_buf$collations))
  }
  invisible(NULL)
}
