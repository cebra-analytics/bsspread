# Helpers for step-by-step path calculation equivalence checks.
#
# Normalises old nested-list storage and new data.table storage into one
# canonical long format, then compares checkpoints captured across branches.

pe_git_info <- function() {
  sha <- tryCatch(
    system("git rev-parse HEAD", intern = TRUE),
    error = function(e) NA_character_
  )
  branch <- tryCatch(
    system("git rev-parse --abbrev-ref HEAD", intern = TRUE),
    error = function(e) NA_character_
  )
  list(sha = sha[[1L]], branch = branch[[1L]])
}

pe_find_test_inputs <- function() {
  candidates <- c(
    file.path(getwd(), "tests/testthat/test_inputs"),
    file.path(getwd(), "tests", "testthat", "test_inputs"),
    Sys.getenv("BSspread_TEST_INPUTS", unset = "")
  )
  candidates <- candidates[nzchar(candidates)]
  for (p in candidates) {
    if (dir.exists(p)) {
      return(normalizePath(p, winslash = "/", mustWork = FALSE))
    }
  }
  stop("Could not find tests/testthat/test_inputs.", call. = FALSE)
}

pe_as_int <- function(x) {
  as.integer(x)
}

pe_origin_keys <- function(x) {
  sort(unique(as.integer(x)))
}

pe_tier_names <- function(flat = FALSE) {
  if (flat) {
    return("cell")
  }
  c("cell", "aggr")
}

pe_empty_canonical <- function() {
  data.frame(
    origin = integer(0L),
    tier = character(0L),
    dest = integer(0L),
    distance = integer(0L),
    direction = integer(0L),
    stringsAsFactors = FALSE
  )
}

pe_add_perm_columns <- function(df, perm_values, perm_id) {
  if (is.null(perm_values) || !length(perm_values)) {
    return(df)
  }
  col <- if (perm_id <= 1L) "perm_dist" else paste0("perm_dist", perm_id)
  df[[col]] <- pe_as_int(perm_values)
  df
}

pe_canonical_from_legacy_paths <- function(paths_obj, flat = FALSE) {
  if (is.null(paths_obj$idx) || !length(paths_obj$idx)) {
    return(pe_empty_canonical())
  }

  origins <- pe_origin_keys(names(paths_obj$idx))
  rows <- list()
  n_perms <- if (is.list(paths_obj$perms)) length(paths_obj$perms) else 0L
  if (!n_perms && is.list(paths_obj$perm_dist) && length(paths_obj$perm_dist)) {
    sample <- paths_obj$perm_dist[[1L]]
    if (is.list(sample) && is.list(sample$cell)) {
      n_perms <- length(sample$cell)
    } else if (is.list(sample)) {
      n_perms <- length(sample)
    } else {
      n_perms <- 1L
    }
  }

  for (origin in origins) {
    key <- as.character(origin)
    idx_entry <- paths_obj$idx[[key]]
    dist_entry <- paths_obj$distances[[key]]
    dir_entry <- if (is.list(paths_obj$directions)) paths_obj$directions[[key]] else NULL
    perm_entry <- if (is.list(paths_obj$perm_dist)) paths_obj$perm_dist[[key]] else NULL

    tiers <- pe_tier_names(flat = flat)
    for (tier in tiers) {
      dest <- if (flat) {
        idx_entry
      } else {
        idx_entry[[tier]]
      }
      if (is.null(dest) || !length(dest)) {
        next
      }
      distance <- if (flat) dist_entry else dist_entry[[tier]]
      direction <- if (is.null(dir_entry)) {
        rep(NA_integer_, length(dest))
      } else if (flat) {
        dir_entry
      } else {
        dir_entry[[tier]]
      }
      chunk <- data.frame(
        origin = rep.int(origin, length(dest)),
        tier = rep_len(tier, length(dest)),
        dest = pe_as_int(dest),
        distance = pe_as_int(distance),
        direction = pe_as_int(direction),
        stringsAsFactors = FALSE
      )
      if (!is.null(perm_entry)) {
        if (flat) {
          if (is.list(perm_entry)) {
            for (perm_id in seq_along(perm_entry)) {
              chunk <- pe_add_perm_columns(chunk, perm_entry[[perm_id]], perm_id)
            }
          } else {
            chunk <- pe_add_perm_columns(chunk, perm_entry, 1L)
          }
        } else if (is.list(perm_entry$cell)) {
          for (perm_id in seq_len(n_perms)) {
            perm_values <- perm_entry[[tier]][[perm_id]]
            chunk <- pe_add_perm_columns(chunk, perm_values, perm_id)
          }
        } else if (!is.null(perm_entry[[tier]])) {
          chunk <- pe_add_perm_columns(chunk, perm_entry[[tier]], 1L)
        }
      }
      rows[[length(rows) + 1L]] <- chunk
    }
  }

  if (!length(rows)) {
    return(pe_empty_canonical())
  }
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out[order(out$origin, out$tier, out$dest), , drop = FALSE]
}

pe_canonical_from_get_paths <- function(paths, origins, flat = FALSE) {
  if (is.null(paths$idx) || !length(paths$idx)) {
    return(pe_empty_canonical())
  }
  storage_like <- list(
    idx = paths$idx,
    distances = paths$distances,
    directions = paths$directions,
    perm_dist = paths$perm_dist
  )
  canonical <- pe_canonical_from_legacy_paths(storage_like, flat = flat)
  canonical[canonical$origin %in% pe_as_int(origins), , drop = FALSE]
}

pe_canonical_from_snapshot <- function(snapshot) {
  if (is.null(snapshot) || is.null(snapshot$storage)) {
    return(pe_empty_canonical())
  }
  flat <- identical(snapshot$region_type, "patch")
  if (identical(snapshot$implementation, "data.table")) {
    dt <- snapshot$storage
    if (!data.table::is.data.table(dt)) {
      dt <- data.table::as.data.table(dt)
    }
    if (nrow(dt) == 0L) {
      return(pe_empty_canonical())
    }
    out <- as.data.frame(dt)
    if (!"direction" %in% names(out)) {
      out$direction <- NA_integer_
    }
    perm_cols <- grep("^perm_dist", names(out), value = TRUE)
    keep <- c("origin", "tier", "dest", "distance", "direction", perm_cols)
    out <- out[, keep, drop = FALSE]
    out <- out[order(out$origin, out$tier, out$dest), , drop = FALSE]
    rownames(out) <- NULL
    return(out)
  }
  if (identical(snapshot$implementation, "legacy_list")) {
    return(pe_canonical_from_legacy_paths(snapshot$storage, flat = flat))
  }
  stop("Unknown snapshot implementation: ", snapshot$implementation, call. = FALSE)
}

pe_compare_vectors <- function(a, b, label, tolerance = 0L) {
  a <- pe_as_int(a)
  b <- pe_as_int(b)
  same_length <- length(a) == length(b)
  max_abs_diff <- if (same_length && length(a)) {
    diffs <- abs(a - b)
    if (all(is.na(diffs))) {
      0
    } else {
      max(diffs, na.rm = TRUE)
    }
  } else {
    NA_real_
  }
  exact <- same_length &&
    (length(a) == 0L || all(a == b | (is.na(a) & is.na(b))))
  within_tol <- same_length &&
    (length(a) == 0L || max_abs_diff <= tolerance)
  list(
    label = label,
    n_a = length(a),
    n_b = length(b),
    exact = exact,
    within_tolerance = within_tol,
    max_abs_diff = max_abs_diff
  )
}

pe_compare_canonical <- function(left, right, left_label = "left",
                                 right_label = "right",
                                 tolerance = 0L) {
  key_cols <- c("origin", "tier", "dest")
  all_cols <- union(names(left), names(right))
  for (col in setdiff(all_cols, names(left))) {
    left[[col]] <- NA
  }
  for (col in setdiff(all_cols, names(right))) {
    right[[col]] <- NA
  }

  left_keys <- paste(left$origin, left$tier, left$dest, sep = "|")
  right_keys <- paste(right$origin, right$tier, right$dest, sep = "|")
  only_left <- setdiff(left_keys, right_keys)
  only_right <- setdiff(right_keys, left_keys)
  common <- intersect(left_keys, right_keys)

  value_cols <- setdiff(all_cols, key_cols)
  col_results <- lapply(value_cols, function(col) {
    l_idx <- match(common, left_keys)
    r_idx <- match(common, right_keys)
    pe_compare_vectors(left[[col]][l_idx], right[[col]][r_idx], col, tolerance)
  })
  names(col_results) <- value_cols

  list(
    left_label = left_label,
    right_label = right_label,
    n_rows_left = nrow(left),
    n_rows_right = nrow(right),
    n_common = length(common),
    n_only_left = length(only_left),
    n_only_right = length(only_right),
    only_left = head(only_left, 10L),
    only_right = head(only_right, 10L),
    columns = col_results,
    equivalent = length(only_left) == 0L &&
      length(only_right) == 0L &&
      all(vapply(col_results, function(x) x$within_tolerance, logical(1L)))
  )
}

pe_compare_get_paths <- function(left_paths, right_paths, origins, flat = FALSE,
                                 left_label = "left", right_label = "right",
                                 tolerance = 0L) {
  left <- pe_canonical_from_get_paths(left_paths, origins, flat = flat)
  right <- pe_canonical_from_get_paths(right_paths, origins, flat = flat)
  pe_compare_canonical(left, right, left_label, right_label, tolerance)
}

pe_compare_snapshots <- function(left_snapshot, right_snapshot,
                                 left_label = "left", right_label = "right",
                                 tolerance = 0L) {
  left <- pe_canonical_from_snapshot(left_snapshot)
  right <- pe_canonical_from_snapshot(right_snapshot)
  pe_compare_canonical(left, right, left_label, right_label, tolerance)
}

pe_compare_graph_summary <- function(left_graph, right_graph) {
  left_graph <- left_graph %||% list()
  right_graph <- right_graph %||% list()
  keys <- union(names(left_graph), names(right_graph))
  diffs <- list()
  for (k in keys) {
    lv <- left_graph[[k]]
    rv <- right_graph[[k]]
    if (!identical(lv, rv)) {
      diffs[[k]] <- list(left = lv, right = rv)
    }
  }
  list(equivalent = length(diffs) == 0L, diffs = diffs)
}

`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}

pe_default_queries <- function(scenario_name) {
  switch(
    scenario_name,
    grid_max_distance = list(
      list(name = "configured_max", max_distance = 20000),
      list(name = "directions", directions = TRUE, max_distance = 20000)
    ),
    grid_two_tier = list(
      list(name = "tiered", directions = TRUE, max_distance = 30000),
      list(name = "tiered_perm", max_distance = 40000, perm_id = 1)
    ),
    grid_permeability = list(
      list(name = "base", max_distance = 40000),
      list(name = "perm", max_distance = 40000, perm_id = 1)
    ),
    patch_cities = list(
      list(name = "configured_max", directions = TRUE, max_distance = 150000)
    ),
    patch_permeability = list(
      list(name = "base", max_distance = 400000),
      list(name = "perm", max_distance = 400000, perm_id = 1)
    ),
    list(
      list(name = "full", directions = FALSE, max_distance = NULL),
      list(name = "directions", directions = TRUE, max_distance = NULL)
    )
  )
}

pe_run_get_paths_query <- function(region, cells, query, flat = FALSE) {
  args <- list(
    cells = cells,
    directions = isTRUE(query$directions),
    max_distance = query$max_distance %||% NULL,
    perm_id = query$perm_id %||% NULL
  )
  do.call(region$get_paths, args)
}

pe_run_get_paths_fast_queries <- function(region, cells, query) {
  if (!is.function(region$get_paths_fast)) {
    return(NULL)
  }
  lapply(cells, function(cell) {
    do.call(region$get_paths_fast, list(
      cell = cell,
      directions = isTRUE(query$directions),
      max_distance = query$max_distance %||% NULL,
      perm_id = query$perm_id %||% NULL
    ))
  })
}

pe_capture_checkpoint <- function(region, step_name, cells_requested,
                                  computed_origins, queries = NULL,
                                  flat = FALSE) {
  cells_requested <- pe_as_int(cells_requested)
  computed_origins <- pe_as_int(computed_origins)
  if (is.null(queries)) {
    queries <- list(list(name = "full"))
  }

  api <- list()
  for (query in queries) {
    qname <- query$name
    get_paths <- pe_run_get_paths_query(region, cells_requested, query, flat)
    api[[qname]] <- list(
      get_paths = get_paths,
      get_paths_canonical = pe_canonical_from_get_paths(
        get_paths,
        cells_requested,
        flat = flat
      )
    )
    fast <- pe_run_get_paths_fast_queries(region, cells_requested, query)
    if (!is.null(fast)) {
      names(fast) <- as.character(cells_requested)
      api[[qname]]$get_paths_fast = fast
    }
  }

  internal <- NULL
  if (is.function(region$paths_debug_snapshot)) {
    internal <- region$paths_debug_snapshot()
  }

  list(
    step = step_name,
    cells_requested = cells_requested,
    computed_origins = computed_origins,
    api = api,
    internal = internal
  )
}

pe_scenarios <- function(test_inputs_dir) {
  greater_melb <- file.path(test_inputs_dir, "greater_melb.tif")
  vic_cities <- file.path(test_inputs_dir, "vic_cities.csv")

  list(
    grid_single_origin = list(
      name = "grid_single_origin",
      flat = FALSE,
      build = function() {
        region <- Region(terra::rast(greater_melb))
        region$configure_paths(directions = TRUE)
        region
      },
      steps = list(
        list(name = "after_first_origin", calculate = 5922L)
      )
    ),
    grid_max_distance = list(
      name = "grid_max_distance",
      flat = FALSE,
      build = function() {
        region <- Region(terra::rast(greater_melb))
        region$configure_paths(max_distance = 20000)
        region
      },
      steps = list(
        list(name = "after_first_origin", calculate = 5922L)
      )
    ),
    grid_two_tier = list(
      name = "grid_two_tier",
      flat = FALSE,
      build = function() {
        region <- Region(terra::rast(greater_melb))
        region$set_aggr(aggr_factor = 5, inner_radius = 10000)
        region$configure_paths(directions = TRUE)
        region
      },
      steps = list(
        list(name = "after_first_origin", calculate = 5922L)
      )
    ),
    grid_permeability = list(
      name = "grid_permeability",
      flat = FALSE,
      build = function() {
        region <- Region(terra::rast(greater_melb))
        region$set_aggr(aggr_factor = 5, inner_radius = 10000)
        template <- region$get_template()
        perm_rast <- template
        perm_rast[region$get_indices()] <- round(region$get_indices()/1560)/10
        perm <- Permeability(perm_rast, region)
        region$configure_paths(permeability = perm)
        region
      },
      steps = list(
        list(name = "after_first_origin", calculate = 5922L)
      )
    ),
    grid_incremental = list(
      name = "grid_incremental",
      flat = FALSE,
      build = function() {
        region <- Region(terra::rast(greater_melb))
        region$configure_paths(directions = TRUE, max_distance = 30000)
        region
      },
      steps = list(
        list(name = "after_origin_5922", calculate = 5922L),
        list(name = "after_origin_5923", calculate = 5923L)
      )
    ),
    patch_cities = list(
      name = "patch_cities",
      flat = TRUE,
      build = function() {
        locations <- utils::read.csv(vic_cities)
        region <- Region(locations)
        region$configure_paths(directions = TRUE, max_distance = 150000)
        region
      },
      steps = list(
        list(name = "after_first_patch", calculate = 1L)
      )
    ),
    patch_permeability = list(
      name = "patch_permeability",
      flat = TRUE,
      build = function() {
        locations <- utils::read.csv(vic_cities)
        region <- Region(locations)
        perm_data <- matrix(c(
          1,  2, 0.8,
          1,  3, 0.7,
          1,  4, 0.7,
          1,  6, 0.6,
          1, 10, 0.8,
          2,  5, 0.6,
          3, 12, 0.6,
          3, 14, 0.5,
          4, 13, 0.5,
          5,  8, 0.6,
          6,  9, 0.6,
          7, 10, 0.5,
          10, 11, 0.8), ncol = 3, byrow = TRUE)
        colnames(perm_data) <- c("i", "j", "weight")
        perm <- Permeability(perm_data, region)
        region$configure_paths(permeability = perm)
        region
      },
      steps = list(
        list(name = "after_patches_1_2", calculate = 1:2)
      )
    )
  )
}

pe_capture_scenario <- function(scenario, run_label, git_info = pe_git_info()) {
  region <- scenario$build()
  computed <- integer(0L)
  checkpoints <- list()

  for (step in scenario$steps) {
    cells <- pe_as_int(step$calculate)
    region$calculate_paths(cells)
    computed <- unique(c(computed, cells))
    queries <- pe_default_queries(scenario$name)
    checkpoints[[length(checkpoints) + 1L]] <- pe_capture_checkpoint(
      region = region,
      step_name = step$name,
      cells_requested = cells,
      computed_origins = computed,
      queries = queries,
      flat = isTRUE(scenario$flat)
    )
  }

  list(
    meta = list(
      run_label = run_label,
      git_sha = git_info$sha,
      git_branch = git_info$branch,
      captured_at = format(Sys.time(), tz = "UTC", usetz = TRUE),
      package_version = utils::packageVersion("bsspread")
    ),
    scenario = scenario$name,
    flat = isTRUE(scenario$flat),
    checkpoints = checkpoints
  )
}

pe_capture_all <- function(run_label, out_dir = "output/paths_equivalence",
                           scenarios = NULL, test_inputs_dir = NULL) {
  if (!dir.exists(out_dir)) {
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  }
  if (is.null(test_inputs_dir)) {
    test_inputs_dir <- pe_find_test_inputs()
  }
  if (is.null(scenarios)) {
    scenarios <- pe_scenarios(test_inputs_dir)
  }

  paths <- list()
  for (scenario_name in names(scenarios)) {
    capture <- pe_capture_scenario(scenarios[[scenario_name]], run_label)
    out_path <- file.path(out_dir, paste0(run_label, "_", scenario_name, ".rds"))
    saveRDS(capture, out_path)
    paths[[scenario_name]] <- out_path
  }
  invisible(paths)
}

pe_print_compare_result <- function(result, title) {
  cat("\n== ", title, " ==\n", sep = "")
  if (isTRUE(result$equivalent)) {
    cat("PASS\n")
    return(invisible(result))
  }
  cat("FAIL\n")
  cat("Rows:", result$n_rows_left, "vs", result$n_rows_right,
      "| common:", result$n_common,
      "| only left:", result$n_only_left,
      "| only right:", result$n_only_right, "\n")
  if (length(result$only_left)) {
    cat("Only left (first 10):", paste(result$only_left, collapse = ", "), "\n")
  }
  if (length(result$only_right)) {
    cat("Only right (first 10):", paste(result$only_right, collapse = ", "), "\n")
  }
  for (col in names(result$columns)) {
    col_res <- result$columns[[col]]
    if (!col_res$within_tolerance) {
      cat("Column", col, ": max abs diff =", col_res$max_abs_diff, "\n")
    }
  }
  invisible(result)
}
