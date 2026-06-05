#!/usr/bin/env Rscript

# Compare path-calculation checkpoints from two branch captures.
#
# Usage:
#   Rscript tools/paths_equivalence_compare.R \
#     <left_dir> <left_label> <right_label> [out_report_csv] [scenario_name]
#
# Example:
#   Rscript tools/paths_equivalence_compare.R \
#     output/paths_equivalence main optimised \
#     output/paths_equivalence/paths_compare.csv

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3L) {
  stop(
    paste(
      "Usage:",
      "Rscript tools/paths_equivalence_compare.R",
      "<capture_dir> <left_label> <right_label> [out_report_csv] [scenario_name]"),
    call. = FALSE
  )
}

capture_dir <- args[[1L]]
left_label <- args[[2L]]
right_label <- args[[3L]]
out_report <- if (length(args) >= 4L) args[[4L]] else file.path(
  capture_dir,
  paste0("paths_compare_", left_label, "_vs_", right_label, ".csv")
)
scenario_filter <- if (length(args) >= 5L) args[[5L]] else NULL

script_args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", script_args, value = TRUE)
script_path <- if (length(file_arg)) {
  sub("^--file=", "", file_arg[[1L]])
} else {
  "tools/paths_equivalence_compare.R"
}
repo_root <- normalizePath(
  file.path(dirname(script_path), ".."),
  winslash = "/",
  mustWork = FALSE
)
if (dir.exists(file.path(repo_root, "R"))) {
  setwd(repo_root)
}

helpers <- file.path(getwd(), "tools/paths_equivalence_helpers.R")
if (!file.exists(helpers)) {
  stop("Missing tools/paths_equivalence_helpers.R", call. = FALSE)
}
source(helpers)

scenario_from_path <- function(path, label) {
  nm <- basename(path)
  prefix <- paste0(label, "_")
  sub("\\.rds$", "", substr(nm, nchar(prefix) + 1L, nchar(nm)))
}

left_files <- list.files(
  capture_dir,
  pattern = paste0("^", left_label, "_.*\\.rds$"),
  full.names = TRUE
)
right_files <- list.files(
  capture_dir,
  pattern = paste0("^", right_label, "_.*\\.rds$"),
  full.names = TRUE
)

left_map <- stats::setNames(
  left_files,
  vapply(left_files, scenario_from_path, "", label = left_label)
)
right_map <- stats::setNames(
  right_files,
  vapply(right_files, scenario_from_path, "", label = right_label)
)
common_scenarios <- intersect(names(left_map), names(right_map))

if (!is.null(scenario_filter)) {
  common_scenarios <- intersect(common_scenarios, scenario_filter)
}

if (!length(common_scenarios)) {
  stop("No common scenario captures found in ", capture_dir, call. = FALSE)
}

pe_compare_fast_vs_get <- function(get_paths, get_paths_fast, cells, flat) {
  if (is.null(get_paths_fast)) {
    return(list(equivalent = NA, reason = "get_paths_fast unavailable"))
  }
  for (cell in cells) {
    key <- as.character(cell)
    fast <- get_paths_fast[[key]]
    if (is.null(fast)) {
      return(list(equivalent = FALSE, reason = paste("missing fast path for", key)))
    }

    if (flat) {
      fast_idx <- fast$idx
      fast_dist <- fast$distances
      fast_dir <- fast$directions
      fast_perm <- fast$perm_dist
      get_idx <- get_paths$idx[[key]]
      get_dist <- get_paths$distances[[key]]
      get_dir <- if (!is.null(get_paths$directions)) get_paths$directions[[key]] else NULL
      get_perm <- if (!is.null(get_paths$perm_dist)) get_paths$perm_dist[[key]] else NULL
    } else {
      fast_idx <- fast$idx
      fast_dist <- fast$distances
      fast_dir <- fast$directions
      fast_perm <- fast$perm_dist
      get_idx <- get_paths$idx[[key]]
      get_dist <- get_paths$distances[[key]]
      get_dir <- if (!is.null(get_paths$directions)) get_paths$directions[[key]] else NULL
      get_perm <- if (!is.null(get_paths$perm_dist)) get_paths$perm_dist[[key]] else NULL
    }

    fast_canonical <- pe_canonical_from_legacy_paths(
      list(
        idx = setNames(list(fast_idx), key),
        distances = setNames(list(fast_dist), key),
        directions = if (!is.null(fast_dir)) setNames(list(fast_dir), key) else NULL,
        perm_dist = if (!is.null(fast_perm)) setNames(list(fast_perm), key) else NULL
      ),
      flat = flat
    )
    get_canonical <- pe_canonical_from_legacy_paths(
      list(
        idx = setNames(list(get_idx), key),
        distances = setNames(list(get_dist), key),
        directions = if (!is.null(get_dir)) setNames(list(get_dir), key) else NULL,
        perm_dist = if (!is.null(get_perm)) setNames(list(get_perm), key) else NULL
      ),
      flat = flat
    )
    cmp <- pe_compare_canonical(fast_canonical, get_canonical, "get_paths_fast", "get_paths")
    if (!cmp$equivalent) {
      return(list(equivalent = FALSE, cell = cell, compare = cmp))
    }
  }
  list(equivalent = TRUE)
}

rows <- list()
n_fail <- 0L
n_pass <- 0L

record_row <- function(scenario, step, comparison, pass, n_only_left = NA_integer_,
                       n_only_right = NA_integer_) {
  if (isTRUE(pass)) {
    n_pass <<- n_pass + 1L
  } else if (!is.na(pass)) {
    n_fail <<- n_fail + 1L
  }
  rows[[length(rows) + 1L]] <<- data.frame(
    scenario = scenario,
    step = step,
    comparison = comparison,
    equivalent = pass,
    n_only_left = n_only_left,
    n_only_right = n_only_right,
    stringsAsFactors = FALSE
  )
}

for (scenario in common_scenarios) {
  left_capture <- readRDS(left_map[[scenario]])
  right_capture <- readRDS(right_map[[scenario]])
  flat <- isTRUE(left_capture$flat)

  n_steps <- min(length(left_capture$checkpoints), length(right_capture$checkpoints))
  for (i in seq_len(n_steps)) {
    left_cp <- left_capture$checkpoints[[i]]
    right_cp <- right_capture$checkpoints[[i]]
    step <- left_cp$step

    query_names <- intersect(names(left_cp$api), names(right_cp$api))
    for (query_name in query_names) {
      cmp_api <- pe_compare_get_paths(
        left_cp$api[[query_name]]$get_paths,
        right_cp$api[[query_name]]$get_paths,
        left_cp$cells_requested,
        flat = flat,
        left_label = left_label,
        right_label = right_label
      )
      pass <- isTRUE(cmp_api$equivalent)
      record_row(
        scenario, step, paste0("cross_branch:get_paths:", query_name), pass,
        cmp_api$n_only_left, cmp_api$n_only_right
      )
      pe_print_compare_result(
        cmp_api,
        paste(scenario, step, "cross-branch API", query_name, sep = " | ")
      )
    }

    if (!is.null(left_cp$internal) && !is.null(right_cp$internal)) {
      cmp_internal <- pe_compare_snapshots(
        left_cp$internal,
        right_cp$internal,
        left_label,
        right_label
      )
      pass <- isTRUE(cmp_internal$equivalent)
      record_row(
        scenario, step, "cross_branch:internal_storage", pass,
        cmp_internal$n_only_left, cmp_internal$n_only_right
      )
      pe_print_compare_result(
        cmp_internal,
        paste(scenario, step, "cross-branch internal storage", sep = " | ")
      )

      cmp_graph <- pe_compare_graph_summary(
        left_cp$internal$graphs,
        right_cp$internal$graphs
      )
      pass <- isTRUE(cmp_graph$equivalent)
      record_row(scenario, step, "cross_branch:graph_summary", pass)
      if (!pass) {
        cat("\n== ", scenario, " | ", step, " | graph summary ==\n", sep = "")
        cat("FAIL\n")
        print(cmp_graph$diffs)
      }
    }

    for (query_name in names(right_cp$api)) {
      fast_cmp <- pe_compare_fast_vs_get(
        right_cp$api[[query_name]]$get_paths,
        right_cp$api[[query_name]]$get_paths_fast,
        right_cp$cells_requested,
        flat = flat
      )
      pass <- fast_cmp$equivalent
      record_row(
        scenario, step,
        paste0("right_branch:get_paths_fast_vs_get_paths:", query_name),
        pass
      )
      if (isFALSE(pass)) {
        cat("\n== ", scenario, " | ", step,
            " | get_paths_fast vs get_paths (", query_name, ") ==\n", sep = "")
        cat("FAIL\n")
        if (!is.null(fast_cmp$compare)) {
          print(fast_cmp$compare)
        } else {
          print(fast_cmp)
        }
      }
    }
  }
}

report <- do.call(rbind, rows)
dir.create(dirname(out_report), recursive = TRUE, showWarnings = FALSE)
utils::write.csv(report, out_report, row.names = FALSE)

cat("\nWrote report:", out_report, "\n")
cat("Checks passed:", n_pass, "\n")
cat("Checks failed:", n_fail, "\n")
if (n_fail > 0L) {
  quit(status = 1L)
}
