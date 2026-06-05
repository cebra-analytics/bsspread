#!/usr/bin/env Rscript

# Capture path-calculation checkpoints for one branch/installation.
#
# Usage:
#   Rscript tools/paths_equivalence_capture.R <run_label> [out_dir] [scenario_name]
#
# Examples:
#   Rscript tools/paths_equivalence_capture.R main output/paths_equivalence
#   Rscript tools/paths_equivalence_capture.R optimised output/paths_equivalence grid_two_tier

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1L) {
  stop(
    paste(
      "Usage:",
      "Rscript tools/paths_equivalence_capture.R",
      "<run_label> [out_dir] [scenario_name]"),
    call. = FALSE
  )
}

run_label <- args[[1L]]
out_dir <- if (length(args) >= 2L) args[[2L]] else "output/paths_equivalence"
scenario_name <- if (length(args) >= 3L) args[[3L]] else NULL

script_args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", script_args, value = TRUE)
script_path <- if (length(file_arg)) {
  sub("^--file=", "", file_arg[[1L]])
} else {
  "tools/paths_equivalence_capture.R"
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

if (!requireNamespace("bsspread", quietly = TRUE)) {
  if (requireNamespace("devtools", quietly = TRUE)) {
    devtools::load_all(quiet = TRUE)
  } else {
    stop("Install or load bsspread before capturing.", call. = FALSE)
  }
} else {
  library(bsspread)
}

test_inputs_dir <- pe_find_test_inputs()
all_scenarios <- pe_scenarios(test_inputs_dir)

if (!is.null(scenario_name)) {
  if (!scenario_name %in% names(all_scenarios)) {
    stop(
      "Unknown scenario: ", scenario_name,
      ". Available: ", paste(names(all_scenarios), collapse = ", "),
      call. = FALSE
    )
  }
  all_scenarios <- all_scenarios[scenario_name]
}

paths <- pe_capture_all(
  run_label = run_label,
  out_dir = out_dir,
  scenarios = all_scenarios,
  test_inputs_dir = test_inputs_dir
)

cat("Captured path equivalence snapshots:\n")
for (nm in names(paths)) {
  cat(" -", nm, "->", paths[[nm]], "\n")
}
