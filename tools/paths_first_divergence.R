#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1L) {
  stop(
    paste(
      "Usage:",
      "Rscript tools/paths_first_divergence.R",
      "<equivalence_comparison.csv> [occupied_delta_threshold] [out_csv]"),
    call. = FALSE
  )
}

cmp_path <- args[[1L]]
threshold <- if (length(args) >= 2L) as.numeric(args[[2L]]) else 0
out_csv <- if (length(args) >= 3L) args[[3L]] else {
  sub("\\.csv$", "_first_divergence.csv", cmp_path)
}

cmp <- utils::read.csv(cmp_path, stringsAsFactors = FALSE)
req <- c("replicate", "timestep", "occupied_n_delta")
if (!all(req %in% names(cmp))) {
  stop(
    "Expected columns in comparison CSV: ",
    paste(req, collapse = ", "),
    call. = FALSE)
}

cmp <- cmp[order(cmp$replicate, cmp$timestep), , drop = FALSE]
diverged <- abs(cmp$occupied_n_delta) > threshold
if ("jaccard_vs_main" %in% names(cmp)) {
  diverged <- diverged | (!is.na(cmp$jaccard_vs_main) & cmp$jaccard_vs_main < 1)
}

first <- do.call(
  rbind,
  lapply(split(cmp, cmp$replicate), function(df) {
    hit <- which(diverged[seq_len(nrow(df))])
    if (!length(hit)) {
      return(NULL)
    }
    i <- hit[[1L]]
    data.frame(
      replicate = df$replicate[[i]],
      first_divergent_timestep = df$timestep[[i]],
      occupied_n_delta = df$occupied_n_delta[[i]],
      new_occupied_n_delta = if ("new_occupied_n_delta" %in% names(df)) {
        df$new_occupied_n_delta[[i]]
      } else {
        NA_real_
      },
      total_pop_delta = if ("total_pop_delta" %in% names(df)) {
        df$total_pop_delta[[i]]
      } else {
        NA_real_
      },
      jaccard_vs_main = if ("jaccard_vs_main" %in% names(df)) {
        df$jaccard_vs_main[[i]]
      } else {
        NA_real_
      },
      stringsAsFactors = FALSE)
  }))

if (is.null(first) || !nrow(first)) {
  cat("No divergence found at threshold", threshold, "\n")
} else {
  utils::write.csv(first, out_csv, row.names = FALSE)
  cat("Wrote first divergence per replicate to:", out_csv, "\n")
  print(first)
}
