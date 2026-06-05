#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3L) {
  stop(
    paste(
      "Usage:",
      "Rscript tools/equivalence_compare.R",
      "<main_metrics.csv> <optimised_metrics.csv> <out_csv>",
      "[main_occupancy.rds optimised_occupancy.rds]"),
    call. = FALSE
  )
}

main_metrics_path <- args[[1L]]
opt_metrics_path <- args[[2L]]
out_csv <- args[[3L]]

main_df <- utils::read.csv(main_metrics_path, stringsAsFactors = FALSE)
opt_df <- utils::read.csv(opt_metrics_path, stringsAsFactors = FALSE)

key_cols <- c("replicate", "timestep")
join_key <- paste(main_df$replicate, main_df$timestep, sep = "_")
opt_key <- paste(opt_df$replicate, opt_df$timestep, sep = "_")

main_df$key <- join_key
opt_df$key <- opt_key

common <- intersect(main_df$key, opt_df$key)
main_df <- main_df[main_df$key %in% common, , drop = FALSE]
opt_df <- opt_df[opt_df$key %in% common, , drop = FALSE]
main_df <- main_df[order(main_df$replicate, main_df$timestep), , drop = FALSE]
opt_df <- opt_df[order(opt_df$replicate, opt_df$timestep), , drop = FALSE]

cmp <- data.frame(
  replicate = main_df$replicate,
  timestep = main_df$timestep,
  occupied_n_main = main_df$occupied_n,
  occupied_n_opt = opt_df$occupied_n,
  occupied_n_delta = opt_df$occupied_n - main_df$occupied_n,
  new_occupied_n_main = main_df$new_occupied_n,
  new_occupied_n_opt = opt_df$new_occupied_n,
  new_occupied_n_delta = opt_df$new_occupied_n - main_df$new_occupied_n,
  total_pop_main = main_df$total_pop,
  total_pop_opt = opt_df$total_pop,
  total_pop_delta = opt_df$total_pop - main_df$total_pop,
  total_s_main = main_df$total_s,
  total_s_opt = opt_df$total_s,
  total_s_ratio = opt_df$total_s / pmax(main_df$total_s, .Machine$double.eps),
  dispersal_s_main = main_df$dispersal_s,
  dispersal_s_opt = opt_df$dispersal_s,
  dispersal_s_ratio = opt_df$dispersal_s / pmax(main_df$dispersal_s, .Machine$double.eps),
  rng_check_main = main_df$rng_check,
  rng_check_opt = opt_df$rng_check,
  stringsAsFactors = FALSE
)

jaccard <- rep(NA_real_, nrow(cmp))
if (length(args) >= 5L) {
  main_occ <- readRDS(args[[4L]])
  opt_occ <- readRDS(args[[5L]])
  for (i in seq_len(nrow(cmp))) {
    k <- paste0("r", cmp$replicate[i], "_t", cmp$timestep[i])
    a <- main_occ[[k]]
    b <- opt_occ[[k]]
    if (is.null(a) || is.null(b)) {
      next
    }
    inter <- length(intersect(a, b))
    uni <- length(union(a, b))
    jaccard[i] <- if (uni == 0L) 1 else inter / uni
  }
}
cmp$jaccard_vs_main <- jaccard

dir.create(dirname(out_csv), recursive = TRUE, showWarnings = FALSE)
utils::write.csv(cmp, out_csv, row.names = FALSE)

summ <- function(x) c(
  median = stats::median(x, na.rm = TRUE),
  mean = mean(x, na.rm = TRUE),
  p05 = as.numeric(stats::quantile(x, 0.05, na.rm = TRUE)),
  p95 = as.numeric(stats::quantile(x, 0.95, na.rm = TRUE))
)

cat("Wrote comparison CSV:", out_csv, "\n")
cat("\n== Runtime ratios (optimised/main) ==\n")
print(rbind(
  dispersal_s_ratio = summ(cmp$dispersal_s_ratio),
  total_s_ratio = summ(cmp$total_s_ratio)
))

cat("\n== Outcome deltas (optimised-main) ==\n")
print(rbind(
  occupied_n_delta = summ(cmp$occupied_n_delta),
  new_occupied_n_delta = summ(cmp$new_occupied_n_delta),
  total_pop_delta = summ(cmp$total_pop_delta)
))

if (all(is.na(cmp$jaccard_vs_main))) {
  cat("\nJaccard not computed (occupancy RDS not provided).\n")
} else {
  cat("\n== Jaccard vs main ==\n")
  print(summ(cmp$jaccard_vs_main))
}
