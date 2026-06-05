#!/usr/bin/env Rscript

# Compare PATHS_EQ_DUMP simulation snapshots (two RDS files from main vs optimised).
#
# Usage:
#   Rscript tools/paths_simulation_compare.R \
#     <main_snapshot.rds> <optimised_snapshot.rds> [out_detail_csv]

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2L) {
  stop(
    paste(
      "Usage:",
      "Rscript tools/paths_simulation_compare.R",
      "<main_snapshot.rds> <optimised_snapshot.rds> [out_detail_csv]"),
    call. = FALSE
  )
}

main_path <- args[[1L]]
opt_path <- args[[2L]]
out_detail <- if (length(args) >= 3L) args[[3L]] else {
  file.path(
    dirname(main_path),
    "paths_equivalence_detail.csv")
}

read_snapshot <- function(path) {
  if (!file.exists(path)) {
    stop("Snapshot not found: ", path, call. = FALSE)
  }
  readRDS(path)
}

stack_tables <- function(snap, side) {
  add_side <- function(df) {
    if (!nrow(df)) {
      return(df)
    }
    df$side <- side
    df
  }
  rbind(
    add_side(snap$tables$cell),
    add_side(snap$tables$aggr))
}

rel_diff <- function(a, b) {
  denom <- pmax(abs(a), abs(b), .Machine$double.eps)
  abs(a - b) / denom
}

main_snap <- read_snapshot(main_path)
opt_snap <- read_snapshot(opt_path)

main_long <- stack_tables(main_snap, "main")
opt_long <- stack_tables(opt_snap, "optimised")

key_cols <- c("tier", "dest_idx")
main_key <- paste(main_long$tier, main_long$dest_idx, sep = ":")
opt_key <- paste(opt_long$tier, opt_long$dest_idx, sep = ":")

all_keys <- unique(c(main_key, opt_key))
detail <- data.frame(
  key = all_keys,
  tier = sub(":.*", "", all_keys),
  dest_idx = as.integer(sub(".*:", "", all_keys)),
  stringsAsFactors = FALSE)

idx_main <- match(detail$key, main_key)
idx_opt <- match(detail$key, opt_key)

for (col in c("distance", "perm_dist", "direction", "dest_p")) {
  detail[[paste0(col, "_main")]] <- if (col %in% names(main_long)) {
    main_long[[col]][idx_main]
  } else {
    NA_real_
  }
  detail[[paste0(col, "_opt")]] <- if (col %in% names(opt_long)) {
    opt_long[[col]][idx_opt]
  } else {
    NA_real_
  }
  detail[[paste0(col, "_delta")]] <- detail[[paste0(col, "_opt")]] -
    detail[[paste0(col, "_main")]]
  detail[[paste0(col, "_rel_diff")]] <- rel_diff(
    detail[[paste0(col, "_main")]],
    detail[[paste0(col, "_opt")]])
}

detail$in_main <- !is.na(idx_main)
detail$in_opt <- !is.na(idx_opt)
detail$match_status <- ifelse(
  detail$in_main & detail$in_opt,
  "both",
  ifelse(detail$in_main, "main_only", "optimised_only"))

num_close <- function(x, tol = 1e-9) {
  is.finite(x) & abs(x) <= tol
}

detail$dest_p_match <- detail$match_status == "both" &
  num_close(detail$dest_p_delta)
detail$distance_match <- detail$match_status == "both" &
  num_close(detail$distance_delta)
detail$perm_dist_match <- detail$match_status == "both" &
  (is.na(detail$perm_dist_main) & is.na(detail$perm_dist_opt) |
     num_close(detail$perm_dist_delta))

utils::write.csv(detail, out_detail, row.names = FALSE)

jaccard_dest <- function(a, b) {
  a <- a[!is.na(a)]
  b <- b[!is.na(b)]
  if (!length(a) && !length(b)) {
    return(1)
  }
  length(intersect(a, b)) / length(union(a, b))
}

cell_main <- main_long[main_long$tier == "cell", "dest_idx", drop = TRUE]
cell_opt <- opt_long[opt_long$tier == "cell", "dest_idx", drop = TRUE]

summary_lines <- c(
  "Paths equivalence comparison",
  paste("main:", main_path),
  paste("optimised:", opt_path),
  paste("detail_csv:", out_detail),
  "",
  "Meta (main vs optimised):",
  paste(
    "  replicate:",
    main_snap$meta$replicate,
    opt_snap$meta$replicate),
  paste(
    "  timestep:",
    main_snap$meta$timestep,
    opt_snap$meta$timestep),
  paste(
    "  origin:",
    main_snap$meta$origin,
    opt_snap$meta$origin),
  paste(
    "  branch:",
    main_snap$meta$branch,
    opt_snap$meta$branch),
  paste(
    "  used_fast_paths:",
    main_snap$meta$used_fast_paths,
    opt_snap$meta$used_fast_paths),
  "",
  "Destination sets (cell tier):",
  paste("  n_cell main/opt:", nrow(main_snap$tables$cell), nrow(opt_snap$tables$cell)),
  paste("  n_aggr main/opt:", nrow(main_snap$tables$aggr), nrow(opt_snap$tables$aggr)),
  paste(
    "  main_only:",
    sum(detail$match_status == "main_only"),
    "optimised_only:",
    sum(detail$match_status == "optimised_only")),
  paste(
    "  dest_p mismatches:",
    sum(detail$match_status == "both" & !detail$dest_p_match)),
  paste(
    "  distance mismatches:",
    sum(detail$match_status == "both" & !detail$distance_match)),
  paste(
    "  perm_dist mismatches:",
    sum(detail$match_status == "both" & !detail$perm_dist_match)),
  paste(
    "  max |dest_p_delta| (matched):",
    if (any(detail$match_status == "both")) {
      max(abs(detail$dest_p_delta[detail$match_status == "both"]), na.rm = TRUE)
    } else {
      NA_real_
    }),
  paste(
    "  max |distance_delta| (matched):",
    if (any(detail$match_status == "both")) {
      max(abs(detail$distance_delta[detail$match_status == "both"]), na.rm = TRUE)
    } else {
      NA_real_
    }),
  paste(
    "  Jaccard(cell dest_idx):",
    jaccard_dest(cell_main, cell_opt)),
  paste(
    "  sum(dest_p) cell main/opt:",
    main_snap$summary$sum_dest_p_cell,
    opt_snap$summary$sum_dest_p_cell))

cat(paste(summary_lines, collapse = "\n"), "\n")
