# RNG stream probes for baseline vs optimised equivalence bisect.
# Non-mutating: snapshots digest(.Random.seed), does not advance the stream.
# Enable from bsspread-compare.R:
#   RNG_PROBE=1 TASK_CPUS=1 Rscript bsspread-compare.R
# Optional filters:
#   RNG_PROBE_TM=4,5          only these timesteps (default: all)
#   RNG_PROBE_ORIGIN_MAX=5    per-origin probes for origins 1..N in dispersal_ready
#   RNG_PROBE_ORIGIN_DUMP=8   non-RNG key/value dump for origin(s) 8,16,...
#   RNG_PROBE_LABEL=dm1       prefix set by compare script per dispersal model

rng_probe_enabled <- function() {
  identical(Sys.getenv("RNG_PROBE"), "1")
}

rng_probe_tm_active <- function(tm) {
  if (!rng_probe_enabled()) {
    return(FALSE)
  }
  raw <- Sys.getenv("RNG_PROBE_TM", unset = "")
  if (!nzchar(raw)) {
    return(TRUE)
  }
  tm %in% as.integer(strsplit(raw, "[, ]+")[[1]])
}

rng_probe_origin_active <- function(origin_i) {
  if (!rng_probe_enabled()) {
    return(FALSE)
  }
  max_o <- as.integer(Sys.getenv("RNG_PROBE_ORIGIN_MAX", unset = "0"))
  max_o > 0L && origin_i <= max_o
}

rng_probe_origin_dump_targets <- function() {
  raw <- Sys.getenv("RNG_PROBE_ORIGIN_DUMP", unset = "")
  if (!nzchar(raw)) {
    return(integer(0))
  }
  as.integer(strsplit(raw, "[, ]+")[[1]])
}

rng_probe_origin_dump_active <- function(tm, origin_i) {
  if (!rng_probe_enabled() || !rng_probe_tm_active(tm)) {
    return(FALSE)
  }
  origin_i %in% rng_probe_origin_dump_targets()
}

rng_probe_origin_dump <- function(phase, tm, origin_i, loc_i = NA_integer_, ...) {
  if (!rng_probe_origin_dump_active(tm, origin_i)) {
    return(invisible(NULL))
  }
  prefix <- Sys.getenv("RNG_PROBE_LABEL", unset = "disperse")
  parts <- vapply(list(...), function(x) {
    if (is.logical(x)) {
      paste0(as.integer(x))
    } else if (is.numeric(x) && length(x) > 1L) {
      paste(format(x, digits = 15, trim = TRUE), collapse = ",")
    } else if (length(x) == 0L) {
      ""
    } else {
      as.character(x)
    }
  }, character(1))
  kv <- paste(names(parts), parts, sep = "=", collapse = " ")
  loc_part <- if (is.na(loc_i)) "NA" else as.character(loc_i)
  message(sprintf(
    "    rng dump %s %s tm=%d o=%d loc=%s %s",
    prefix, phase, tm, origin_i, loc_part, kv))
  invisible(NULL)
}

#' Digest of the current R RNG state (does not advance the stream)
#' @export
rng_state_digest <- function() {
  digest::digest(.Random.seed, algo = "xxhash64")
}

#' Timestep RNG probe (bsspread-compare harness)
#' @export
rng_probe <- function(label, tm, r = NA_integer_) {
  if (!rng_probe_tm_active(tm)) {
    return(invisible(NULL))
  }
  r_part <- if (is.na(r)) "" else sprintf(" r=%d", r)
  message(sprintf(
    "  rng probe %s tm=%d%s: %s",
    label, tm, r_part, rng_state_digest()))
  invisible(NULL)
}

rng_probe_disperse <- function(label, tm, origin_i = NA_integer_) {
  if (!rng_probe_tm_active(tm)) {
    return(invisible(NULL))
  }
  prefix <- Sys.getenv("RNG_PROBE_LABEL", unset = "disperse")
  tag <- rng_state_digest()
  if (is.na(origin_i)) {
    message(sprintf(
      "    rng probe %s %s tm=%d: %s",
      prefix, label, tm, tag))
  } else {
    message(sprintf(
      "    rng probe %s %s tm=%d o=%d: %s",
      prefix, label, tm, origin_i, tag))
  }
  invisible(NULL)
}

rng_probe_origin <- function(label, tm, origin_i) {
  if (rng_probe_origin_active(origin_i)) {
    rng_probe_disperse(label, tm, origin_i)
  }
  invisible(NULL)
}
