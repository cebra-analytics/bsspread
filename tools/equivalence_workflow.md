# Behavioural Equivalence Workflow (manual branch switching)

This workflow assumes you manually:

1. check out/install the branch version into your Docker container, then
2. run the simulation script.

The helper files are:

- `tools/equivalence_metrics.R` (collect metrics during a run)
- `tools/equivalence_compare.R` (compare outputs from two runs)

## 1) Collect metrics for each branch

Source `tools/equivalence_metrics.R` in your run script and create a collector:

```r
eq <- eq_init(
  run_label = "main_seedset_a",
  git_sha = "edc76da...",
  branch = "main",
  scenario = "single_initial_location",
  out_dir = "/workdir/output/equivalence",
  save_occupancy = TRUE
)
```

Inside your timestep loop, after each timestep completes, append:

```r
eq <- eq_capture_step(
  state = eq,
  n = n,
  timestep = tm,
  replicate = r,
  grow_s = elapsed(t0, t1),
  dispersal_s = elapsed(t1, t2),
  impacts_s = elapsed(t2, t3),
  actions_s = elapsed(t3, t4),
  other_s = elapsed(t4, t5),
  total_s = elapsed(t0, t5),
  rng_check = runif(1)
)
```

After the run:

```r
eq_write_outputs(eq)
```

Repeat with the optimised branch, using a different `run_label`.

## 2) Compare main vs optimised outputs

```bash
Rscript tools/equivalence_compare.R \
  /path/main_seedset_a_metrics.csv \
  /path/optimised_seedset_a_metrics.csv \
  /path/equivalence_comparison.csv \
  /path/main_seedset_a_occupancy_sets.rds \
  /path/optimised_seedset_a_occupancy_sets.rds
```

This writes:

- a per-replicate/per-timestep comparison CSV with runtime ratios and outcome deltas
- console summary for:
  - runtime ratios (`optimised/main`)
  - behavioural deltas (`optimised - main`)
  - Jaccard overlap vs main (if occupancy RDS provided)

## Profile only `calculate_dispersals`

Set env vars before running `temp-bsspread.R`:

```bash
export PROFVIS_CALCULATE_DISPERSALS=1
export PROFVIS_TM=25
export PROFVIS_R=1
export PROFVIS_OUTPUT=/workdir/output/profvis_calculate_dispersals_r1_tm25.html
```

Run with serial cores (`cpus=1` or no parallel) for readable profiles.  
`region$calculate_paths()` still runs once at the start of `disperse()`, but the
profvis sample window covers only the serial `calculate_dispersals` loop.

## 3) Path-level dump at one (replicate, timestep, origin)

Use this after `equivalence_compare.R` to inspect the first divergent point.

Find the earliest divergent timestep per replicate:

```bash
Rscript tools/paths_first_divergence.R \
  /path/equivalence_comparison.csv \
  0 \
  /path/first_divergence.csv
```

Pick `(replicate, timestep, origin)` from that output (origin is a **region cell
index**, not the packed occupied-slot index). In your replicate loop, set:

```bash
export PATHS_EQ_DUMP=1
export PATHS_EQ_R=1
export PATHS_EQ_REPLICATE=1
export PATHS_EQ_TM=7
export PATHS_EQ_ORIGIN=5922
export PATHS_EQ_RUN_LABEL=main
export PATHS_EQ_OUT_DIR=/workdir/output/paths_equivalence
export PATHS_EQ_RNG_CHECK=$(python3 -c 'import random; print(random.random())')
```

Run once on `main`, once on the optimised branch (change `PATHS_EQ_RUN_LABEL`).
Reinstall the package after checkout so `Dispersal.R` hooks are active.

Compare the two snapshots:

```bash
Rscript tools/paths_simulation_compare.R \
  /workdir/output/paths_equivalence/paths_eq_main_r1_tm7_origin5922.rds \
  /workdir/output/paths_equivalence/paths_eq_optimised_r1_tm7_origin5922.rds \
  /workdir/output/paths_equivalence/paths_detail_r1_tm7_origin5922.csv
```

For step-by-step path-cache comparison (test scenarios), see
`tools/paths_equivalence_workflow.md` and `tools/paths_equivalence_compare.R`.

The detail CSV is one row per destination (`tier` + `dest_idx`) with side-by-side
`distance`, `perm_dist`, `direction`, and `dest_p` (relative probabilities before
sampling). Console output summarises set mismatches and max deltas.
