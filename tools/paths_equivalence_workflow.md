# Path Calculation Equivalence Workflow

This workflow compares **old vs new path storage** at each step of
`configure_paths()` → `calculate_paths()` → `get_paths()`, for fixed test
scenarios. It complements the simulation-level checks in
`tools/equivalence_workflow.md`.

## What gets compared

For each scenario step (e.g. after `calculate_paths(5922)`):

| Checkpoint | What it checks |
|------------|----------------|
| `cross_branch:get_paths:*` | Public API output on main vs optimised branch |
| `cross_branch:internal_storage` | Normalised path cache (nested lists vs `data.table`) |
| `cross_branch:graph_summary` | Graph vertex/edge counts for permeability |
| `right_branch:get_paths_fast_vs_get_paths:*` | New branch only: fast single-origin lookup vs legacy adapter |

All path rows are normalised to a **canonical long format**:

`origin | tier | dest | distance | direction | perm_dist*`

## Files

- `tools/paths_equivalence_helpers.R` — normalisation and comparison utilities
- `tools/paths_equivalence_capture.R` — run scenarios and save RDS checkpoints
- `tools/paths_equivalence_compare.R` — compare two capture sets

Built-in scenarios mirror `tests/testthat/test_region.R`:

- `grid_single_origin`, `grid_max_distance`, `grid_two_tier`, `grid_permeability`
- `grid_incremental` (5922 then 5923)
- `patch_cities`, `patch_permeability`

## 1) Add debug snapshot on the old branch (one-time)

The new branch exposes `region$paths_debug_snapshot()`. For internal storage
comparison on **main**, add this temporary method to both `Region.SpatRaster`
and `Region.data.frame` before `return(self)`:

```r
self$paths_debug_snapshot <- function() {
  list(
    implementation = "legacy_list",
    region_type = "grid", # or "patch" in Region.data.frame
    two_tier = is.list(aggr), # FALSE for patch region
    store_directions = is.list(paths$directions),
    max_distance = paths$max_distance,
    n_perms = if (is.list(paths$perms)) length(paths$perms) else 0L,
    computed_origins = NULL,
    storage = list(
      idx = paths$idx,
      distances = paths$distances,
      directions = paths$directions,
      perm_dist = paths$perm_dist,
      perms = paths$perms
    ),
    graphs = NULL
  )
}
```

Adjust `region_type` / `two_tier` for the patch variant. Without this, cross-branch
**API** comparisons still work; internal storage comparisons are skipped.

## 2) Capture checkpoints on each branch

From the repo root, with `bsspread` installed or loaded:

```bash
# On main (old list storage)
Rscript tools/paths_equivalence_capture.R main output/paths_equivalence

# On optimised branch (data.table storage)
Rscript tools/paths_equivalence_capture.R optimised output/paths_equivalence
```

Optional: one scenario only

```bash
Rscript tools/paths_equivalence_capture.R main output/paths_equivalence grid_two_tier
```

Each run writes RDS files like `output/paths_equivalence/main_grid_two_tier.rds`.

From R instead:

```r
devtools::load_all()
source("tools/paths_equivalence_helpers.R")
pe_capture_all("optimised", out_dir = "output/paths_equivalence")
```

## 3) Compare captures

```bash
Rscript tools/paths_equivalence_compare.R \
  output/paths_equivalence \
  main \
  optimised \
  output/paths_equivalence/paths_compare_main_vs_optimised.csv
```

Console output shows PASS/FAIL per checkpoint. The CSV is a summary table.
Exit code is non-zero if any check fails.

## Interpreting failures

- **`n_only_left` / `n_only_right` > 0** — destination sets differ for some
  `(origin, tier, dest)` keys. Inspect the first 10 keys printed.
- **Column diffs** — same destinations but different distances/directions/perms.
  Integer columns should match exactly unless you set a tolerance in helpers.
- **`get_paths_fast` mismatch** — bug in the legacy adapter or single-origin
  query path on the new branch (does not necessarily mean simulation divergence).
- **Graph summary mismatch** — permeability graph construction changed; follow up
  with permeability-specific scenarios (`grid_permeability`, `patch_permeability`).

## Extending scenarios

Add entries to `pe_scenarios()` in `paths_equivalence_helpers.R`. Each scenario
needs:

- `build()` — create and configure a `Region`
- `steps` — list of `list(name = "...", calculate = cells)`
- optional `flat = TRUE` for patch regions

Query variants come from `pe_default_queries()`; override per scenario there if
needed.
