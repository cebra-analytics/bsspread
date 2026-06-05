# Core scaling benchmark (`core_scaling.csv`)

One row per timed repetition of a **single benchmark timestep**, resumed from a fixed checkpoint.

## Workflow

1. Create a **master** checkpoint at end of tm 29, replicate 2 (`checkpoint_tm29_r2.rds`). By default `create-ckpt` runs **replicate 1 through full `time_steps`** (warms the shared `region` path cache), then replicate 2 through tm 29 and saves. Override warmup length with `CHECKPOINT_WARMUP_TIME_END`. The RDS stores `r` and `tm_next`; no extra “replicate” flag is needed on resume.
2. For each `cpus` in the sweep, `CHECKPOINT_RESUME=1` on that file and run **only** tm 30 (`BENCHMARK_TIME_END=30`), 10 times. The compare script runs **only the replicate in the checkpoint** (r=2), so replicate 1 is not re-simulated each rep.
3. Plot median `dispersal_s` (or `total_s`) vs `cpus`.

Use `benchmarks/run_core_scaling.sh` (see script header for env vars).

**Default `BENCHMARK_SCALING_REPS=0`**: one R process per rep (fresh `doParallel` cluster; best match for per-rep `dispersal_s` on 2+ cores).

Optional `BENCHMARK_SCALING_REPS=10`: one R session per CPU count — scenario built once, frozen checkpoint in memory, **paths-only copy** per rep (no `readRDS`). Faster total wall time; per-rep times can be higher than mode 0. Set `BENCHMARK_SCALING_GC=1` to run `gc()` between reps if RSS grows.

## Columns

| Column | Type | Description |
|--------|------|-------------|
| `recorded_at` | ISO 8601 | Wall time when the row was appended |
| `run_started_at` | ISO 8601 | UTC time when `run_core_scaling.sh run` started (same for all rows in that invocation) |
| `run_id` | string | Short label for the run (default `YYYYMMDD_HHMMSS` UTC; override with `SCALING_RUN_ID`) |
| `git_sha` | string | `git rev-parse HEAD` at run time |
| `hostname` | string | Machine identifier |
| `cpus` | int | `TASK_CPUS` / parallel cores for this rep |
| `rep_idx` | int | 1-based repeat index within `(cpus, …)` group |
| `warmup_discarded` | 0/1 | `1` if this row was the first rep and excluded from medians |
| `checkpoint_file` | path | Master checkpoint RDS (read-only for reps) |
| `checkpoint_saved_tm` | int | Timestep saved in checkpoint (e.g. 29) |
| `replicate` | int | Population replicate (default 2) |
| `benchmark_timestep` | int | Timestep executed (default 30) |
| `path_origins` | int | From checkpoint load message (`N path origins`) |
| `grow_s` | float | Growth phase seconds |
| `dispersal_s` | float | **Primary metric** for scaling plot |
| `dispersal_model_1_s` … `_3_s` | float | Per-model wall times (may be blank if log order ambiguous) |
| `impacts_s` | float | Impacts phase seconds |
| `actions_s` | float | Actions phase seconds |
| `other_s` | float | User function / collate / tail |
| `total_s` | float | Full timestep wall clock |
| `gc_s` | float | GC time attributed to this timestep |
| `rng_check` | float | Per-tm `runif(1)` sanity value (not equivalence golden) |
| `log_file` | path | Captured stdout for this rep |
| `notes` | string | Free text (outliers, GC spikes, etc.) |

## Analysis

- Filter by `run_id` or `run_started_at` to keep one sweep together.
- Group by `cpus`; summarise `dispersal_s` with **median** and IQR over reps where `warmup_discarded=0`.
- Optional reference curve: `dispersal_s[cpus=1] / cpus` (and same for `total_s`) for ideal linear speedup. See `benchmarks/plot_core_scaling.R`.
- Keep `path_origins` constant across rows; if it varies, the checkpoint was mutated.

## Requirements

- `params.json` must not force `cpus=1` unless `BENCHMARK_ALLOW_PARALLEL=1` is set (see `bsspread-compare.R`).
- Do **not** overwrite the master checkpoint; each rep uses `CHECKPOINT_RESUME=1` on the same file.
- Default checkpoint path is `platform:env` → `outputdir` in `params.json` (not necessarily `BENCHMARK_WORKDIR`). Override with `MASTER_CKPT` or `BENCHMARK_OUTPUTDIR`.
