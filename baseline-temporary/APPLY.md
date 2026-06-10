# Applying RNG probes to baseline package on VM

Copy these four files into your baseline `bsspread` package `R/` directory (and root for DESCRIPTION/NAMESPACE):

| File here | Destination on VM |
|-----------|---------------------|
| `RngProbe.R` | `R/RngProbe.R` |
| `Dispersal.R` | `R/Dispersal.R` |
| `DESCRIPTION` | `DESCRIPTION` |
| `NAMESPACE` | `NAMESPACE` |

`bsspread-compare.R` stays as in the repo (current version with `RNG_PROBE` support).

Then reinstall:

```bash
cd /path/to/bsspread
Rscript -e 'devtools::install(quick=TRUE)'
# or: R CMD INSTALL .
```

Run baseline probe capture:

```bash
RNG_PROBE=1 RNG_PROBE_TM=4,5 RNG_PROBE_ORIGIN_MAX=10 \
TASK_CPUS=1 BENCHMARK_REPLICATE=1 \
Rscript bsspread-compare.R 2>&1 | tee rng_probe_baseline.log
```

Use the same command on the optimised package install for `rng_probe_optimised.log`.

Probes log `digest::digest(.Random.seed, algo = "xxhash64")` and do not advance the RNG stream.
End-of-timestep summary lines use `rng state=<hash>` (not `runif`).

Origin dump (no RNG consumption), e.g. for tm=3 origin 8:

```bash
RNG_PROBE=1 RNG_PROBE_TM=3 RNG_PROBE_ORIGIN_DUMP=8 \
TASK_CPUS=1 BENCHMARK_REPLICATE=1 Rscript bsspread-compare.R
```
