#!/usr/bin/env bash
# benchmarks/run_core_scaling.sh
#
# Core scaling study: fixed checkpoint (rep=2, end of tm=29) -> repeated tm=30.
#
# Usage:
#   ./benchmarks/run_core_scaling.sh create-ckpt   # once: run to tm=29 and save master ckpt
#   ./benchmarks/run_core_scaling.sh run           # sweep cpus × reps, append core_scaling.csv
#
# Required env (or edit defaults below):
#   BENCHMARK_WORKDIR   Directory containing params.json (scenario run dir)
#   BENCHMARK_SCRIPT    Path to bsspread-compare.R (default: repo root)
#
# Optional:
#   BENCHMARK_OUTPUTDIR Scenario output dir (default: platform:env$outputdir from params.json)
#   MASTER_CKPT         Read-only checkpoint (default: $BENCHMARK_OUTPUTDIR/checkpoint_tm29_r2.rds)
#   CORE_SCALING_CSV    Output CSV (default: repo benchmarks/core_scaling.csv)
#   CORES               Space-separated list (default: "1 2 4 8 16")
#   REPS_PER_CPU        Repetitions per cpu setting (default: 10)
#   DISCARD_WARMUP      1 = mark rep_idx 1 as warmup_discarded (default: 1)
#   CHECKPOINT_REPLICATE Replicate to save at CHECKPOINT_SAVED_TM (default: 2)
#   CHECKPOINT_WARMUP_REPS Comma-separated reps run first in same R session to warm
#                         region path cache (default: 1). Set empty to disable.
#   CHECKPOINT_WARMUP_TIME_END Last tm for warmup reps (default: unset = full time_steps)
#   CHECKPOINT_SAVED_TM Saved tm for checkpoint replicate (default: 29)
#   BENCHMARK_TIMESTEP  Timestep to run on resume (default: 30)
#   SCALING_RUN_ID        Optional label for this run (default: UTC timestamp slug YYYYMMDD_HHMMSS)
#   BENCHMARK_SCALING_REPS Reps per cpu in one R session (default: REPS_PER_CPU; set 0 for legacy 1 R/rep)
#   BENCHMARK_REPLICATE Optional override; if unset, resume uses r from the checkpoint file

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"

BENCHMARK_WORKDIR="${BENCHMARK_WORKDIR:-}"
BENCHMARK_SCRIPT="${BENCHMARK_SCRIPT:-${REPO_ROOT}/bsspread-compare.R}"
CORE_SCALING_CSV="${CORE_SCALING_CSV:-${SCRIPT_DIR}/core_scaling.csv}"
CORES="${CORES:-1 2 4 8 16}"
REPS_PER_CPU="${REPS_PER_CPU:-10}"
DISCARD_WARMUP="${DISCARD_WARMUP:-1}"
CHECKPOINT_REPLICATE="${CHECKPOINT_REPLICATE:-2}"
CHECKPOINT_WARMUP_REPS="${CHECKPOINT_WARMUP_REPS:-1}"
CHECKPOINT_SAVED_TM="${CHECKPOINT_SAVED_TM:-29}"
BENCHMARK_TIMESTEP="${BENCHMARK_TIMESTEP:-30}"
LOG_DIR="${LOG_DIR:-${SCRIPT_DIR}/scaling_logs}"

MASTER_CKPT="${MASTER_CKPT:-}"

usage() {
  sed -n '2,20p' "$0" | sed 's/^# \{0,1\}//'
  echo ""
  echo "Commands: create-ckpt | run"
  exit 1
}

die() { echo "run_core_scaling.sh: $*" >&2; exit 1; }

# Expand leading ~ in paths
expand_path() {
  local p="$1"
  if [[ "${p}" == "~"* ]]; then
    echo "${p/#\~/$HOME}"
  else
    echo "${p}"
  fi
}

# R saves checkpoints under params[['platform:env']]$outputdir, not BENCHMARK_WORKDIR
resolve_benchmark_outputdir() {
  if [[ -n "${BENCHMARK_OUTPUTDIR:-}" ]]; then
    expand_path "${BENCHMARK_OUTPUTDIR}"
    return 0
  fi
  local params="${BENCHMARK_WORKDIR}/params.json"
  local out=""
  if command -v python3 >/dev/null 2>&1; then
    out="$(python3 - "${params}" <<'PY' 2>/dev/null || true
import json, sys
with open(sys.argv[1]) as f:
    p = json.load(f)
env = p.get("platform:env") or {}
print(env.get("outputdir", ""))
PY
)"
  fi
  if [[ -z "${out}" ]]; then
    die "Set BENCHMARK_OUTPUTDIR or fix params.json platform:env\$outputdir"
  fi
  expand_path "${out}"
}

default_master_ckpt() {
  echo "$(resolve_benchmark_outputdir)/checkpoint_tm${CHECKPOINT_SAVED_TM}_r${CHECKPOINT_REPLICATE}.rds"
}

# Parse "Wrote checkpoint ... -> /path/checkpoint_tm29_r2.rds" from R log
ckpt_path_from_log() {
  local logfile="$1"
  local line path
  line="$(grep -E 'Wrote checkpoint \(tm_next=' "${logfile}" | tail -1 || true)"
  [[ -n "${line}" ]] || return 1
  path="$(echo "${line}" | sed -n 's/.*-> \(.*\)$/\1/p')"
  [[ -n "${path}" ]] || return 1
  expand_path "${path}"
}

[[ $# -ge 1 ]] || usage
CMD="$1"

if [[ -z "${BENCHMARK_WORKDIR}" ]]; then
  die "Set BENCHMARK_WORKDIR to the scenario directory (params.json)."
fi
if [[ ! -f "${BENCHMARK_WORKDIR}/params.json" ]]; then
  die "params.json not found in BENCHMARK_WORKDIR=${BENCHMARK_WORKDIR}"
fi
if [[ ! -f "${BENCHMARK_SCRIPT}" ]]; then
  die "BENCHMARK_SCRIPT not found: ${BENCHMARK_SCRIPT}"
fi

BENCHMARK_OUTPUTDIR="$(resolve_benchmark_outputdir)"
MASTER_CKPT="${MASTER_CKPT:-$(default_master_ckpt)}"
MASTER_CKPT="$(expand_path "${MASTER_CKPT}")"

git_sha() {
  (cd "${REPO_ROOT}" && git rev-parse HEAD 2>/dev/null) || echo "unknown"
}

HOSTNAME_SHORT="$(hostname -s 2>/dev/null || hostname)"

# Reps per cpu in one R process (0 = one Rscript per rep; best per-rep timings on 2+ cores)
BENCHMARK_SCALING_REPS="${BENCHMARK_SCALING_REPS:-0}"

# Shared R / benchmark env
export BENCHMARK_ALLOW_PARALLEL=1
export BENCHMARK_TIME_END="${BENCHMARK_TIMESTEP}"
# Do not set BENCHMARK_REPLICATE on resume: R uses checkpoint_data$r (e.g. 2 from checkpoint_tm29_r2.rds)

run_compare() {
  # shellcheck disable=SC2086
  (cd "${BENCHMARK_WORKDIR}" && Rscript "${BENCHMARK_SCRIPT}" "$@")
}

# Parse one summary line (BENCHMARK_REP=N tm=... or legacy tm=...)
parse_summary_line() {
  local summary="$1"
  [[ -n "${summary}" ]] || return 1

  PARSED_REP_IDX="$(echo "${summary}" | sed -n 's/^BENCHMARK_REP=\([0-9]*\).*/\1/p')"
  [[ -n "${PARSED_REP_IDX}" ]] || PARSED_REP_IDX="1"

  PARSED_GROW="$(echo "${summary}" | sed -n 's/.*grow=\([0-9.]*\)s.*/\1/p')"
  PARSED_DISPERSAL="$(echo "${summary}" | sed -n 's/.*dispersal=\([0-9.]*\)s.*/\1/p')"
  PARSED_IMPACTS="$(echo "${summary}" | sed -n 's/.*impacts=\([0-9.]*\)s.*/\1/p')"
  PARSED_ACTIONS="$(echo "${summary}" | sed -n 's/.*actions=\([0-9.]*\)s.*/\1/p')"
  PARSED_OTHER="$(echo "${summary}" | sed -n 's/.*other=\([0-9.]*\)s.*/\1/p')"
  PARSED_TOTAL="$(echo "${summary}" | sed -n 's/.*total=\([0-9.]*\)s.*/\1/p')"
  PARSED_GC="$(echo "${summary}" | sed -n 's/.*gc=\([0-9.]*\)s.*/\1/p')"
  PARSED_RNG="$(echo "${summary}" | sed -n 's/.*rng check=\([0-9.]*\).*/\1/p')"
  return 0
}

parse_log_metrics() {
  local logfile="$1"
  local summary
  summary="$(grep -E "tm=${BENCHMARK_TIMESTEP} r=${CHECKPOINT_REPLICATE} \\|" "${logfile}" | tail -1 || true)"
  [[ -n "${summary}" ]] || return 1

  local path_origins
  path_origins="$(grep -E 'Loaded checkpoint' "${logfile}" | tail -1 | sed -n 's/.*, \([0-9][0-9]*\) path origins).*/\1/p')"
  [[ -n "${path_origins}" ]] || path_origins=""

  local m1 m2 m3
  m1="$(grep -E '^    dispersal model 1:' "${logfile}" | tail -1 | sed -n 's/.*: \([0-9.]*\)s/\1/p')"
  m2="$(grep -E '^    dispersal model 2:' "${logfile}" | tail -1 | sed -n 's/.*: \([0-9.]*\)s/\1/p')"
  m3="$(grep -E '^    dispersal model 3:' "${logfile}" | tail -1 | sed -n 's/.*: \([0-9.]*\)s/\1/p')"

  parse_summary_line "${summary}" || return 1
  PARSED_PATH_ORIGINS="${path_origins}"
  PARSED_M1="${m1}"
  PARSED_M2="${m2}"
  PARSED_M3="${m3}"
  return 0
}

append_scaling_log_rows() {
  local logfile="$1" cpus="$2"
  local path_origins
  path_origins="$(grep -E 'Loaded checkpoint' "${logfile}" | head -1 | sed -n 's/.*, \([0-9][0-9]*\) path origins).*/\1/p')"
  [[ -n "${path_origins}" ]] || path_origins=""

  local pattern="^(BENCHMARK_REP=[0-9]+ )?tm=${BENCHMARK_TIMESTEP} r=${CHECKPOINT_REPLICATE} \\|"
  local summaries=()
  while IFS= read -r line; do
    summaries+=("${line}")
  done < <(grep -E "${pattern}" "${logfile}" || true)

  if [[ ${#summaries[@]} -eq 0 ]]; then
    return 1
  fi

  local summary rep_idx warmup
  for summary in "${summaries[@]}"; do
    parse_summary_line "${summary}" || continue
    rep_idx="${PARSED_REP_IDX}"
    warmup=0
    if [[ "${DISCARD_WARMUP}" == "1" && "${rep_idx}" -eq 1 ]]; then
      warmup=1
    fi
    PARSED_PATH_ORIGINS="${path_origins}"
    PARSED_M1=""
    PARSED_M2=""
    PARSED_M3=""
    append_csv_row "${cpus}" "${rep_idx}" "${warmup}" "${logfile}" ""
    echo "    rep ${rep_idx}: dispersal=${PARSED_DISPERSAL}s total=${PARSED_TOTAL}s"
  done
  return 0
}

CORE_SCALING_HEADER="recorded_at,run_started_at,run_id,git_sha,hostname,cpus,rep_idx,warmup_discarded,checkpoint_file,checkpoint_saved_tm,replicate,benchmark_timestep,path_origins,grow_s,dispersal_s,dispersal_model_1_s,dispersal_model_2_s,dispersal_model_3_s,impacts_s,actions_s,other_s,total_s,gc_s,rng_check,log_file,notes"

ensure_core_scaling_csv() {
  if [[ ! -f "${CORE_SCALING_CSV}" ]] || [[ ! -s "${CORE_SCALING_CSV}" ]]; then
    echo "${CORE_SCALING_HEADER}" >> "${CORE_SCALING_CSV}"
    return 0
  fi
  if ! head -1 "${CORE_SCALING_CSV}" | grep -q 'run_started_at'; then
    die "core_scaling.csv uses an old schema (no run_started_at). Back up the file and start fresh, or add run_started_at,run_id after recorded_at in the header."
  fi
}

append_csv_row() {
  local cpus="$1" rep_idx="$2" warmup="$3" logfile="$4" notes="$5"
  local ts sha
  ts="$(date -u +"%Y-%m-%dT%H:%M:%SZ")"
  sha="$(git_sha)"

  ensure_core_scaling_csv

  # Escape commas in notes
  notes="${notes//,/;}"

  printf '%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n' \
    "${ts}" "${SCALING_RUN_STARTED_AT}" "${SCALING_RUN_ID}" "${sha}" "${HOSTNAME_SHORT}" \
    "${cpus}" "${rep_idx}" "${warmup}" \
    "${MASTER_CKPT}" "${CHECKPOINT_SAVED_TM}" "${CHECKPOINT_REPLICATE}" "${BENCHMARK_TIMESTEP}" \
    "${PARSED_PATH_ORIGINS}" \
    "${PARSED_GROW}" "${PARSED_DISPERSAL}" "${PARSED_M1}" "${PARSED_M2}" "${PARSED_M3}" \
    "${PARSED_IMPACTS}" "${PARSED_ACTIONS}" "${PARSED_OTHER}" "${PARSED_TOTAL}" "${PARSED_GC}" "${PARSED_RNG}" \
    "${logfile}" "${notes}" >> "${CORE_SCALING_CSV}"
}

create_ckpt() {
  echo "Creating master checkpoint: ${MASTER_CKPT}"
  if [[ -n "${CHECKPOINT_WARMUP_REPS}" ]]; then
    if [[ -n "${CHECKPOINT_WARMUP_TIME_END:-}" ]]; then
      echo "  Warmup: replicate(s) ${CHECKPOINT_WARMUP_REPS} through tm=${CHECKPOINT_WARMUP_TIME_END}"
    else
      echo "  Warmup: replicate(s) ${CHECKPOINT_WARMUP_REPS} through full time_steps (params.json)"
    fi
    echo "  Then: replicate ${CHECKPOINT_REPLICATE} through tm=${CHECKPOINT_SAVED_TM}, save and exit"
  else
    echo "  (runs replicate ${CHECKPOINT_REPLICATE} through tm=${CHECKPOINT_SAVED_TM}, then exits)"
  fi
  mkdir -p "${LOG_DIR}"
  local logfile="${LOG_DIR}/create_ckpt_$(date +%Y%m%d_%H%M%S).log"

  export CHECKPOINT_TM="${CHECKPOINT_SAVED_TM}"
  export CHECKPOINT_REPLICATE="${CHECKPOINT_REPLICATE}"
  export CHECKPOINT_STOP=1
  unset CHECKPOINT_RESUME CHECKPOINT_FILE_IN CHECKPOINT_FILE BENCHMARK_REPLICATE
  export CHECKPOINT_WARMUP_REPS
  # Unset = warmup runs full time_steps; set CHECKPOINT_WARMUP_TIME_END to cap warmup only
  if [[ -n "${CHECKPOINT_WARMUP_TIME_END:-}" ]]; then
    export CHECKPOINT_WARMUP_TIME_END
  else
    unset CHECKPOINT_WARMUP_TIME_END
  fi
  # Cap only the checkpoint replicate (r=2); not the warmup replicate
  export BENCHMARK_TIME_END="${CHECKPOINT_SAVED_TM}"

  if [[ -n "${CHECKPOINT_FILE_IN:-}" ]]; then
    export CHECKPOINT_RESUME=1
    export CHECKPOINT_FILE_IN="${CHECKPOINT_FILE_IN}"
    echo "  Resume replicate ${CHECKPOINT_REPLICATE} from ${CHECKPOINT_FILE_IN}"
    echo "  (warmup reps still run from tm=1 unless BENCHMARK_REPLICATE limits loop)"
  elif [[ -n "${CHECKPOINT_WARMUP_REPS}" ]]; then
    echo "  One R session: warmup then checkpoint replicate (set TASK_CPUS for parallel)"
  else
    echo "  Full run from tm=1 for replicate ${CHECKPOINT_REPLICATE} only."
  fi

  # Optional: force write path (default in R is input.env$outputdir, same as MASTER_CKPT)
  export CHECKPOINT_FILE_OUT="${MASTER_CKPT}"

  run_compare 2>&1 | tee "${logfile}"

  if [[ ! -f "${MASTER_CKPT}" ]]; then
    local from_log
    from_log="$(ckpt_path_from_log "${logfile}" || true)"
    if [[ -n "${from_log}" && -f "${from_log}" ]]; then
      echo "Note: checkpoint found at ${from_log} (script expected ${MASTER_CKPT})" >&2
      MASTER_CKPT="${from_log}"
    else
      die "Checkpoint not found at ${MASTER_CKPT} (see ${logfile})"
    fi
  fi
  echo "Wrote ${MASTER_CKPT}"
}

run_scaling() {
  [[ -f "${MASTER_CKPT}" ]] || die "Master checkpoint missing: ${MASTER_CKPT} (run create-ckpt first)"

  SCALING_RUN_STARTED_AT="$(date -u +"%Y-%m-%dT%H:%M:%SZ")"
  SCALING_RUN_ID="${SCALING_RUN_ID:-$(date -u +"%Y%m%d_%H%M%S")}"
  export SCALING_RUN_STARTED_AT SCALING_RUN_ID

  mkdir -p "${LOG_DIR}"
  echo "Master checkpoint: ${MASTER_CKPT}"
  echo "Run id: ${SCALING_RUN_ID} (started ${SCALING_RUN_STARTED_AT})"
  echo "Appending results to: ${CORE_SCALING_CSV}"
  echo "Cores: ${CORES} | Reps per cpu: ${REPS_PER_CPU}"
  if [[ "${BENCHMARK_SCALING_REPS:-0}" -gt 0 ]]; then
    echo "Mode: one R session per cpu (${BENCHMARK_SCALING_REPS} reps; frozen ckpt, paths copy/rep)"
  else
    echo "Mode: one R session per rep (default; fresh parallel cluster; best for scaling curve)"
  fi

  for cpus in ${CORES}; do
    echo "--- cpus=${cpus} ---"
    export TASK_CPUS="${cpus}"
    export CHECKPOINT_RESUME=1
    export CHECKPOINT_FILE_IN="${MASTER_CKPT}"
    unset CHECKPOINT_STOP CHECKPOINT_TM CHECKPOINT_FILE_OUT BENCHMARK_REPLICATE

    if [[ "${BENCHMARK_SCALING_REPS:-0}" -gt 0 ]]; then
      export BENCHMARK_SCALING_REPS
      local logfile="${LOG_DIR}/cpus${cpus}_$(date +%Y%m%d_%H%M%S).log"
      echo "  ${BENCHMARK_SCALING_REPS} reps -> ${logfile}"
      if run_compare > "${logfile}" 2>&1; then
        :
      else
        echo "  WARNING: Rscript exited non-zero; see ${logfile}" >&2
      fi
      if append_scaling_log_rows "${logfile}" "${cpus}"; then
        :
      else
        echo "  WARNING: could not parse scaling rep lines in ${logfile}" >&2
      fi
    else
      export BENCHMARK_SCALING_REPS=0
      for ((rep_idx = 1; rep_idx <= REPS_PER_CPU; rep_idx++)); do
        local warmup=0
        if [[ "${DISCARD_WARMUP}" == "1" && "${rep_idx}" -eq 1 ]]; then
          warmup=1
        fi

        local logfile="${LOG_DIR}/cpus${cpus}_rep${rep_idx}_$(date +%Y%m%d_%H%M%S).log"

        echo "  rep ${rep_idx}/${REPS_PER_CPU} -> ${logfile}"
        if run_compare > "${logfile}" 2>&1; then
          :
        else
          echo "  WARNING: Rscript exited non-zero; see ${logfile}" >&2
        fi

        if parse_log_metrics "${logfile}"; then
          append_csv_row "${cpus}" "${rep_idx}" "${warmup}" "${logfile}" ""
          echo "    dispersal=${PARSED_DISPERSAL}s total=${PARSED_TOTAL}s path_origins=${PARSED_PATH_ORIGINS}"
        else
          append_csv_row "${cpus}" "${rep_idx}" "${warmup}" "${logfile}" "parse_failed"
          echo "    WARNING: could not parse tm=${BENCHMARK_TIMESTEP} line" >&2
        fi
      done
    fi
  done

  echo "Done. CSV: ${CORE_SCALING_CSV}"
}

case "${CMD}" in
  create-ckpt) create_ckpt ;;
  run) run_scaling ;;
  *) usage ;;
esac
