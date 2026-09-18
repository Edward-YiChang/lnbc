#!/usr/bin/env bash
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "$0")" && pwd)"
SLURM_DIR="$REPO_ROOT/slurm"
mkdir -p "$SLURM_DIR"
SCRIPT_TMP_DIR="${SLURM_SCRIPT_TMPDIR:-${TMPDIR:-/tmp}}"
mkdir -p "$SCRIPT_TMP_DIR"

ACCOUNT="${SLURM_ACCOUNT:-${CC_ACCOUNT:-}}"
TIME_LIMIT="2-00:00:00"
MEMORY="16G"
NODES="8"
CPUS="10"
R_MODULE_CMD='module load StdEnv/2023; module load r/4.5.0; R_USER_LIB=$(Rscript -e '"'"'cat(Sys.getenv("R_LIBS_USER", unset = file.path(Sys.getenv("HOME"), "R", "library")))'"'"'); export R_LIBS_USER="$R_USER_LIB"; mkdir -p "$R_LIBS_USER"'

B="500"
MODE="final"
N_REPS="1000"
REP_START="1"
REP_END=""
SEED="123"
FPR_POINTS="101"
MC_DRAWS="10000"
DRY_RUN=false
LOCAL=false
WORKER=false
OUT_ROOT="$REPO_ROOT/result"

METHODS=(EMBC EMBC_fixk NP-inc LBNP LBNP-log)
DISTS=(lognorm weibull weibull2)
N_VALS=(100 300 500)
J_KEYS=(3 5)
PI_KEYS=(75 9)

usage() {
  cat <<'EOF'
Usage:
  bash run_all.sh [--all|options] --account def-xxx

Runs the unified bootstrap workflow. Each outer replication generates one
dataset per scenario and writes point-estimate, CI coverage, joint-region, and
ROC-grid rows under result/.

Options:
  --all            Use the full grid. Default for this script.
  --pilot          Use B=5 and mode=pilot.
  --final          Use B=500 and mode=final. Default.
  --B N            Override bootstrap resamples.
  --mode NAME      Override output mode label. Mostly for internal workers.
  --reps N         Number of outer replications/jobs. Default: 1000.
  --rep-start N    First replication index to submit. Default: 1.
  --rep-end N      Last replication index to submit. Default: --reps.
  --rep-id N       Run/submit a single replication index.
  --methods CSV    Methods, e.g. EMBC,EMBC_fixk,NP,NP-inc,LBNP,LBNP-log.
  --method M       Single method.
  --dists CSV      Distributions.
  --dist D         Single distribution.
  --n-vals CSV     Sample sizes.
  --n N            Single sample size.
  --J-keys CSV     J keys.
  --J K            Single J key.
  --pi-keys CSV    Pi keys.
  --pi K           Single pi key.
  --seed N         Base seed for bootstrap resampling. Default: 123.
  --fpr-points N   Number of ROC grid points. Default: 101.
  --mc-draws N     Monte Carlo draws for joint area. Default: 10000.
  --out-dir DIR    Output directory root. Default: result.
  --local          Run directly with local Rscript instead of sbatch.
  --account A      SLURM account. Can also use SLURM_ACCOUNT or CC_ACCOUNT.
  --time T         Time per job. Default: 2-00:00:00
  --mem M          Memory per job. Default: 16G
  --nodes N        Nodes per submitted replication job. Default: 8
  --cpus N         Cores per node/task. Default: 10
  --r-module CMD   Module/setup command run in each SLURM job.
  --dry-run        Print generated job scripts without submitting.
  --help           Show this help.
EOF
}

split_csv() {
  local raw="$1"
  local target="$2"
  local values
  local IFS=','
  read -r -a values <<< "$raw"
  eval "$target=(\"\${values[@]}\")"
}

join_csv() {
  local IFS=','
  echo "$*"
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --all)
      METHODS=(EMBC EMBC_fixk NP-inc LBNP LBNP-log)
      DISTS=(lognorm weibull weibull2)
      N_VALS=(100 300 500)
      J_KEYS=(3 5)
      PI_KEYS=(75 9)
      shift
      ;;
    --pilot) B="5"; MODE="pilot"; shift ;;
    --final) B="500"; MODE="final"; shift ;;
    --B) B="$2"; MODE="custom"; shift 2 ;;
    --mode) MODE="$2"; shift 2 ;;
    --reps) N_REPS="$2"; shift 2 ;;
    --rep-start) REP_START="$2"; shift 2 ;;
    --rep-end) REP_END="$2"; shift 2 ;;
    --rep-id|--rep) REP_START="$2"; REP_END="$2"; shift 2 ;;
    --methods) split_csv "$2" METHODS; shift 2 ;;
    --method) METHODS=("$2"); shift 2 ;;
    --dists) split_csv "$2" DISTS; shift 2 ;;
    --dist) DISTS=("$2"); shift 2 ;;
    --n-vals) split_csv "$2" N_VALS; shift 2 ;;
    --n) N_VALS=("$2"); shift 2 ;;
    --J-keys) split_csv "$2" J_KEYS; shift 2 ;;
    --J) J_KEYS=("$2"); shift 2 ;;
    --pi-keys) split_csv "$2" PI_KEYS; shift 2 ;;
    --pi) PI_KEYS=("$2"); shift 2 ;;
    --seed) SEED="$2"; shift 2 ;;
    --fpr-points) FPR_POINTS="$2"; shift 2 ;;
    --mc-draws) MC_DRAWS="$2"; shift 2 ;;
    --out-dir) OUT_ROOT="$2"; shift 2 ;;
    --local) LOCAL=true; shift ;;
    --account|--sbatch-account) ACCOUNT="$2"; shift 2 ;;
    --time|--sbatch-time) TIME_LIMIT="$2"; shift 2 ;;
    --mem|--sbatch-mem) MEMORY="$2"; shift 2 ;;
    --nodes|--sbatch-nodes) NODES="$2"; shift 2 ;;
    --cpus|--sbatch-cpus) CPUS="$2"; shift 2 ;;
    --r-module) R_MODULE_CMD="$2"; shift 2 ;;
    --dry-run) DRY_RUN=true; shift ;;
    --worker) WORKER=true; shift ;;
    --help|-h) usage; exit 0 ;;
    *) echo "ERROR: unknown option: $1" >&2; usage; exit 1 ;;
  esac
done

if [[ -z "$REP_END" ]]; then
  REP_END="$N_REPS"
fi

if ! $LOCAL && ! $DRY_RUN && ! $WORKER && [[ -z "$ACCOUNT" ]]; then
  echo "ERROR: --account is required unless SLURM_ACCOUNT or CC_ACCOUNT is set." >&2
  exit 1
fi

METHODS_CSV="$(join_csv "${METHODS[@]}")"
DISTS_CSV="$(join_csv "${DISTS[@]}")"
N_VALS_CSV="$(join_csv "${N_VALS[@]}")"
J_KEYS_CSV="$(join_csv "${J_KEYS[@]}")"
PI_KEYS_CSV="$(join_csv "${PI_KEYS[@]}")"

run_worker() {
  local rep_id="$REP_START"
  local cores

  cores="${SLURM_CPUS_PER_TASK:-$CPUS}"

  export OMP_NUM_THREADS=1
  export OPENBLAS_NUM_THREADS=1
  export MKL_NUM_THREADS=1
  export VECLIB_MAXIMUM_THREADS=1

  Rscript - <<RSCRIPT
repo_root <- "$REPO_ROOT"
source(file.path(repo_root, "bootstrap_core.R"))

cores <- suppressWarnings(as.integer("$cores"))
if (is.na(cores) || cores < 1) cores <- 1

invisible(lnbc_run_unified_replication_grid(
  rep_id = as.integer("$rep_id"),
  methods = strsplit("$METHODS_CSV", ",", fixed = TRUE)[[1]],
  dists = strsplit("$DISTS_CSV", ",", fixed = TRUE)[[1]],
  n_vals = as.integer(strsplit("$N_VALS_CSV", ",", fixed = TRUE)[[1]]),
  J_keys = strsplit("$J_KEYS_CSV", ",", fixed = TRUE)[[1]],
  pi_keys = strsplit("$PI_KEYS_CSV", ",", fixed = TRUE)[[1]],
  B = as.integer("$B"),
  total_reps = as.integer("$N_REPS"),
  seed = as.integer("$SEED"),
  fpr_grid = lnbc_default_roc_grid(as.integer("$FPR_POINTS")),
  mc_draws = as.integer("$MC_DRAWS"),
  repo_root = repo_root,
  out_dir = "$OUT_ROOT",
  mode = "$MODE",
  cores = cores,
  append = TRUE,
  verbose = TRUE
))
RSCRIPT
}

if $WORKER; then
  run_worker
  exit 0
fi

submit_rep_job() {
  local rep_id="$1"
  local rep_pad job_name sh_file job_setup

  if $LOCAL; then
    echo "running local LNBC_UNIFIED_REP_$(printf "%04d" "$rep_id")_${MODE}_B${B}"
    REP_START="$rep_id"
    run_worker
    return
  fi

  rep_pad="$(printf "%04d" "$rep_id")"
  job_name="LNBC_UNIFIED_REP_${rep_pad}_${MODE}_B${B}"
  sh_file="$(mktemp "$SCRIPT_TMP_DIR/${job_name}_sh_XXXXXX")"
  job_setup="$R_MODULE_CMD"

  cat > "$sh_file" <<SHSCRIPT
#!/bin/bash
#SBATCH --job-name=$job_name
#SBATCH --output=$SLURM_DIR/${job_name}-%j.out
#SBATCH --account=$ACCOUNT
#SBATCH --time=$TIME_LIMIT
#SBATCH --mem=$MEMORY
#SBATCH --nodes=$NODES
#SBATCH --ntasks=$NODES
#SBATCH --cpus-per-task=$CPUS
#SBATCH --chdir=$REPO_ROOT

set -euo pipefail
$job_setup
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1

IFS=',' read -r -a METHODS_ARR <<< "$METHODS_CSV"
IFS=',' read -r -a DISTS_ARR <<< "$DISTS_CSV"
IFS=',' read -r -a N_VALS_ARR <<< "$N_VALS_CSV"
IFS=',' read -r -a J_KEYS_ARR <<< "$J_KEYS_CSV"
IFS=',' read -r -a PI_KEYS_ARR <<< "$PI_KEYS_CSV"

for method in "\${METHODS_ARR[@]}"; do
  for dist in "\${DISTS_ARR[@]}"; do
    for n_val in "\${N_VALS_ARR[@]}"; do
      for J_key in "\${J_KEYS_ARR[@]}"; do
        for pi_key in "\${PI_KEYS_ARR[@]}"; do
          srun --exclusive --nodes=1 --ntasks=1 --cpus-per-task="\${SLURM_CPUS_PER_TASK:-$CPUS}" \\
            bash "$REPO_ROOT/run_all.sh" \\
            --worker --rep-id "$rep_id" --reps "$N_REPS" --B "$B" --mode "$MODE" \\
            --method "\$method" --dist "\$dist" --n "\$n_val" --J "\$J_key" --pi "\$pi_key" \\
            --seed "$SEED" --fpr-points "$FPR_POINTS" --mc-draws "$MC_DRAWS" \\
            --out-dir "$OUT_ROOT" --cpus "\${SLURM_CPUS_PER_TASK:-$CPUS}" &
        done
      done
    done
  done
done
wait
SHSCRIPT

  if $DRY_RUN; then
    echo "----- $job_name -----"
    cat "$sh_file"
    echo
    rm -f "$sh_file"
  else
    sbatch "$sh_file" >/dev/null
    rm -f "$sh_file"
    echo "submitted $job_name"
  fi
}

{
  echo "LNBC unified bootstrap submissions"
  echo "mode: $MODE"
  echo "B: $B"
  echo "replications/jobs: $REP_START-$REP_END of $N_REPS"
  echo "seed: $SEED"
  echo "fpr points: $FPR_POINTS"
  echo "joint area MC draws: $MC_DRAWS"
  echo "out: $OUT_ROOT"
  echo "temporary submit scripts: $SCRIPT_TMP_DIR"
  echo "nodes per job: $NODES"
  echo "cores per node/task: $CPUS"
  echo "methods: ${METHODS[*]}"
  echo "distributions: ${DISTS[*]}"
} | tee "$SLURM_DIR/submit_unified.log"

for rep_id in $(seq "$REP_START" "$REP_END"); do
  submit_rep_job "$rep_id"
done

echo "Done."
