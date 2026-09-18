#!/usr/bin/env bash
set -euo pipefail

# Submit from code/: sbatch --array=1-1000 --account=ACCOUNT \
#   --cpus-per-task=10 --mem=8G --time=1-12:00:00 gof/run_gof.slurm lognormal
cd "${SLURM_SUBMIT_DIR:-$PWD}"
[[ -f bootstrap_core.R ]] || { echo "Submit from the code directory." >&2; exit 1; }
case "${1:-}" in
  lognormal) script=gof/gof_symm.R ;;
  weibull) script=gof/gof_symm_weibull.R ;;
  *) echo "Expected lognormal or weibull." >&2; exit 1 ;;
esac
module load StdEnv/2023
module load r/4.5.0
export R_LIBS_USER="${R_LIBS_USER:-$HOME/R/library}"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1
Rscript "$script" --rep-id "${SLURM_ARRAY_TASK_ID:?Submit with --array}" \
  --cores "${SLURM_CPUS_PER_TASK:-1}" --B 500 --scenario both
