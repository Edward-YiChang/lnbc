#!/usr/bin/env bash
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "$0")" && pwd)"

R_MODULE_CMD="${R_MODULE_CMD:-}"
if [[ -n "$R_MODULE_CMD" ]]; then
  # shellcheck disable=SC1090
  eval "$R_MODULE_CMD"
fi

if ! command -v Rscript >/dev/null 2>&1; then
  if command -v module >/dev/null 2>&1; then
    module load StdEnv/2023
    module load r/4.5.0
  elif [[ -f /etc/profile.d/modules.sh ]]; then
    # shellcheck disable=SC1091
    source /etc/profile.d/modules.sh
    module load StdEnv/2023
    module load r/4.5.0
  fi
fi

if ! command -v Rscript >/dev/null 2>&1; then
  cat >&2 <<'EOF'
ERROR: Rscript was not found.

Load R before running this script, for example:
  module load StdEnv/2023
  module load r/4.5.0
  bash setup_server.sh

Or pass your server's module command:
  R_MODULE_CMD='module load StdEnv/2023; module load r/4.5.0' bash setup_server.sh
EOF
  exit 1
fi

if [[ -z "${R_LIBS_USER:-}" ]]; then
  R_DEFAULT_LIB="$(Rscript -e 'cat(Sys.getenv("R_LIBS_USER", unset = file.path(Sys.getenv("HOME"), "R", "library")))')"
  export R_LIBS_USER="$R_DEFAULT_LIB"
else
  export R_LIBS_USER
fi
export CRAN_REPO="${CRAN_REPO:-https://cloud.r-project.org}"
export MAKEFLAGS="${MAKEFLAGS:--j1}"

if [[ -z "${R_MAKEVARS_USER:-}" ]]; then
  MAKEVARS_FILE="$REPO_ROOT/.R_Makevars_install"
  cat > "$MAKEVARS_FILE" <<'EOF'
CXXFLAGS=-O1
CXX11FLAGS=-O1
CXX14FLAGS=-O1
CXX17FLAGS=-O1
EOF
  export R_MAKEVARS_USER="$MAKEVARS_FILE"
fi

mkdir -p "$R_LIBS_USER"
Rscript "$REPO_ROOT/install_packages.R"
