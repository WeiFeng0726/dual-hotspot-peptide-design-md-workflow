#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT_LOCAL="$(cd "${SCRIPT_DIR}/../.." && pwd)"
ANALYSIS_DIR="${REPO_ROOT_LOCAL}/analysis/peptide_design"

RUN_ROOT="/mnt/d/Project/cw/ly/Proteina-Complexa-runs/ly_b404_b405_pep8_12_bulk64_af2_rescore"
LOG_FILE="${RUN_ROOT}/rescore.log"
PID_FILE="${RUN_ROOT}/rescore.pid"

mkdir -p "${RUN_ROOT}"
cd "${RUN_ROOT}"

source /home/fatcat/Proteina-Complexa/.venv/bin/activate
export PYTHONPATH=/home/fatcat/Proteina-Complexa/src:${PYTHONPATH:-}

nohup python "${ANALYSIS_DIR}/rescore_existing_af2.py" \
  > "${LOG_FILE}" 2>&1 < /dev/null &

pid=$!
echo "${pid}" > "${PID_FILE}"
echo "PID=${pid}"
echo "LOG=${LOG_FILE}"
