#!/usr/bin/env bash
set -euo pipefail

RUN_ROOT="/mnt/d/Project/cw/ly/Proteina-Complexa-runs/dual_hotspot_peptide_campaign"
LOG_FILE="${RUN_ROOT}/campaign.log"
PID_FILE="${RUN_ROOT}/campaign.pid"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
WORKER_SCRIPT="${SCRIPT_DIR}/run_dual_hotspot_peptide_campaign.sh"

mkdir -p "${RUN_ROOT}"
cd "${RUN_ROOT}"

nohup "${WORKER_SCRIPT}" > "${LOG_FILE}" 2>&1 < /dev/null &

pid=$!
echo "${pid}" > "${PID_FILE}"
echo "PID=${pid}"
echo "RUN_ROOT=${RUN_ROOT}"
echo "LOG=${LOG_FILE}"
