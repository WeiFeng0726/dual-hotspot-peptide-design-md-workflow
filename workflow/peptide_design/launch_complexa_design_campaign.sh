#!/usr/bin/env bash
set -euo pipefail

RUN_NAME="${RUN_NAME:?RUN_NAME is required}"
RUN_ROOT="/mnt/d/Project/cw/ly/Proteina-Complexa-runs/${RUN_NAME}"
LOG_FILE="${RUN_ROOT}/pipeline.log"
PID_FILE="${RUN_ROOT}/pipeline.pid"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
WORKER_SCRIPT="${SCRIPT_DIR}/run_complexa_design_campaign.sh"

mkdir -p "${RUN_ROOT}"
cd "${RUN_ROOT}"

nohup "${WORKER_SCRIPT}" > "${LOG_FILE}" 2>&1 < /dev/null &

pid=$!
echo "${pid}" > "${PID_FILE}"
echo "PID=${pid}"
echo "RUN_NAME=${RUN_NAME}"
echo "RUN_ROOT=${RUN_ROOT}"
echo "LOG=${LOG_FILE}"
