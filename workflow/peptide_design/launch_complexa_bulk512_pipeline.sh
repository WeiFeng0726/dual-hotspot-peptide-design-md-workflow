#!/usr/bin/env bash
set -euo pipefail

RUN_NAME="ly_b404_b405_pep8_12_bulk512_noreward"
RUN_ROOT="/mnt/d/Project/cw/ly/Proteina-Complexa-runs/${RUN_NAME}"
RESCORE_ROOT="/mnt/d/Project/cw/ly/Proteina-Complexa-runs/${RUN_NAME}_af2_rescore"
SCREEN_ROOT="/mnt/d/Project/cw/ly/Complexa_B404_B405_secondary_screening_bulk512"
LOG_FILE="${RUN_ROOT}/pipeline.log"
PID_FILE="${RUN_ROOT}/pipeline.pid"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
WORKER_SCRIPT="${SCRIPT_DIR}/run_complexa_bulk512_pipeline.sh"

mkdir -p /mnt/d/Project/cw/ly/Proteina-Complexa-storage/data/target_data/custom
mkdir -p "${RUN_ROOT}"
cd "${RUN_ROOT}"

nohup "${WORKER_SCRIPT}" > "${LOG_FILE}" 2>&1 < /dev/null &

pid=$!
echo "${pid}" > "${PID_FILE}"
echo "PID=${pid}"
echo "RUN_NAME=${RUN_NAME}"
echo "RUN_ROOT=${RUN_ROOT}"
echo "RESCORE_ROOT=${RESCORE_ROOT}"
echo "SCREEN_ROOT=${SCREEN_ROOT}"
echo "LOG=${LOG_FILE}"
