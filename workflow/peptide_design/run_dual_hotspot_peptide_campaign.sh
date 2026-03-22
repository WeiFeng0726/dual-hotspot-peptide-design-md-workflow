#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
WORKER_SCRIPT="${SCRIPT_DIR}/run_complexa_design_campaign.sh"

echo "[campaign] start 404_405"
RUN_NAME="ly_b404_b405_pep8_12_bulk3072_noreward" \
TASK_NAME="LY_B404_B405_PEP8_12" \
HOTSPOT_SOURCE_IDS="404,405" \
REQUIRED_HOTSPOT_COVERAGE="2" \
TARGET_LEN="156" \
TARGET_INPUT_START="265" \
NSAMPLES="3072" \
BATCH_SIZE="8" \
SEARCH_ALGO="single-pass" \
ENABLE_REWARD_MODEL="0" \
ENABLE_REFINEMENT="0" \
SCREEN_ROOT="/mnt/d/Project/cw/ly/Complexa_B404_B405_secondary_screening_bulk3072" \
FINAL_SELECTION_LIMIT="10" \
"${WORKER_SCRIPT}"

echo "[campaign] start 348"
RUN_NAME="ly_b348_pep8_12_bulk3072_noreward" \
TASK_NAME="LY_B348_PEP8_12" \
HOTSPOT_SOURCE_IDS="348" \
REQUIRED_HOTSPOT_COVERAGE="1" \
TARGET_LEN="156" \
TARGET_INPUT_START="265" \
NSAMPLES="3072" \
BATCH_SIZE="8" \
SEARCH_ALGO="single-pass" \
ENABLE_REWARD_MODEL="0" \
ENABLE_REFINEMENT="0" \
SCREEN_ROOT="/mnt/d/Project/cw/ly/Complexa_B348_secondary_screening_bulk3072" \
FINAL_SELECTION_LIMIT="10" \
"${WORKER_SCRIPT}"

echo "[campaign] done"
