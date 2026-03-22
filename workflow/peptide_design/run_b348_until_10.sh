#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT_LOCAL="$(cd "${SCRIPT_DIR}/../.." && pwd)"
WORKER_SCRIPT="${SCRIPT_DIR}/run_complexa_design_campaign.sh"
AGG_SCRIPT="${REPO_ROOT_LOCAL}/analysis/peptide_design/aggregate_b348_candidates.py"
PYTHON_BIN="/home/fatcat/Proteina-Complexa/.venv/bin/python"
ROOT="/mnt/d/Project/cw/ly"
TARGET_TOP=10
NSAMPLES="${NSAMPLES:-3072}"
BATCH_SIZE="${BATCH_SIZE:-8}"
MAX_BATCHES="${MAX_BATCHES:-10}"
SEEDS=(5 11 17 23 29 37 43 53 61 71)

aggregate_count() {
  if [[ -d "${ROOT}/Complexa_B348_aggregate_top10" ]]; then
    local summary="${ROOT}/Complexa_B348_aggregate_top10/summary.txt"
    if [[ -f "${summary}" ]]; then
      awk -F': ' '/Top N available/ {print $2}' "${summary}"
      return 0
    fi
  fi
  echo 0
}

run_aggregate() {
  ROOT="${ROOT}" OUTPUT_DIR="${ROOT}/Complexa_B348_aggregate_top10" TOP_N="${TARGET_TOP}" \
    "${PYTHON_BIN}" "${AGG_SCRIPT}"
}

echo "[b348] start aggregate-driven expansion"
run_aggregate || true

for idx in "${!SEEDS[@]}"; do
  batch_num=$((idx + 1))
  current_count="$(aggregate_count)"
  echo "[b348] current_strict_hits=${current_count}"
  if (( current_count >= TARGET_TOP )); then
    echo "[b348] target reached before batch ${batch_num}"
    break
  fi
  if (( batch_num > MAX_BATCHES )); then
    echo "[b348] reached MAX_BATCHES=${MAX_BATCHES}"
    break
  fi

  seed="${SEEDS[$idx]}"
  run_name="ly_b348_pep8_12_bulk3072_seed${seed}"
  screen_root="${ROOT}/Complexa_B348_secondary_screening_seed${seed}"
  if [[ "${seed}" == "5" ]]; then
    run_name="ly_b348_pep8_12_bulk3072_noreward"
    screen_root="${ROOT}/Complexa_B348_secondary_screening_bulk3072"
  fi

  if [[ -f "${screen_root}/selection_summary.txt" ]]; then
    echo "[b348] batch already completed for seed=${seed}, aggregating"
    run_aggregate
    continue
  fi

  if [[ "${seed}" == "5" ]]; then
    progress_csv="${ROOT}/Proteina-Complexa-runs/${run_name}_af2_rescore/af2_rescore_progress.csv"
    if [[ -f "${progress_csv}" ]]; then
      echo "[b348] seed=${seed} batch appears to be in progress already; waiting for completion"
      while [[ ! -f "${screen_root}/selection_summary.txt" ]]; do
        if [[ -f "${progress_csv}" ]]; then
          ts="$(date '+%F %T')"
          nrows="$(python3 - <<PY
import csv
from pathlib import Path
p = Path("${progress_csv}")
with p.open() as f:
    print(sum(1 for _ in csv.DictReader(f)))
PY
)"
          echo "[b348] ${ts} progress rows=${nrows}"
        fi
        sleep 120
      done
      echo "[b348] existing seed=${seed} batch finished"
      run_aggregate
      continue
    fi
  fi

  echo "[b348] start batch seed=${seed} run_name=${run_name}"
  RUN_NAME="${run_name}" \
  TASK_NAME="LY_B348_PEP8_12" \
  HOTSPOT_SOURCE_IDS="348" \
  REQUIRED_HOTSPOT_COVERAGE="1" \
  TARGET_LEN="156" \
  TARGET_INPUT_START="265" \
  NSAMPLES="${NSAMPLES}" \
  BATCH_SIZE="${BATCH_SIZE}" \
  SEARCH_ALGO="single-pass" \
  ENABLE_REWARD_MODEL="0" \
  ENABLE_REFINEMENT="0" \
  SCREEN_ROOT="${screen_root}" \
  FINAL_SELECTION_LIMIT="${TARGET_TOP}" \
  SEED="${seed}" \
  "${WORKER_SCRIPT}"

  echo "[b348] batch done seed=${seed}"
  run_aggregate
done

echo "[b348] finished"
