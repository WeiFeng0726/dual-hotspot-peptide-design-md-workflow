#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT_LOCAL="$(cd "${SCRIPT_DIR}/../.." && pwd)"
ANALYSIS_DIR="${REPO_ROOT_LOCAL}/analysis/peptide_design"

RUN_NAME="ly_b404_b405_pep8_12_bulk512_noreward"
RUN_ROOT="/mnt/d/Project/cw/ly/Proteina-Complexa-runs/${RUN_NAME}"
REPO_ROOT="/home/fatcat/Proteina-Complexa"
SEARCH_OUTPUT_ROOT="${RUN_ROOT}/inference/search_binder_local_pipeline_LY_B404_B405_PEP8_12_${RUN_NAME}"
RESCORE_ROOT="/mnt/d/Project/cw/ly/Proteina-Complexa-runs/${RUN_NAME}_af2_rescore"
SCREEN_ROOT="/mnt/d/Project/cw/ly/Complexa_B404_B405_secondary_screening_bulk512"

mkdir -p /mnt/d/Project/cw/ly/Proteina-Complexa-storage/data/target_data/custom
mkdir -p "${RUN_ROOT}"
cd "${RUN_ROOT}"

source "${REPO_ROOT}/.venv/bin/activate"
source "${REPO_ROOT}/env.sh"
export PYTHONPATH="${REPO_ROOT}/src:${PYTHONPATH:-}"

echo "[stage] generate start"
complexa generate "${REPO_ROOT}/configs/search_binder_local_pipeline.yaml" \
  ++run_name="${RUN_NAME}" \
  ++ckpt_path="${REPO_ROOT}/ckpts" \
  ++autoencoder_ckpt_path="${REPO_ROOT}/ckpts/complexa_ae.ckpt" \
  ++generation.task_name=LY_B404_B405_PEP8_12 \
  ++generation.dataloader.batch_size=8 \
  ++generation.dataloader.dataset.nres.nsamples=512 \
  ++generation.search.algorithm=single-pass \
  ++generation.reward_model=null \
  ++gen_njobs=1

echo "[stage] af2 rescore start"
INPUT_ROOT="${SEARCH_OUTPUT_ROOT}" \
OUTPUT_ROOT="${RESCORE_ROOT}" \
TARGET_LEN=156 \
AF2_DIR="${REPO_ROOT}/community_models/ckpts/AF2" \
python "${ANALYSIS_DIR}/rescore_existing_af2.py"

echo "[stage] secondary screen start"
SOURCE_CSV="${RESCORE_ROOT}/af2_rescore_ranked.csv" \
OUTPUT_DIR="${SCREEN_ROOT}" \
TARGET_INPUT_START=265 \
HOTSPOT_SOURCE_IDS=404,405 \
python "${ANALYSIS_DIR}/secondary_screen_complexa.py"

echo "[stage] pipeline done"
