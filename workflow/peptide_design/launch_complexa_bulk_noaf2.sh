#!/usr/bin/env bash
set -euo pipefail

RUN_NAME="ly_b404_b405_pep8_12_bulk64_noreward"
RUN_ROOT="/mnt/d/Project/cw/ly/Proteina-Complexa-runs/${RUN_NAME}"
REPO_ROOT="/home/fatcat/Proteina-Complexa"
LOG_FILE="${RUN_ROOT}/generate.log"
PID_FILE="${RUN_ROOT}/generate.pid"

mkdir -p /mnt/d/Project/cw/ly/Proteina-Complexa-storage/data/target_data/custom
mkdir -p "${RUN_ROOT}"
cd "${RUN_ROOT}"

source "${REPO_ROOT}/.venv/bin/activate"
source "${REPO_ROOT}/env.sh"

nohup complexa generate "${REPO_ROOT}/configs/search_binder_local_pipeline.yaml" \
  ++run_name="${RUN_NAME}" \
  ++ckpt_path="${REPO_ROOT}/ckpts" \
  ++autoencoder_ckpt_path="${REPO_ROOT}/ckpts/complexa_ae.ckpt" \
  ++generation.task_name=LY_B404_B405_PEP8_12 \
  ++generation.dataloader.batch_size=8 \
  ++generation.dataloader.dataset.nres.nsamples=64 \
  ++generation.search.algorithm=single-pass \
  ++generation.reward_model=null \
  ++gen_njobs=1 \
  > "${LOG_FILE}" 2>&1 < /dev/null &

pid=$!
echo "${pid}" > "${PID_FILE}"
echo "PID=${pid}"
echo "RUN_NAME=${RUN_NAME}"
echo "RUN_ROOT=${RUN_ROOT}"
echo "LOG=${LOG_FILE}"
