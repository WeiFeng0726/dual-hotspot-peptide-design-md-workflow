#!/usr/bin/env bash
set -euo pipefail

RUN_ROOT="/mnt/d/Project/cw/ly/Proteina-Complexa-runs/ly_b404_b405_pep8_12_probe"
REPO_ROOT="/home/fatcat/Proteina-Complexa"
LOG_FILE="${RUN_ROOT}/generate_probe.log"
PID_FILE="${RUN_ROOT}/generate_probe.pid"

mkdir -p /mnt/d/Project/cw/ly/Proteina-Complexa-storage/data/target_data
mkdir -p "${RUN_ROOT}"
cd "${RUN_ROOT}"

source "${REPO_ROOT}/.venv/bin/activate"
source "${REPO_ROOT}/env.sh"

nohup complexa generate "${REPO_ROOT}/configs/search_binder_local_pipeline.yaml" \
  ++run_name=ly_b404_b405_pep8_12_probe \
  ++ckpt_path="${REPO_ROOT}/ckpts" \
  ++autoencoder_ckpt_path="${REPO_ROOT}/ckpts/complexa_ae.ckpt" \
  ++generation.task_name=LY_B404_B405_PEP8_12 \
  ++generation.dataloader.batch_size=2 \
  ++generation.dataloader.dataset.nres.nsamples=2 \
  ++generation.search.algorithm=single-pass \
  ++gen_njobs=1 \
  > "${LOG_FILE}" 2>&1 < /dev/null &

pid=$!
echo "${pid}" > "${PID_FILE}"
echo "PID=${pid}"
echo "LOG=${LOG_FILE}"
