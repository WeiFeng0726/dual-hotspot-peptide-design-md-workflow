#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT_LOCAL="$(cd "${SCRIPT_DIR}/../.." && pwd)"
ANALYSIS_DIR="${REPO_ROOT_LOCAL}/analysis/peptide_design"

REPO_ROOT="${REPO_ROOT:-/home/fatcat/Proteina-Complexa}"
RUN_NAME="${RUN_NAME:?RUN_NAME is required}"
TASK_NAME="${TASK_NAME:?TASK_NAME is required}"
HOTSPOT_SOURCE_IDS="${HOTSPOT_SOURCE_IDS:?HOTSPOT_SOURCE_IDS is required}"
TARGET_LEN="${TARGET_LEN:-156}"
TARGET_INPUT_START="${TARGET_INPUT_START:-265}"
NSAMPLES="${NSAMPLES:-1024}"
BATCH_SIZE="${BATCH_SIZE:-4}"
SEARCH_ALGO="${SEARCH_ALGO:-best-of-n}"
BEST_OF_N_REPLICAS="${BEST_OF_N_REPLICAS:-4}"
ENABLE_REWARD_MODEL="${ENABLE_REWARD_MODEL:-1}"
ENABLE_REFINEMENT="${ENABLE_REFINEMENT:-1}"
REFINEMENT_ALGO="${REFINEMENT_ALGO:-sequence_hallucination}"
GEN_NJOBS="${GEN_NJOBS:-1}"
FILTER_LIMIT="${FILTER_LIMIT:-$NSAMPLES}"
SCREEN_ROOT="${SCREEN_ROOT:-/mnt/d/Project/cw/ly/${RUN_NAME}_secondary_screening}"
SEED="${SEED:-5}"
AF2_WEIGHT_IPAE="${AF2_WEIGHT_IPAE:--0.60}"
AF2_WEIGHT_PLDDT="${AF2_WEIGHT_PLDDT:--0.25}"
AF2_WEIGHT_IPTM="${AF2_WEIGHT_IPTM:--0.25}"
AF2_WEIGHT_ICON="${AF2_WEIGHT_ICON:--0.10}"
REFINE_SOFT="${REFINE_SOFT:-false}"
REFINE_GREEDY="${REFINE_GREEDY:-true}"
REFINE_TEMP_ITERS="${REFINE_TEMP_ITERS:-45}"
REFINE_HARD_ITERS="${REFINE_HARD_ITERS:-5}"
REFINE_GREEDY_ITERS="${REFINE_GREEDY_ITERS:-20}"
REFINE_LOSS_HELIX="${REFINE_LOSS_HELIX:-0.0}"
REFINE_LOSS_RG="${REFINE_LOSS_RG:-0.20}"
REFINE_LOSS_IPTM="${REFINE_LOSS_IPTM:-0.10}"
MIN_CONTACT_BINDER_RESIDUES="${MIN_CONTACT_BINDER_RESIDUES:-2}"
MAX_HOTSPOT_DISTANCE="${MAX_HOTSPOT_DISTANCE:-5.0}"
MIN_BINDER_PLDDT="${MIN_BINDER_PLDDT:-0.58}"
MIN_I_PTM="${MIN_I_PTM:-0.20}"
MAX_MIN_IPAE="${MAX_MIN_IPAE:-0.45}"
FINAL_SELECTION_LIMIT="${FINAL_SELECTION_LIMIT:-20}"

if [[ -z "${REQUIRED_HOTSPOT_COVERAGE:-}" ]]; then
  IFS=',' read -r -a HOTSPOT_LIST <<< "${HOTSPOT_SOURCE_IDS}"
  REQUIRED_HOTSPOT_COVERAGE="${#HOTSPOT_LIST[@]}"
fi

RUN_ROOT="/mnt/d/Project/cw/ly/Proteina-Complexa-runs/${RUN_NAME}"
SEARCH_OUTPUT_ROOT="${RUN_ROOT}/inference/search_binder_local_pipeline_${TASK_NAME}_${RUN_NAME}"
RESCORE_ROOT="/mnt/d/Project/cw/ly/Proteina-Complexa-runs/${RUN_NAME}_af2_rescore"

mkdir -p "${RUN_ROOT}"
cd "${RUN_ROOT}"

source "${REPO_ROOT}/.venv/bin/activate"
source "${REPO_ROOT}/env.sh"
export PYTHONPATH="${REPO_ROOT}/src:${PYTHONPATH:-}"

echo "[config] run_name=${RUN_NAME}"
echo "[config] task_name=${TASK_NAME}"
echo "[config] hotspots=${HOTSPOT_SOURCE_IDS}"
echo "[config] nsamples=${NSAMPLES}"
echo "[config] search_algo=${SEARCH_ALGO}"
echo "[config] best_of_n_replicas=${BEST_OF_N_REPLICAS}"
echo "[config] refinement=${ENABLE_REFINEMENT}"
echo "[config] seed=${SEED}"

generate_cmd=(
  complexa generate "${REPO_ROOT}/configs/search_binder_local_pipeline.yaml"
  "++run_name=${RUN_NAME}"
  "++ckpt_path=${REPO_ROOT}/ckpts"
  "++autoencoder_ckpt_path=${REPO_ROOT}/ckpts/complexa_ae.ckpt"
  "++generation.task_name=${TASK_NAME}"
  "++generation.dataloader.batch_size=${BATCH_SIZE}"
  "++generation.dataloader.dataset.nres.nsamples=${NSAMPLES}"
  "++generation.search.algorithm=${SEARCH_ALGO}"
  "++generation.search.max_batch_size=${BATCH_SIZE}"
  "++generation.search.best_of_n.replicas=${BEST_OF_N_REPLICAS}"
  "++generation.filter.filter_samples_limit=${FILTER_LIMIT}"
  "++gen_njobs=${GEN_NJOBS}"
  "++seed=${SEED}"
)

if [[ "${ENABLE_REWARD_MODEL}" == "0" ]]; then
  generate_cmd+=("++generation.reward_model=null")
else
  generate_cmd+=(
    "++generation.reward_model.reward_models.af2folding.reward_weights.i_pae=${AF2_WEIGHT_IPAE}"
    "++generation.reward_model.reward_models.af2folding.reward_weights.plddt=${AF2_WEIGHT_PLDDT}"
    "++generation.reward_model.reward_models.af2folding.reward_weights.i_ptm=${AF2_WEIGHT_IPTM}"
    "++generation.reward_model.reward_models.af2folding.reward_weights.i_con=${AF2_WEIGHT_ICON}"
  )
fi

if [[ "${ENABLE_REFINEMENT}" == "0" ]]; then
  generate_cmd+=("++generation.refinement.algorithm=null")
else
  generate_cmd+=(
    "++generation.refinement.algorithm=${REFINEMENT_ALGO}"
    "++generation.refinement.enable_soft_optimization=${REFINE_SOFT}"
    "++generation.refinement.enable_greedy_optimization=${REFINE_GREEDY}"
    "++generation.refinement.n_temp_iters=${REFINE_TEMP_ITERS}"
    "++generation.refinement.n_hard_iters=${REFINE_HARD_ITERS}"
    "++generation.refinement.n_greedy_iters=${REFINE_GREEDY_ITERS}"
    "++generation.refinement.loss_weights.helix_binder=${REFINE_LOSS_HELIX}"
    "++generation.refinement.loss_weights.rg=${REFINE_LOSS_RG}"
    "++generation.refinement.loss_weights.i_ptm=${REFINE_LOSS_IPTM}"
  )
fi

echo "[stage] generate start"
"${generate_cmd[@]}"

echo "[stage] af2 rescore start"
INPUT_ROOT="${SEARCH_OUTPUT_ROOT}" \
OUTPUT_ROOT="${RESCORE_ROOT}" \
TARGET_LEN="${TARGET_LEN}" \
AF2_DIR="${REPO_ROOT}/community_models/ckpts/AF2" \
python "${ANALYSIS_DIR}/rescore_existing_af2.py"

echo "[stage] secondary screen start"
SOURCE_CSV="${RESCORE_ROOT}/af2_rescore_ranked.csv" \
OUTPUT_DIR="${SCREEN_ROOT}" \
TARGET_INPUT_START="${TARGET_INPUT_START}" \
HOTSPOT_SOURCE_IDS="${HOTSPOT_SOURCE_IDS}" \
REQUIRED_HOTSPOT_COVERAGE="${REQUIRED_HOTSPOT_COVERAGE}" \
MIN_CONTACT_BINDER_RESIDUES="${MIN_CONTACT_BINDER_RESIDUES}" \
MAX_HOTSPOT_DISTANCE="${MAX_HOTSPOT_DISTANCE}" \
MIN_BINDER_PLDDT="${MIN_BINDER_PLDDT}" \
MIN_I_PTM="${MIN_I_PTM}" \
MAX_MIN_IPAE="${MAX_MIN_IPAE}" \
FINAL_SELECTION_LIMIT="${FINAL_SELECTION_LIMIT}" \
python "${ANALYSIS_DIR}/secondary_screen_complexa.py"

echo "[stage] pipeline done"
