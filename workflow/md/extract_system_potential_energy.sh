#!/usr/bin/env bash
set -euo pipefail

GMX_BIN="${GMX_BIN:-/usr/local/gmx2025.4/bin/gmx}"

log() {
    printf '[%s] %s\n' "$(date '+%F %T')" "$*"
}

require_file() {
    [[ -f "$1" ]] || {
        printf 'Missing file: %s\n' "$1" >&2
        exit 1
    }
}

if [[ $# -lt 1 ]]; then
    printf 'Usage: %s RUN_DIR [RUN_DIR ...]\n' "$0" >&2
    exit 1
fi

for run_dir in "$@"; do
    prod_dir="${run_dir}/07_md_50ns"
    analysis_dir="${run_dir}/10_system_energy"
    energy_file="${prod_dir}/md_50ns.edr"
    output_xvg="${analysis_dir}/system_potential_energy.xvg"

    require_file "$energy_file"
    mkdir -p "$analysis_dir"

    log "Extracting whole-system potential energy for ${run_dir}"
    printf 'Potential\n0\n' | "$GMX_BIN" energy -f "$energy_file" -o "$output_xvg" >/dev/null
done
