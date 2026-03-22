#!/usr/bin/env bash
set -euo pipefail

GMX_BIN="${GMX_BIN:-/usr/local/gmx2025.4/bin/gmx}"
OMP_THREADS="${OMP_THREADS:-8}"

log() {
    printf '[%s] %s\n' "$(date '+%F %T')" "$*"
}

require_file() {
    [[ -f "$1" ]] || {
        printf 'Missing file: %s\n' "$1" >&2
        exit 1
    }
}

write_rerun_mdp() {
    local mdp_path="$1"
    cat > "$mdp_path" <<'EOF'
integrator              = md
nsteps                  = 1
dt                      = 0.002
nstxout                 = 0
nstvout                 = 0
nstfout                 = 0
nstxout-compressed      = 0
nstlog                  = 0
nstenergy               = 1
nstcalcenergy           = 1
continuation            = yes
constraint-algorithm    = lincs
constraints             = h-bonds
lincs_iter              = 1
lincs_order             = 4
cutoff-scheme           = Verlet
nstlist                 = 20
ns_type                 = grid
coulombtype             = PME
rcoulomb                = 1.0
rvdw                    = 1.0
vdwtype                 = Cut-off
tcoupl                  = no
pcoupl                  = no
pbc                     = xyz
DispCorr                = no
energygrps              = Protein LIG
EOF
}

extract_term() {
    local edr_path="$1"
    local primary_name="$2"
    local alternate_name="$3"
    local output_path="$4"

    if printf '%s\n0\n' "$primary_name" | "$GMX_BIN" energy -f "$edr_path" -o "$output_path" >/dev/null 2>&1; then
        return 0
    fi

    printf '%s\n0\n' "$alternate_name" | "$GMX_BIN" energy -f "$edr_path" -o "$output_path" >/dev/null
}

if [[ $# -lt 1 ]]; then
    printf 'Usage: %s RUN_DIR [RUN_DIR ...]\n' "$0" >&2
    exit 1
fi

for run_dir in "$@"; do
    prod_dir="${run_dir}/07_md_50ns"
    analysis_dir="${run_dir}/09_interaction_energy"
    mkdir -p "$analysis_dir"

    prod_tpr="${prod_dir}/md_50ns.tpr"
    prod_xtc="${prod_dir}/md_50ns.xtc"
    topol_top="${prod_dir}/topol.top"
    ndx_path="${analysis_dir}/energy_groups.ndx"
    mdp_path="${analysis_dir}/rerun_energy.mdp"
    tpr_path="${analysis_dir}/rerun_energy.tpr"
    deffnm="${analysis_dir}/rerun_energy"
    edr_path="${analysis_dir}/rerun_energy.edr"

    require_file "$prod_tpr"
    require_file "$prod_xtc"
    require_file "$topol_top"

    log "Preparing interaction-energy rerun for ${run_dir}"
    printf '13\nname 13 LIG\nq\n' | "$GMX_BIN" make_ndx -f "$prod_tpr" -o "$ndx_path" >/dev/null

    write_rerun_mdp "$mdp_path"

    "$GMX_BIN" grompp \
        -f "$mdp_path" \
        -c "$prod_dir/md_50ns.gro" \
        -p "$topol_top" \
        -n "$ndx_path" \
        -o "$tpr_path" \
        -maxwarn 1 >/dev/null

    "$GMX_BIN" mdrun \
        -s "$tpr_path" \
        -deffnm "$deffnm" \
        -rerun "$prod_xtc" \
        -ntmpi 1 \
        -ntomp "$OMP_THREADS" \
        -nb cpu \
        -pme cpu >/dev/null

    extract_term "$edr_path" 'LJ-SR:Protein-LIG' 'LJ-SR:LIG-Protein' "${analysis_dir}/lj_sr_protein_lig.xvg"
    extract_term "$edr_path" 'Coul-SR:Protein-LIG' 'Coul-SR:LIG-Protein' "${analysis_dir}/coul_sr_protein_lig.xvg"

    log "Finished interaction-energy rerun for ${run_dir}"
done
