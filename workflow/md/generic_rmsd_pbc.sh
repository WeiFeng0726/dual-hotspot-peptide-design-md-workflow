#!/usr/bin/env bash
set -euo pipefail

GMX_BIN="${GMX_BIN:-/usr/local/gmx2025.4/bin/gmx}"

Tpr=""
Xtc=""
OutDir=""

while [[ $# -gt 0 ]]; do
    case "$1" in
        --tpr) Tpr="$2"; shift 2 ;;
        --xtc) Xtc="$2"; shift 2 ;;
        --out-dir) OutDir="$2"; shift 2 ;;
        *)
            printf 'Unknown argument: %s\n' "$1" >&2
            exit 1
            ;;
    esac
done

[[ -f "$Tpr" ]] || { printf 'Missing TPR: %s\n' "$Tpr" >&2; exit 1; }
[[ -f "$Xtc" ]] || { printf 'Missing XTC: %s\n' "$Xtc" >&2; exit 1; }
[[ -n "$OutDir" ]] || { printf '--out-dir is required\n' >&2; exit 1; }

mkdir -p "$OutDir"
cd "$OutDir"

printf '12 & ! a H*\nq\n' | "$GMX_BIN" make_ndx -f "$Tpr" -o rmsd_analysis.ndx >/dev/null
printf '0\n' | "$GMX_BIN" trjconv -s "$Tpr" -f "$Xtc" -o md_nojump.xtc -pbc nojump >/dev/null
printf '1\n0\n' | "$GMX_BIN" trjconv -s "$Tpr" -f md_nojump.xtc -o md_pbcfix.xtc -pbc mol -center -ur compact >/dev/null
printf 'Backbone\nBackbone\n' | "$GMX_BIN" rms -s "$Tpr" -f md_pbcfix.xtc -n rmsd_analysis.ndx -o rmsd_backbone_pbcfix.xvg -tu ns >/dev/null
printf 'Backbone\nOther_&_!H*\n' | "$GMX_BIN" rms -s "$Tpr" -f md_pbcfix.xtc -n rmsd_analysis.ndx -o rmsd_ligand_heavy_fit_backbone_pbcfix.xvg -tu ns >/dev/null
