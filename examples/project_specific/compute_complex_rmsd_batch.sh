#!/usr/bin/env bash
set -euo pipefail

GMX_BIN="${GMX_BIN:-/usr/local/gmx2025.4/bin/gmx}"

run_one() {
  local tpr="$1"
  local xtc="$2"
  local ndx="$3"
  local out="$4"
  printf 'Complex\nComplex\n' | "$GMX_BIN" rms -s "$tpr" -f "$xtc" -n "$ndx" -o "$out" -tu ns >/dev/null
}

run_one '/mnt/d/Project/cw/lf/MD/s_cdtc-2_c_model/gmx_amber99sb_ildn_gaff2_50ns/07_md_50ns/md_50ns.tpr' '/mnt/d/Project/cw/lf/MD/s_cdtc-2_c_model/gmx_amber99sb_ildn_gaff2_50ns/08_rmsd/md_pbcfix.xtc' '/mnt/d/Project/cw/lf/MD/s_cdtc-2_c_model/gmx_amber99sb_ildn_gaff2_50ns/11_binding_energy/index_energy.ndx' '/mnt/d/Project/cw/lf/MD/s_cdtc-2_c_model/gmx_amber99sb_ildn_gaff2_50ns/08_rmsd/rmsd_complex_allatom_pbcfix.xvg'
run_one '/mnt/d/Project/cw/lf/MD/s_cdtc-2_2c_model/md_0_10.tpr' '/mnt/d/Project/cw/lf/MD/s_cdtc-2_2c_model/08_rmsd/md_pbcfix.xtc' '/mnt/d/Project/cw/lf/MD/s_cdtc-2_2c_model/11_binding_energy/index_energy.ndx' '/mnt/d/Project/cw/lf/MD/s_cdtc-2_2c_model/08_rmsd/rmsd_complex_allatom_pbcfix.xvg'
run_one '/mnt/d/Project/cw/lf/MD/s_cdtc-2_3c-2_model/gmx_amber99sb_ildn_gaff2_50ns/07_md_50ns/md_50ns.tpr' '/mnt/d/Project/cw/lf/MD/s_cdtc-2_3c-2_model/gmx_amber99sb_ildn_gaff2_50ns/08_rmsd/md_pbcfix.xtc' '/mnt/d/Project/cw/lf/MD/s_cdtc-2_3c-2_model/gmx_amber99sb_ildn_gaff2_50ns/11_binding_energy/index_energy.ndx' '/mnt/d/Project/cw/lf/MD/s_cdtc-2_3c-2_model/gmx_amber99sb_ildn_gaff2_50ns/08_rmsd/rmsd_complex_allatom_pbcfix.xvg'
run_one '/mnt/d/Project/cw/lf/MD/s_cdtc-2_214a_c_model/md_0_10.tpr' '/mnt/d/Project/cw/lf/MD/s_cdtc-2_214a_c_model/08_rmsd/md_pbcfix.xtc' '/mnt/d/Project/cw/lf/MD/s_cdtc-2_214a_c_model/11_binding_energy/index_energy.ndx' '/mnt/d/Project/cw/lf/MD/s_cdtc-2_214a_c_model/08_rmsd/rmsd_complex_allatom_pbcfix.xvg'
run_one '/mnt/d/Project/cw/lf/MD/s_cdtc-2_214a_2c_model/gmx_amber99sb_ildn_gaff2_50ns/07_md_50ns/md_50ns.tpr' '/mnt/d/Project/cw/lf/MD/s_cdtc-2_214a_2c_model/gmx_amber99sb_ildn_gaff2_50ns/08_rmsd/md_pbcfix.xtc' '/mnt/d/Project/cw/lf/MD/s_cdtc-2_214a_2c_model/gmx_amber99sb_ildn_gaff2_50ns/11_binding_energy/index_energy.ndx' '/mnt/d/Project/cw/lf/MD/s_cdtc-2_214a_2c_model/gmx_amber99sb_ildn_gaff2_50ns/08_rmsd/rmsd_complex_allatom_pbcfix.xvg'
run_one '/mnt/d/Project/cw/lf/MD/s_cdtc-2_214g_c_model/gmx_amber99sb_ildn_gaff2_50ns/07_md_50ns/md_50ns.tpr' '/mnt/d/Project/cw/lf/MD/s_cdtc-2_214g_c_model/gmx_amber99sb_ildn_gaff2_50ns/08_rmsd/md_pbcfix.xtc' '/mnt/d/Project/cw/lf/MD/s_cdtc-2_214g_c_model/gmx_amber99sb_ildn_gaff2_50ns/11_binding_energy/index_energy.ndx' '/mnt/d/Project/cw/lf/MD/s_cdtc-2_214g_c_model/gmx_amber99sb_ildn_gaff2_50ns/08_rmsd/rmsd_complex_allatom_pbcfix.xvg'
run_one '/mnt/d/Project/cw/lf/MD/s_cdtc-2_214g_2c_model/gmx_amber99sb_ildn_gaff2_50ns/07_md_50ns/md_50ns.tpr' '/mnt/d/Project/cw/lf/MD/s_cdtc-2_214g_2c_model/gmx_amber99sb_ildn_gaff2_50ns/08_rmsd/md_pbcfix.xtc' '/mnt/d/Project/cw/lf/MD/s_cdtc-2_214g_2c_model/gmx_amber99sb_ildn_gaff2_50ns/11_binding_energy/index_energy.ndx' '/mnt/d/Project/cw/lf/MD/s_cdtc-2_214g_2c_model/gmx_amber99sb_ildn_gaff2_50ns/08_rmsd/rmsd_complex_allatom_pbcfix.xvg'
run_one '/mnt/d/Project/cw/lf/MD/s_cdtc-2_214g_3c-2_model/gmx_amber99sb_ildn_gaff2_50ns/07_md_50ns/md_50ns.tpr' '/mnt/d/Project/cw/lf/MD/s_cdtc-2_214g_3c-2_model/gmx_amber99sb_ildn_gaff2_50ns/08_rmsd/md_pbcfix.xtc' '/mnt/d/Project/cw/lf/MD/s_cdtc-2_214g_3c-2_model/gmx_amber99sb_ildn_gaff2_50ns/11_binding_energy/index_energy.ndx' '/mnt/d/Project/cw/lf/MD/s_cdtc-2_214g_3c-2_model/gmx_amber99sb_ildn_gaff2_50ns/08_rmsd/rmsd_complex_allatom_pbcfix.xvg'
