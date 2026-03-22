#!/usr/bin/env bash
set -euo pipefail

if [[ $# -lt 1 ]]; then
  echo "Usage: $0 /mnt/d/.../system_dir [more system dirs...]" >&2
  exit 1
fi

extract_group_atoms() {
  local ndx_file="$1"
  local group_name="$2"
  python3 - "$ndx_file" "$group_name" <<'PY'
from pathlib import Path
import sys

ndx_path = Path(sys.argv[1])
target = sys.argv[2]
current = None
atoms = []
capturing = False

for raw_line in ndx_path.read_text().splitlines():
    line = raw_line.strip()
    if not line:
        continue
    if line.startswith("[") and line.endswith("]"):
        if capturing and atoms:
            break
        current = line[1:-1].strip()
        capturing = current == target
        continue
    if capturing:
        atoms.extend(line.split())

if not atoms:
    raise SystemExit(1)

print(" ".join(atoms))
PY
}

write_index_group() {
  local handle="$1"
  local group_name="$2"
  local atoms="$3"
  printf '[ %s ]\n' "$group_name" >>"$handle"
  read -r -a atom_array <<<"$atoms"
  local idx=0
  local line=""
  for atom in "${atom_array[@]}"; do
    line+=$(printf '%6s' "$atom")
    idx=$((idx + 1))
    if (( idx % 15 == 0 )); then
      printf '%s\n' "$line" >>"$handle"
      line=""
    fi
  done
  if [[ -n "$line" ]]; then
    printf '%s\n' "$line" >>"$handle"
  fi
  printf '\n' >>"$handle"
}

detect_ligand_group_name() {
  local ndx_file="$1"
  local molecule_name="$2"
  python3 - "$ndx_file" "$molecule_name" <<'PY'
from pathlib import Path
import sys

ndx_path = Path(sys.argv[1])
molecule_name = sys.argv[2]
group_names = []
for raw_line in ndx_path.read_text().splitlines():
    line = raw_line.strip()
    if line.startswith("[") and line.endswith("]"):
        group_names.append(line[1:-1].strip())

excluded = {
    "System",
    "Protein",
    "Protein-H",
    "C-alpha",
    "Backbone",
    "MainChain",
    "MainChain+Cb",
    "MainChain+H",
    "SideChain",
    "SideChain-H",
    "Prot-Masses",
    "non-Protein",
    "Other",
    "Water",
    "SOL",
    "non-Water",
    "Ion",
    "Water_and_ions",
    "NA",
    "CL",
    "K",
    "MG",
    "CA",
    "ZN",
}

if molecule_name in group_names:
    print(molecule_name)
    raise SystemExit

candidates = [name for name in group_names if name not in excluded]
if len(candidates) == 1:
    print(candidates[0])
    raise SystemExit

for preferred in ("LIG", "LGB", "UNL"):
    if preferred in candidates:
        print(preferred)
        raise SystemExit

raise SystemExit(1)
PY
}

for system_root in "$@"; do
  if [[ -d "$system_root/gmx_amber99sb_ildn_gaff2_50ns/07_md_50ns" && -f "$system_root/gmx_amber99sb_ildn_gaff2_50ns/03_build/topol.top" ]]; then
    run_base="$system_root/gmx_amber99sb_ildn_gaff2_50ns"
    run_dir="$run_base/07_md_50ns"
    top_dir="$run_base/03_build"
    traj_name="md_50ns"
    if [[ -f "$run_base/mdp/md_50ns.mdp" ]]; then
      mdp_src="$run_base/mdp/md_50ns.mdp"
    else
      mdp_src="$run_base/mdp/md_production.mdp"
    fi
  else
    run_base="$system_root"
    run_dir="$system_root"
    top_dir="$system_root"
    mdp_src="$system_root/md.mdp"
    traj_name="md_0_10"
  fi

  tpr="$run_dir/${traj_name}.tpr"
  xtc="$run_dir/${traj_name}.xtc"
  topol="$top_dir/topol.top"
  complex_gro="$top_dir/complex.gro"
  out_bind="$run_base/11_binding_energy"
  out_complex="$run_base/12_complex_only_energy"
  mkdir -p "$out_bind" "$out_complex"
  if compgen -G "$top_dir/*.itp" >/dev/null; then
    cp "$top_dir"/*.itp "$out_complex"/
  fi
  for ff_dir in "$top_dir"/*.ff; do
    if [[ -d "$ff_dir" ]]; then
      cp -r "$ff_dir" "$out_complex"/
    fi
  done

  ligand_name=$(
    awk '
      BEGIN { in_molecules = 0; count = 0 }
      /^\[ *molecules *\]/ { in_molecules = 1; next }
      in_molecules {
        if ($0 ~ /^[[:space:]]*;/ || $0 ~ /^[[:space:]]*$/) next
        count++
        if (count == 2) {
          print $1
          exit
        }
      }
    ' "$topol"
  )
  if [[ -z "$ligand_name" ]]; then
    echo "Could not determine ligand name from $topol" >&2
    exit 1
  fi

  tmp_index="$out_bind/index_tmp.ndx"
  custom_index="$out_bind/index_energy.ndx"
  printf 'r %s\nq\n' "$ligand_name" | gmx make_ndx -f "$tpr" -o "$tmp_index" >/dev/null 2>&1
  ligand_group_name=$(detect_ligand_group_name "$tmp_index" "$ligand_name")

  protein_atoms=$(extract_group_atoms "$tmp_index" "Protein")
  ligand_atoms=$(extract_group_atoms "$tmp_index" "$ligand_group_name")
  system_atoms=$(extract_group_atoms "$tmp_index" "System")
  : >"$custom_index"
  write_index_group "$custom_index" "System" "$system_atoms"
  write_index_group "$custom_index" "Protein" "$protein_atoms"
  write_index_group "$custom_index" "Ligand" "$ligand_atoms"
  write_index_group "$custom_index" "Complex" "$protein_atoms $ligand_atoms"

  bind_mdp="$out_bind/rerun_binding.mdp"
  cp "$mdp_src" "$bind_mdp"
  {
    printf '\n; Added for rerun interaction energy analysis\n'
    printf 'energygrps = Protein Ligand\n'
  } >>"$bind_mdp"

  bind_tpr="$out_bind/rerun_binding.tpr"
  bind_edr="$out_bind/rerun_binding.edr"
  gmx grompp -f "$bind_mdp" -c "$run_dir/${traj_name}.gro" -p "$topol" -n "$custom_index" -o "$bind_tpr" -maxwarn 1 >/dev/null 2>&1
  gmx mdrun -s "$bind_tpr" -rerun "$xtc" -deffnm "$out_bind/rerun_binding" >/dev/null 2>&1

  printf 'Coul-SR:Protein-Ligand\nLJ-SR:Protein-Ligand\n0\n' | gmx energy -f "$bind_edr" -o "$out_bind/binding_components.xvg" >/dev/null 2>&1

  python3 - "$out_bind/binding_components.xvg" "$out_bind/binding_energy_total.tsv" <<'PY'
from pathlib import Path
import sys

xvg_path = Path(sys.argv[1])
tsv_path = Path(sys.argv[2])
rows = []
for line in xvg_path.read_text().splitlines():
    if not line or line[0] in "#@":
        continue
    parts = line.split()
    time_ps = float(parts[0])
    coul = float(parts[1])
    lj = float(parts[2])
    rows.append((time_ps, coul, lj, coul + lj))

with tsv_path.open("w", encoding="utf-8") as handle:
    handle.write("time_ps\tcoul_sr_kj_per_mol\tlj_sr_kj_per_mol\tbinding_energy_kj_per_mol\n")
    for time_ps, coul, lj, total in rows:
        handle.write(f"{time_ps:.3f}\t{coul:.6f}\t{lj:.6f}\t{total:.6f}\n")
PY

  python3 - "$out_bind/binding_energy_total.tsv" "$out_bind/binding_energy_total.xvg" <<'PY'
from pathlib import Path
import sys

tsv_path = Path(sys.argv[1])
xvg_path = Path(sys.argv[2])
lines = tsv_path.read_text().splitlines()[1:]
with xvg_path.open("w", encoding="utf-8") as handle:
    handle.write('@    title "Protein-Ligand Binding Energy"\n')
    handle.write('@    xaxis  label "Time (ps)"\n')
    handle.write('@    yaxis  label "Binding Energy (kJ/mol)"\n')
    for line in lines:
        time_ps, _coul, _lj, total = line.split("\t")
        handle.write(f"{time_ps} {total}\n")
PY

  complex_top="$out_complex/topol_complex_only.top"
  awk '
    BEGIN { in_molecules = 0; kept = 0 }
    /^\[ *molecules *\]/ {
      in_molecules = 1
      print
      next
    }
    in_molecules {
      if ($0 ~ /^[[:space:]]*;/ || $0 ~ /^[[:space:]]*$/) {
        print
        next
      }
      kept++
      if (kept <= 2) {
        print
      }
      next
    }
    { print }
  ' "$topol" >"$complex_top"

  complex_mdp="$out_complex/rerun_complex_only.mdp"
  cp "$mdp_src" "$complex_mdp"
  complex_tpr="$out_complex/rerun_complex_only.tpr"
  complex_xtc="$out_complex/complex_only.xtc"
  complex_edr="$out_complex/rerun_complex_only.edr"

  printf '3\n' | gmx trjconv -s "$tpr" -f "$xtc" -n "$custom_index" -o "$complex_xtc" -pbc mol -ur compact >/dev/null 2>&1
  gmx grompp -f "$complex_mdp" -c "$complex_gro" -p "$complex_top" -o "$complex_tpr" -maxwarn 1 >/dev/null 2>&1
  gmx mdrun -s "$complex_tpr" -rerun "$complex_xtc" -deffnm "$out_complex/rerun_complex_only" >/dev/null 2>&1
  printf 'Potential\n0\n' | gmx energy -f "$complex_edr" -o "$out_complex/complex_only_potential.xvg" >/dev/null 2>&1

  rm -f "$tmp_index"
done
