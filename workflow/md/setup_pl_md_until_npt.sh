#!/usr/bin/env bash
set -euo pipefail

GMX_BIN="/usr/local/gmx2025.4/bin/gmx"
OBABEL_BIN="/home/fatcat/miniconda3/bin/obabel"
ACPYPE_BIN="/home/fatcat/miniconda3/bin/acpype"

INPUT_DIR="/mnt/d/Project/cw/lf/MD/s_cdtc-2_214a_2c_model"
RUN_DIR="${INPUT_DIR}/gmx_amber99sb_ildn_gaff2_50ns"

PROTEIN_PDB="${INPUT_DIR}/protein.pdb"
LIGAND_MOL2="${INPUT_DIR}/ligand_B_fixed.mol2"

INPUT_STAGE="${RUN_DIR}/00_input"
LIG_STAGE="${RUN_DIR}/01_ligand_param"
TOP_STAGE="${RUN_DIR}/02_topology"
BUILD_STAGE="${RUN_DIR}/03_build"
EM_STAGE="${RUN_DIR}/04_em"
NVT_STAGE="${RUN_DIR}/05_nvt"
NPT_STAGE="${RUN_DIR}/06_npt"
MDP_STAGE="${RUN_DIR}/mdp"

mkdir -p "${INPUT_STAGE}" "${LIG_STAGE}" "${TOP_STAGE}" "${BUILD_STAGE}" "${EM_STAGE}" "${NVT_STAGE}" "${NPT_STAGE}" "${MDP_STAGE}"

cp "${PROTEIN_PDB}" "${INPUT_STAGE}/protein.pdb"
cp "${LIGAND_MOL2}" "${INPUT_STAGE}/ligand_B_fixed.mol2"

cat > "${MDP_STAGE}/ions.mdp" <<'EOF'
integrator      = steep
nsteps          = 500
emtol           = 1000.0
emstep          = 0.01
nstlist         = 20
cutoff-scheme   = Verlet
ns_type         = grid
coulombtype     = PME
rcoulomb        = 1.0
rvdw            = 1.0
pbc             = xyz
EOF

cat > "${MDP_STAGE}/em.mdp" <<'EOF'
integrator      = steep
nsteps          = 50000
emtol           = 1000.0
emstep          = 0.01
nstenergy       = 500
nstlog          = 500
cutoff-scheme   = Verlet
nstlist         = 20
ns_type         = grid
coulombtype     = PME
rcoulomb        = 1.0
rvdw            = 1.0
pbc             = xyz
constraints     = h-bonds
constraint-algorithm = lincs
EOF

cat > "${MDP_STAGE}/nvt.mdp" <<'EOF'
define                  = -DPOSRES -DPOSRES_LIG
integrator              = md
nsteps                  = 50000
dt                      = 0.002
nstxout-compressed      = 5000
nstenergy               = 1000
nstlog                  = 1000
continuation            = no
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
tcoupl                  = V-rescale
tc-grps                 = System
tau_t                   = 0.1
ref_t                   = 300
pcoupl                  = no
pbc                     = xyz
gen_vel                 = yes
gen_temp                = 300
gen_seed                = -1
DispCorr                = EnerPres
EOF

cat > "${MDP_STAGE}/npt.mdp" <<'EOF'
define                  = -DPOSRES -DPOSRES_LIG
integrator              = md
nsteps                  = 50000
dt                      = 0.002
nstxout-compressed      = 5000
nstenergy               = 1000
nstlog                  = 1000
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
tcoupl                  = V-rescale
tc-grps                 = System
tau_t                   = 0.1
ref_t                   = 300
pcoupl                  = C-rescale
pcoupltype              = isotropic
tau_p                   = 2.0
ref_p                   = 1.0
compressibility         = 4.5e-5
refcoord_scaling        = com
pbc                     = xyz
gen_vel                 = no
DispCorr                = EnerPres
EOF

cat > "${MDP_STAGE}/md_50ns.mdp" <<'EOF'
integrator              = md
nsteps                  = 25000000
dt                      = 0.002
nstxout-compressed      = 5000
nstenergy               = 5000
nstlog                  = 5000
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
tcoupl                  = V-rescale
tc-grps                 = System
tau_t                   = 0.1
ref_t                   = 300
pcoupl                  = Parrinello-Rahman
pcoupltype              = isotropic
tau_p                   = 2.0
ref_p                   = 1.0
compressibility         = 4.5e-5
pbc                     = xyz
gen_vel                 = no
DispCorr                = EnerPres
EOF

cd "${LIG_STAGE}"
if [[ ! -f "LGB.acpype/LGB_GMX.itp" ]]; then
    cp "${INPUT_STAGE}/ligand_B_fixed.mol2" .
    "${OBABEL_BIN}" ligand_B_fixed.mol2 -O ligand_B_H.mol2 -h
    "${ACPYPE_BIN}" -i ligand_B_H.mol2 -b LGB -a gaff2 -c bcc -n 0
fi

cd "${TOP_STAGE}"
"${GMX_BIN}" pdb2gmx \
    -f "${INPUT_STAGE}/protein.pdb" \
    -o protein_processed.gro \
    -p topol.top \
    -i posre.itp \
    -ff amber99sb-ildn \
    -water tip3p \
    -ignh

cp "${LIG_STAGE}/LGB.acpype/LGB_GMX.itp" "${TOP_STAGE}/LGB_GMX.itp"
cp "${LIG_STAGE}/LGB.acpype/LGB_GMX.gro" "${TOP_STAGE}/LGB_GMX.gro"
cp "${LIG_STAGE}/LGB.acpype/posre_LGB.itp" "${TOP_STAGE}/posre_LGB.itp"

python3 - <<'PY'
from pathlib import Path

top_path = Path("/mnt/d/Project/cw/lf/MD/s_cdtc-2_214a_2c_model/gmx_amber99sb_ildn_gaff2_50ns/02_topology/topol.top")
lines = top_path.read_text().splitlines()
forcefield_idx = next(i for i, line in enumerate(lines) if line.strip() == '#include "amber99sb-ildn.ff/forcefield.itp"')
protein_start = next(i for i, line in enumerate(lines) if line.strip() == "[ moleculetype ]")
protein_posre_idx = next(i for i, line in enumerate(lines) if line.strip() == '; Include Position restraint file')
water_idx = next(i for i, line in enumerate(lines) if line.strip() == '; Include water topology')
system_idx = next(i for i, line in enumerate(lines) if line.strip() == "[ system ]")

protein_itp = Path("/mnt/d/Project/cw/lf/MD/s_cdtc-2_214a_2c_model/gmx_amber99sb_ildn_gaff2_50ns/02_topology/topol_Protein_chain_A.itp")
protein_itp.write_text("\n".join(lines[protein_start:protein_posre_idx]) + "\n")

water_and_ions_block = lines[water_idx:system_idx]

new_topology = [
    "; Combined protein-ligand topology",
    "; Protein generated by pdb2gmx, ligand generated by acpype/antechamber",
    "",
    "; Include forcefield parameters",
    lines[forcefield_idx],
    '#include "LGB_GMX.itp"',
    '#ifdef POSRES_LIG',
    '#include "posre_LGB.itp"',
    '#endif',
    '#include "topol_Protein_chain_A.itp"',
    '#ifdef POSRES',
    '#include "posre.itp"',
    '#endif',
    "",
]
new_topology.extend(water_and_ions_block)
new_topology.extend(
    [
        "",
        "[ system ]",
        "; Name",
        "Protein-Ligand complex",
        "",
        "[ molecules ]",
        "; Compound        #mols",
        "Protein_chain_A     1",
        "LGB                 1",
    ]
)

top_path.write_text("\n".join(new_topology) + "\n")
PY

python3 - <<'PY'
from pathlib import Path

protein = Path("/mnt/d/Project/cw/lf/MD/s_cdtc-2_214a_2c_model/gmx_amber99sb_ildn_gaff2_50ns/02_topology/protein_processed.gro")
ligand = Path("/mnt/d/Project/cw/lf/MD/s_cdtc-2_214a_2c_model/gmx_amber99sb_ildn_gaff2_50ns/02_topology/LGB_GMX.gro")
out = Path("/mnt/d/Project/cw/lf/MD/s_cdtc-2_214a_2c_model/gmx_amber99sb_ildn_gaff2_50ns/03_build/complex.gro")

def read_gro(path):
    lines = path.read_text().splitlines()
    natoms = int(lines[1].strip())
    atoms = lines[2:2 + natoms]
    box = lines[2 + natoms]
    return lines[0], atoms, box

_, protein_atoms, protein_box = read_gro(protein)
_, ligand_atoms, _ = read_gro(ligand)

out.write_text(
    "Protein-Ligand complex\n"
    f"{len(protein_atoms) + len(ligand_atoms):5d}\n"
    + "\n".join(protein_atoms + ligand_atoms)
    + "\n"
    + protein_box
    + "\n"
)
PY

cp "${TOP_STAGE}/topol.top" "${BUILD_STAGE}/topol.top"
cp "${TOP_STAGE}/topol_Protein_chain_A.itp" "${BUILD_STAGE}/topol_Protein_chain_A.itp"
cp "${TOP_STAGE}/posre.itp" "${BUILD_STAGE}/posre.itp"
cp "${TOP_STAGE}/LGB_GMX.itp" "${BUILD_STAGE}/LGB_GMX.itp"
cp "${TOP_STAGE}/posre_LGB.itp" "${BUILD_STAGE}/posre_LGB.itp"

cd "${BUILD_STAGE}"
"${GMX_BIN}" editconf -f complex.gro -o boxed.gro -c -d 1.0 -bt dodecahedron
"${GMX_BIN}" solvate -cp boxed.gro -cs spc216.gro -o solvated.gro -p topol.top
"${GMX_BIN}" grompp -f "${MDP_STAGE}/ions.mdp" -c solvated.gro -p topol.top -o ions.tpr -maxwarn 1
printf "SOL\n" | "${GMX_BIN}" genion -s ions.tpr -o solv_ions.gro -p topol.top -pname NA -nname CL -neutral -conc 0.15

cd "${EM_STAGE}"
cp "${BUILD_STAGE}/solv_ions.gro" .
cp "${BUILD_STAGE}/topol.top" .
cp "${BUILD_STAGE}/topol_Protein_chain_A.itp" .
cp "${BUILD_STAGE}/posre.itp" .
cp "${BUILD_STAGE}/LGB_GMX.itp" .
cp "${BUILD_STAGE}/posre_LGB.itp" .
"${GMX_BIN}" grompp -f "${MDP_STAGE}/em.mdp" -c solv_ions.gro -p topol.top -o em.tpr
"${GMX_BIN}" mdrun -deffnm em

cd "${NVT_STAGE}"
cp "${EM_STAGE}/em.gro" .
cp "${EM_STAGE}/topol.top" .
cp "${EM_STAGE}/topol_Protein_chain_A.itp" .
cp "${EM_STAGE}/posre.itp" .
cp "${EM_STAGE}/LGB_GMX.itp" .
cp "${EM_STAGE}/posre_LGB.itp" .
"${GMX_BIN}" grompp -f "${MDP_STAGE}/nvt.mdp" -c em.gro -r em.gro -p topol.top -o nvt.tpr
"${GMX_BIN}" mdrun -deffnm nvt

cd "${NPT_STAGE}"
cp "${NVT_STAGE}/nvt.gro" .
cp "${NVT_STAGE}/nvt.cpt" .
cp "${NVT_STAGE}/topol.top" .
cp "${NVT_STAGE}/topol_Protein_chain_A.itp" .
cp "${NVT_STAGE}/posre.itp" .
cp "${NVT_STAGE}/LGB_GMX.itp" .
cp "${NVT_STAGE}/posre_LGB.itp" .
"${GMX_BIN}" grompp -f "${MDP_STAGE}/npt.mdp" -c nvt.gro -r nvt.gro -t nvt.cpt -p topol.top -o npt.tpr
"${GMX_BIN}" mdrun -deffnm npt

printf "Temperature\n0\n" | "${GMX_BIN}" energy -f "${NVT_STAGE}/nvt.edr" -o "${NVT_STAGE}/temperature.xvg" > "${NVT_STAGE}/temperature_summary.txt"
printf "Pressure\n0\n" | "${GMX_BIN}" energy -f "${NPT_STAGE}/npt.edr" -o "${NPT_STAGE}/pressure.xvg" > "${NPT_STAGE}/pressure_summary.txt"
printf "Density\n0\n" | "${GMX_BIN}" energy -f "${NPT_STAGE}/npt.edr" -o "${NPT_STAGE}/density.xvg" > "${NPT_STAGE}/density_summary.txt"
printf "Potential\n0\n" | "${GMX_BIN}" energy -f "${EM_STAGE}/em.edr" -o "${EM_STAGE}/potential.xvg" > "${EM_STAGE}/potential_summary.txt"

python3 - <<'PY'
from pathlib import Path

run_dir = Path("/mnt/d/Project/cw/lf/MD/s_cdtc-2_214a_2c_model/gmx_amber99sb_ildn_gaff2_50ns")

def parse_xvg(path):
    xs = []
    ys = []
    for line in path.read_text().splitlines():
        if not line or line[0] in "#@":
            continue
        parts = line.split()
        xs.append(float(parts[0]))
        ys.append(float(parts[1]))
    return xs, ys

def summarize(path):
    xs, ys = parse_xvg(path)
    n = len(ys)
    tail = max(1, n // 5)
    return {
        "points": n,
        "start": ys[0],
        "end": ys[-1],
        "avg": sum(ys) / n,
        "tail_avg": sum(ys[-tail:]) / tail,
    }

temp = summarize(run_dir / "05_nvt" / "temperature.xvg")
press = summarize(run_dir / "06_npt" / "pressure.xvg")
density = summarize(run_dir / "06_npt" / "density.xvg")
potential = summarize(run_dir / "04_em" / "potential.xvg")

summary = [
    "Protein-ligand equilibration summary",
    "",
    f"Run directory: {run_dir}",
    "",
    "Assumptions:",
    "- Protein force field: amber99sb-ildn",
    "- Ligand force field: GAFF2 with AM1-BCC charges via acpype/antechamber",
    "- Water model: TIP3P",
    "- Box type: dodecahedron, 1.0 nm padding",
    "- Equilibration: 100 ps NVT + 100 ps NPT at 300 K",
    "",
    f"EM potential energy: start {potential['start']:.3f}, end {potential['end']:.3f}, last-window avg {potential['tail_avg']:.3f}",
    f"NVT temperature (K): start {temp['start']:.3f}, end {temp['end']:.3f}, avg {temp['avg']:.3f}, last-window avg {temp['tail_avg']:.3f}",
    f"NPT pressure (bar): start {press['start']:.3f}, end {press['end']:.3f}, avg {press['avg']:.3f}, last-window avg {press['tail_avg']:.3f}",
    f"NPT density (kg/m^3): start {density['start']:.3f}, end {density['end']:.3f}, avg {density['avg']:.3f}, last-window avg {density['tail_avg']:.3f}",
    "",
    "The production 50 ns mdp file is prepared in mdp/md_50ns.mdp but production has not been started.",
]

(run_dir / "equilibration_summary.txt").write_text("\n".join(summary) + "\n")
PY

touch "${RUN_DIR}/READY_FOR_REVIEW"
echo "Equilibration finished. Review files under ${RUN_DIR} before starting production MD."
