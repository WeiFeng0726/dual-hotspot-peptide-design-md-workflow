# MD Workflow

This part of the repository supports MD-based interpretation of a protein interface before peptide design.

## Purpose

The MD workflow is used to answer the following questions:

- Is the complex stable enough for interface interpretation?
- Which residues remain involved in persistent contacts or hydrogen bonds?
- Can the equilibrium window support hotspot selection?

## Main Shell Scripts

### `workflow/md/setup_pl_md_until_npt.sh`

Project setup through NPT for a protein-ligand style workflow.

### `workflow/md/extract_system_potential_energy.sh`

Extract whole-system potential energy from a finished GROMACS run.

### `workflow/md/generic_extract_system_energy.sh`

Generic helper for energy extraction from `.edr` files.

### `workflow/md/generic_rmsd_pbc.sh`

Trajectory correction helper for RMSD workflows after PBC handling.

### `workflow/md/compute_complex_and_binding_energy.sh`

Extract complex-only potential and interaction-energy terms from a finished run.

### `workflow/md/md_interaction_energy_pipeline.sh`

Convenience pipeline for interaction-energy reruns on one or more systems.

## Main Plotting Scripts

### `analysis/md/plot_rmsd_analysis.py`

Summarizes RMSD behavior and estimates whether the late trajectory window looks near-stationary.

### `analysis/md/plot_system_energy.py`

Creates system-level energy overlays for multiple runs.

### `analysis/md/plot_interaction_energy.py`

Builds interaction-energy overlays and exports summary tables.

### `analysis/md/plot_complex_rmsd_grid.py`

Creates a panel comparison of complex RMSD across multiple systems.

## Expected Inputs

- `.tpr`
- `.xtc`
- `.edr`
- `.gro`
- `.top`
- optional index files

## Typical Outcome

The MD workflow should produce:

- RMSD plots
- energy tables
- interaction-energy plots
- a justified equilibrium window
- hotspot candidates supported by dynamic evidence
