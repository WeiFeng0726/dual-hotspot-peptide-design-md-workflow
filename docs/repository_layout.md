# Repository Layout

This repository contains both reusable workflow code and project-specific execution wrappers.

## Design Principles

- Keep MD utilities separate from peptide-design utilities.
- Keep analysis scripts separate from launch scripts.
- Preserve the original project wrappers instead of rewriting the entire production workflow.
- Document which scripts are generic and which are tied to the original filesystem layout.

## Structure

### `workflow/md/`

Primary Bash entry points for GROMACS-related tasks.

Use these when you need to:

- prepare or continue an MD workflow
- extract energies
- perform trajectory-level processing in shell form

### `workflow/peptide_design/`

Primary Bash entry points for Proteina-Complexa campaigns.

These scripts orchestrate:

- candidate generation
- AF2 rescoring
- secondary filtering
- repeated sampling
- background launch behavior

### `analysis/md/`

Python plotting scripts used after MD calculations finish.

These scripts are focused on:

- RMSD summaries
- energy summaries
- publication-ready figures

### `analysis/peptide_design/`

Python scripts used after candidate generation finishes.

These scripts handle:

- rescoring
- hotspot-aware filtering
- campaign aggregation
- pairwise peptide-combination ranking

### `reporting/`

PPT generation scripts that convert final project outputs into a presentation deck.

### `examples/project_specific/`

Files kept for traceability but not recommended as the first entry points for new work.

## Generic vs Project-Specific

### More reusable

- `workflow/md/*`
- `analysis/md/*`
- `analysis/peptide_design/evaluate_joint_peptide_combos.py`
- `analysis/peptide_design/secondary_screen_complexa.py`

### More project-specific

- the `launch_*` wrappers
- batch-expansion control scripts
- PPT generators
- example helper files in `examples/project_specific/`

## What Is Not Included

- raw trajectories
- model checkpoints
- large result tables
- exported final figures
- final PPT files
