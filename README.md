# Dual-Hotspot Peptide Design and MD Workflow

This repository organizes the molecular dynamics re-analysis, hotspot-guided peptide design, pairwise peptide-combination screening, and report-generation scripts used in a dual-hotspot peptide design project.

The project logic is:

1. Re-analyze MD trajectories to identify stable interface hotspots.
2. Run hotspot-guided 8-12 aa peptide design campaigns with Proteina-Complexa.
3. Rescore and filter peptide candidates with AF2-derived metrics and hotspot-contact criteria.
4. Expand a secondary hotspot campaign until enough candidates are collected.
5. Screen all peptide pairs for potential joint use.
6. Generate a presentation-ready PPT report from the final results.

The repository name is intended to match that purpose directly: a reproducible workflow for **dual-hotspot peptide design plus MD-based support analysis**.

## Repository Layout

```text
.
|-- workflow/
|   |-- md/
|   `-- peptide_design/
|-- analysis/
|   |-- md/
|   `-- peptide_design/
|-- reporting/
|-- docs/
|-- examples/
|   `-- project_specific/
|-- requirements.txt
|-- .gitignore
`-- README.md
```

## Directory Summary

### `workflow/md/`

Bash workflows for GROMACS-based post-processing and interaction-energy extraction.

- `compute_complex_and_binding_energy.sh`
- `extract_system_potential_energy.sh`
- `generic_extract_system_energy.sh`
- `generic_rmsd_pbc.sh`
- `md_interaction_energy_pipeline.sh`
- `setup_pl_md_until_npt.sh`

### `workflow/peptide_design/`

Bash workflows and launch wrappers for Proteina-Complexa peptide-design campaigns.

- large-scale hotspot campaigns
- AF2 rescoring pipelines
- B348 aggregate-until-target automation
- dual-hotspot campaign launchers

### `analysis/md/`

Plotting scripts for RMSD and energy analysis.

### `analysis/peptide_design/`

Python analysis scripts for:

- AF2 rescoring
- hotspot-aware secondary screening
- aggregate ranking across repeated campaigns
- pairwise peptide-combination screening
- automated handoff from final B348 expansion to joint screening

### `reporting/`

PPT-building scripts used to generate the project presentation deck from the MD and design outputs.

### `examples/project_specific/`

Small helper files that are useful for the original project context but are not meant to be the primary reusable entry points.

## Main Workflows

### 1. MD Re-analysis

Use the scripts in `workflow/md/` and `analysis/md/` to:

- extract system potential energy
- compute complex-only and interaction-energy terms
- perform PBC-corrected RMSD processing
- create publication-style RMSD and energy plots

### 2. Single-Hotspot Peptide Design

Use `workflow/peptide_design/run_complexa_design_campaign.sh` to:

- generate peptide candidates with Proteina-Complexa
- rescore candidates with AF2
- apply secondary filtering based on hotspot coverage, distance, and confidence metrics

### 3. Secondary-Hotspot Expansion

Use `workflow/peptide_design/run_b348_until_10.sh` together with:

- `analysis/peptide_design/aggregate_b348_candidates.py`
- `analysis/peptide_design/finalize_b348_round_then_joint.py`

This workflow keeps sampling additional batches until enough B348 candidates are collected, or until the final stopping logic is triggered.

### 4. Joint Peptide-Pair Screening

Use `analysis/peptide_design/evaluate_joint_peptide_combos.py` to:

- merge peptide complexes into a common target frame
- measure peptide-peptide clashes
- preserve hotspot-specific contact logic
- rank all peptide pairs by static compatibility

### 5. Reporting

Use the scripts in `reporting/` to generate PPT reports that summarize:

- MD stability
- hotspot evidence
- peptide family structure
- filtered candidate sets
- joint-pair rankings

## Software Requirements

### Core external tools

- WSL or Linux shell
- GROMACS
- Proteina-Complexa
- Python 3.10+

### Python packages for plotting/reporting

Install:

```bash
pip install -r requirements.txt
```

The peptide-design rescoring scripts additionally require a working Proteina-Complexa environment because they import `proteinfoundation`.

## Important Reproducibility Notes

- Several peptide-design workflows are still intentionally tied to the original project path conventions under `D:\Project\cw\ly` and `/mnt/d/Project/cw/ly`.
- WSL-facing scripts assume a Proteina-Complexa installation under `/home/fatcat/Proteina-Complexa` unless overridden.
- The launch wrappers are preserved because they were part of the actual production workflow, but they should be treated as **project-specific execution helpers**, not generic package-style interfaces.

## Recommended Reading Order

1. `docs/repository_layout.md`
2. `docs/md_workflow.md`
3. `docs/peptide_design_workflow.md`
4. `docs/reporting_workflow.md`

## Current Scope

This repository is focused on **code organization and workflow reproducibility**. It does not bundle large trajectory files, model checkpoints, or result directories.
