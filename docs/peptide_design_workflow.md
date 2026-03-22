# Peptide Design Workflow

This workflow uses Proteina-Complexa to generate and filter short peptide candidates around MD-supported hotspot residues.

## Purpose

The peptide-design pipeline answers:

- Which 8-12 aa peptides can plausibly target a chosen hotspot?
- Which candidates survive stricter AF2-based rescoring?
- Can peptides from two separate hotspots be combined without obvious geometric conflicts?

## Workflow Stages

### 1. Candidate generation

Entry script:

- `workflow/peptide_design/run_complexa_design_campaign.sh`

This script runs a Proteina-Complexa generation campaign for one hotspot definition.

### 2. AF2 rescoring

Primary analysis script:

- `analysis/peptide_design/rescore_existing_af2.py`

This script rescans generated structures and exports:

- `af2_rescore_progress.csv`
- `af2_rescore_ranked.csv`
- optional top-N summaries

### 3. Secondary filtering

Primary analysis script:

- `analysis/peptide_design/secondary_screen_complexa.py`

This script combines:

- hotspot coverage
- hotspot distance
- binder pLDDT
- iPTM
- minimum interface PAE

to define a stricter shortlist.

### 4. Repeated expansion for a secondary hotspot

Primary workflow script:

- `workflow/peptide_design/run_b348_until_10.sh`

Supporting analysis scripts:

- `analysis/peptide_design/aggregate_b348_candidates.py`
- `analysis/peptide_design/finalize_b348_round_then_joint.py`

This part of the workflow is designed to keep sampling until enough B348 candidates are available for downstream pairing.

### 5. Joint-use peptide pairing

Primary analysis script:

- `analysis/peptide_design/evaluate_joint_peptide_combos.py`

This script:

- aligns complexes to a shared target frame
- merges peptide candidates
- measures peptide-peptide distance and clashes
- evaluates hotspot retention
- ranks peptide pairs by static compatibility

## Launch Wrappers

The `launch_*` scripts in `workflow/peptide_design/` are preserved because they were used in production, but they are convenience wrappers rather than portable package entry points.

## Key Filtering Logic

The peptide-design analysis treats a candidate as stronger when it shows:

- better hotspot retention
- smaller hotspot distance
- stronger AF2-based structural confidence
- stronger interface confidence
- cleaner peptide-peptide geometry in the joint-use stage

## Practical Note

The current scripts intentionally preserve the original path conventions from the working project. If you adapt this repository for a new system, the first step should be to externalize those paths into environment variables or configuration files.
