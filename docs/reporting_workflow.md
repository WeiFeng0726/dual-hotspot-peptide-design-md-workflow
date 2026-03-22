# Reporting Workflow

This repository includes scripts used to transform the final MD and peptide-design outputs into a PPT presentation.

## Included Scripts

### `reporting/build_ly_md_complexa_report.py`

The earlier report builder used during the first round of analysis.

### `reporting/build_ly_md_complexa_report_20260322.py`

The updated report builder that incorporates:

- expanded 404/405 peptide results
- expanded B348 peptide results
- dual-peptide combination screening
- fully Chinese presentation output for the project meeting

## Why Reporting Scripts Are Kept

These files are not just cosmetic. They encode:

- how the final figures were assembled
- which summary statistics were emphasized
- how the project story was structured for communication

## Dependencies

The reporting scripts require:

- `matplotlib`
- `numpy`
- `pandas`
- `seaborn`
- `biopython`
- `python-pptx`

## Output Philosophy

The reporting workflow was designed to:

- preserve the logical chain from MD evidence to hotspot choice
- separate the importance of the primary and secondary hotspots
- explain why the selected metrics support prioritizing some candidates over others
- present single-peptide and paired-peptide outcomes in one coherent narrative
