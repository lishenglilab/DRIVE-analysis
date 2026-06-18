# DRIVE-analysis

`DRIVE-analysis` contains the downstream R and Python scripts used to summarize DRIVE prediction results, generate manuscript figures, and process validation analyses.

The DRIVE model framework itself is maintained separately:

- https://github.com/lishenglilab/DRIVE

## Overview

This repository is currently organized in a flat, figure-oriented style rather than a nested `scripts/` layout.

- Main-figure scripts such as `fig1a.R`, `fig2.R`, `fig4f.R`, `fig5c.R` are placed in the repository root.
- Supplementary scripts such as `supfig_s7.R`, `supfig_s9.py`, `supfig_s12.R` are also placed in the repository root.
- Input tables are stored alongside the scripts unless noted otherwise in the script header.

## Scope

This repository covers downstream analysis only, including:

- benchmark comparison across mixed, cell-blind, and drug-blind settings;
- ensemble analysis;
- residual correlation analysis;
- natural product prediction summaries;
- CRC candidate prioritization;
- experimental validation plots;
- main and supplementary figure generation.

