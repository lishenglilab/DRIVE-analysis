# DRIVE-analysis

This repository contains the downstream analysis scripts used to reproduce the statistical analyses, result summaries, and figures reported in the DRIVE study.

The source code of the DRIVE prediction framework itself is maintained separately at:

- https://github.com/lishenglilab/DRIVE

## Repository scope

`DRIVE-analysis` is intended for analyses performed after model training and prediction, including:

- benchmark result collection and summary;
- performance comparison across mixed, cell-blind, and drug-blind settings;
- ensemble performance analysis;
- residual correlation and marginal ensemble gain analysis;
- natural product prediction result summary;
- CRC candidate compound prioritization analysis;
- experimental result processing;
- generation of main and supplementary figures.

This repository does not contain the core implementation of the DRIVE prediction workflow. The model execution pipeline, prediction scripts, and workflow-related files are provided in the `DRIVE` repository.

## Suggested repository structure

```text
DRIVE-analysis/
├── README.md
├── data/
│   ├── processed/
│   └── example/
├── scripts/
│   ├── 01_benchmark_summary/
│   ├── 02_split_comparison/
│   ├── 03_ensemble_analysis/
│   ├── 04_residual_correlation/
│   ├── 05_natural_product_analysis/
│   ├── 06_crc_prioritization/
│   └── 07_experimental_analysis/
├── results/
│   ├── tables/
│   └── figures/
└── environment/
    └── requirements.txt
```

The exact folder names may differ from the manuscript figure numbering, but each script is named to indicate the corresponding analysis module.

## Data availability

Processed input tables required to run the analysis scripts are provided when redistribution is permitted.

Large files, third-party datasets, and raw pharmacogenomic resources are not directly redistributed in this repository when restricted by the original data providers. In these cases, the scripts expect preprocessed tables generated from the sources described in the manuscript.

Example input files are provided in `data/example/` to demonstrate the expected data format.

## Main analysis modules

### 1. Benchmark performance summary

Scripts in `scripts/01_benchmark_summary/` summarize model performance across evaluation settings.

Typical inputs include tables containing:

- model name;
- evaluation setting;
- cross-validation fold;
- RMSE;
- R²;
- PCC;
- SCC;
- NDCG;
- NWPC.

Typical outputs include summary tables and figures comparing model performance across mixed, cell-blind, and drug-blind settings.

### 2. Split setting comparison

Scripts in `scripts/02_split_comparison/` compare model behavior under different generalization scenarios:

- mixed setting;
- cell-blind setting;
- drug-blind setting.

These analyses were used to support the comparison of interpolation-like prediction, unseen-cell-line prediction, and unseen-drug prediction.

### 3. Ensemble analysis

Scripts in `scripts/03_ensemble_analysis/` analyze the performance of different ensemble strategies and model combinations.

Typical analyses include:

- comparison between single models and ensemble models;
- comparison between random forest meta-learner and simpler ensemble strategies;
- relationship between ensemble size and RMSE;
- selection of the final DRIVE ensemble configuration.

### 4. Residual correlation analysis

Scripts in `scripts/04_residual_correlation/` evaluate the correlation structure of residuals among base models and examine whether lower residual correlation is associated with larger marginal ensemble gain.

Typical outputs include residual correlation matrices, marginal gain summaries, and related statistical plots.

### 5. Natural product analysis

Scripts in `scripts/05_natural_product_analysis/` summarize DRIVE predictions for the natural product library.

Typical analyses include:

- predicted sensitivity distribution;
- drug class or source category summaries;
- candidate compound ranking;
- filtering of CRC-related candidate compounds.

### 6. CRC candidate prioritization

Scripts in `scripts/06_crc_prioritization/` document the downstream prioritization of candidate compounds for colorectal cancer validation.

The prioritization may include:

- DRIVE-predicted sensitivity;
- consistency across CRC cell lines;
- ranking stability;
- compound availability;
- structural representativeness;
- prior biological evidence when applicable.

### 7. Experimental result analysis

Scripts in `scripts/07_experimental_analysis/` process experimental validation results, including cell viability assays and in vivo tumor growth analysis.

Typical outputs include:

- dose-response curves;
- fitted IC50 or log-transformed IC50 values;
- tumor volume curves;
- endpoint tumor weight summaries;
- statistical comparison tables.

## Running the analysis

After preparing the required input tables, the analysis scripts can be run according to the corresponding module.

For example:

```bash
Rscript scripts/01_benchmark_summary/benchmark_summary.R
Rscript scripts/03_ensemble_analysis/ensemble_size_analysis.R
Rscript scripts/07_experimental_analysis/cck8_dose_response_fit.R
```

Python scripts, if used, can be run as:

```bash
python scripts/04_residual_correlation/residual_correlation.py
```

Please check the header of each script for required input files and output paths.

## Software environment

The analysis scripts were written mainly in R and Python.

A typical environment includes:

- R;
- Python;
- tidyverse or data.table for tabular processing;
- ggplot2 for visualization;
- pandas and numpy for numerical analysis;
- scipy or scikit-learn for statistical calculations when required.

Detailed package versions are provided in the `environment/` folder where applicable.

## Output files

Generated tables are saved under:

```text
results/tables/
```

Generated figures are saved under:

```text
results/figures/
```

The output file names are designed to correspond to the relevant analysis modules and manuscript figures.

## Relationship to the DRIVE framework repository

This repository is complementary to the DRIVE framework repository.

- `DRIVE`: source code and workflow for model prediction.
- `DRIVE-analysis`: downstream analysis scripts for summarizing prediction results, generating figures, and processing validation data.

Users who only want to run DRIVE predictions should refer to the `DRIVE` repository. Users who want to reproduce the manuscript analyses and figures should use this repository.

## Notes

Because several original datasets are subject to redistribution restrictions, this repository focuses on reproducible scripts, processed tables that can be shared, and example input files. When raw data are not included, users should obtain them from the original resources described in the manuscript.

## Contact

For questions related to this repository, please contact the corresponding authors listed in the manuscript.
