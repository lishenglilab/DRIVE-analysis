#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import warnings
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

warnings.filterwarnings("ignore")

WORK_DIR = "path/to/your/project"
os.chdir(WORK_DIR)

EXTERNAL_TRUTH_FILE = "gdsc_ic50.csv"
CELL_MAP_FILE = "cellline2.csv"

MODELS_TO_ANALYZE = {
    "DeepTTA": "DeepTTC_predictions",
    "DIPK": "DIPK_predictions",
    "GPDRP_GAT": "GPDRP_predictions_GAT",
    "GPDRP_GCN": "GPDRP_predictions_GCN",
    "Precily": "Precily_predictions",
    "BANDRP": "bandrp_predictions_part_1",
    "paccmann": "paccmann_predictions",
    "GraphDRP_GATNet": "GraphDRP_predictions_GATNet",
    "GraphDRP_GAT_GCN": "GraphDRP_predictions_GAT_GCN",
}

PARSING_RULES = {
    "DeepTTA": {"drug": "DrugName", "cell": "COSMIC_ID", "pred": "Predicted_LN_IC50"},
    "DIPK": {"drug": "resolved_drug_name", "cell": "cell_line_name", "pred": "predicted_ic50"},
    "GPDRP_GAT": {"drug": "drug_name", "cell": "cell_line_name", "pred": "predicted_ic50_original"},
    "GPDRP_GCN": {"drug": "drug_name", "cell": "cell_line_name", "pred": "predicted_ic50_original"},
    "Precily": {"drug": "drug_name", "cell": "cell_line_name", "pred": "predicted_ic50"},
    "BANDRP": {"drug": "DrugName", "cell": "CellLineID", "pred": "PredictedValue"},
    "paccmann": {"drug": "drug_name", "cell": "cell_line_name", "pred": "predicted_value_denormalized"},
    "GraphDRP_GATNet": {"drug": "drug_name", "cell": "cell_line_name", "pred": "IC50_original"},
    "GraphDRP_GAT_GCN": {"drug": "drug_name", "cell": "cell_line_name", "pred": "IC50_original"},
}

def create_cell_map(filepath):
    cell_df = pd.read_csv(filepath, index_col=0)
    cell_map = {}
    for _, row in cell_df.iterrows():
        canonical_id = str(row["ID"]).strip()
        for col in ["ID", "cell.names", "cosmic.id", "cell.names_nerd", "cell.names.1"]:
            if col in cell_df.columns and pd.notna(row[col]):
                for alias in str(row[col]).split("|"):
                    cell_map[alias.strip().lower()] = canonical_id
    return cell_map

def load_and_sync_data():
    cell_map = create_cell_map(CELL_MAP_FILE)
    raw_truth = pd.read_csv(EXTERNAL_TRUTH_FILE)
    id_col = raw_truth.columns[0]
    truth_df = raw_truth.melt(id_vars=id_col, var_name="DrugName", value_name="TrueValue").dropna()
    truth_df["DrugName"] = truth_df["DrugName"].astype(str).str.lower().str.strip()
    truth_df["CellLineID"] = truth_df[id_col].astype(str).str.lower().str.strip().map(cell_map)
    truth_df = truth_df.dropna(subset=["CellLineID"])
    if truth_df["TrueValue"].mean() > 20:
        truth_df["TrueValue"] = np.log(truth_df["TrueValue"].clip(lower=1e-5))
    final_df = truth_df[["DrugName", "CellLineID", "TrueValue"]].copy()

    for nickname, file_base in MODELS_TO_ANALYZE.items():
        fname = f"{file_base}.csv"
        if not os.path.exists(fname):
            raise FileNotFoundError(f"Missing prediction file: {fname}")
        rule = PARSING_RULES[nickname]
        pred_df = pd.read_csv(fname, low_memory=False)
        sub = pred_df[[rule["drug"], rule["cell"], rule["pred"]]].copy()
        sub.columns = ["DrugName", "cell_identifier", nickname]
        sub["DrugName"] = sub["DrugName"].astype(str).str.lower().str.strip()
        sub["CellLineID"] = sub["cell_identifier"].astype(str).str.lower().str.strip().map(cell_map)
        sub = (
            sub.dropna(subset=["CellLineID", nickname])
            .groupby(["DrugName", "CellLineID"], as_index=False)[nickname]
            .mean()
        )
        final_df = pd.merge(final_df, sub, on=["DrugName", "CellLineID"], how="inner")
    return final_df

def plot_residual_correlation(df):
    model_names = list(MODELS_TO_ANALYZE.keys())
    error_df = pd.DataFrame()
    for model in model_names:
        if model in df.columns:
            error_df[model] = df[model] - df["TrueValue"]

    corr_matrix = error_df.corr(method="spearman")
    corr_matrix.to_csv("SupFig_S7_residual_correlation_matrix.csv")

    plt.rcParams["font.family"] = "serif"
    plt.rcParams["font.serif"] = "Times New Roman"

    plt.figure(figsize=(11, 9))
    mask = np.triu(np.ones_like(corr_matrix, dtype=bool), k=1)
    sns.heatmap(
        corr_matrix, mask=mask, annot=True, fmt=".2f", cmap="YlOrRd",
        vmin=0.2, vmax=1.0, square=True, linewidths=0.8, linecolor="white",
        annot_kws={"size": 10},
        cbar_kws={"shrink": 0.7, "label": "Spearman Correlation of Residuals"},
    )
    plt.title("Supplementary Fig. S10 | Residual Correlation among Base Models",
              fontsize=15, pad=20, fontweight="bold")
    plt.xticks(rotation=45, ha="right", fontsize=11)
    plt.yticks(rotation=0, fontsize=11)
    plt.xlabel("")
    plt.ylabel("")
    plt.tight_layout()
    plt.savefig("SupFig_S10_residual_correlation.pdf", bbox_inches="tight")
    plt.close()

if __name__ == "__main__":
    final_eval_df = load_and_sync_data()
    if final_eval_df.empty:
        raise ValueError("Merged dataset is empty. Please check DrugName and CellLineID matching.")
    plot_residual_correlation(final_eval_df)
