#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Purpose:
    Generate Fig. 3D directional SHAP importance from the final ensemble model.

Input files:
    1. best_model_RandomForest_rmse.pkl
    2. gdsc_ic50.csv
    3. cellline2.csv
    4. Prediction files listed in FILENAME_FINDER

Output files:
    1. Fig3D_directional_shap_importance.csv
    2. Fig3D_input_merge_summary.csv
    3. Fig3D_directional_shap_importance.pdf
"""

import os
import pickle
import warnings
from functools import reduce
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import pearsonr

try:
    import shap
    HAS_SHAP = True
except ImportError:
    shap = None
    HAS_SHAP = False

warnings.filterwarnings("ignore")

WORK_DIR = Path(__file__).resolve().parent
os.chdir(WORK_DIR)

MODEL_PKL_FILE = "best_model_RandomForest_rmse.pkl"
EXTERNAL_TRUTH_FILE = "gdsc_ic50.csv"
CELL_MAP_FILE = "cellline2.csv"

ELITE_MODEL_POOL = [
    "DeepTTA", "DIPK", "GPDRP_GAT", "GPDRP_GCN", "Precily",
    "BANDRP", "paccmann", "GraphDRP_GATNet", "GraphDRP_GAT_GCN"
]

FILENAME_FINDER = {
    "BANDRP": "bandrp_predictions_part_1",
    "DeepTTA": "DeepTTC_predictions",
    "DIPK": "DIPK_predictions",
    "GPDRP_GAT": "GPDRP_predictions_GAT",
    "GPDRP_GCN": "GPDRP_predictions_GCN",
    "GraphDRP_GATNet": "GraphDRP_predictions_GATNet",
    "GraphDRP_GAT_GCN": "GraphDRP_predictions_GAT_GCN",
    "paccmann": "paccmann_predictions",
    "Precily": "Precily_predictions",
}

PARSING_RULES = {
    "BANDRP": {"drug": "DrugName", "cell": "CellLineID", "pred": "PredictedValue"},
    "DeepTTA": {"drug": "DrugName", "cell": "COSMIC_ID", "pred": "Predicted_LN_IC50"},
    "DIPK": {"drug": "resolved_drug_name", "cell": "cell_line_name", "pred": "predicted_ic50"},
    "GPDRP_GAT": {"drug": "drug_name", "cell": "cell_line_name", "pred": "predicted_ic50_original"},
    "GPDRP_GCN": {"drug": "drug_name", "cell": "cell_line_name", "pred": "predicted_ic50_original"},
    "GraphDRP_GATNet": {"drug": "drug_name", "cell": "cell_line_name", "pred": "IC50_original"},
    "GraphDRP_GAT_GCN": {"drug": "drug_name", "cell": "cell_line_name", "pred": "IC50_original"},
    "paccmann": {"drug": "drug_name", "cell": "cell_line_name", "pred": "predicted_value_denormalized"},
    "Precily": {"drug": "drug_name", "cell": "cell_line_name", "pred": "predicted_ic50"},
}

def create_cell_map(filepath):
    cell_df = pd.read_csv(filepath, index_col=0)
    cell_map = {}
    for _, row in cell_df.iterrows():
        canonical_id = str(row["ID"]).strip()
        for col in ["ID", "cell.names", "cosmic.id", "cell.names_nerd", "cell.names.1"]:
            if col in cell_df.columns and pd.notna(row[col]):
                for alias in str(row[col]).split("|"):
                    alias_clean = alias.strip().lower()
                    if alias_clean:
                        cell_map[alias_clean] = canonical_id
    return cell_map

def load_truth(filepath, cell_map):
    raw = pd.read_csv(filepath)
    id_col = raw.columns[0]
    truth = raw.melt(id_vars=id_col, var_name="DrugName", value_name="TrueValue").dropna()
    truth.rename(columns={id_col: "CellLineID_original"}, inplace=True)
    truth["DrugName"] = truth["DrugName"].astype(str).str.lower().str.strip()
    truth["CellLineID"] = truth["CellLineID_original"].astype(str).str.lower().str.strip().map(cell_map)
    truth = truth.dropna(subset=["CellLineID"])
    if truth["TrueValue"].mean() > 20:
        truth["TrueValue"] = np.log(truth["TrueValue"].clip(lower=1e-5))
    return truth[["DrugName", "CellLineID", "TrueValue"]]

def load_predictions(cell_map):
    all_preds = []
    for model_name in ELITE_MODEL_POOL:
        filepath = f"{FILENAME_FINDER[model_name]}.csv"
        if not os.path.exists(filepath):
            raise FileNotFoundError(f"Missing prediction file: {filepath}")
        rule = PARSING_RULES[model_name]
        df = pd.read_csv(filepath, low_memory=False)
        sub = df[[rule["drug"], rule["cell"], rule["pred"]]].copy()
        sub.columns = ["DrugName", "cell_identifier", model_name]
        sub["DrugName"] = sub["DrugName"].astype(str).str.lower().str.strip()
        sub["CellLineID"] = sub["cell_identifier"].astype(str).str.lower().str.strip().map(cell_map)
        sub = (
            sub.dropna(subset=["CellLineID", model_name])
            .groupby(["DrugName", "CellLineID"], as_index=False)[model_name]
            .mean()
        )
        all_preds.append(sub[["DrugName", "CellLineID", model_name]])
    return reduce(lambda left, right: pd.merge(left, right, on=["DrugName", "CellLineID"], how="inner"), all_preds)

def main():
    cell_map = create_cell_map(CELL_MAP_FILE)
    truth_df = load_truth(EXTERNAL_TRUTH_FILE, cell_map)
    pred_df = load_predictions(cell_map)
    final_df = pd.merge(pred_df, truth_df, on=["DrugName", "CellLineID"], how="inner")

    if final_df.empty:
        raise ValueError("Merged dataset is empty. Please check DrugName and CellLineID matching.")

    with open(MODEL_PKL_FILE, "rb") as f:
        saved = pickle.load(f)

    model = saved["model"]
    features = [f.replace("-", "_") for f in saved["features"]]
    final_df.columns = [c.replace("-", "_") for c in final_df.columns]

    missing_features = [feature for feature in features if feature not in final_df.columns]
    if missing_features:
        raise ValueError(
            "Merged dataset is missing model features required by the saved meta-learner: "
            + ", ".join(missing_features)
        )

    X_explain = final_df[features].fillna(final_df[features].median())
    if HAS_SHAP:
        explainer = shap.Explainer(model, X_explain)
        shap_values = explainer(X_explain)
        importance_values = np.abs(shap_values.values).mean(axis=0)
        direction_reference = shap_values.values
        importance_label = "Mean_Absolute_SHAP_Value"
        title_text = "Directional SHAP Importance for RandomForest Meta-learner"
        source_label = "SHAP"
    else:
        importance_values = getattr(model, "feature_importances_", np.zeros(len(features)))
        direction_reference = np.column_stack([model.predict(X_explain)] * len(features))
        importance_label = "Fallback_Feature_Importance"
        title_text = "Directional Feature Importance for RandomForest Meta-learner"
        source_label = "RandomForest feature_importances_"

    signs = []
    for i, feature in enumerate(features):
        corr, _ = pearsonr(X_explain[feature].values.flatten(), direction_reference[:, i].flatten())
        signs.append("+" if np.isnan(corr) or corr >= 0 else "-")

    summary = pd.DataFrame({
        "Base_Model_Feature": features,
        importance_label: importance_values,
        "Direction_Sign": signs,
        "Importance_Source": source_label
    }).sort_values(importance_label, ascending=False)

    summary.to_csv("Fig3D_directional_shap_importance.csv", index=False)
    pd.DataFrame(
        [
            {"dataset": "truth", "rows": len(truth_df), "unique_drugs": truth_df["DrugName"].nunique(), "unique_cells": truth_df["CellLineID"].nunique()},
            {"dataset": "predictions", "rows": len(pred_df), "unique_drugs": pred_df["DrugName"].nunique(), "unique_cells": pred_df["CellLineID"].nunique()},
            {"dataset": "merged", "rows": len(final_df), "unique_drugs": final_df["DrugName"].nunique(), "unique_cells": final_df["CellLineID"].nunique()},
        ]
    ).to_csv("Fig3D_input_merge_summary.csv", index=False)

    sns.set_theme(style="whitegrid")
    plt.rcParams["font.family"] = "serif"
    plt.rcParams["font.serif"] = "Times New Roman"

    fig, ax = plt.subplots(figsize=(12, 8))
    sns.barplot(x=importance_label, y="Base_Model_Feature", data=summary, color="#d62728", ax=ax)

    for y_pos, (_, row) in enumerate(summary.iterrows()):
        value = row[importance_label]
        sign = row["Direction_Sign"]
        ax.text(value + 0.01, y_pos, f" {sign}{value:.3f}", va="center", fontsize=12, color="#d62728")

    ax.set_title(title_text, fontsize=16, pad=16)
    ax.set_xlabel(importance_label.replace("_", " "))
    ax.set_ylabel("Base model feature")
    ax.set_xlim(0, summary[importance_label].max() * 1.18)

    plt.tight_layout()
    plt.savefig("Fig3D_directional_shap_importance.pdf", bbox_inches="tight")
    plt.close()

if __name__ == "__main__":
    main()
