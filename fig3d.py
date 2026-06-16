#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import pickle
import warnings
from functools import reduce

import numpy as np
import pandas as pd
import shap
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import pearsonr

warnings.filterwarnings("ignore")

WORK_DIR = "path/to/your/project"
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

    with open(MODEL_PKL_FILE, "rb") as f:
        saved = pickle.load(f)

    model = saved["model"]
    features = [f.replace("-", "_") for f in saved["features"]]
    final_df.columns = [c.replace("-", "_") for c in final_df.columns]

    X_explain = final_df[features].fillna(final_df[features].median())
    explainer = shap.Explainer(model, X_explain)
    shap_values = explainer(X_explain)

    mean_abs_shap = np.abs(shap_values.values).mean(axis=0)
    signs = []
    for i, feature in enumerate(features):
        corr, _ = pearsonr(X_explain[feature].values.flatten(), shap_values.values[:, i].flatten())
        signs.append("+" if corr > 0 else "-")

    summary = pd.DataFrame({
        "Base_Model_Feature": features,
        "Mean_Absolute_SHAP_Value": mean_abs_shap,
        "Direction_Sign": signs
    }).sort_values("Mean_Absolute_SHAP_Value", ascending=False)

    summary.to_csv("Fig3D_directional_shap_importance.csv", index=False)

    sns.set_theme(style="whitegrid")
    plt.rcParams["font.family"] = "serif"
    plt.rcParams["font.serif"] = "Times New Roman"

    fig, ax = plt.subplots(figsize=(12, 8))
    sns.barplot(x="Mean_Absolute_SHAP_Value", y="Base_Model_Feature", data=summary, color="#d62728", ax=ax)

    for y_pos, (_, row) in enumerate(summary.iterrows()):
        value = row["Mean_Absolute_SHAP_Value"]
        sign = row["Direction_Sign"]
        ax.text(value + 0.01, y_pos, f" {sign}{value:.3f}", va="center", fontsize=12, color="#d62728")

    ax.set_title("Directional SHAP Importance for RandomForest Meta-learner", fontsize=16, pad=16)
    ax.set_xlabel("mean(|SHAP value|)")
    ax.set_ylabel("Base model feature")
    ax.set_xlim(0, summary["Mean_Absolute_SHAP_Value"].max() * 1.18)

    plt.tight_layout()
    plt.savefig("Fig3D_directional_shap_importance.pdf", bbox_inches="tight")
    plt.close()

if __name__ == "__main__":
    main()
