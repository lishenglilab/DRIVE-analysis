#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import warnings
from functools import reduce

import numpy as np
import pandas as pd
import shap
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import spearmanr
from sklearn.model_selection import KFold
from sklearn.metrics import mean_squared_error
from sklearn.ensemble import RandomForestRegressor
from tqdm import tqdm

warnings.filterwarnings("ignore")

WORK_DIR = "path/to/your/project"
os.chdir(WORK_DIR)

CELL_MAP_FILE = "cellline2.csv"
TRAINING_TRUTH_FILE = "ic50.csv"

ELITE_MODEL_POOL = [
    "DeepTTA", "DIPK", "GPDRP-GAT", "GPDRP-GCN", "Precily", "BANDRP",
    "paccmann", "GraphDRP-GATNet", "GraphDRP-GAT_GCN", "DeepAEG", "DeepCDR",
    "GADRP", "GPDRP-GIN", "GPDRP-GINTransformer", "GraphDRP-GCNNet",
    "GraphDRP-GINConvNet", "NERD", "DeepCCDS"
]
FEATURES = [f.replace("-", "_") for f in ELITE_MODEL_POOL]

FULL_MODEL_POOL = {
    "bandrp_predictions_part_1": "BANDRP",
    "DeepAEG_predictions": "DeepAEG",
    "DeepCDR_predictions": "DeepCDR",
    "DeepTTC_predictions": "DeepTTA",
    "DIPK_predictions": "DIPK",
    "GADRP_predictions": "GADRP",
    "GPDRP_predictions_GAT": "GPDRP-GAT",
    "GPDRP_predictions_GCN": "GPDRP-GCN",
    "GPDRP_predictions_GIN": "GPDRP-GIN",
    "GPDRP_predictions_GINTransformer": "GPDRP-GINTransformer",
    "GraphDRP_predictions_GATNet": "GraphDRP-GATNet",
    "GraphDRP_predictions_GAT_GCN": "GraphDRP-GAT_GCN",
    "GraphDRP_predictions_GCNNet": "GraphDRP-GCNNet",
    "GraphDRP_predictions_GINConvNet": "GraphDRP-GINConvNet",
    "NERD_predictions": "NERD",
    "paccmann_predictions": "paccmann",
    "Precily_predictions": "Precily",
    "DeepCCDS_predictions": "DeepCCDS",
}

PARSING_RULES = {
    "BANDRP": {"drug": "DrugName", "cell": "CellLineID", "pred": "PredictedValue", "format": "long"},
    "DeepAEG": {"drug": "Drug_ID", "cell": None, "pred": None, "format": "wide"},
    "DeepCDR": {"drug": "drug_name", "cell": "cell_line_id", "pred": "predicted_ln_IC50", "format": "long"},
    "DeepTTA": {"drug": "DrugName", "cell": "COSMIC_ID", "pred": "Predicted_LN_IC50", "format": "long"},
    "DIPK": {"drug": "resolved_drug_name", "cell": "cell_line_name", "pred": "predicted_ic50", "format": "long"},
    "GADRP": {"drug": "DrugName", "cell": "CellLineName", "pred": "Predicted_IC50", "format": "long"},
    "GPDRP_GAT": {"drug": "drug_name", "cell": "cell_line_name", "pred": "predicted_ic50_original", "format": "long"},
    "GPDRP_GCN": {"drug": "drug_name", "cell": "cell_line_name", "pred": "predicted_ic50_original", "format": "long"},
    "GPDRP_GIN": {"drug": "drug_name", "cell": "cell_line_name", "pred": "predicted_ic50_original", "format": "long"},
    "GPDRP_GINTransformer": {"drug": "drug_name", "cell": "cell_line_name", "pred": "predicted_ic50_original", "format": "long"},
    "GraphDRP_GATNet": {"drug": "drug_name", "cell": "cell_line_name", "pred": "IC50_original", "format": "long"},
    "GraphDRP_GAT_GCN": {"drug": "drug_name", "cell": "cell_line_name", "pred": "IC50_original", "format": "long"},
    "GraphDRP_GCNNet": {"drug": "drug_name", "cell": "cell_line_name", "pred": "IC50_original", "format": "long"},
    "GraphDRP_GINConvNet": {"drug": "drug_name", "cell": "cell_line_name", "pred": "IC50_original", "format": "long"},
    "NERD": {"drug": "drug_name", "cell": "cell_line_name", "pred": "predicted_ic50_unscaled", "format": "long"},
    "paccmann": {"drug": "drug_name", "cell": "cell_line_name", "pred": "predicted_value_denormalized", "format": "long"},
    "Precily": {"drug": "drug_name", "cell": "cell_line_name", "pred": "predicted_ic50", "format": "long"},
    "DeepCCDS": {"drug": "drug_name", "cell": "cell_name", "pred": "predicted_IC50", "format": "long"},
}

def create_cell_map(filepath):
    cell_df = pd.read_csv(filepath, index_col=0)
    cell_map = {}
    for _, row in cell_df.iterrows():
        canonical_id = str(row["ID"]).strip()
        for col in ["ID", "cell.names", "cosmic.id", "cell.names_nerd", "cell.names.1"]:
            if col in cell_df.columns and pd.notna(row[col]):
                for syn in str(row[col]).split("|"):
                    cell_map[syn.strip().lower()] = canonical_id
    return cell_map

def load_prediction_data(pred_dir, cell_map):
    all_predictions = []
    for filename_base, method_name in tqdm(FULL_MODEL_POOL.items(), desc="Loading predictions"):
        rule = PARSING_RULES.get(method_name.replace("-", "_"))
        filepath = os.path.join(pred_dir, f"{filename_base}.csv")
        if not rule or not os.path.exists(filepath):
            continue
        df = pd.read_csv(filepath, low_memory=False)
        if rule["format"] == "wide":
            long = pd.melt(df, id_vars=[rule["drug"]], var_name="cell_identifier", value_name=method_name)
            long.rename(columns={rule["drug"]: "DrugName"}, inplace=True)
        else:
            long = df[[rule["drug"], rule["cell"], rule["pred"]]].copy()
            long.columns = ["DrugName", "cell_identifier", method_name]
        long["CellLineID"] = long["cell_identifier"].astype(str).str.strip().str.lower().map(cell_map)
        long.dropna(subset=["CellLineID", method_name], inplace=True)
        all_predictions.append(long[["DrugName", "CellLineID", method_name]])
    if len(all_predictions) == 0:
        raise ValueError("No prediction files were loaded.")
    return reduce(lambda left, right: pd.merge(left, right, on=["DrugName", "CellLineID"], how="outer"), all_predictions)

def run_shap_cv_analysis():
    cell_map = create_cell_map(CELL_MAP_FILE)
    truth_wide = pd.read_csv(TRAINING_TRUTH_FILE)
    truth_wide.rename(columns={truth_wide.columns[0]: "CellLineID"}, inplace=True)
    truth_long = truth_wide.melt(id_vars=["CellLineID"], var_name="DrugName", value_name="TrueValue_lnIC50")
    preds = load_prediction_data("./", cell_map)
    df = pd.merge(preds, truth_long, on=["DrugName", "CellLineID"], how="inner")
    df.columns = [c.replace("-", "_") for c in df.columns]
    current_features = [f for f in FEATURES if f in df.columns]
    X = df[current_features].fillna(df[current_features].median())
    y = df["TrueValue_lnIC50"]

    kf = KFold(n_splits=5, shuffle=True, random_state=42)
    shap_records, perf_records = [], []

    for fold, (train_idx, test_idx) in enumerate(kf.split(X), start=1):
        X_train, X_test = X.iloc[train_idx], X.iloc[test_idx]
        y_train, y_test = y.iloc[train_idx], y.iloc[test_idx]
        model = RandomForestRegressor(
            n_estimators=50, max_depth=10, random_state=42, n_jobs=-1,
            min_samples_leaf=5, max_features="sqrt"
        )
        model.fit(X_train, y_train)
        explainer = shap.TreeExplainer(model)
        shap_values = explainer.shap_values(X_test)
        mas = np.abs(shap_values).mean(axis=0)
        for i, feature in enumerate(current_features):
            shap_records.append({"Fold": fold, "Model": feature, "SHAP_Importance": mas[i]})
            perf_records.append({
                "Fold": fold,
                "Model": feature,
                "Standalone_RMSE": np.sqrt(mean_squared_error(y_test, X_test[feature])),
            })
    return pd.DataFrame(shap_records), pd.DataFrame(perf_records)

def plot_shap_stability(shap_df, perf_df):
    sns.set_theme(style="white")
    plt.rcParams["font.family"] = "serif"
    plt.rcParams["font.serif"] = "Times New Roman"

    fig = plt.figure(figsize=(18, 12))
    gs = plt.GridSpec(2, 2, height_ratios=[1, 1])

    ax1 = fig.add_subplot(gs[0, :])
    order = shap_df.groupby("Model")["SHAP_Importance"].mean().sort_values(ascending=False).index
    sns.boxplot(x="Model", y="SHAP_Importance", data=shap_df, order=order, ax=ax1, palette="vlag", fliersize=0)
    sns.stripplot(x="Model", y="SHAP_Importance", data=shap_df, order=order, ax=ax1, color=".3", size=4, jitter=True)
    ax1.set_title("A. SHAP Importance Distribution across 5 Folds", fontsize=16, loc="left", pad=15)
    ax1.set_xticklabels(ax1.get_xticklabels(), rotation=45, ha="right")
    ax1.set_ylabel("Mean Absolute SHAP Value")

    ax2 = fig.add_subplot(gs[1, 0])
    shap_df = shap_df.copy()
    shap_df["Rank"] = shap_df.groupby("Fold")["SHAP_Importance"].rank(ascending=False)
    rank_pivot = shap_df.pivot(index="Model", columns="Fold", values="Rank").loc[order]
    sns.heatmap(rank_pivot, annot=True, cmap="YlGnBu_r", ax=ax2, cbar_kws={"label": "Rank (1=Best)"})
    ax2.set_title("B. Feature Rank Stability across Folds", fontsize=16, loc="left", pad=15)

    ax3 = fig.add_subplot(gs[1, 1])
    merged = pd.merge(shap_df, perf_df, on=["Fold", "Model"])
    sns.regplot(
        x="Standalone_RMSE", y="SHAP_Importance", data=merged, ax=ax3,
        scatter_kws={"alpha": 0.5, "color": "#d62728"}, line_kws={"color": "black"}
    )
    r, p = spearmanr(merged["Standalone_RMSE"], merged["SHAP_Importance"])
    ax3.text(0.05, 0.95, f"Spearman R: {r:.3f}\np-value: {p:.2e}",
             transform=ax3.transAxes, va="top",
             bbox=dict(boxstyle="round", facecolor="white", alpha=0.5))
    ax3.set_title("C. SHAP Contribution vs Standalone RMSE", fontsize=16, loc="left", pad=15)
    ax3.set_xlabel("Standalone Model RMSE (Lower is Better)")
    ax3.set_ylabel("SHAP Meta-Model Contribution")

    plt.tight_layout()
    plt.savefig("SupFig_S6_SHAP_stability.pdf", dpi=300, bbox_inches="tight")
    plt.close()

if __name__ == "__main__":
    shap_results, perf_results = run_shap_cv_analysis()
    shap_results.to_csv("SupFig_S6_SHAP_importance_by_fold.csv", index=False)
    perf_results.to_csv("SupFig_S6_standalone_rmse_by_fold.csv", index=False)
    plot_shap_stability(shap_results, perf_results)
