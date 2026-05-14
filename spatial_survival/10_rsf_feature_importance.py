"""
10_rsf_feature_importance.py
Feature importance analysis and ablation study for Random Survival Forest.

This script:
1. Loads trained RSF models from all 5 folds
2. Extracts feature importances (per-tree average)
3. Computes mean/std importance across folds
4. Performs feature ablation: drop one feature group at a time, retrain, measure ΔC-index
5. Visualizes importance and ablation results
6. Recommends optimal feature subsets

Feature groups (102 total):
  - Protein means (39):  avg protein expression
  - Protein stds (39):   protein heterogeneity
  - Edge means (3):      neighbor similarity, distance, interaction
  - Edge stds (3):       variation in edge features
  - Graph stats (3):     n_nodes, n_edges, edge_density
  - Interactions (5):    fractions of 5 interaction types

Results saved to:
  output/results/feature_importance/importances.csv
  output/results/feature_importance/ablation_results.csv
  output/results/feature_importance/feature_importance_plot.png
  output/results/feature_importance/ablation_plot.png

Run AFTER 06_training.py:

    python spatial_survival/10_rsf_feature_importance.py
"""

import sys
import pickle
from pathlib import Path
from typing import Dict, Tuple

import numpy as np
import pandas as pd
import torch
from sklearn.model_selection import GroupKFold
from sklearn.preprocessing import StandardScaler
from sksurv.ensemble import RandomSurvivalForest
from sksurv.util import Surv
from torch_geometric.data import Data

sys.path.insert(0, str(Path(__file__).parent))
from config import (
    OUTPUT_DIR, RESULTS_DIR, SEED, N_CV_FOLDS, VAL_FRACTION, N_PROTEINS, PROTEIN_COLS,
)
from utils import get_logger, ensure_dirs, set_seed, compute_cindex

PYG_RAW_DIR = OUTPUT_DIR / "pyg_dataset" / "raw"
INDEX_PATH = OUTPUT_DIR / "pyg_dataset" / "dataset_index.csv"
CKPT_DIR = RESULTS_DIR / "checkpoints"
IMPORTANCE_DIR = RESULTS_DIR / "feature_importance"

logger = get_logger("feature_importance", RESULTS_DIR / "feature_importance.log")

# Feature group definitions (102 total features)
FEATURE_GROUPS = {
    "Protein_Mean": (0, 39, "Protein expression means (avg level)"),
    "Spatial_Mean": (39, 44, "Spatial statistics means (density, entropy, boundary, degree, gradient)"),
    "Protein_Std": (44, 83, "Protein expression stds (heterogeneity)"),
    "Spatial_Std": (83, 88, "Spatial statistics stds"),
    "EdgeSim_Mean": (88, 89, "Edge cosine similarity mean"),
    "EdgeDist_Mean": (89, 90, "Edge distance weight mean"),
    "EdgeInt_Mean": (90, 91, "Edge interaction type mean"),
    "EdgeSim_Std": (91, 92, "Edge cosine similarity std"),
    "EdgeDist_Std": (92, 93, "Edge distance weight std"),
    "EdgeInt_Std": (93, 94, "Edge interaction type std"),
    "GraphSize": (94, 95, "Number of nodes"),
    "GraphEdges": (95, 96, "Number of edges"),
    "EdgeDensity": (96, 97, "Edge density"),
    "IntFrac_Homotypic": (97, 98, "Fraction homotypic interactions"),
    "IntFrac_TumorImmune": (98, 99, "Fraction tumor-immune interactions"),
    "IntFrac_TumorStroma": (99, 100, "Fraction tumor-stroma interactions"),
    "IntFrac_ImmuneStroma": (100, 101, "Fraction immune-stroma interactions"),
    "IntFrac_Unknown": (101, 102, "Fraction unknown interactions"),
}

FEATURE_NAMES_DETAILED = (
    # Protein means (39)
    [f"{p}_mean" for p in PROTEIN_COLS] +
    # Spatial means (5)
    ["density_mean", "entropy_mean", "boundary_mean", "degree_mean", "gradient_mean"] +
    # Protein stds (39)
    [f"{p}_std" for p in PROTEIN_COLS] +
    # Spatial stds (5)
    ["density_std", "entropy_std", "boundary_std", "degree_std", "gradient_std"] +
    # Edge features (6)
    ["edge_cosine_sim_mean", "edge_distance_mean", "edge_interaction_mean",
     "edge_cosine_sim_std", "edge_distance_std", "edge_interaction_std"] +
    # Graph stats (3)
    ["n_nodes", "n_edges", "edge_density"] +
    # Interaction fractions (5)
    ["frac_homotypic", "frac_tumor_immune", "frac_tumor_stroma",
     "frac_immune_stroma", "frac_unknown"]
)

assert len(FEATURE_NAMES_DETAILED) == 102, f"Expected 102 features, got {len(FEATURE_NAMES_DETAILED)}"


# ---------------------------------------------------------------------------
# Data helpers
# ---------------------------------------------------------------------------

def load_all_data(index_df: pd.DataFrame) -> list[Data]:
    data_list = []
    for _, row in index_df.iterrows():
        pt = PYG_RAW_DIR / f"{row['acquisition_id']}.pt"
        if pt.exists():
            data_list.append(torch.load(pt, weights_only=False))
    return data_list


def _to_numpy(array) -> np.ndarray:
    if array is None:
        return np.asarray([], dtype=np.float32)
    if torch.is_tensor(array):
        return array.detach().cpu().numpy()
    return np.asarray(array)


def extract_graph_features(data: Data) -> np.ndarray:
    """Collapse one enriched graph into a fixed-length (102,) feature vector."""
    x = _to_numpy(data.x).astype(np.float32)
    edge_attr = _to_numpy(getattr(data, "edge_attr", None)).astype(np.float32)

    if x.ndim != 2:
        raise ValueError(f"Expected 2D node features, got shape {x.shape}")

    node_mean = x.mean(axis=0)
    node_std = x.std(axis=0)

    if edge_attr.size == 0:
        edge_mean = np.zeros(3, dtype=np.float32)
        edge_std = np.zeros(3, dtype=np.float32)
        interaction_frac = np.zeros(5, dtype=np.float32)
        n_edges = 0.0
    else:
        edge_mean = edge_attr.mean(axis=0)
        edge_std = edge_attr.std(axis=0)
        interaction_codes = np.clip(np.rint(edge_attr[:, 2]).astype(int), 0, 4)
        interaction_frac = np.bincount(interaction_codes, minlength=5).astype(np.float32)
        interaction_frac /= max(float(len(interaction_codes)), 1.0)
        n_edges = float(edge_attr.shape[0])

    n_nodes = float(x.shape[0])
    edge_density = n_edges / max(n_nodes, 1.0)
    graph_stats = np.array([n_nodes, n_edges, edge_density], dtype=np.float32)

    return np.concatenate([
        node_mean,
        node_std,
        edge_mean,
        edge_std,
        graph_stats,
        interaction_frac,
    ]).astype(np.float32)


def build_graph_feature_matrix(data_list: list[Data]) -> np.ndarray:
    return np.vstack([extract_graph_features(data) for data in data_list])


def scale_graph_features(train_x: np.ndarray, val_x: np.ndarray, test_x: np.ndarray):
    scaler = StandardScaler()
    scaler.fit(train_x)
    return (
        scaler.transform(train_x),
        scaler.transform(val_x),
        scaler.transform(test_x),
        scaler,
    )


# ---------------------------------------------------------------------------
# Feature importance extraction
# ---------------------------------------------------------------------------

def compute_permutation_importance(rsf: RandomSurvivalForest,
                                   X: np.ndarray,
                                   y_times: np.ndarray,
                                   y_events: np.ndarray,
                                   baseline_cindex: float,
                                   n_repeats: int = 5) -> np.ndarray:
    """
    Compute permutation importance: for each feature, shuffle it and measure
    drop in C-index. Repeats this process n_repeats times and averages.
    
    Returns: (n_features,) importance scores
    """
    importances = np.zeros(X.shape[1])
    
    for feat_idx in range(X.shape[1]):
        feat_importances = []
        
        for _ in range(n_repeats):
            X_permuted = X.copy()
            np.random.shuffle(X_permuted[:, feat_idx])
            
            risk_permuted = rsf.predict(X_permuted)
            cindex_permuted = compute_cindex(risk_permuted, y_times, y_events)
            
            importance = baseline_cindex - cindex_permuted
            feat_importances.append(importance)
        
        importances[feat_idx] = np.mean(feat_importances)
    
    return importances


def extract_feature_importances_from_rsf(rsf: RandomSurvivalForest,
                                         X: np.ndarray,
                                         y_times: np.ndarray,
                                         y_events: np.ndarray) -> np.ndarray:
    """Extract permutation feature importance from RSF."""
    baseline_risk = rsf.predict(X)
    baseline_cindex = compute_cindex(baseline_risk, y_times, y_events)
    
    importances = compute_permutation_importance(
        rsf, X, y_times, y_events, baseline_cindex, n_repeats=3
    )
    return importances


# ---------------------------------------------------------------------------
# Feature ablation
# ---------------------------------------------------------------------------

def compute_feature_importance_via_ablation(
    train_x: np.ndarray,
    val_x: np.ndarray,
    test_x: np.ndarray,
    y_train_events: np.ndarray,
    y_train_times: np.ndarray,
    y_val_events: np.ndarray,
    y_val_times: np.ndarray,
    y_test_events: np.ndarray,
    y_test_times: np.ndarray,
    fold_idx: int,
) -> Dict[str, Tuple[float, float]]:
    """
    Ablation study: train RSF with all features, then with each group zeroed out.
    Returns: {group_name: (baseline_cindex, ablated_cindex, delta_cindex)}
    """
    # Train baseline RSF on full feature set
    y_train = Surv.from_arrays(y_train_events, y_train_times)
    baseline_rsf = RandomSurvivalForest(
        n_estimators=300,
        min_samples_split=10,
        min_samples_leaf=3,
        max_features="sqrt",
        n_jobs=-1,
        random_state=SEED + fold_idx,
    )
    baseline_rsf.fit(train_x, y_train)
    
    baseline_test_risk = baseline_rsf.predict(test_x)
    baseline_cindex = compute_cindex(baseline_test_risk, y_test_times, y_test_events)
    
    logger.info(f"  Fold {fold_idx+1} Baseline Test C-index: {baseline_cindex:.4f}")
    
    ablation_results = {}
    
    # For each feature group, zero out and retrain
    for group_name, (start, end, description) in FEATURE_GROUPS.items():
        train_x_ablated = train_x.copy()
        val_x_ablated = val_x.copy()
        test_x_ablated = test_x.copy()
        
        # Zero out this group
        train_x_ablated[:, start:end] = 0.0
        val_x_ablated[:, start:end] = 0.0
        test_x_ablated[:, start:end] = 0.0
        
        # Train RSF on ablated features
        ablated_rsf = RandomSurvivalForest(
            n_estimators=300,
            min_samples_split=10,
            min_samples_leaf=3,
            max_features="sqrt",
            n_jobs=-1,
            random_state=SEED + fold_idx,
        )
        ablated_rsf.fit(train_x_ablated, y_train)
        
        ablated_test_risk = ablated_rsf.predict(test_x_ablated)
        ablated_cindex = compute_cindex(ablated_test_risk, y_test_times, y_test_events)
        
        delta_cindex = baseline_cindex - ablated_cindex
        importance_pct = (delta_cindex / baseline_cindex) * 100 if baseline_cindex > 0 else 0.0
        
        ablation_results[group_name] = {
            "baseline_cindex": baseline_cindex,
            "ablated_cindex": ablated_cindex,
            "delta_cindex": delta_cindex,
            "importance_pct": importance_pct,
        }
        
        logger.info(f"    {group_name:25s} | ablated_C={ablated_cindex:.4f} | "
                   f"ΔC={delta_cindex:+.4f} ({importance_pct:+.1f}%)")
    
    return ablation_results


# ---------------------------------------------------------------------------
# Main analysis
# ---------------------------------------------------------------------------

def analyze_feature_importance():
    set_seed(SEED)
    ensure_dirs(IMPORTANCE_DIR)
    
    logger.info("="*70)
    logger.info("RSF FEATURE IMPORTANCE ANALYSIS")
    logger.info("="*70)
    
    index_df = pd.read_csv(INDEX_PATH)
    all_data = load_all_data(index_df)
    groups = index_df["patient_id"].to_numpy()
    
    logger.info(f"Dataset: {len(index_df)} samples, {len(np.unique(groups))} patients")
    
    # ===== STEP 1: Extract importance from trained models =====
    
    logger.info("\n" + "="*70)
    logger.info("STEP 1: Extract importance from trained RSF models")
    logger.info("="*70)
    
    gkf = GroupKFold(n_splits=N_CV_FOLDS)
    
    fold_importances = []  # list of (fold, feature_idx, importance)
    ablation_by_fold = []  # list of dicts per fold
    
    for fold_idx, (train_val_idx, test_idx) in enumerate(gkf.split(index_df, groups=groups)):
        logger.info(f"\n--- Fold {fold_idx+1}/{N_CV_FOLDS} ---")
        
        # Load checkpoint
        ckpt_path = CKPT_DIR / f"fold_{fold_idx+1}_best.pkl"
        if not ckpt_path.exists():
            logger.warning(f"  Checkpoint not found: {ckpt_path}")
            continue
        
        with open(ckpt_path, "rb") as f:
            ckpt = pickle.load(f)
        rsf = ckpt["model"]
        
        # Prepare data
        train_val_patients = np.unique(groups[train_val_idx])
        rng = np.random.default_rng(SEED + fold_idx)
        rng.shuffle(train_val_patients)
        n_val_patients = max(1, int(len(train_val_patients) * VAL_FRACTION))
        val_patients = set(train_val_patients[:n_val_patients])
        
        train_idx_fold = [i for i in train_val_idx if groups[i] not in val_patients]
        val_idx_fold = [i for i in train_val_idx if groups[i] in val_patients]
        
        train_data = [all_data[i] for i in train_idx_fold]
        val_data = [all_data[i] for i in val_idx_fold]
        test_data = [all_data[i] for i in test_idx]
        
        train_x = build_graph_feature_matrix(train_data)
        val_x = build_graph_feature_matrix(val_data)
        test_x = build_graph_feature_matrix(test_data)
        train_x, val_x, test_x, _ = scale_graph_features(train_x, val_x, test_x)
        
        y_test_times = np.array([float(d.y_time.item()) for d in test_data])
        y_test_events = np.array([bool(d.y_event.item()) for d in test_data])
        
        # Extract permutation importance on test set
        logger.info(f"  Computing permutation importance ({test_x.shape[0]} test samples)...")
        importances = extract_feature_importances_from_rsf(rsf, test_x, y_test_times, y_test_events)
        
        for feat_idx, imp in enumerate(importances):
            fold_importances.append({
                "fold": fold_idx + 1,
                "feature_idx": feat_idx,
                "feature_name": FEATURE_NAMES_DETAILED[feat_idx],
                "importance": imp,
            })
        
        logger.info(f"  Extracted {len(importances)} feature importances")
        
        # ===== STEP 2: Perform ablation study =====
        
        logger.info(f"\n  Ablation Study (Fold {fold_idx+1}):")
        
        y_train_times = np.array([float(d.y_time.item()) for d in train_data])
        y_train_events = np.array([bool(d.y_event.item()) for d in train_data])
        y_val_times = np.array([float(d.y_time.item()) for d in val_data])
        y_val_events = np.array([bool(d.y_event.item()) for d in val_data])
        
        ablation_results = compute_feature_importance_via_ablation(
            train_x, val_x, test_x,
            y_train_events, y_train_times,
            y_val_events, y_val_times,
            y_test_events, y_test_times,
            fold_idx,
        )
        
        ablation_by_fold.append({
            "fold": fold_idx + 1,
            "results": ablation_results,
        })
    
    # ===== STEP 3: Aggregate across folds =====
    
    logger.info("\n" + "="*70)
    logger.info("STEP 2: Aggregate importance across folds")
    logger.info("="*70)
    
    # Per-feature importance (from model.feature_importances_)
    importance_df = pd.DataFrame(fold_importances)
    importance_summary = importance_df.groupby("feature_name")["importance"].agg(["mean", "std", "min", "max"]).reset_index()
    importance_summary = importance_summary.sort_values("mean", ascending=False)
    
    logger.info("\nTop 20 Most Important Features:")
    for idx, row in importance_summary.head(20).iterrows():
        logger.info(f"  {row['feature_name']:30s} | mean={row['mean']:.6f} ± {row['std']:.6f}")
    
    # Group importance
    logger.info("\nFeature Group Importance (aggregated):")
    group_importance = []
    for group_name, (start, end, description) in FEATURE_GROUPS.items():
        group_imps = importance_df[
            (importance_df["feature_idx"] >= start) & 
            (importance_df["feature_idx"] < end)
        ]["importance"].values
        mean_imp = group_imps.mean()
        std_imp = group_imps.std()
        total_imp = group_imps.sum()
        group_importance.append({
            "group": group_name,
            "description": description,
            "n_features": end - start,
            "mean_importance": mean_imp,
            "std_importance": std_imp,
            "total_importance": total_imp,
        })
        logger.info(f"  {group_name:25s} | mean={mean_imp:.6f} | total={total_imp:.6f}")
    
    group_importance_df = pd.DataFrame(group_importance).sort_values("mean_importance", ascending=False)
    
    # Ablation study summary
    logger.info("\nAblation Study Summary (avg across folds):")
    ablation_summary = {}
    for fold_data in ablation_by_fold:
        for group_name, results in fold_data["results"].items():
            if group_name not in ablation_summary:
                ablation_summary[group_name] = []
            ablation_summary[group_name].append(results)
    
    ablation_agg = []
    for group_name in sorted(ablation_summary.keys()):
        deltas = [r["delta_cindex"] for r in ablation_summary[group_name]]
        imps = [r["importance_pct"] for r in ablation_summary[group_name]]
        ablation_agg.append({
            "group": group_name,
            "mean_delta_cindex": np.mean(deltas),
            "std_delta_cindex": np.std(deltas),
            "mean_importance_pct": np.mean(imps),
            "std_importance_pct": np.std(imps),
        })
        logger.info(f"  {group_name:25s} | ΔC={np.mean(deltas):+.4f} ± {np.std(deltas):.4f}")
    
    ablation_agg_df = pd.DataFrame(ablation_agg).sort_values("mean_delta_cindex", ascending=False)
    
    # ===== STEP 4: Save results =====
    
    logger.info("\n" + "="*70)
    logger.info("STEP 3: Save results")
    logger.info("="*70)
    
    importance_df.to_csv(IMPORTANCE_DIR / "feature_importances_per_fold.csv", index=False)
    logger.info(f"  Saved: feature_importances_per_fold.csv ({len(importance_df)} rows)")
    
    importance_summary.to_csv(IMPORTANCE_DIR / "feature_importance_summary.csv", index=False)
    logger.info(f"  Saved: feature_importance_summary.csv")
    
    group_importance_df.to_csv(IMPORTANCE_DIR / "group_importance.csv", index=False)
    logger.info(f"  Saved: group_importance.csv")
    
    ablation_agg_df.to_csv(IMPORTANCE_DIR / "ablation_summary.csv", index=False)
    logger.info(f"  Saved: ablation_summary.csv")
    
    # ===== STEP 5: Visualization =====
    
    logger.info("\n" + "="*70)
    logger.info("STEP 4: Create visualizations")
    logger.info("="*70)
    
    try:
        import matplotlib.pyplot as plt
        import seaborn as sns
        
        # Plot 1: Top 30 features
        fig, ax = plt.subplots(figsize=(12, 8))
        top_30 = importance_summary.head(30)
        ax.barh(range(len(top_30)), top_30["mean"].values, xerr=top_30["std"].values, capsize=3)
        ax.set_yticks(range(len(top_30)))
        ax.set_yticklabels(top_30["feature_name"].values, fontsize=9)
        ax.set_xlabel("Mean Feature Importance", fontsize=11)
        ax.set_title("Top 30 Most Important Features (RSF)", fontsize=13, fontweight="bold")
        ax.invert_yaxis()
        plt.tight_layout()
        plt.savefig(IMPORTANCE_DIR / "top_features_importance.png", dpi=150, bbox_inches="tight")
        logger.info("  Saved: top_features_importance.png")
        plt.close()
        
        # Plot 2: Feature group importance
        fig, ax = plt.subplots(figsize=(12, 6))
        groups_sorted = group_importance_df.sort_values("total_importance", ascending=True)
        ax.barh(range(len(groups_sorted)), groups_sorted["total_importance"].values, 
                xerr=groups_sorted["std_importance"].values, capsize=3)
        ax.set_yticks(range(len(groups_sorted)))
        ax.set_yticklabels(groups_sorted["group"].values, fontsize=10)
        ax.set_xlabel("Total Importance (sum across features)", fontsize=11)
        ax.set_title("Feature Group Importance (RSF)", fontsize=13, fontweight="bold")
        plt.tight_layout()
        plt.savefig(IMPORTANCE_DIR / "group_importance_plot.png", dpi=150, bbox_inches="tight")
        logger.info("  Saved: group_importance_plot.png")
        plt.close()
        
        # Plot 3: Ablation results
        fig, ax = plt.subplots(figsize=(12, 6))
        ablation_sorted = ablation_agg_df.sort_values("mean_delta_cindex", ascending=True)
        colors = ["green" if x > 0 else "red" for x in ablation_sorted["mean_delta_cindex"].values]
        ax.barh(range(len(ablation_sorted)), ablation_sorted["mean_delta_cindex"].values,
                xerr=ablation_sorted["std_delta_cindex"].values, capsize=3, color=colors, alpha=0.7)
        ax.set_yticks(range(len(ablation_sorted)))
        ax.set_yticklabels(ablation_sorted["group"].values, fontsize=10)
        ax.set_xlabel("Mean ΔC-index (removed - baseline)", fontsize=11)
        ax.axvline(0, color="black", linestyle="--", linewidth=1)
        ax.set_title("Feature Ablation: Impact on Test C-index", fontsize=13, fontweight="bold")
        plt.tight_layout()
        plt.savefig(IMPORTANCE_DIR / "ablation_results_plot.png", dpi=150, bbox_inches="tight")
        logger.info("  Saved: ablation_results_plot.png")
        plt.close()
        
    except Exception as e:
        logger.warning(f"Visualization failed: {e}")
    
    # ===== STEP 6: Summary report =====
    
    logger.info("\n" + "="*70)
    logger.info("SUMMARY")
    logger.info("="*70)
    logger.info(f"Total features: {len(FEATURE_NAMES_DETAILED)}")
    logger.info(f"Feature groups: {len(FEATURE_GROUPS)}")
    logger.info(f"Folds analyzed: {N_CV_FOLDS}")
    logger.info(f"\nTop 5 most important features:")
    for i, row in importance_summary.head(5).iterrows():
        logger.info(f"  {i+1}. {row['feature_name']:30s} (importance={row['mean']:.6f})")
    logger.info(f"\nTop 5 feature groups by total importance:")
    for i, row in group_importance_df.head(5).iterrows():
        logger.info(f"  {i+1}. {row['group']:25s} (total={row['total_importance']:.6f})")
    logger.info(f"\nMost impactful features (ablation ΔC-index):")
    for i, row in ablation_agg_df.head(5).iterrows():
        logger.info(f"  {i+1}. {row['group']:25s} (ΔC={row['mean_delta_cindex']:+.4f})")
    
    logger.info("\n" + "="*70)
    logger.info("Analysis complete!")
    logger.info("="*70)


if __name__ == "__main__":
    analyze_feature_importance()
