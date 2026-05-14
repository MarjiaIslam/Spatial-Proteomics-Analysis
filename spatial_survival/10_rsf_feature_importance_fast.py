"""
10_rsf_feature_importance_fast.py
Fast feature importance analysis via group-level ablation (no per-feature permutation).

This script:
1. For each feature group: zero out and retrain RSF
2. Measure ΔC-index = baseline - ablated
3. Compute mean/std across folds
4. Visualize results

Much faster than per-feature permutation importance.

Results saved to:
  output/results/feature_importance/ablation_results.csv
  output/results/feature_importance/ablation_summary.csv
  output/results/feature_importance/ablation_plot.png

Run AFTER 06_training.py:

    python spatial_survival/10_rsf_feature_importance_fast.py
"""

import sys
import pickle
from pathlib import Path

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

logger = get_logger("feature_importance_fast", RESULTS_DIR / "feature_importance_fast.log")

# Feature group definitions (102 total)
FEATURE_GROUPS = {
    "Protein_Mean": (0, 39, "Protein expression means (avg level)"),
    "Spatial_Mean": (39, 44, "Spatial statistics means"),
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
# Feature ablation
# ---------------------------------------------------------------------------

def compute_group_ablation_importance(
    train_x: np.ndarray,
    val_x: np.ndarray,
    test_x: np.ndarray,
    y_train_events: np.ndarray,
    y_train_times: np.ndarray,
    y_test_events: np.ndarray,
    y_test_times: np.ndarray,
    fold_idx: int,
) -> dict:
    """
    Train baseline RSF, then ablate each group and measure ΔC-index.
    
    Returns: {group_name: {baseline_C, ablated_C, delta_C, importance_pct}}
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
    
    logger.info(f"  Baseline test C-index: {baseline_cindex:.4f}")
    
    results = {}
    
    # Ablate each group
    for group_name, (start, end, description) in FEATURE_GROUPS.items():
        train_x_ablated = train_x.copy()
        test_x_ablated = test_x.copy()
        
        # Zero out this group
        train_x_ablated[:, start:end] = 0.0
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
        importance_pct = (delta_cindex / baseline_cindex * 100) if baseline_cindex > 0 else 0.0
        
        results[group_name] = {
            "baseline_cindex": baseline_cindex,
            "ablated_cindex": ablated_cindex,
            "delta_cindex": delta_cindex,
            "importance_pct": importance_pct,
        }
        
        logger.info(f"    {group_name:25s} | ablated_C={ablated_cindex:.4f} | "
                   f"ΔC={delta_cindex:+.4f} ({importance_pct:+.1f}%)")
    
    return results


# ---------------------------------------------------------------------------
# Main analysis
# ---------------------------------------------------------------------------

def analyze_feature_importance():
    set_seed(SEED)
    ensure_dirs(IMPORTANCE_DIR)
    
    logger.info("="*70)
    logger.info("RSF FEATURE IMPORTANCE ANALYSIS (GROUP-LEVEL ABLATION)")
    logger.info("="*70)
    
    index_df = pd.read_csv(INDEX_PATH)
    all_data = load_all_data(index_df)
    groups = index_df["patient_id"].to_numpy()
    
    logger.info(f"Dataset: {len(index_df)} samples, {len(np.unique(groups))} patients\n")
    
    gkf = GroupKFold(n_splits=N_CV_FOLDS)
    ablation_by_fold = []
    
    for fold_idx, (train_val_idx, test_idx) in enumerate(gkf.split(index_df, groups=groups)):
        logger.info(f"\n{'='*70}")
        logger.info(f"Fold {fold_idx+1}/{N_CV_FOLDS}")
        logger.info(f"{'='*70}")
        
        # Prepare data
        train_val_patients = np.unique(groups[train_val_idx])
        rng = np.random.default_rng(SEED + fold_idx)
        rng.shuffle(train_val_patients)
        n_val_patients = max(1, int(len(train_val_patients) * VAL_FRACTION))
        val_patients = set(train_val_patients[:n_val_patients])
        
        train_idx_fold = [i for i in train_val_idx if groups[i] not in val_patients]
        test_idx_fold = [i for i in test_idx]
        
        train_data = [all_data[i] for i in train_idx_fold]
        test_data = [all_data[i] for i in test_idx_fold]
        
        train_x = build_graph_feature_matrix(train_data)
        test_x = build_graph_feature_matrix(test_data)
        
        # Use dummy val_x for scaling
        val_x = test_x[:5] if len(test_x) > 5 else test_x
        train_x, val_x, test_x, _ = scale_graph_features(train_x, val_x, test_x)
        
        y_train_times = np.array([float(d.y_time.item()) for d in train_data])
        y_train_events = np.array([bool(d.y_event.item()) for d in train_data])
        y_test_times = np.array([float(d.y_time.item()) for d in test_data])
        y_test_events = np.array([bool(d.y_event.item()) for d in test_data])
        
        logger.info(f"Train: {len(train_data)} | Test: {len(test_data)}\n")
        logger.info("Group-level ablation study:")
        
        # Ablation study
        ablation_results = compute_group_ablation_importance(
            train_x, val_x, test_x,
            y_train_events, y_train_times,
            y_test_events, y_test_times,
            fold_idx,
        )
        
        ablation_by_fold.append({
            "fold": fold_idx + 1,
            "results": ablation_results,
        })
    
    # ===== Aggregate across folds =====
    
    logger.info("\n" + "="*70)
    logger.info("AGGREGATED RESULTS (across 5 folds)")
    logger.info("="*70)
    
    ablation_agg = []
    for group_name in sorted(FEATURE_GROUPS.keys()):
        deltas = [fold_data["results"][group_name]["delta_cindex"] for fold_data in ablation_by_fold]
        imps = [fold_data["results"][group_name]["importance_pct"] for fold_data in ablation_by_fold]
        
        ablation_agg.append({
            "group": group_name,
            "n_folds": len(deltas),
            "mean_delta_cindex": np.mean(deltas),
            "std_delta_cindex": np.std(deltas),
            "min_delta": np.min(deltas),
            "max_delta": np.max(deltas),
            "mean_importance_pct": np.mean(imps),
            "std_importance_pct": np.std(imps),
        })
    
    ablation_agg_df = pd.DataFrame(ablation_agg).sort_values("mean_delta_cindex", ascending=False)
    
    logger.info("\nFeature Group Importance Ranking:")
    logger.info("-" * 100)
    for idx, row in ablation_agg_df.iterrows():
        logger.info(f"  {row['group']:25s} | ΔC={row['mean_delta_cindex']:+.4f} ± {row['std_delta_cindex']:.4f} | "
                   f"Impact={row['mean_importance_pct']:+.1f}% | Range=[{row['min_delta']:+.3f}, {row['max_delta']:+.3f}]")
    
    # ===== Save results =====
    
    logger.info("\n" + "="*70)
    logger.info("SAVING RESULTS")
    logger.info("="*70)
    
    ablation_agg_df.to_csv(IMPORTANCE_DIR / "ablation_summary.csv", index=False)
    logger.info(f"  Saved: ablation_summary.csv")
    
    # Save per-fold results
    fold_results_list = []
    for fold_data in ablation_by_fold:
        for group_name, res in fold_data["results"].items():
            fold_results_list.append({
                "fold": fold_data["fold"],
                "group": group_name,
                "baseline_cindex": res["baseline_cindex"],
                "ablated_cindex": res["ablated_cindex"],
                "delta_cindex": res["delta_cindex"],
                "importance_pct": res["importance_pct"],
            })
    fold_results_df = pd.DataFrame(fold_results_list)
    fold_results_df.to_csv(IMPORTANCE_DIR / "ablation_results.csv", index=False)
    logger.info(f"  Saved: ablation_results.csv")
    
    # ===== Visualization =====
    
    logger.info("\n" + "="*70)
    logger.info("CREATING VISUALIZATIONS")
    logger.info("="*70)
    
    try:
        import matplotlib.pyplot as plt
        
        # Plot: Ablation results sorted by importance
        fig, ax = plt.subplots(figsize=(12, 8))
        df_sorted = ablation_agg_df.sort_values("mean_delta_cindex", ascending=True)
        colors = ["green" if x > 0 else "red" for x in df_sorted["mean_delta_cindex"].values]
        
        y_pos = np.arange(len(df_sorted))
        ax.barh(y_pos, df_sorted["mean_delta_cindex"].values,
               xerr=df_sorted["std_delta_cindex"].values, capsize=3, color=colors, alpha=0.7)
        ax.set_yticks(y_pos)
        ax.set_yticklabels(df_sorted["group"].values, fontsize=10)
        ax.set_xlabel("Mean ΔC-index (baseline - ablated)", fontsize=11, fontweight="bold")
        ax.axvline(0, color="black", linestyle="--", linewidth=1.5)
        ax.set_title("Feature Group Importance: Ablation Study Results\n(Green = Important, Red = Not Important)", 
                    fontsize=13, fontweight="bold")
        ax.grid(axis="x", alpha=0.3)
        
        plt.tight_layout()
        plt.savefig(IMPORTANCE_DIR / "ablation_plot.png", dpi=150, bbox_inches="tight")
        logger.info("  Saved: ablation_plot.png")
        plt.close()
        
        # Summary plot: bars grouped by importance level
        fig, ax = plt.subplots(figsize=(10, 6))
        
        critical = ablation_agg_df[ablation_agg_df["mean_delta_cindex"] >= 0.03]
        moderate = ablation_agg_df[(ablation_agg_df["mean_delta_cindex"] >= 0.01) & 
                                   (ablation_agg_df["mean_delta_cindex"] < 0.03)]
        minor = ablation_agg_df[ablation_agg_df["mean_delta_cindex"] < 0.01]
        
        sizes = [len(critical), len(moderate), len(minor)]
        labels = [f"Critical\n(ΔC ≥ 0.03)\n{len(critical)} groups", 
                 f"Moderate\n(0.01 ≤ ΔC < 0.03)\n{len(moderate)} groups",
                 f"Minor\n(ΔC < 0.01)\n{len(minor)} groups"]
        colors_pie = ["#2ecc71", "#f39c12", "#e74c3c"]
        
        ax.pie(sizes, labels=labels, colors=colors_pie, autopct="%1.0f%%", startangle=90,
              textprops={"fontsize": 11, "fontweight": "bold"})
        ax.set_title("Feature Group Classification by Importance", fontsize=13, fontweight="bold")
        
        plt.tight_layout()
        plt.savefig(IMPORTANCE_DIR / "importance_classification.png", dpi=150, bbox_inches="tight")
        logger.info("  Saved: importance_classification.png")
        plt.close()
        
    except Exception as e:
        logger.warning(f"Visualization failed: {e}")
    
    # ===== Summary Report =====
    
    logger.info("\n" + "="*70)
    logger.info("SUMMARY REPORT")
    logger.info("="*70)
    
    logger.info(f"\nTotal feature groups analyzed: {len(ablation_agg_df)}")
    logger.info(f"Total folds: {N_CV_FOLDS}")
    
    logger.info(f"\n✅ CRITICAL FEATURES (ΔC-index ≥ 0.03):")
    for _, row in ablation_agg_df[ablation_agg_df["mean_delta_cindex"] >= 0.03].iterrows():
        logger.info(f"   • {row['group']:25s} (ΔC = {row['mean_delta_cindex']:+.4f})")
    
    logger.info(f"\n⚠️  MODERATE FEATURES (0.01 ≤ ΔC < 0.03):")
    for _, row in ablation_agg_df[(ablation_agg_df["mean_delta_cindex"] >= 0.01) & 
                                   (ablation_agg_df["mean_delta_cindex"] < 0.03)].iterrows():
        logger.info(f"   • {row['group']:25s} (ΔC = {row['mean_delta_cindex']:+.4f})")
    
    logger.info(f"\n⚪ MINOR FEATURES (ΔC < 0.01):")
    for _, row in ablation_agg_df[ablation_agg_df["mean_delta_cindex"] < 0.01].iterrows():
        logger.info(f"   • {row['group']:25s} (ΔC = {row['mean_delta_cindex']:+.4f})")
    
    logger.info("\n" + "="*70)
    logger.info("Analysis complete!")
    logger.info("="*70)


if __name__ == "__main__":
    analyze_feature_importance()
