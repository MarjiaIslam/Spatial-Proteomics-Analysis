"""
06_training.py
Patient-level GroupKFold cross-validation for an enriched Random Survival Forest.

Each graph is collapsed into a fixed-length feature vector built from the
enriched graph representation:
    - node feature mean, std
    - edge feature mean, std
    - basic graph size statistics
    - interaction-type fractions

For each fold:
    - Train set: samples from (N_CV_FOLDS - 1) patient groups
    - Val   set: 10 % of train patients held out for reporting
    - Test  set: samples from the held-out patient group
    - StandardScaler fit on train graph features only (no leakage)
    - RandomSurvivalForest fit on train features

Results saved to:
    output/results/training_logs/fold_<k>_log.csv
    output/results/checkpoints/fold_<k>_best.pkl
    output/results/cv_results.csv

Run:
        python spatial_survival/06_training.py
"""

import sys
import copy
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

import importlib

sys.path.insert(0, str(Path(__file__).parent))
from config import (
    OUTPUT_DIR, RESULTS_DIR, SEED, N_CV_FOLDS, VAL_FRACTION, N_PROTEINS,
)
from utils import get_logger, ensure_dirs, set_seed, compute_cindex

_model_module = importlib.import_module("05_graphsage_model")
GraphSAGESurvival = _model_module.GraphSAGESurvival

PYG_RAW_DIR  = OUTPUT_DIR / "pyg_dataset" / "raw"
INDEX_PATH   = OUTPUT_DIR / "pyg_dataset" / "dataset_index.csv"
CKPT_DIR     = RESULTS_DIR / "checkpoints"
LOG_DIR      = RESULTS_DIR / "training_logs"

logger = get_logger("training", RESULTS_DIR / "training.log")


# ---------------------------------------------------------------------------
# Dataset helpers
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
    """Collapse one enriched graph into a fixed-length feature vector."""
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


def scale_graph_features(train_x: np.ndarray,
                         val_x: np.ndarray,
                         test_x: np.ndarray):
    scaler = StandardScaler()
    scaler.fit(train_x)
    return (
        scaler.transform(train_x),
        scaler.transform(val_x),
        scaler.transform(test_x),
        scaler,
    )


def scale_node_features(train_data: list[Data],
                         val_data:   list[Data],
                         test_data:  list[Data]):
    """
    Fit StandardScaler on train node features, apply to all splits in-place.
    Only the first N_PROTEINS features are continuous and scaled;
    the 6 spatial features (indices 39–44) are also scaled.
    """
    scaler = StandardScaler()

    # Fit on train
    train_x = np.vstack([d.x.numpy() for d in train_data])
    scaler.fit(train_x)

    def _apply(data_list):
        for d in data_list:
            x_np  = d.x.numpy()
            x_sc  = scaler.transform(x_np)
            d.x   = torch.tensor(x_sc, dtype=torch.float32)

    _apply(train_data)
    _apply(val_data)
    _apply(test_data)
    return scaler


# ---------------------------------------------------------------------------
# One training epoch
# ---------------------------------------------------------------------------

def train_epoch(model, loader, optimizer, criterion, device):
    model.train()
    total_loss = 0.0
    n_batches  = 0
    for batch in loader:
        batch = batch.to(device)
        optimizer.zero_grad()
        risk   = model(batch).squeeze(1)
        times  = batch.y_time.squeeze()
        events = batch.y_event.squeeze()
        loss   = criterion(risk, times, events)
        if torch.isnan(loss):
            continue
        loss.backward()
        torch.nn.utils.clip_grad_norm_(model.parameters(), max_norm=1.0)
        optimizer.step()
        total_loss += loss.item()
        n_batches  += 1
    return total_loss / max(n_batches, 1)


# ---------------------------------------------------------------------------
# Evaluation
# ---------------------------------------------------------------------------

@torch.no_grad()
def evaluate(model, loader, device):
    model.eval()
    risks_all  = []
    times_all  = []
    events_all = []
    for batch in loader:
        batch = batch.to(device)
        risk   = model(batch).squeeze(1).cpu().numpy()
        times  = batch.y_time.squeeze().cpu().numpy()
        events = batch.y_event.squeeze().cpu().numpy()
        risks_all.append(risk)
        times_all.append(times)
        events_all.append(events)

    risks  = np.concatenate(risks_all)
    times  = np.concatenate(times_all)
    events = np.concatenate(events_all)
    ci     = compute_cindex(risks, times, events)
    return ci, risks, times, events


# ---------------------------------------------------------------------------
# Main training loop
# ---------------------------------------------------------------------------

def train():
    set_seed(SEED)
    ensure_dirs(CKPT_DIR, LOG_DIR, RESULTS_DIR)

    index_df = pd.read_csv(INDEX_PATH)
    logger.info(f"Dataset: {len(index_df)} samples, {index_df['patient_id'].nunique()} patients")

    all_data = load_all_data(index_df)
    groups   = index_df["patient_id"].to_numpy()

    gkf = GroupKFold(n_splits=N_CV_FOLDS)
    fold_results = []

    for fold_idx, (train_val_idx, test_idx) in enumerate(gkf.split(index_df, groups=groups)):
        logger.info(f"\n{'='*60}")
        logger.info(f"Fold {fold_idx+1}/{N_CV_FOLDS}")

        train_val_patients = np.unique(groups[train_val_idx])
        rng = np.random.default_rng(SEED + fold_idx)
        rng.shuffle(train_val_patients)
        n_val_patients = max(1, int(len(train_val_patients) * VAL_FRACTION))
        val_patients = set(train_val_patients[:n_val_patients])

        train_idx = [i for i in train_val_idx if groups[i] not in val_patients]
        val_idx   = [i for i in train_val_idx if groups[i] in val_patients]

        train_data = [all_data[i] for i in train_idx]
        val_data   = [all_data[i] for i in val_idx]
        test_data  = [all_data[i] for i in test_idx]

        logger.info(f"  Train: {len(train_data)} | Val: {len(val_data)} | Test: {len(test_data)}")

        train_x = build_graph_feature_matrix(train_data)
        val_x   = build_graph_feature_matrix(val_data)
        test_x  = build_graph_feature_matrix(test_data)
        train_x, val_x, test_x, scaler = scale_graph_features(train_x, val_x, test_x)

        y_train_times  = np.array([float(d.y_time.item()) for d in train_data])
        y_train_events = np.array([bool(d.y_event.item()) for d in train_data])
        y_val_times    = np.array([float(d.y_time.item()) for d in val_data])
        y_val_events   = np.array([bool(d.y_event.item()) for d in val_data])
        y_test_times   = np.array([float(d.y_time.item()) for d in test_data])
        y_test_events  = np.array([bool(d.y_event.item()) for d in test_data])

        y_train = Surv.from_arrays(y_train_events, y_train_times)

        rsf = RandomSurvivalForest(
            n_estimators=300,
            min_samples_split=10,
            min_samples_leaf=3,
            max_features="sqrt",
            n_jobs=-1,
            random_state=SEED,
        )
        rsf.fit(train_x, y_train)

        val_risk  = rsf.predict(val_x)
        test_risk = rsf.predict(test_x)
        val_ci    = compute_cindex(val_risk, y_val_times, y_val_events)
        test_ci   = compute_cindex(test_risk, y_test_times, y_test_events)

        logger.info(f"  Fold {fold_idx+1} | val_C={val_ci:.4f} | test_C={test_ci:.4f}")

        ckpt_path = CKPT_DIR / f"fold_{fold_idx+1}_best.pkl"
        with open(ckpt_path, "wb") as handle:
            pickle.dump({
                "model": rsf,
                "scaler": scaler,
                "feature_dim": train_x.shape[1],
                "fold": fold_idx + 1,
            }, handle)

        log_df = pd.DataFrame([{
            "fold": fold_idx + 1,
            "n_train": len(train_data),
            "n_val": len(val_data),
            "n_test": len(test_data),
            "val_cindex": val_ci,
            "test_cindex": test_ci,
        }])
        log_df.to_csv(LOG_DIR / f"fold_{fold_idx+1}_log.csv", index=False)

        pred_df = pd.DataFrame({
            "acquisition_id": [index_df.iloc[i]["acquisition_id"] for i in test_idx],
            "risk_score": test_risk,
            "y_time": y_test_times,
            "y_event": y_test_events.astype(int),
            "fold": fold_idx + 1,
        })
        pred_df.to_csv(RESULTS_DIR / f"fold_{fold_idx+1}_predictions.csv", index=False)

        fold_results.append({
            "fold": fold_idx + 1,
            "n_train": len(train_data),
            "n_val": len(val_data),
            "n_test": len(test_data),
            "val_ci": val_ci,
            "test_ci": test_ci,
        })

    results_df = pd.DataFrame(fold_results)
    results_df.to_csv(RESULTS_DIR / "cv_results.csv", index=False)

    logger.info("\n=== RSF Cross-Validation Summary ===")
    logger.info(f"\n{results_df.to_string(index=False)}")
    logger.info(f"\nMean test C-index: {results_df['test_ci'].mean():.4f} "
                f"± {results_df['test_ci'].std():.4f}")


if __name__ == "__main__":
    train()
