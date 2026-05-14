"""
11_mlp_enriched.py
Train an MLP on enriched graph-level features (same 102-dim vector as RSF)
using patient-level GroupKFold CV, then compare to RSF results.

Run:
    python spatial_survival/11_mlp_enriched.py
"""

import sys
import copy
from pathlib import Path

import numpy as np
import pandas as pd
import torch
import torch.nn as nn
from torch.optim import Adam
from torch.optim.lr_scheduler import ReduceLROnPlateau
from sklearn.model_selection import GroupKFold
from sklearn.preprocessing import StandardScaler
from torch_geometric.data import Data

sys.path.insert(0, str(Path(__file__).parent))
from config import (
    OUTPUT_DIR, RESULTS_DIR, SEED, N_CV_FOLDS, VAL_FRACTION,
    LEARNING_RATE, WEIGHT_DECAY, MAX_EPOCHS, PATIENCE, LR_PATIENCE,
)
from utils import get_logger, ensure_dirs, set_seed, compute_cindex, CoxPHLoss

PYG_RAW_DIR = OUTPUT_DIR / "pyg_dataset" / "raw"
INDEX_PATH = OUTPUT_DIR / "pyg_dataset" / "dataset_index.csv"
MLP_DIR = RESULTS_DIR / "mlp_enriched"

logger = get_logger("mlp_enriched", RESULTS_DIR / "mlp_enriched.log")


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


# ---------------------------------------------------------------------------
# Model
# ---------------------------------------------------------------------------

class MLPEnriched(nn.Module):
    def __init__(self, in_dim: int, hidden_dim: int = 128, dropout: float = 0.3):
        super().__init__()
        self.net = nn.Sequential(
            nn.Linear(in_dim, hidden_dim),
            nn.ReLU(),
            nn.Dropout(dropout),
            nn.Linear(hidden_dim, hidden_dim // 2),
            nn.ReLU(),
            nn.Dropout(dropout),
            nn.Linear(hidden_dim // 2, 1),
        )

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        return self.net(x)


# ---------------------------------------------------------------------------
# Training helpers
# ---------------------------------------------------------------------------

def train_full_batch(model, optimizer, criterion, x, times, events):
    model.train()
    optimizer.zero_grad()
    risk = model(x).squeeze(1)
    loss = criterion(risk, times, events)
    if torch.isnan(loss):
        return float("nan")
    loss.backward()
    torch.nn.utils.clip_grad_norm_(model.parameters(), 1.0)
    optimizer.step()
    return float(loss.item())


@torch.no_grad()
def eval_full_batch(model, x, times, events):
    model.eval()
    risk = model(x).squeeze(1).cpu().numpy()
    t = times.cpu().numpy()
    e = events.cpu().numpy()
    ci = compute_cindex(risk, t, e)
    return ci, risk, t, e


# ---------------------------------------------------------------------------
# Main training loop
# ---------------------------------------------------------------------------

def train():
    set_seed(SEED)
    ensure_dirs(MLP_DIR, RESULTS_DIR)

    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")

    index_df = pd.read_csv(INDEX_PATH)
    logger.info(
        f"Dataset: {len(index_df)} samples, {index_df['patient_id'].nunique()} patients"
    )

    all_data = load_all_data(index_df)
    groups = index_df["patient_id"].to_numpy()

    logger.info("Building graph-level feature matrix (one-time) ...")
    features_all = build_graph_feature_matrix(all_data)
    y_times_all = np.array([float(d.y_time.item()) for d in all_data])
    y_events_all = np.array([float(d.y_event.item()) for d in all_data])

    gkf = GroupKFold(n_splits=N_CV_FOLDS)
    fold_results = []

    for fold_idx, (train_val_idx, test_idx) in enumerate(
        gkf.split(index_df, groups=groups)
    ):
        logger.info("=" * 60)
        logger.info(f"Fold {fold_idx + 1}/{N_CV_FOLDS}")

        train_val_patients = np.unique(groups[train_val_idx])
        rng = np.random.default_rng(SEED + fold_idx)
        rng.shuffle(train_val_patients)
        n_val_patients = max(1, int(len(train_val_patients) * VAL_FRACTION))
        val_patients = set(train_val_patients[:n_val_patients])

        train_idx = [i for i in train_val_idx if groups[i] not in val_patients]
        val_idx = [i for i in train_val_idx if groups[i] in val_patients]

        train_data = [all_data[i] for i in train_idx]
        val_data = [all_data[i] for i in val_idx]
        test_data = [all_data[i] for i in test_idx]

        logger.info(
            f"  Train: {len(train_data)} | Val: {len(val_data)} | Test: {len(test_data)}"
        )

        train_x = features_all[train_idx]
        val_x = features_all[val_idx]
        test_x = features_all[test_idx]
        train_x, val_x, test_x, scaler = scale_graph_features(train_x, val_x, test_x)

        y_train_times = y_times_all[train_idx]
        y_train_events = y_events_all[train_idx]
        y_val_times = y_times_all[val_idx]
        y_val_events = y_events_all[val_idx]
        y_test_times = y_times_all[test_idx]
        y_test_events = y_events_all[test_idx]

        x_train_t = torch.tensor(train_x, dtype=torch.float32, device=device)
        x_val_t = torch.tensor(val_x, dtype=torch.float32, device=device)
        x_test_t = torch.tensor(test_x, dtype=torch.float32, device=device)

        t_train_t = torch.tensor(y_train_times, dtype=torch.float32, device=device)
        e_train_t = torch.tensor(y_train_events, dtype=torch.float32, device=device)
        t_val_t = torch.tensor(y_val_times, dtype=torch.float32, device=device)
        e_val_t = torch.tensor(y_val_events, dtype=torch.float32, device=device)
        t_test_t = torch.tensor(y_test_times, dtype=torch.float32, device=device)
        e_test_t = torch.tensor(y_test_events, dtype=torch.float32, device=device)

        set_seed(SEED)
        model = MLPEnriched(in_dim=train_x.shape[1]).to(device)
        criterion = CoxPHLoss()
        optimizer = Adam(model.parameters(), lr=LEARNING_RATE, weight_decay=WEIGHT_DECAY)
        scheduler = ReduceLROnPlateau(
            optimizer, mode="max", patience=LR_PATIENCE, factor=0.5, min_lr=1e-6
        )

        best_val_ci = -1.0
        best_state = None
        no_improve = 0

        for epoch in range(1, MAX_EPOCHS + 1):
            train_loss = train_full_batch(
                model, optimizer, criterion, x_train_t, t_train_t, e_train_t
            )

            val_ci, _, _, _ = eval_full_batch(model, x_val_t, t_val_t, e_val_t)
            if not np.isnan(val_ci):
                scheduler.step(val_ci)

            if not np.isnan(val_ci) and val_ci > best_val_ci:
                best_val_ci = val_ci
                no_improve = 0
                best_state = copy.deepcopy(model.state_dict())
            else:
                no_improve += 1

            if epoch == 1 or epoch % 25 == 0:
                logger.info(
                    f"  Epoch {epoch:3d} | loss={train_loss:.4f} | val_C={val_ci:.4f}"
                )

            if no_improve >= PATIENCE:
                break

        if best_state is not None:
            model.load_state_dict(best_state)

        test_ci, risks, times, events = eval_full_batch(
            model, x_test_t, t_test_t, e_test_t
        )

        logger.info(
            f"  Fold {fold_idx + 1} | val_C={best_val_ci:.4f} | test_C={test_ci:.4f}"
        )

        pred_df = pd.DataFrame({
            "acquisition_id": [index_df.iloc[i]["acquisition_id"] for i in test_idx],
            "risk_score": risks,
            "y_time": times,
            "y_event": events.astype(int),
            "fold": fold_idx + 1,
        })
        pred_df.to_csv(MLP_DIR / f"fold_{fold_idx + 1}_predictions.csv", index=False)

        fold_results.append({
            "fold": fold_idx + 1,
            "n_train": len(train_data),
            "n_val": len(val_data),
            "n_test": len(test_data),
            "val_ci": best_val_ci,
            "test_ci": test_ci,
        })

    results_df = pd.DataFrame(fold_results)
    results_df.to_csv(MLP_DIR / "cv_results.csv", index=False)

    logger.info("=== MLP Enriched Cross-Validation Summary ===")
    logger.info(f"\n{results_df.to_string(index=False)}")
    logger.info(
        f"\nMean test C-index: {results_df['test_ci'].mean():.4f} "
        f"+- {results_df['test_ci'].std():.4f}"
    )

    # Comparison to RSF if available
    comparison_rows = []
    rsf_path = RESULTS_DIR / "cv_results.csv"
    if rsf_path.exists():
        rsf_df = pd.read_csv(rsf_path)
        comparison_rows.append({
            "model": "EnrichedRSF",
            "mean_test_cindex": rsf_df["test_ci"].mean(),
            "std_test_cindex": rsf_df["test_ci"].std(),
        })

    comparison_rows.append({
        "model": "EnrichedMLP",
        "mean_test_cindex": results_df["test_ci"].mean(),
        "std_test_cindex": results_df["test_ci"].std(),
    })

    comparison_df = pd.DataFrame(comparison_rows)
    comparison_df.to_csv(MLP_DIR / "comparison.csv", index=False)
    logger.info("\n=== Model Comparison (mean test C-index) ===")
    logger.info(f"\n{comparison_df.to_string(index=False)}")


if __name__ == "__main__":
    train()
