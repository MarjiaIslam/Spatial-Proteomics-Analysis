"""
MIBI Protein → Gene-Expression Dataset Conversion Pipeline
============================================================
Reads arcsinh-normalized MIBI-SoC protein data, maps markers to HUGO gene
symbols using UniProt-verified mappings, and produces a NEW gene-expression
dataset CSV with gene names as columns.

Input:  data/dataset_info/labeled_arcsinh_norm_data.csv  (already normalized)
        data/dataset_info/cell_locations_and_labels.csv
        data/dataset_info/sample_metadata.csv

Primary Output:
  output/gene_expression_dataset.csv          ← new dataset (genes as columns)
  output/uniprot_protein_to_gene_mapping.csv  ← mapping provenance table

Secondary Outputs:
  output/mibi_as_gene_expression.h5ad         ← AnnData for Scanpy workflows
  output/cell_metadata.csv                    ← cell annotations + UMAP + spatial
"""

import os
import warnings
import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
from scipy import stats
from pathlib import Path

from protein_to_gene_map import (
    UNIPROT_MAP, PROTEIN_TO_GENE, GENE_COLUMN_NAMES,
    ISOFORM_CONFLICTS, AMBIGUOUS_MARKERS, save_mapping_csv,
)

warnings.filterwarnings("ignore", category=FutureWarning)

# ──────────────────────────────────────────────────────────
# 0.  PATHS & UniProt MAPPING EXPORT
# ──────────────────────────────────────────────────────────
DATA_DIR  = Path("data/dataset_info")
OUT_DIR   = Path("output")
OUT_DIR.mkdir(exist_ok=True)

EXPR_FILE     = DATA_DIR / "labeled_arcsinh_norm_data.csv"
SPATIAL_FILE  = DATA_DIR / "cell_locations_and_labels.csv"
META_FILE     = DATA_DIR / "sample_metadata.csv"

# Save the full UniProt mapping provenance table first
print("── Step 0: Exporting UniProt protein-to-gene mapping ──")
mapping_df = save_mapping_csv(OUT_DIR / "uniprot_protein_to_gene_mapping.csv")


# ──────────────────────────────────────────────────────────
# 1.  LOAD DATA
# ──────────────────────────────────────────────────────────
print("── Step 1: Loading data ──")

expr_df = pd.read_csv(EXPR_FILE)
print(f"  Expression matrix: {expr_df.shape[0]:,} cells × {expr_df.shape[1]} columns")

spatial_df = pd.read_csv(SPATIAL_FILE)
print(f"  Spatial info:      {spatial_df.shape[0]:,} cells")

meta_df = pd.read_csv(META_FILE)
print(f"  Sample metadata:   {meta_df.shape[0]} samples")


# ──────────────────────────────────────────────────────────
# 2.  SEPARATE EXPRESSION FROM METADATA COLUMNS
# ──────────────────────────────────────────────────────────
print("\n── Step 2: Separating expression from annotations ──")

# Protein marker columns = all columns in PROTEIN_TO_GENE keys
protein_cols = [c for c in expr_df.columns if c in PROTEIN_TO_GENE]
meta_cols_in_expr = [c for c in expr_df.columns if c not in PROTEIN_TO_GENE]

print(f"  Protein columns found: {len(protein_cols)}/40")
print(f"  Annotation columns:    {meta_cols_in_expr}")

# Expression matrix (cells × proteins)
X_protein = expr_df[protein_cols].values.astype(np.float32)

# Map column names: protein → gene symbol (using UniProt-verified GENE_COLUMN_NAMES)
gene_names = [GENE_COLUMN_NAMES[p] for p in protein_cols]

print(f"  Gene symbols mapped:   {len(gene_names)} (via UniProt Swiss-Prot)")
print(f"  UniProt mapping source: output/uniprot_protein_to_gene_mapping.csv")


# ──────────────────────────────────────────────────────────
# 3.  HANDLE PTPRC ISOFORMS  (CD45 / CD45RA / CD45RO)
# ──────────────────────────────────────────────────────────
print("\n── Step 3: Handling PTPRC isoform markers ──")

# Strategy: Keep all three as separate features with suffixed names
# (PTPRC, PTPRC_RA, PTPRC_RO) so they remain distinguishable.
# For pathway analysis, a single "PTPRC" summary can be derived later.
#
# This avoids information loss: CD45RA marks naive T cells,
# CD45RO marks memory T cells — collapsing them loses biology.

ptprc_indices = [i for i, g in enumerate(gene_names) if "PTPRC" in g]
print(f"  PTPRC-related features kept separate: {[gene_names[i] for i in ptprc_indices]}")
print("  → For gene-set analysis, use PTPRC only (total CD45)")


# ──────────────────────────────────────────────────────────
# 4.  BUILD AnnData OBJECT
# ──────────────────────────────────────────────────────────
print("\n── Step 4: Building AnnData object ──")

adata = ad.AnnData(
    X=X_protein,
    dtype=np.float32,
)
adata.var_names = gene_names
adata.var["protein_name"]       = protein_cols
adata.var["gene_symbol"]        = [PROTEIN_TO_GENE[p] for p in protein_cols]
adata.var["uniprot_accession"]  = [UNIPROT_MAP[p][1] for p in protein_cols]
adata.var["uniprot_protein"]    = [UNIPROT_MAP[p][2] for p in protein_cols]

# Mark ambiguous and isoform genes
adata.var["is_ambiguous"] = [p in AMBIGUOUS_MARKERS for p in protein_cols]
adata.var["is_isoform"]   = ["PTPRC" in g for g in gene_names]

# Store raw expression as a layer
adata.layers["arcsinh_norm"] = X_protein.copy()

# Cell annotations from expression file
if "sample_id" in expr_df.columns:
    adata.obs["sample_id"] = expr_df["sample_id"].values
if "cell_id" in expr_df.columns:
    adata.obs["cell_id"] = expr_df["cell_id"].values.astype(str)
if "cluster" in expr_df.columns:
    adata.obs["cluster_id"] = expr_df["cluster"].values
if "cluster_label" in expr_df.columns:
    adata.obs["cell_type"] = expr_df["cluster_label"].values

# Generate unique obs index
adata.obs_names = [
    f"{sid}_{cid}"
    for sid, cid in zip(adata.obs["sample_id"], adata.obs["cell_id"])
]

print(f"  AnnData: {adata.shape[0]:,} cells × {adata.shape[1]} genes")


# ──────────────────────────────────────────────────────────
# 5.  ADD SPATIAL COORDINATES
# ──────────────────────────────────────────────────────────
print("\n── Step 5: Adding spatial coordinates ──")

# Build a lookup key that matches between the two dataframes
spatial_df["_key"] = (
    spatial_df["ACQUISITION_ID"].astype(str) + "_" +
    spatial_df["CELL_ID"].astype(str)
)
spatial_lookup = spatial_df.set_index("_key")[["X", "Y"]].to_dict("index")

coords = np.zeros((adata.n_obs, 2), dtype=np.float64)
matched = 0
for i, obs_name in enumerate(adata.obs_names):
    if obs_name in spatial_lookup:
        coords[i, 0] = spatial_lookup[obs_name]["X"]
        coords[i, 1] = spatial_lookup[obs_name]["Y"]
        matched += 1

adata.obsm["spatial"] = coords
print(f"  Spatial coords matched: {matched:,}/{adata.n_obs:,} cells")


# ──────────────────────────────────────────────────────────
# 6.  ADD SAMPLE-LEVEL CLINICAL METADATA
# ──────────────────────────────────────────────────────────
print("\n── Step 6: Merging clinical metadata ──")

meta_indexed = meta_df.set_index("acquisition_id")
clinical_cols = [
    "patient_id", "tissue_type", "primarysite", "larynxpharynx",
    "hpvstatus_new", "survival_status", "survival_day", "recurred",
]

for col in clinical_cols:
    if col in meta_indexed.columns:
        mapping = meta_indexed[col].to_dict()
        adata.obs[col] = adata.obs["sample_id"].map(mapping)
        non_null = adata.obs[col].notna().sum()
        print(f"  {col}: {non_null:,} cells annotated")


# ──────────────────────────────────────────────────────────
# 7.  TECHNOLOGY-SPECIFIC PREPROCESSING
# ──────────────────────────────────────────────────────────
print("\n── Step 7: Technology-specific preprocessing ──")

# 7a. Data is already arcsinh-normalized — NO further normalization needed.
#     arcsinh(x/cofactor) is the standard for mass-spec proteomics (CyTOF/MIBI).
#     This is analogous to log1p in scRNA-seq.
print("  ✓ Data already arcsinh-normalized (equivalent to log1p for RNA)")

# 7b. Clip extreme outliers (>99.9th percentile per gene)
#     MIBI can have hot pixels / segmentation artifacts
clip_threshold = 0.999
for j in range(adata.n_vars):
    col = adata.X[:, j]
    upper = np.quantile(col, clip_threshold)
    n_clipped = (col > upper).sum()
    adata.X[col > upper, j] = upper
    if n_clipped > 0:
        pass  # silently clip

print(f"  ✓ Clipped values above {clip_threshold*100}th percentile per gene")

# 7c. Min-max scale to [0, 1] range per gene (optional, helps comparability)
#     Store this as a separate layer — keep arcsinh_norm as ground truth
from sklearn.preprocessing import MinMaxScaler
scaler = MinMaxScaler()
adata.layers["scaled_0_1"] = scaler.fit_transform(adata.X).astype(np.float32)
print("  ✓ Min-max scaled layer added ('scaled_0_1')")

# 7d. Z-score layer (useful for heatmaps, gene-set scoring)
from scipy.stats import zscore
adata.layers["zscore"] = zscore(adata.X, axis=0).astype(np.float32)
print("  ✓ Z-score layer added ('zscore')")


# ──────────────────────────────────────────────────────────
# 8.  SAVE — PRIMARY OUTPUT: new gene-expression dataset CSV
#     (saved BEFORE UMAP so large datasets don't block output)
# ──────────────────────────────────────────────────────────
print("\n── Step 8: Saving primary outputs ──")

# ─── A. PRIMARY OUTPUT: Combined gene-expression dataset CSV ───
gene_expr_df = pd.DataFrame(
    adata.X, columns=adata.var_names, index=adata.obs_names
)

# Attach cell-level metadata columns
gene_expr_df.insert(0, "sample_id",    adata.obs["sample_id"].values)
gene_expr_df.insert(1, "cell_id",      adata.obs["cell_id"].values)
gene_expr_df.insert(2, "cell_type",    adata.obs["cell_type"].values)
gene_expr_df.insert(3, "cluster_id",   adata.obs["cluster_id"].values)
gene_expr_df.insert(4, "spatial_x",    adata.obsm["spatial"][:, 0])
gene_expr_df.insert(5, "spatial_y",    adata.obsm["spatial"][:, 1])

# Attach clinical metadata
for col in ["patient_id", "tissue_type", "primarysite", "hpvstatus_new",
            "survival_status", "survival_day", "recurred"]:
    if col in adata.obs.columns:
        gene_expr_df.insert(len([c for c in gene_expr_df.columns
                                 if c not in list(adata.var_names)]),
                            col, adata.obs[col].values)

dataset_path = OUT_DIR / "gene_expression_dataset.csv"
gene_expr_df.to_csv(dataset_path, index_label="cell_barcode")
n_meta_cols = len([c for c in gene_expr_df.columns if c not in list(adata.var_names)])
n_gene_cols = len(adata.var_names)
print(f"  ★ PRIMARY OUTPUT: {dataset_path}")
print(f"    {gene_expr_df.shape[0]:,} cells × "
      f"({n_meta_cols} metadata + {n_gene_cols} gene) columns")

# ─── B. Gene-only matrix (no metadata) ───
pure_expr = pd.DataFrame(adata.X, columns=adata.var_names, index=adata.obs_names)
pure_expr.to_csv(OUT_DIR / "gene_expression_matrix.csv")
print(f"  Pure matrix: {OUT_DIR / 'gene_expression_matrix.csv'}")

# ─── C. AnnData (without UMAP — added in step 9) ───
out_path = OUT_DIR / "mibi_as_gene_expression.h5ad"


# ──────────────────────────────────────────────────────────
# 9.  OPTIONAL: SCANPY EMBEDDING (PCA → UMAP → Leiden)
#     For datasets >500K cells this can take 30-60 minutes.
#     Set RUN_EMBEDDING=False to skip and still get all CSVs.
# ──────────────────────────────────────────────────────────
RUN_EMBEDDING = True

if RUN_EMBEDDING:
    print("\n── Step 9: Running Scanpy embedding pipeline ──")
    print("  (PCA + UMAP on 2M cells — this may take 20-60 min)")

    # PCA
    n_components = min(30, adata.n_vars - 1)
    sc.pp.pca(adata, n_comps=n_components)
    print(f"  PCA: {n_components} components")

    # Neighbors (kNN graph)
    sc.pp.neighbors(adata, n_neighbors=15, n_pcs=n_components)
    print("  Neighbors computed")

    # UMAP
    sc.tl.umap(adata)
    print("  UMAP computed")

    # Leiden clustering
    sc.tl.leiden(adata, resolution=0.8, key_added="leiden_0.8")
    print("  Leiden clustering computed")

    # ─── Save metadata with UMAP coords ───
    obs_export = adata.obs.copy()
    obs_export["UMAP_1"] = adata.obsm["X_umap"][:, 0]
    obs_export["UMAP_2"] = adata.obsm["X_umap"][:, 1]
    obs_export["spatial_x"] = adata.obsm["spatial"][:, 0]
    obs_export["spatial_y"] = adata.obsm["spatial"][:, 1]
    obs_export.to_csv(OUT_DIR / "cell_metadata.csv")
    print(f"  Metadata: {OUT_DIR / 'cell_metadata.csv'}")
else:
    n_components = min(30, adata.n_vars - 1)
    print("\n── Step 9: Skipped (RUN_EMBEDDING=False) ──")
    obs_export = adata.obs.copy()
    obs_export["spatial_x"] = adata.obsm["spatial"][:, 0]
    obs_export["spatial_y"] = adata.obsm["spatial"][:, 1]
    obs_export.to_csv(OUT_DIR / "cell_metadata.csv")
    print(f"  Metadata (no UMAP): {OUT_DIR / 'cell_metadata.csv'}")

# ─── Save AnnData ───
adata.write(out_path)
print(f"  AnnData:  {out_path}")

# ─── Summary ───
print(f"\n{'='*60}")
print(f"  PIPELINE COMPLETE")
print(f"{'='*60}")
print(f"  Cells:            {adata.n_obs:,}")
print(f"  Gene features:    {adata.n_vars} (from 39 MIBI protein markers)")
print(f"  Mapping source:   UniProt Swiss-Prot (reviewed human entries)")
print(f"  Normalization:    arcsinh (pre-applied, untouched)")
print(f"  Extra layers:     arcsinh_norm, scaled_0_1, zscore")
umap_status = f"PCA ({n_components}D), UMAP (2D), Leiden" if RUN_EMBEDDING else "skipped"
print(f"  Embeddings:       {umap_status}")
print(f"")
print(f"  Output files:")
print(f"    {dataset_path}  ← USE THIS (new gene-expression dataset)")
print(f"    {OUT_DIR / 'uniprot_protein_to_gene_mapping.csv'}  ← mapping provenance")
print(f"    {out_path}  ← for Scanpy/Python workflows")
print(f"    {OUT_DIR / 'cell_metadata.csv'}  ← for R interop")
print(f"{'='*60}")
