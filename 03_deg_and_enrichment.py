"""
DEG Analysis and Gene Set Enrichment on Protein-as-Gene Data
=============================================================
Performs differential expression and enrichment analysis using
the converted MIBI protein data.

Expects: output/mibi_as_gene_expression.h5ad
"""

import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
from pathlib import Path

OUT_DIR = Path("output")
FIG_DIR = OUT_DIR / "deg_gsea_results"
FIG_DIR.mkdir(exist_ok=True, parents=True)

from protein_to_gene_map import PROTEIN_TO_GENE

# ──────────────────────────────────────────────────────────
# Load
# ──────────────────────────────────────────────────────────
adata = sc.read_h5ad(OUT_DIR / "mibi_as_gene_expression.h5ad")
print(f"Loaded: {adata.shape[0]:,} cells × {adata.shape[1]} genes")


# ──────────────────────────────────────────────────────────
# 1.  DEG: Cell type vs rest (Wilcoxon)
# ──────────────────────────────────────────────────────────
print("\n── DEG: each cell type vs rest ──")
sc.tl.rank_genes_groups(adata, groupby="cell_type", method="wilcoxon",
                        pts=True, key_added="deg_celltype")

# Export
all_degs = []
for ct in adata.obs["cell_type"].cat.categories:
    df = sc.get.rank_genes_groups_df(adata, group=ct, key="deg_celltype")
    df["cell_type"] = ct
    all_degs.append(df)
deg_df = pd.concat(all_degs)
deg_df.to_csv(FIG_DIR / "deg_celltype_vs_rest.csv", index=False)
print(f"  Saved {len(deg_df)} DEG results → deg_celltype_vs_rest.csv")

# Dotplot
sc.pl.rank_genes_groups_dotplot(
    adata, n_genes=5, key="deg_celltype", show=False,
    standard_scale="var"
)
plt.savefig(FIG_DIR / "deg_dotplot.png", dpi=150, bbox_inches="tight")
plt.close()


# ──────────────────────────────────────────────────────────
# 2.  DEG: HPV+ vs HPV- (within tumor cells)
# ──────────────────────────────────────────────────────────
print("\n── DEG: HPV+ vs HPV- in tumor cells ──")

tumor_mask = adata.obs["cell_type"].str.startswith("Tumor")
adata_tumor = adata[tumor_mask].copy()

if "hpvstatus_new" in adata_tumor.obs.columns:
    adata_tumor.obs["hpv_group"] = adata_tumor.obs["hpvstatus_new"].map(
        {0: "HPV-", 0.0: "HPV-", 1: "HPV+", 1.0: "HPV+"}
    ).astype("category")

    if adata_tumor.obs["hpv_group"].nunique() >= 2:
        sc.tl.rank_genes_groups(adata_tumor, groupby="hpv_group",
                                method="wilcoxon", key_added="deg_hpv")
        hpv_degs = sc.get.rank_genes_groups_df(adata_tumor, group="HPV+",
                                                key="deg_hpv")
        hpv_degs.to_csv(FIG_DIR / "deg_hpv_pos_vs_neg_tumor.csv", index=False)
        print(f"  Saved HPV+ vs HPV- DEGs ({len(hpv_degs)} genes)")


# ──────────────────────────────────────────────────────────
# 3.  Gene Set Scoring (manual curated sets)
# ──────────────────────────────────────────────────────────
print("\n── Gene Set Scoring ──")

# Since we have only 40 genes, formal GSEA databases won't match well.
# Instead, define biologically meaningful gene sets from our panel.
GENE_SETS = {
    "T_cell_activation": ["CD3E", "CD4", "CD8A", "CD69", "GZMB", "ICOS"],
    "Immune_checkpoint": ["PDCD1", "CD274", "CTLA4", "ICOS"],
    "Antigen_presentation": ["HLA-DRA", "ITGAX", "ITGAM", "CD68"],
    "Tumor_markers": ["KRT18", "ANO1", "MKI67", "CDKN2A"],
    "Immune_infiltration": ["PTPRC", "CD3E", "MS4A1", "CD68", "ITGAM"],
    "Regulatory_T_cell": ["CD4", "FOXP3", "CTLA4", "ICOS"],
    "Cytotoxicity": ["CD8A", "GZMB", "NCAM1", "B3GAT1"],
    "Stromal": ["VIM", "ACTA2", "COL4A1", "PECAM1"],
}

for gs_name, genes in GENE_SETS.items():
    available = [g for g in genes if g in adata.var_names]
    if len(available) < 2:
        continue
    sc.tl.score_genes(adata, gene_list=available, score_name=gs_name)
    print(f"  {gs_name}: {len(available)}/{len(genes)} genes used")

# Plot gene set scores on UMAP
score_cols = [gs for gs in GENE_SETS.keys() if gs in adata.obs.columns]
if score_cols:
    sc.pl.umap(adata, color=score_cols, ncols=3, frameon=False, show=False)
    plt.savefig(FIG_DIR / "geneset_scores_umap.png", dpi=150, bbox_inches="tight")
    plt.close()
    print(f"  Saved: geneset_scores_umap.png")

# Violin plots per cell type
for gs in score_cols:
    sc.pl.violin(adata, keys=gs, groupby="cell_type", rotation=90, show=False)
    plt.savefig(FIG_DIR / f"geneset_violin_{gs}.png", dpi=150, bbox_inches="tight")
    plt.close()


# ──────────────────────────────────────────────────────────
# 4.  DEG: Recurrence vs No Recurrence
# ──────────────────────────────────────────────────────────
print("\n── DEG: Recurred vs Not Recurred ──")

if "recurred" in adata.obs.columns:
    adata.obs["recurrence_group"] = adata.obs["recurred"].map(
        {0: "No recurrence", 0.0: "No recurrence",
         1: "Recurred", 1.0: "Recurred"}
    ).astype("category")

    if adata.obs["recurrence_group"].nunique() >= 2:
        sc.tl.rank_genes_groups(adata, groupby="recurrence_group",
                                method="wilcoxon", key_added="deg_recurrence")
        rec_degs = sc.get.rank_genes_groups_df(
            adata, group="Recurred", key="deg_recurrence"
        )
        rec_degs.to_csv(FIG_DIR / "deg_recurred_vs_not.csv", index=False)
        print(f"  Saved recurrence DEGs ({len(rec_degs)} genes)")


# ──────────────────────────────────────────────────────────
# 5.  Save updated AnnData with scores
# ──────────────────────────────────────────────────────────
adata.write(OUT_DIR / "mibi_as_gene_expression.h5ad")
print(f"\n✅ DEG/GSEA pipeline complete. Results in: {FIG_DIR}")
