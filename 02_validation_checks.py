"""
Validation: Does the converted protein data behave like real gene expression?
=============================================================================
Runs sanity checks and generates diagnostic plots.

Expects: output/mibi_as_gene_expression.h5ad (from 01_protein_to_gene_pipeline.py)
"""

import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
from scipy import stats
from pathlib import Path

OUT_DIR = Path("output")
FIG_DIR = OUT_DIR / "validation_figures"
FIG_DIR.mkdir(exist_ok=True, parents=True)


# ──────────────────────────────────────────────────────────
# Load
# ──────────────────────────────────────────────────────────
print("Loading converted data...")
adata = sc.read_h5ad(OUT_DIR / "mibi_as_gene_expression.h5ad")
print(f"  {adata.shape[0]:,} cells × {adata.shape[1]} genes\n")


# ══════════════════════════════════════════════════════════
# CHECK 1: Distribution shape — gene expression is roughly
#           log-normal; arcsinh-transformed protein data
#           should look similar to log1p-transformed RNA.
# ══════════════════════════════════════════════════════════
print("── Check 1: Per-gene distribution shape ──")

fig, axes = plt.subplots(5, 8, figsize=(24, 15))
axes = axes.flatten()
for i, gene in enumerate(adata.var_names):
    ax = axes[i] if i < len(axes) else None
    if ax is None:
        break
    vals = adata.X[:, i]
    ax.hist(vals, bins=60, color="steelblue", alpha=0.7, density=True)
    ax.set_title(gene, fontsize=8)
    ax.tick_params(labelsize=6)
    # Shapiro on subsample
    sub = np.random.choice(vals, min(5000, len(vals)), replace=False)
    _, p_normal = stats.normaltest(sub)
    ax.text(0.95, 0.95, f"p={p_normal:.1e}", transform=ax.transAxes,
            fontsize=5, ha="right", va="top")

for j in range(i + 1, len(axes)):
    axes[j].set_visible(False)

plt.suptitle("Per-gene value distributions (arcsinh-normalized)\n"
             "Should resemble log-normal RNA distributions", fontsize=12)
plt.tight_layout()
plt.savefig(FIG_DIR / "check1_per_gene_distributions.png", dpi=150)
plt.close()
print("  Saved: check1_per_gene_distributions.png")


# ══════════════════════════════════════════════════════════
# CHECK 2: Marker-cell-type concordance — known markers
#           should be high in their expected cell types.
# ══════════════════════════════════════════════════════════
print("\n── Check 2: Marker–cell type concordance ──")

EXPECTED_MARKERS = {
    "CD3E":   ["CD4 T cell", "CD8 T cell"],
    "CD4":    ["CD4 T cell"],
    "CD8A":   ["CD8 T cell"],
    "MS4A1":  ["B cell"],                       # CD20
    "CD68":   ["Macrophage"],
    "PECAM1": ["Vessel", "Lymph vessel"],        # CD31
    "VIM":    ["Stromal / Fibroblast"],
    "MKI67":  ["Tumor (Ki67+)"],
    "PDPN":   ["Tumor (Podo+)", "Lymph vessel"], # Podoplanin
    "KRT18":  ["Tumor", "Tumor (CD15+)", "Tumor (CD20+)",
               "Tumor (CD21+)", "Tumor (Ki67+)", "Tumor (Podo+)"],
    "ITGAM":  ["Macrophage", "Granulocyte", "APC"],  # CD11b
    "FOXP3":  ["CD4 T cell"],                    # Tregs subset
    "GZMB":   ["CD8 T cell"],
}

concordance_results = []
for gene, expected_types in EXPECTED_MARKERS.items():
    if gene not in adata.var_names:
        continue
    gene_idx = list(adata.var_names).index(gene)
    vals = adata.X[:, gene_idx]

    mask_expected = adata.obs["cell_type"].isin(expected_types)
    mean_in  = vals[mask_expected].mean()
    mean_out = vals[~mask_expected].mean()
    fc = mean_in / max(mean_out, 1e-6)

    # Mann-Whitney U test
    u_stat, u_pval = stats.mannwhitneyu(
        vals[mask_expected], vals[~mask_expected], alternative="greater"
    )

    status = "✓ PASS" if (fc > 1.2 and u_pval < 0.01) else "✗ FAIL"
    concordance_results.append({
        "gene": gene,
        "expected_types": ", ".join(expected_types),
        "mean_in": round(mean_in, 3),
        "mean_out": round(mean_out, 3),
        "fold_change": round(fc, 2),
        "p_value": f"{u_pval:.2e}",
        "status": status,
    })
    print(f"  {status}  {gene:8s}  FC={fc:.2f}  p={u_pval:.2e}  "
          f"(expected in: {', '.join(expected_types)})")

conc_df = pd.DataFrame(concordance_results)
conc_df.to_csv(FIG_DIR / "check2_marker_concordance.csv", index=False)


# ══════════════════════════════════════════════════════════
# CHECK 3: Gene-gene correlation structure — should see
#           known co-expression (e.g., CD3E with CD4/CD8A)
# ══════════════════════════════════════════════════════════
print("\n── Check 3: Gene-gene correlation matrix ──")

corr = np.corrcoef(adata.X.T)
corr_df = pd.DataFrame(corr, index=adata.var_names, columns=adata.var_names)

fig, ax = plt.subplots(figsize=(14, 12))
im = ax.imshow(corr, cmap="RdBu_r", vmin=-0.5, vmax=0.5)
ax.set_xticks(range(len(adata.var_names)))
ax.set_yticks(range(len(adata.var_names)))
ax.set_xticklabels(adata.var_names, rotation=90, fontsize=6)
ax.set_yticklabels(adata.var_names, fontsize=6)
plt.colorbar(im, label="Pearson r")
plt.title("Gene-gene correlation (arcsinh-normalized protein)")
plt.tight_layout()
plt.savefig(FIG_DIR / "check3_gene_correlation_matrix.png", dpi=150)
plt.close()
corr_df.to_csv(FIG_DIR / "check3_correlation_matrix.csv")
print("  Saved: check3_gene_correlation_matrix.png")

# Spot-check expected correlations
expected_pos = [("CD3E", "CD4"), ("CD3E", "CD8A"), ("ITGAM", "CD68")]
expected_neg = [("CD3E", "KRT18"), ("MS4A1", "CD3E")]

print("  Expected positive correlations:")
for g1, g2 in expected_pos:
    if g1 in corr_df.index and g2 in corr_df.columns:
        r = corr_df.loc[g1, g2]
        status = "✓" if r > 0.05 else "✗"
        print(f"    {status} {g1} vs {g2}: r={r:.3f}")

print("  Expected negative/low correlations:")
for g1, g2 in expected_neg:
    if g1 in corr_df.index and g2 in corr_df.columns:
        r = corr_df.loc[g1, g2]
        status = "✓" if r < 0.3 else "✗"
        print(f"    {status} {g1} vs {g2}: r={r:.3f}")


# ══════════════════════════════════════════════════════════
# CHECK 4: UMAP colored by cell type — clusters should
#           separate similarly to scRNA-seq.
# ══════════════════════════════════════════════════════════
print("\n── Check 4: UMAP visualization ──")

sc.pl.umap(adata, color=["cell_type"], frameon=False, show=False,
           title="Cell types (protein-derived)")
plt.savefig(FIG_DIR / "check4_umap_celltypes.png", dpi=150, bbox_inches="tight")
plt.close()

# Also plot key markers on UMAP
key_markers = ["CD3E", "CD8A", "MS4A1", "CD68", "KRT18", "VIM"]
available = [m for m in key_markers if m in adata.var_names]
sc.pl.umap(adata, color=available, frameon=False, show=False, ncols=3)
plt.savefig(FIG_DIR / "check4_umap_key_markers.png", dpi=150, bbox_inches="tight")
plt.close()
print("  Saved: check4_umap_celltypes.png, check4_umap_key_markers.png")


# ══════════════════════════════════════════════════════════
# CHECK 5: Rank genes (DEG) per cell type — top markers
#           should match known biology.
# ══════════════════════════════════════════════════════════
print("\n── Check 5: DEG analysis per cell type ──")

sc.tl.rank_genes_groups(adata, groupby="cell_type", method="wilcoxon")
sc.pl.rank_genes_groups(adata, n_genes=10, show=False)
plt.savefig(FIG_DIR / "check5_deg_per_celltype.png", dpi=150, bbox_inches="tight")
plt.close()

# Extract top 5 markers per cell type
deg_results = []
for ct in adata.obs["cell_type"].unique():
    try:
        degs = sc.get.rank_genes_groups_df(adata, group=ct)
        top5 = degs.head(5)
        for _, row in top5.iterrows():
            deg_results.append({
                "cell_type": ct,
                "gene": row["names"],
                "score": round(row["scores"], 2),
                "logfoldchange": round(row["logfoldchanges"], 2),
                "pval_adj": f"{row['pvals_adj']:.2e}",
            })
    except Exception:
        pass

deg_df = pd.DataFrame(deg_results)
deg_df.to_csv(FIG_DIR / "check5_top_degs_per_celltype.csv", index=False)
print("  Saved: check5_deg_per_celltype.png")


# ══════════════════════════════════════════════════════════
# CHECK 6: Dynamic range comparison — protein vs typical RNA
# ══════════════════════════════════════════════════════════
print("\n── Check 6: Dynamic range summary ──")

gene_stats = pd.DataFrame({
    "gene": adata.var_names,
    "mean":   adata.X.mean(axis=0),
    "std":    adata.X.std(axis=0),
    "min":    adata.X.min(axis=0),
    "max":    adata.X.max(axis=0),
    "median": np.median(adata.X, axis=0),
    "pct_zero": (adata.X == 0).sum(axis=0) / adata.n_obs * 100,
})
gene_stats.to_csv(FIG_DIR / "check6_gene_stats.csv", index=False)

print(f"  Mean expression range:  {gene_stats['mean'].min():.3f} – {gene_stats['mean'].max():.3f}")
print(f"  Max value range:        {gene_stats['max'].min():.3f} – {gene_stats['max'].max():.3f}")
print(f"  Sparsity (% zeros):     {gene_stats['pct_zero'].min():.1f}% – {gene_stats['pct_zero'].max():.1f}%")
print("  NOTE: Protein data is much LESS sparse than scRNA-seq (typically 80-95% zeros)")
print("        Low sparsity is expected and correct for MIBI data.")


# ══════════════════════════════════════════════════════════
# SUMMARY
# ══════════════════════════════════════════════════════════
print("\n" + "="*60)
print("VALIDATION SUMMARY")
print("="*60)
n_pass = sum(1 for r in concordance_results if "PASS" in r["status"])
n_total = len(concordance_results)
print(f"  Marker concordance:     {n_pass}/{n_total} markers pass")
print(f"  Distribution plots:     see check1_per_gene_distributions.png")
print(f"  Correlation structure:  see check3_gene_correlation_matrix.png")
print(f"  UMAP separation:        see check4_umap_celltypes.png")
print(f"  DEG analysis:           see check5_deg_per_celltype.png")
print(f"\n  All outputs in: {FIG_DIR.resolve()}")
print("="*60)
