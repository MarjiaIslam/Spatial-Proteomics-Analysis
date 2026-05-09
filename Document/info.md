
## 11. File Inventory <a name="11-files"></a>

### Code files
| File | Language | Purpose |
|---|---|---|
| `protein_to_gene_map.py` | Python | UniProt-verified mapping (40 markers → gene symbols + accessions) |
| `01_protein_to_gene_pipeline.py` | Python | Main conversion pipeline: protein data → gene-expression dataset |
| `02_validation_checks.py` | Python | 6 automated validation checks with diagnostic plots |
| `03_deg_and_enrichment.py` | Python | DEG analysis and custom gene-set scoring |
| `04_cellchat_analysis.R` | R | CellChat cell-cell communication |

### Output files (after running pipeline)
| File | Description |
|---|---|
| `output/gene_expression_dataset.csv` | **PRIMARY OUTPUT** — new gene-expression dataset |
| `output/uniprot_protein_to_gene_mapping.csv` | Full mapping provenance table |
| `output/mibi_as_gene_expression.h5ad` | AnnData object for Scanpy |
| `output/gene_expression_matrix.csv` | Pure expression matrix (no metadata) |
| `output/cell_metadata.csv` | Cell annotations + UMAP + spatial coords |
| `output/validation_figures/` | Diagnostic plots from validation checks |
| `output/deg_gsea_results/` | DEG tables and gene-set score plots |

---

## 12. How to Run <a name="12-how-to-run"></a>

### Prerequisites
```bash
pip install numpy pandas scanpy anndata scikit-learn matplotlib scipy
```

### Execution order
```bash
# 1. Verify mapping (optional, prints table)
python protein_to_gene_map.py

# 2. Run main conversion pipeline → produces gene_expression_dataset.csv
python 01_protein_to_gene_pipeline.py

# 3. Run validation checks → produces diagnostic figures
python 02_validation_checks.py

# 4. Run DEG / gene set analysis (optional)
python 03_deg_and_enrichment.py

# 5. Run CellChat in R (optional)
Rscript 04_cellchat_analysis.R
```

### Quick verification after Step 2
```python
import pandas as pd
df = pd.read_csv("output/gene_expression_dataset.csv", index_col="cell_barcode", nrows=5)
print(df.columns.tolist())  # Should show metadata + gene names
print(df.shape)              # Should be (5, ~52 columns)
```

---
'01_protein_to_gene_pipeline.py' is the main script that performs the conversion from protein expression data to a gene-expression-like dataset. It follows these steps:
```
Step 0  ─ Export UniProt mapping to CSV (provenance)
Step 1  ─ Load expression matrix, spatial coords, clinical metadata
Step 2  ─ Identify protein columns, rename them → gene symbols
Step 3  ─ Handle PTPRC isoforms (CD45/RA/RO → separate columns)
Step 4  ─ Build AnnData (X = expression, var = gene metadata w/ UniProt IDs)
Step 5  ─ Match & attach spatial X,Y coordinates to each cell
Step 6  ─ Join clinical metadata (HPV, survival, recurrence) per cell
Step 7  ─ Clip 99.9th percentile outliers; add scaled_0_1 and zscore layers
Step 8  ─ Run PCA → UMAP → Leiden clustering
Step 9  ─ Save 4 outputs:
            ★ gene_expression_dataset.csv  ← primary new dataset
              uniprot_protein_to_gene_mapping.csv
              mibi_as_gene_expression.h5ad
              cell_metadata.csv
````
02_validation_checks.py performs 6 validation checks to ensure the conversion is biologically meaningful and technically sound. It generates diagnostic plots for each check, saved in `output/validation_figures/`. The checks include:
```
Check 1 ─ Histogram per gene → should be unimodal/bimodal, not flat
Check 2 ─ Marker-cell concordance → CD3E should be high in T cells, etc.
           Reports fold-change + Mann-Whitney p-value per known marker
Check 3 ─ Gene-gene correlation matrix → CD3E-CD4 should be positive,
           CD3E-KRT18 should be low/negative
Check 4 ─ UMAP colored by cell type → cell types should form distinct islands
Check 5 ─ Wilcoxon DEGs per cell type → top genes should match known biology
Check 6 ─ Dynamic range stats → continuous [0,~3], near-zero sparsity
```

03_deg_and_enrichment.py performs differential expression analysis using Wilcoxon rank-sum tests to identify marker genes for each cell type. It also computes custom gene-set scores (e.g., T cell activation) based on the 40 measured genes, and generates plots of these scores across cell types.

```Section 1 ─ DEG: every cell type vs all others (Wilcoxon)
Section 2 ─ DEG: HPV+ vs HPV− in tumor cells only
Section 3 ─ Gene set scoring: 8 curated sets (T cell activation, checkpoint,
             cytotoxicity, stromal, etc.) scored per cell → UMAP + violin plots
Section 4 ─ DEG: Recurred vs not recurred
Section 5 ─ Re-saves .h5ad with gene set scores attached to cells
```

04_cellchat_analysis.R uses the CellChat R package to analyze cell-cell communication based on the converted gene-expression dataset. It identifies significant ligand-receptor interactions between cell types and visualizes them with circle plots and heatmaps.

```
Step 1 ─ Load CSVs, transpose to genes × cells format CellChat expects
Step 2 ─ Create CellChat object; check overlap of our 39 genes with LR database
Step 3 ─ Identify over-expressed genes/interactions per cell type
Step 4 ─ Compute communication probabilities (permutation test)
Step 5 ─ Export: network plot, heatmap, signaling roles, LR pair table (.csv)
```

