# MIBI Protein → Gene Expression Conversion: Complete Documentation
## Dayao et al. HNSCC MIBI-SoC Dataset (40 markers, 570K+ cells)

**Date prepared**: April 2026  
**Mapping source**: UniProt Swiss-Prot (reviewed human entries, organism_id:9606)

---

## Table of Contents

1. [Overview — What Was Done](#1-overview)
2. [Input Data Summary](#2-input-data)
3. [Protein-to-Gene Mapping (UniProt)](#3-mapping)
4. [Pipeline Steps in Detail](#4-pipeline-steps)
5. [Output Dataset Description](#5-output-dataset)
6. [Technology-Specific Issues (MIBI)](#6-mibi-issues)
7. [Transformation & Normalization Decisions](#7-normalization)
8. [Assumptions & Limitations](#8-assumptions)
9. [Validation Checks](#9-validation)
10. [Downstream Analyses](#10-downstream)
11. [File Inventory](#11-files)
12. [How to Run](#12-how-to-run)

---

## 1. Overview — What Was Done <a name="1-overview"></a>

**Goal**: Convert single-cell protein expression data (MIBI-SoC, 40 antibody markers) into a gene-expression-like dataset where columns are HUGO gene symbols instead of protein/antibody names.

**What was produced**: A new CSV dataset (`gene_expression_dataset.csv`) containing:
- ~570,000 rows (one per cell)
- 40 gene-symbol columns (expression values, arcsinh-normalized)
- Cell annotations: cell type, sample ID, spatial coordinates
- Clinical metadata: patient ID, HPV status, survival, recurrence

**Why**: To enable analysis with standard scRNA-seq computational tools (Scanpy, Seurat) and biological interpretation via gene symbols (DEG, gene set scoring, CellChat).

---

## 2. Input Data Summary <a name="2-input-data"></a>

| Property | Value |
|---|---|
| **Source publication** | Dayao et al., "Deriving spatial features from in situ proteomics imaging to enhance cancer survival analysis" |
| **Technology** | MIBI-SoC (Multiplexed Ion Beam Imaging – Staining on Chip) |
| **Cancer type** | Head & Neck Squamous Cell Carcinoma (HNSCC) |
| **Patients** | 7 (UPMC_c001 – UPMC_c007) |
| **Tissue regions** | 378 (passed QC) |
| **Total cells** | ~570,000+ |
| **Protein markers** | 40 antibodies |
| **Cell types** | 16 categories (assigned by kNN classification) |
| **Normalization** | arcsinh-transformed (pre-applied by authors) |
| **Survival data** | Yes (survival_status, survival_day) |

### Input files used:
| File | Description |
|---|---|
| `data/dataset_info/labeled_arcsinh_norm_data.csv` | Arcsinh-normalized protein expression per cell + cluster labels |
| `data/dataset_info/cell_locations_and_labels.csv` | Spatial (X,Y) coordinates, cell IDs, cell type labels |
| `data/dataset_info/sample_metadata.csv` | Clinical metadata per sample (patient, HPV, survival, recurrence) |
| `data/dataset_info/marker_names.csv` | List of 40 protein marker names |

---

## 3. Protein-to-Gene Mapping (UniProt) <a name="3-mapping"></a>

### 3.1 Mapping Method

Each of the 40 MIBI antibody markers was mapped to its corresponding HUGO gene symbol using **UniProt Swiss-Prot** (reviewed human entries only). The mapping was verified by querying:

```
https://rest.uniprot.org/uniprotkb/search?query=gene:{GENE}+AND+organism_id:9606+AND+reviewed:true
```

Each mapping records:
- **MIBI marker name** (as in the original data)
- **HUGO gene symbol** (primary gene name from UniProt)
- **UniProt accession** (Swiss-Prot ID)
- **UniProt protein name** (full name from database)

### 3.2 Complete Mapping Table

| MIBI Marker | Gene Symbol | UniProt ID | Output Column | UniProt Protein Name |
|---|---|---|---|---|
| CD31 | PECAM1 | A0A075B728 | PECAM1 | Platelet endothelial cell adhesion molecule |
| CD57 | B3GAT1 | Q9P2W7 | B3GAT1 | Galactosylgalactosylxylosylprotein 3-beta-glucuronosyltransferase 1 |
| CD4 | CD4 | P01730 | CD4 | T-cell surface glycoprotein CD4 |
| CD15 | FUT4 | P22083 | FUT4 | Alpha-(1,3)-fucosyltransferase 4 |
| FoxP3 | FOXP3 | Q9BZS1 | FOXP3 | Forkhead box protein P3 |
| CD16 | FCGR3A | P08637 | FCGR3A | Low affinity immunoglobulin gamma Fc region receptor III-A |
| CD20 | MS4A1 | P11836 | MS4A1 | B-lymphocyte antigen CD20 |
| CD45RO | PTPRC | P08575 | PTPRC_RO | Receptor-type tyrosine-protein phosphatase C (CD45RO isoform) |
| CD38 | CD38 | P28907 | CD38 | ADP-ribosyl cyclase/cyclic ADP-ribose hydrolase 1 |
| CD34 | CD34 | P28906 | CD34 | Hematopoietic progenitor cell antigen CD34 |
| CD11b | ITGAM | P11215 | ITGAM | Integrin alpha-M |
| CD68 | CD68 | P34810 | CD68 | Macrosialin |
| CD134 | TNFRSF4 | P43489 | TNFRSF4 | Tumor necrosis factor receptor superfamily member 4 |
| TMEM16A | ANO1 | Q5XXA6 | ANO1 | Anoctamin-1 |
| PanCK | KRT18 | P05783 | KRT18 | Keratin, type I cytoskeletal 18 (**representative for pan-CK cocktail**) |
| Podoplanin | PDPN | Q86YL7 | PDPN | Podoplanin |
| CD45 | PTPRC | P08575 | PTPRC | Receptor-type tyrosine-protein phosphatase C (total CD45) |
| GranzymeB | GZMB | P10144 | GZMB | Granzyme B |
| CD49f | ITGA6 | P23229 | ITGA6 | Integrin alpha-6 |
| CD11c | ITGAX | P20702 | ITGAX | Integrin alpha-X |
| CD47 | CD47 | Q08722 | CD47 | Leukocyte surface antigen CD47 |
| CD8 | CD8A | P01732 | CD8A | T-cell surface glycoprotein CD8 alpha chain |
| CD117 | KIT | P10721 | KIT | Mast/stem cell growth factor receptor Kit |
| Vimentin | VIM | P08670 | VIM | Vimentin |
| CD69 | CD69 | Q07108 | CD69 | Early activation antigen CD69 |
| aSMA | ACTA2 | P62736 | ACTA2 | Actin, aortic smooth muscle |
| CD14 | CD14 | P08571 | CD14 | Monocyte differentiation antigen CD14 |
| CD21 | CR2 | P20023 | CR2 | Complement receptor type 2 |
| HLA-DR | HLA-DRA | P01903 | HLA-DRA | HLA class II histocompatibility antigen, DR alpha chain |
| PDL1 | CD274 | Q9NZQ7 | CD274 | Programmed cell death 1 ligand 1 (PD-L1) |
| CD56 | NCAM1 | P13591 | NCAM1 | Neural cell adhesion molecule 1 |
| p16 | CDKN2A | P42771 | CDKN2A | Cyclin-dependent kinase inhibitor 2A |
| PD1 | PDCD1 | Q15116 | PDCD1 | Programmed cell death protein 1 |
| CD45RA | PTPRC | P08575 | PTPRC_RA | Receptor-type tyrosine-protein phosphatase C (CD45RA isoform) |
| ICOS | ICOS | Q9Y6W8 | ICOS | Inducible T-cell costimulator |
| CD152 | CTLA4 | P16410 | CTLA4 | Cytotoxic T-lymphocyte protein 4 |
| Ki67 | MKI67 | P46013 | MKI67 | Proliferation marker protein Ki-67 |
| CollagenIV | COL4A1 | P02462 | COL4A1 | Collagen alpha-1(IV) chain |
| CD3e | CD3E | P07766 | CD3E | T-cell surface glycoprotein CD3 epsilon chain |

### 3.3 Special Cases

#### PTPRC Isoforms (CD45 / CD45RA / CD45RO)
All three antibodies detect the **same gene** (*PTPRC*, P08575) but different splice variants:
- **CD45** (total) → column `PTPRC`
- **CD45RA** (naive T cells, exon A present) → column `PTPRC_RA`
- **CD45RO** (memory T cells, all variable exons spliced out) → column `PTPRC_RO`

**Decision**: Keep all three as separate columns with suffixes. Collapsing would lose the naive/memory T cell distinction, which is biologically important. For gene-set overlap, only `PTPRC` (total) is used.

#### PanCK (Pan-Cytokeratin)
The PanCK antibody cocktail detects multiple keratins (KRT5, KRT8, KRT14, KRT18, KRT19). It is impossible to determine which specific keratin is contributing signal.

**Decision**: Mapped to **KRT18** (P05783) as the most commonly referenced representative. Flagged as `is_ambiguous=True` in the mapping table.

#### HLA-DR
The HLA-DR antibody detects a heterodimer (HLA-DRA + HLA-DRB1). 

**Decision**: Mapped to **HLA-DRA** (P01903, alpha chain) as the primary chain. Flagged as ambiguous.

#### CD8
CD8 is a heterodimer (CD8A + CD8B).

**Decision**: Mapped to **CD8A** (P01732, alpha chain) as the more commonly referenced subunit.

---

## 4. Pipeline Steps in Detail <a name="4-pipeline-steps"></a>

### Step 0: Export UniProt mapping
- Saves `uniprot_protein_to_gene_mapping.csv` — the full provenance table
- This file documents exactly how each protein was mapped and to which UniProt entry

### Step 1: Load data
- Reads the arcsinh-normalized expression matrix (~570K cells × 43 columns)
- Reads spatial coordinates (~570K cells with X, Y positions)
- Reads clinical metadata (378 samples with survival, HPV status, etc.)

### Step 2: Separate expression from annotations
- Identifies which columns are protein markers (40 markers)
- Which columns are metadata (sample_id, cell_id, cluster, cluster_label)
- Renames protein columns → gene symbols using the UniProt mapping

### Step 3: Handle PTPRC isoforms
- CD45 → `PTPRC`, CD45RA → `PTPRC_RA`, CD45RO → `PTPRC_RO`
- All linked to same UniProt entry P08575 but kept as separate features

### Step 4: Build AnnData object
- Expression matrix as `adata.X`
- Gene metadata in `adata.var` (protein name, gene symbol, UniProt ID, ambiguity flags)
- Cell annotations in `adata.obs` (sample, cell type, cluster)
- Raw arcsinh values stored in `adata.layers["arcsinh_norm"]`

### Step 5: Add spatial coordinates
- Matches cells between expression and spatial files by `sample_id + cell_id` key
- Stored in `adata.obsm["spatial"]`

### Step 6: Merge clinical metadata
- Maps sample-level clinical data to each cell
- Columns: patient_id, tissue_type, primarysite, hpvstatus_new, survival_status, survival_day, recurred

### Step 7: Technology-specific preprocessing
- **No renormalization** — data is already arcsinh-transformed (see Section 7)
- **99.9th percentile clipping** per gene — removes MIBI hot-pixel artifacts
- **Min-max scaling** layer [0,1] for comparisons
- **Z-score** layer for heatmaps and gene-set scoring

### Step 8: Scanpy embedding
- PCA (up to 30 components on all 40 gene features)
- UMAP (2D, 30 neighbors)
- Leiden clustering (resolution 0.8)

### Step 9: Save outputs
- Primary: `gene_expression_dataset.csv` (new gene-expression dataset with metadata)
- AnnData: `mibi_as_gene_expression.h5ad`
- Cell metadata: `cell_metadata.csv` (includes UMAP coordinates)
- Pure expression: `gene_expression_matrix.csv` (genes only, no metadata)

---

## 5. Output Dataset Description <a name="5-output-dataset"></a>

### Primary output: `output/gene_expression_dataset.csv`

This is the **new gene-expression dataset** — the main deliverable.

| Column type | Columns | Description |
|---|---|---|
| **Index** | `cell_barcode` | Unique cell ID (`{sample_id}_{cell_id}`) |
| **Metadata** | `sample_id` | Tissue region ID (e.g., UPMC_c001_v001_r001_reg001) |
| | `cell_id` | Cell number within sample |
| | `cell_type` | One of 16 cell types (kNN-assigned) |
| | `cluster_id` | Numeric cluster (0–15) |
| | `spatial_x`, `spatial_y` | Spatial coordinates (pixels) |
| | `patient_id` | De-identified patient |
| | `tissue_type` | Primary tumor, etc. |
| | `primarysite` | Anatomical location |
| | `hpvstatus_new` | HPV status (0/1) |
| | `survival_status` | 0=alive, 1=dead |
| | `survival_day` | Days to event/censoring |
| | `recurred` | Recurrence flag (0/1) |
| **Gene expression** | `PECAM1`, `B3GAT1`, `CD4`, `FUT4`, `FOXP3`, `FCGR3A`, `MS4A1`, `PTPRC_RO`, `CD38`, `CD34`, `ITGAM`, `CD68`, `TNFRSF4`, `ANO1`, `KRT18`, `PDPN`, `PTPRC`, `GZMB`, `ITGA6`, `ITGAX`, `CD47`, `CD8A`, `KIT`, `VIM`, `CD69`, `ACTA2`, `CD14`, `CR2`, `HLA-DRA`, `CD274`, `NCAM1`, `CDKN2A`, `PDCD1`, `PTPRC_RA`, `ICOS`, `CTLA4`, `MKI67`, `COL4A1`, `CD3E` | Arcsinh-normalized expression values |

### How to use the new dataset
```python
import pandas as pd

df = pd.read_csv("output/gene_expression_dataset.csv", index_col="cell_barcode")

# Gene expression columns only
gene_cols = df.columns[df.columns.get_loc("PECAM1"):]  # everything after metadata
expr = df[gene_cols]

# Get expression for T cells only
t_cells = df[df["cell_type"].isin(["CD4 T cell", "CD8 T cell"])]
```

---

## 6. Technology-Specific Issues (MIBI) <a name="6-mibi-issues"></a>

### 6.1 Signal Spillover (Ion Beam Crosstalk)
MIBI uses metal-conjugated antibodies measured by mass spectrometry. Adjacent mass channels can bleed signal into each other.
- **Impact**: Low-level false positives for markers on adjacent channels
- **Mitigation**: Authors applied QC; we additionally clip at 99.9th percentile per gene

### 6.2 Segmentation Artifacts
Cell segmentation in tissue images is imperfect. In dense regions, signal from neighboring cells can leak across boundaries.
- **Impact**: Mixed cell-type signals; a tumor cell may show low-level immune markers
- **Mitigation**: The kNN cell classification with probability scores partially addresses this. Low-confidence cells (kNN.prob < 0.5) can be filtered if desired.

### 6.3 No Dropout / Zero Inflation
Unlike scRNA-seq, MIBI measures protein in every cell. There is **no dropout** or zero-inflation.
- **Impact**: ~0% zeros in the data (vs ~90% in scRNA-seq). Tools that model zero-inflation (ZINB, MAST, DCA) are **inappropriate** for this data.
- **Mitigation**: Use Wilcoxon rank-sum or t-tests for DEG. Do not apply scRNA-seq imputation methods (MAGIC, scVI denoising, etc.).

### 6.4 Batch Effects
7 patients, 378 regions — patient-level and slide-level batch effects exist.
- **Impact**: Technical variation may confound biological comparisons across patients
- **Mitigation**: If comparing across patients, consider Harmony, scVI, or ComBat batch correction

### 6.5 Panel Size Limitation
Only 40 proteins measured (vs ~20,000 genes in scRNA-seq).
- **Impact**: Cannot perform genome-wide analyses, pathway enrichment has very low coverage
- **Mitigation**: Use custom gene sets from the panel; do not attempt GO/KEGG enrichment

---

## 7. Transformation & Normalization Decisions <a name="7-normalization"></a>

### Why NO Further Normalization Was Applied

The input data (`labeled_arcsinh_norm_data.csv`) is already **arcsinh-normalized** by the original authors. This is the standard normalization for mass-spectrometry-based proteomics data:

```
arcsinh(x / cofactor)   where cofactor = 5 (typical for MIBI/CyTOF)
```

**Arcsinh is the protein-data equivalent of log1p for RNA-seq:**
- Both are variance-stabilizing transforms
- Both compress high values while preserving low-value resolution
- Both produce roughly log-normal–like distributions
- Scanpy/Seurat treat them identically for PCA, clustering, DEG

**What would go WRONG if we re-normalized:**
- `sc.pp.normalize_total()` → designed for count data; protein data doesn't have library-size variation the same way
- `sc.pp.log1p()` → double-transform (arcsinh already applied); would over-compress
- `scran` normalization → assumes dropout/UMI counts; would introduce artifacts

### Additional layers provided
| Layer | Transform | Use case |
|---|---|---|
| `arcsinh_norm` | Raw arcsinh values (untouched) | Ground truth, DEG, most analyses |
| `scaled_0_1` | Min-max to [0,1] per gene | Cross-gene comparisons, heatmaps |
| `zscore` | Z-score per gene | Gene-set scoring, standardized comparisons |

### Why NO Imputation Was Applied

MIBI data has near-zero sparsity. Imputation methods (MAGIC, scVI, SAVER) were designed to recover dropout in sparse scRNA-seq data. Applying them to MIBI data would:
- Smooth away real biological variation between cells
- Introduce artificial correlations between genes
- Violate the methods' statistical assumptions (they assume dropout noise)

---

## 8. Assumptions & Limitations <a name="8-assumptions"></a>

### What This Conversion CAN Do
| Capability | Status | Notes |
|---|---|---|
| Differential expression (DEG) | ✅ Works | Wilcoxon on arcsinh values |
| Cell clustering & embedding | ✅ Works | PCA + UMAP + Leiden |
| Gene set scoring (custom sets) | ✅ Works | 8 curated gene sets from panel |
| Cell-cell communication (CellChat) | ⚠️ Limited | Only LR pairs within 40 genes |
| Spatial gene expression analyses | ✅ Works | Coordinates preserved |
| Survival modeling | ✅ Works | Clinical metadata attached |

### What This Conversion CANNOT Do
| Limitation | Reason |
|---|---|
| Genome-wide GSEA (GO/KEGG/Reactome) | Only 40 genes — far too few for ranked enrichment |
| Trajectory / RNA velocity | No spliced/unspliced information |
| Gene regulatory network inference | 40 genes cannot capture TF-target relationships |
| Imputation of unmeasured genes | No principled way to predict ~20,000 genes from 40 |
| Direct quantitative RNA-protein comparison | Different scales, different noise models |

### Key Assumptions

1. **Protein → gene is 1:1** — Mostly true for CD markers. Exceptions flagged (PanCK, HLA-DR, PTPRC isoforms).

2. **Protein abundance reflects gene activity** — The correlation between protein and mRNA is ~0.4–0.6 across published studies. Post-translational regulation, protein half-life, and secretion can decouple them. Results reflect **protein-level biology**, not transcription.

3. **arcsinh ≈ log1p for computational tools** — Both are variance-stabilizing transforms. Scanpy, Seurat, and most DE methods treat them identically. This is mathematically sound.

4. **40 genes are sufficient for the intended analyses** — For DEG between cell types and custom gene-set scoring, yes. For genome-wide analyses, no.

---

## 9. Validation Checks <a name="9-validation"></a>

Script `02_validation_checks.py` runs 6 automated checks:

### Check 1: Per-gene distribution shape
- Histograms for all 40 genes
- Protein data should show unimodal or bimodal distributions (NOT uniform)
- arcsinh-transformed values should resemble log-normal (similar to log1p RNA)

### Check 2: Marker–cell type concordance
- For known marker–cell type pairs (e.g., CD3E should be high in T cells), tests whether:
  - Mean expression is higher in expected cell types (fold change > 1.2)
  - Difference is statistically significant (Mann-Whitney p < 0.01)
- Expected results: CD3E high in T cells, MS4A1 high in B cells, CD68 high in macrophages, etc.

### Check 3: Gene-gene correlation structure
- Correlation matrix of all 40 genes
- Expected positive correlations: CD3E-CD4, CD3E-CD8A (T cell co-expression)
- Expected low/negative correlations: CD3E-KRT18 (immune vs tumor)

### Check 4: UMAP visualization
- Cell types should form visually distinct clusters
- Key markers should show expected spatial patterns on UMAP

### Check 5: DEG per cell type
- Wilcoxon rank-sum test: each cell type vs all others
- Top markers for each cell type should match known biology

### Check 6: Dynamic range
- Summary statistics per gene (mean, std, min, max, % zeros)
- Expected: continuous values in [0, ~3], NOT binary, near-zero sparsity

---

## 10. Downstream Analyses <a name="10-downstream"></a>

### 10.1 DEG & Gene Set Scoring (`03_deg_and_enrichment.py`)
- DEGs per cell type (Wilcoxon, all 40 genes)
- DEGs for HPV+ vs HPV- in tumor cells
- DEGs for recurred vs non-recurred
- 8 custom gene sets scored per cell:
  - T cell activation, Immune checkpoint, Antigen presentation, Tumor markers
  - Immune infiltration, Regulatory T cell, Cytotoxicity, Stromal

### 10.2 CellChat (`04_cellchat_analysis.R`)
- Cell-cell communication inference using CellChatDB.human
- Limited to LR pairs where both ligand and receptor are among the 40 genes
- Outputs: interaction network, heatmap, signaling roles, LR pair table

---

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

## Citation Note

When publishing results from this conversion, state:

> "Single-cell protein expression data from MIBI-SoC imaging (40-marker panel, Dayao et al.)
> was mapped to HUGO gene symbols using UniProt Swiss-Prot (reviewed human entries) and
> analyzed using standard single-cell RNA-seq computational frameworks. Expression values
> represent arcsinh-normalized protein abundance, not mRNA. Gene set analyses use curated
> panel-specific gene sets. PTPRC splice variants (CD45/CD45RA/CD45RO) were retained as
> separate features. Pan-cytokeratin was mapped to KRT18 as a representative gene."
