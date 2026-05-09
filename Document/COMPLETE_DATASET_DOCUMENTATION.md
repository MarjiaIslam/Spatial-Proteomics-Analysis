# Complete Dataset Documentation
## Dayao et al. Head & Neck Cancer MIBI-SoC Analysis

**Dataset**: UPMC Head & Neck Cancer Cohort  
**Patients**: 7 (UPMC_c001 to UPMC_c007)  
**Total Samples**: 378 tissue regions  
**Cell Types**: 16 categories  
**Biomarkers**: 40 immunohistochemistry markers  
**Survival Data**: Yes (primary outcome)

---

## Table of Contents

1. [Core Input Files](#core-input-files) - `dataset_info.tar.gz`
2. [Biomarker Expression](#biomarker-expression) - `biomarker_expr_summary.tar.gz`
3. [Neighborhood Matrices](#neighborhood-matrices) - `all_k_neighborhood_mats.tar.gz`
4. [Ripley K Functions](#ripley-k-functions) - `k_fns_norm_by_uw_qc_labeled.csv`
5. [Comparison Methods](#comparison-methods) - `patwa_comparisons.tar.gz`, `denvar.tar.gz`
6. [Output Results](#output-results) - `rsf_risk_scores.tar.gz`

---

# 1. CORE INPUT FILES

## Archive: `dataset_info.tar.gz`

**Contains**: 5 core data files with metadata, markers, and cell-level expression

### A. `marker_names.csv`

**Purpose**: Reference list of all biomarkers measured  
**Format**: Plain text, one marker per line  
**Total**: 40 markers

**Content**:
```
CD31
CD57
CD4
CD15
FoxP3
CD16
CD20
CD45RO
CD38
CD34
CD11b
CD68
CD134
TMEM16A
PanCK
Podoplanin
CD45
GranzymeB
CD49f
CD11c
CD47
CD8
CD117
Vimentin
CD69
aSMA
CD14
CD21
HLA-DR
PDL1
CD56
p16
PD1
CD45RA
ICOS
CD152
Ki67
CollagenIV
CD3e
```

**Grouped by Function**:
| Category | Markers |
|----------|---------|
| T Cell | CD3e, CD4, CD8, CD45, CD45RA, CD45RO, CD69, FoxP3, GranzymeB, ICOS, CD152 |
| B Cell | CD20, CD21, CD38, CD134 |
| Myeloid | CD11b, CD11c, CD14, CD15, CD16, CD68, HLA-DR |
| Checkpoint | PD1, PDL1, CD56, CD57, CD117, p16, Ki67, ICOS, CD152 |
| Structural | CollagenIV, aSMA, Podoplanin, Vimentin, PanCK, CD31, CD34, CD49f, TMEM16A |

---

### B. `qc_acq_ids_labeled.csv`

**Purpose**: List of all samples that passed QC  
**Format**: Plain text, one sample ID per line  
**Total**: 378 samples

**Content** (sample):
```
UPMC_c001_v001_r001_reg001
UPMC_c001_v001_r001_reg002
UPMC_c001_v001_r001_reg004
UPMC_c001_v001_r001_reg005
...
UPMC_c007_v001_r001_reg063
```

**Format Breakdown**: `UPMC_c{patient}_v{visit}_r{replicate}_reg{region}`
- **Patient**: c001–c007 (7 patients)
- **Visit**: v001 (single visit)
- **Replicate**: r001 (single replicate)
- **Region**: reg001–reg066+ (multiple tissue regions per patient)

**Distribution**:
| Patient | Regions | Patient | Regions |
|---------|---------|---------|---------|
| c001 | 50 | c005 | 64 |
| c002 | 36 | c006 | 66 |
| c003 | 34 | c007 | 63 |
| c004 | 65 | TOTAL | 378 |

---

### C. `sample_metadata.csv`

**Purpose**: Clinical metadata and survival outcomes for each sample  
**Format**: CSV with 21 columns  
**Rows**: 378 (one per sample)

**Columns**:
| Column | Description | Example |
|--------|-------------|---------|
| `acquisition_id` | Sample ID | UPMC_c001_v001_r001_reg001 |
| `patient_id` | De-identified patient | 2836 |
| `coverslip_label` | Slide location | L4a2 |
| `tissue_type` | Tissue classification | Primary tumor |
| `tissue_subtype` | Subtype | (empty for most) |
| `primarysite` | Anatomical location | Larynx, Supraglottic, Pharynx |
| `larynxpharynx` | Anatomical region | Larynx |
| `smokinghx` | Smoking history | 0–2 (encoded) |
| `alcoholhx` | Alcohol history | 0–5 (encoded) |
| `hpvstatus` | HPV status | 0–5 (encoded) |
| `hpvstatus_new` | Refined HPV status | 0–1 (binary) |
| `isrecurrencetumor` | Recurrence flag | 0 or 1 |
| `recurredlocal` | Local recurrence | 0–1 |
| `recurredregional` | Regional recurrence | 0–1 |
| `recurreddistantmet` | Distant metastasis | 0–1 |
| `pathorclinstage` | Tumor stage | Path, Clinical |
| `status` | Vital status code | DOC, DOD, NED |
| `recurred` | Any recurrence | 0–1 |
| `survival_status` | Death indicator | 0 (alive), 1 (dead) |
| `survival_day` | Days to death/censoring | Integer (1–2500+) |

**Example Data** (3 rows):

| acquisition_id | patient_id | coverslip_label | tissue_type | primarysite | larynxpharynx | smokinghx | hpvstatus | pathorclinstage | status | survival_status | survival_day | recurred | recurredlocal | recurredregional | recurreddistantmet |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| UPMC_c001_v001_r001_reg001 | 2836 | L4a2 | Primary tumor | Supraglottic | Larynx | 2.0 | 5.0 | Path | DOC | 1.0 | 1001.0 | 0.0 | 0.0 | 0.0 | 0.0 |
| UPMC_c001_v001_r001_reg002 | 2830 | L4a2 | Primary tumor | Larynx | Larynx | 2.0 | 5.0 | Path | DOD | 1.0 | 831.0 | 1.0 | 1.0 | 1.0 | 0.0 |
| UPMC_c001_v001_r001_reg004 | 2823 | L4a2 | Primary tumor | Pyriform Sinus | Larynx | 1.0 | 0.0 | Clinical | DOD | 1.0 | 61.0 | 0.0 | 0.0 | 0.0 | 0.0 |

---

### D. `cell_locations_and_labels.csv`

**Purpose**: Spatial coordinates and cell type classification for every cell  
**Format**: CSV with 10 columns, **570,000+ rows** (one per cell)  
**Size**: >50 MB

**Columns**:
| Column | Type | Description |
|--------|------|-------------|
| `REGION` | int | Tissue region number |
| `TILE_NUM` | int | Tile within region |
| `CELL_ID` | int | Unique cell ID within sample |
| `ACQUISITION_ID` | str | Sample ID (e.g., UPMC_c001_v001_r001_reg001) |
| `X` | float | X-coordinate (pixels) |
| `Y` | float | Y-coordinate (pixels) |
| `SIZE` | int | Cell size (pixels) |
| `CLUSTER_ID` | int | Cell type numerical ID (0–15) |
| `CLUSTER_LABEL` | str | Cell type name (see below) |
| `kNN.prob` | float | k-NN classification probability (0–1) |

**16 Cell Types** (CLUSTER_LABEL values):
| ID | Cell Type | Biological Category |
|----|-----------|------------------|
| 0 | APC | Antigen-Presenting Cell |
| 1 | B cell | B Lymphocyte |
| 2 | CD4 T cell | T Lymphocyte |
| 3 | CD8 T cell | T Lymphocyte |
| 4 | Granulocyte | Myeloid |
| 5 | Lymph vessel | Vascular |
| 6 | Macrophage | Myeloid |
| 7 | Naive immune cell | Immune (unspecified) |
| 8 | Stromal / Fibroblast | Stromal |
| 9 | Tumor | Malignant |
| 10 | Tumor (CD15+) | Malignant (CD15+ marker) |
| 11 | Tumor (CD20+) | Malignant (CD20+ marker) |
| 12 | Tumor (CD21+) | Malignant (CD21+ marker) |
| 13 | Tumor (Ki67+) | Malignant (Ki67+ proliferating) |
| 14 | Tumor (Podo+) | Malignant (Podoplanin+ marker) |
| 15 | Vessel | Blood vessel |

**Example Data** (5 rows):

| kNN.prob | REGION | TILE_NUM | CELL_ID | ACQUISITION_ID | X | Y | SIZE | CLUSTER_ID | CLUSTER_LABEL |
|---|---|---|---|---|---|---|---|---|---|
| 0.65 | 1 | 1 | 1 | UPMC_c001_v001_r001_reg001 | 1618 | 3 | 87 | 14 | Tumor (Podo+) |
| 1.0 | 1 | 1 | 2 | UPMC_c001_v001_r001_reg001 | 2181 | 3 | 137 | 11 | Tumor (CD20+) |
| 0.75 | 1 | 1 | 3 | UPMC_c001_v001_r001_reg001 | 590 | 4 | 194 | 2 | CD4 T cell |
| 1.0 | 1 | 1 | 4 | UPMC_c001_v001_r001_reg001 | 618 | 4 | 196 | 2 | CD4 T cell |
| 0.82 | 1 | 1 | 5 | UPMC_c001_v001_r001_reg001 | 850 | 12 | 156 | 6 | Macrophage |

**Usage**:
- Extract spatial coordinates: (X, Y)
- Get cell type: CLUSTER_LABEL
- Link to expression: Use CELL_ID + ACQUISITION_ID

---

### E. `labeled_arcsinh_norm_data.pkl`

**Purpose**: Cell-level normalized biomarker expression data  
**Format**: Python pickle (pandas DataFrame)  
**Size**: Large (~500 MB+), 570,000+ rows × 43 columns

**Columns** (40 biomarkers + metadata):
```
CD31, CD57, CD4, CD15, FoxP3, CD16, CD20, CD45RO, CD38, CD34,
CD11b, CD68, CD134, TMEM16A, PanCK, Podoplanin, CD45, GranzymeB,
CD49f, CD11c, CD47, CD8, CD117, Vimentin, CD69, aSMA, CD14, CD21,
HLA-DR, PDL1, CD56, p16, PD1, CD45RA, ICOS, CD152, Ki67, CollagenIV, CD3e,
+ metadata: sample_id, cell_id, cluster_id, cluster_label
```

**Data Format**: Arcsinh-normalized (inverse hyperbolic sine) expression values

**Example Structure** (conceptual):
```
   sample_id                  cell_id  CD31   CD57   CD4   ...  cluster_label
0  UPMC_c001_v001_r001_reg001      1   1.23  0.56  2.34  ...  CD4 T cell
1  UPMC_c001_v001_r001_reg001      2   0.08  0.02  0.10  ...  Tumor
2  UPMC_c001_v001_r001_reg001      3   3.45  2.01  0.05  ...  Macrophage
...
570000 UPMC_c007_v001_r001_reg063  12000  1.87  0.34  1.92  ...  Vessel
```

**Loading in Python**:
```python
import pickle
import pandas as pd

with open('labeled_arcsinh_norm_data.pkl', 'rb') as f:
    df = pickle.load(f)

# Get all cells from one sample
sample_data = df[df['sample_id'] == 'UPMC_c001_v001_r001_reg001']
print(sample_data.shape)  # (5000, 43) - 5000 cells in this sample

# Get CD4 expression for all cells
cd4_expr = df['CD4'].values  # Array of 570,000 values

# Get cells of type CD4 T cell
cd4_cells = df[df['cluster_label'] == 'CD4 T cell']
print(cd4_cells.shape)  # Number of CD4 T cells
```

---

# 2. BIOMARKER EXPRESSION

## Archive: `biomarker_expr_summary.tar.gz`

**Creates**: `expression_biomarkers/` directory  
**Contains**: 378 CSV files (one per tissue sample)  
**Purpose**: Sample-level summary of biomarker expression

### Files Naming
```
expression_biomarkers/
├── UPMC_c001_v001_r001_reg001_cell_info.csv
├── UPMC_c001_v001_r001_reg002_cell_info.csv
└── ...
└── UPMC_c007_v001_r001_reg063_cell_info.csv
```

### File Structure

**Format**: 2 rows × 40 columns  
**Rows**: `summed`, `average`

**Example: `UPMC_c001_v001_r001_reg001_cell_info.csv`**

| | CD31 | CD57 | CD4 | CD15 | FoxP3 | CD16 | CD20 | CD45RO | CD38 | CD34 | ... | CD3e |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| summed | 1081720315 | 727836831 | 1397186671 | 3321760577 | 722255116 | 758147326 | ... | 561456885 |
| average | 199.62 | 134.31 | 257.83 | 612.98 | 133.28 | 139.91 | ... | 103.61 |

**Data Interpretation**:
- **`summed` row**: Total pixel intensity across all cells in the sample for each marker
- **`average` row**: Mean intensity per pixel (`summed / number_of_pixels`)

**Real Value Example**:
- CD31 summed: 1,081,720,315 (endothelial marker, cumulative)
- CD31 average: 199.62 (mean per pixel)

**Usage in Feature Engineering**:
```python
import pandas as pd

# Read sample expression
expr = pd.read_csv('expression_biomarkers/UPMC_c001_v001_r001_reg001_cell_info.csv', 
                   index_col=0)

# Get average row (used for Feature 1: Biomarker Region)
average_values = expr.loc['average']  # CD31=199.62, CD57=134.31, ...

# Apply arcsinh normalization
import numpy as np
arcsinh_values = np.arcsinh(average_values)  # Feature matrix input
```

---

# 3. NEIGHBORHOOD MATRICES

## Archive: `all_k_neighborhood_mats.tar.gz`

**Creates**: 20 directories (`k1/` to `k20/`)  
**Contains**: 378 numpy arrays per directory (7,560+ files total)  
**Purpose**: Cell-type spatial neighbor relationships

### Files Organization
```
k1/          ← k=1 nearest neighbor
├── UPMC_c001_v001_r001_reg001.npy
├── UPMC_c001_v001_r001_reg002.npy
└── ...

k2/          ← k=2 nearest neighbors
└── ...

...

k20/         ← k=20 nearest neighbors
├── UPMC_c001_v001_r001_reg001.npy
├── UPMC_c001_v001_r001_reg002.npy
└── ...
```

### File Format

**Data Type**: NumPy array (.npy binary format)  
**Shape**: (16, 16) - cell type interaction matrix  
**Content**: Count of neighbor relationships

**Raw Data Example** (conceptual 16×16 matrix for sample with k=10):

| Cell Type | Tumor | CD4 T | CD8 T | Macro | B cell | Stromal | ... |
|---|---|---|---|---|---|---|---|
| Tumor | 7350 | 2100 | 1400 | 2625 | 875 | 1200 | ... |
| CD4 T cell | 600 | 1500 | 1080 | 480 | 720 | 1100 | ... |
| CD8 T cell | 320 | 600 | 1200 | 400 | 240 | 900 | ... |
| Macrophage | 1500 | 750 | 900 | 2625 | 600 | 1300 | ... |
| B cell | 500 | 800 | 500 | 700 | 1400 | 600 | ... |
| Stromal | 1200 | 1100 | 900 | 1300 | 600 | 2000 | ... |
| ... | ... | ... | ... | ... | ... | ... | ... |

**Interpretation**: 
- Row = cell type A
- Column = cell type B  
- Value = how many times cell type A has a cell type B neighbor (in k nearest neighbors)

**Normalization** (Row normalization):
```
Each row sums to k × (number of cells of that type)
Example: if Tumor has 1750 cells, row sum = 1750 × 10 = 17,500
```

**Loading and Using**:
```python
import numpy as np

# Load k=10 neighborhood matrix for a sample
neigh_mat = np.load('k10/UPMC_c001_v001_r001_reg001.npy')
print(neigh_mat.shape)  # (16, 16)

# Get proportions of neighbor types for Tumor (row 9)
tumor_neighbors = neigh_mat[9] / neigh_mat[9].sum()
# Result: [0.42, 0.12, 0.08, 0.15, ...]
# Means: 42% of Tumor neighbors are Tumors, 12% are CD4T, etc.
```

**Feature 4 Usage** (Neighborhood Matrix Feature):
```python
# Flatten to 256 features for model input
features = neigh_mat.flatten()  # (256,)
# This becomes one row in the feature matrix
```

---

# 4. RIPLEY K FUNCTIONS

## File: `k_fns_norm_by_uw_qc_labeled.csv`

**Purpose**: Spatial clustering patterns for each biomarker  
**Format**: CSV  
**Size**: >50 MB  
**Rows**: 307 (one per sample)  
**Columns**: 80 (40 markers × 2 radii) + acquisition_id

### What is Ripley K?

**Definition**: Measures spatial clustering of cells with high biomarker expression

- **K(r) > 1**: Cells are **clustered** (more neighbors than random)
- **K(r) = 1**: Cells are **randomly distributed**
- **K(r) < 1**: Cells are **dispersed** (fewer neighbors than random)

### File Structure

**Radii Used**: r=30 pixels and r=80 pixels  
**Markers**: All 40 biomarkers  
**Total Features**: 40 × 2 = 80

**Column Names** (pattern):
```
CD31_r30, CD31_r80, CD57_r30, CD57_r80, CD4_r30, CD4_r80, ..., CD3e_r30, CD3e_r80, acquisition_id
```

**Example Data**:

| acquisition_id | CD31_r30 | CD31_r80 | CD57_r30 | CD57_r80 | CD4_r30 | CD4_r80 | ... | CD3e_r30 | CD3e_r80 |
|---|---|---|---|---|---|---|---|---|---|
| UPMC_c001_v001_r001_reg001 | 1.234 | 2.567 | 0.876 | 1.234 | 1.456 | 2.890 | ... | 0.987 | 1.456 |
| UPMC_c001_v001_r001_reg002 | 0.987 | 1.876 | 1.123 | 1.567 | 1.234 | 1.890 | ... | 1.123 | 1.987 |
| UPMC_c001_v001_r001_reg004 | 1.456 | 2.890 | 0.654 | 0.987 | 0.765 | 1.234 | ... | 0.876 | 1.345 |

**Interpretation Examples**:
- CD31_r30 = 1.234 → Endothelial cells are clustered at 30 pixel radius
- CD57_r30 = 0.876 → NK cells (CD57+) are dispersed at 30 pixel radius
- CD4_r80 = 2.890 → CD4 T cells are highly clustered at 80 pixel radius

### Usage
```python
import pandas as pd

# Load K-function features
kfns = pd.read_csv('k_fns_norm_by_uw_qc_labeled.csv', index_col='acquisition_id')

# Get clustering pattern for CD4 at 30-pixel radius
cd4_clustering_r30 = kfns['CD4_r30']

# Merge with survival data for modeling
from base_paper__dataset_details import Feature5_ripley
# This becomes Feature 5 in the model
```

---

# 5. COMPARISON METHODS

## A. Patwa et al. Comparison

### Archive: `patwa_comparisons.tar.gz`

**Creates**: `comparisons/` directory  
**Contains**: 3 CSV files  
**Purpose**: Biomarker-based features for comparison model

### File 1: `binary_biomarker_expr_cells_qc_labeled.csv`

**Purpose**: Binary classification of biomarker expression per cell  
**Format**: CSV, >50 MB  
**Rows**: 570,000+ (one per cell)  
**Columns**: 40 (one per biomarker) + cell/sample metadata

**Content**: Binary (0/1) values
- **1** = cell expression above threshold for that marker
- **0** = cell expression below threshold

**Example**:

| cell_id | acquisition_id | CD31 | CD57 | CD4 | CD15 | ... | CD3e |
|---|---|---|---|---|---|---|---|
| 1 | UPMC_c001_v001_r001_reg001 | 0 | 0 | 1 | 0 | ... | 1 |
| 2 | UPMC_c001_v001_r001_reg001 | 1 | 1 | 0 | 1 | ... | 0 |
| 3 | UPMC_c001_v001_r001_reg001 | 0 | 0 | 1 | 0 | ... | 0 |

---

### File 2: `biomarker_frac_positivity_qc_labeled.csv`

**Purpose**: Fraction of cells positive for each marker per sample  
**Format**: CSV, >50 MB  
**Rows**: 378 (one per sample)  
**Columns**: 40 (one per marker) + acquisition_id

**Content**: Proportions (0–1)
- Value = Number of CD31+ cells / Total cells in sample

**Example**:

| acquisition_id | CD31 | CD57 | CD4 | CD15 | ... | CD3e |
|---|---|---|---|---|---|---|
| UPMC_c001_v001_r001_reg001 | 0.45 | 0.12 | 0.38 | 0.25 | ... | 0.67 |
| UPMC_c001_v001_r001_reg002 | 0.32 | 0.18 | 0.42 | 0.19 | ... | 0.71 |

**Interpretation**:
- CD31: 45% of cells positive (endothelial/vessel cells)
- CD57: 12% of cells positive (NK cells)
- CD4: 38% of cells positive (CD4 T cells)

---

### File 3: `interaction_biomarker_features.csv`

**Purpose**: Pairwise biomarker co-expression counts  
**Format**: CSV  
**Rows**: 378 (one per sample)  
**Columns**: ~800 (all pairwise combinations) + acquisition_id

**Column Naming**: `{marker1}-{marker2}`

**Example Columns**:
```
CD31-CD31, CD31-CD57, CD31-CD4, CD31-CD15, ...,
CD57-CD57, CD57-CD4, CD57-CD15, ...,
CD4-CD4, CD4-CD15, ...,
...
CD3e-CD3e, acquisition_id
```

**Example Data**:

| acquisition_id | CD31-CD31 | CD31-CD57 | CD31-CD4 | ... | CD3e-CD3e |
|---|---|---|---|---|---|
| UPMC_c001_v001_r001_reg001 | 4038.0 | 5656.0 | 9306.0 | ... | 561456885 |
| UPMC_c001_v001_r001_reg002 | 8856.0 | 6929.0 | 18448.0 | ... | 735628901 |

**Interpretation**:
- CD31-CD31: Count of CD31+ cells co-expressing CD31 (self)
- CD31-CD57: Count of cells co-expressing CD31 and CD57
- CD4-CD8: Count of cells co-expressing CD4 and CD8 (rare but possible)

---

## B. DenVar Comparison

### Archive: `denvar.tar.gz`

**Creates**: `denvar/` directory  
**Contains**: 80 files (40 markers × 2 file types) + 3 special files

### File Type 1: `{marker}_DenVar_clusters.csv`

**Purpose**: Density variance cluster assignments per sample  
**Example**: `CD31_DenVar_clusters.csv`

**Format**: CSV  
**Rows**: 378 (one per sample)  
**Columns**: Multiple cluster assignments

**Content**: 
- For each sample, DenVar algorithm assigns cells to clusters
- Result: Which cluster (0, 1, 2, ...) each cell belongs to for this marker

**Example**:

| sample_id | cell_cluster_1 | cell_cluster_2 | cell_cluster_3 | ... |
|---|---|---|---|---|
| UPMC_c001_v001_r001_reg001 | 3 | 1 | 3 | ... |
| UPMC_c001_v001_r001_reg002 | 1 | 2 | 1 | ... |

---

### File Type 2: `{marker}_JSD.RDS`

**Purpose**: Jensen-Shannon Distance matrix between samples  
**Format**: R object (.RDS)  
**Size**: 378 × 378 symmetric matrix

**Content**: 
- Distance values (0–1) between each pair of samples
- Based on distributional similarity of marker expression
- Calculated using Jensen-Shannon divergence

**Example Structure** (conceptual):

| | reg001 | reg002 | reg003 | ... |
|---|---|---|---|---|
| reg001 | 0.0 | 0.234 | 0.198 | ... |
| reg002 | 0.234 | 0.0 | 0.156 | ... |
| reg003 | 0.198 | 0.156 | 0.0 | ... |
| ... | ... | ... | ... | ... |

---

### Special DenVar Files

**1. `Marker_data_TNBC_MIBI.csv`**
- Master marker data in TNBC context
- Raw or processed marker measurements

**2. `qc_norm_marker_data_DenVar.RDS`**
- Quality-controlled normalized marker data
- R format for downstream analysis

**3. `qc_norm_marker_data_range01_DenVar.RDS`**
- Same data but normalized to [0, 1] range
- Useful for certain algorithms

---

# 6. OUTPUT RESULTS

## Archive: `rsf_risk_scores.tar.gz`

**Creates**: `rsf_risk_scores/` directory  
**Contains**: Multiple result files  
**Purpose**: Random Survival Forest model predictions and metrics

### File 1: `sample_risk_info.csv`

**Purpose**: Per-sample risk scores and survival outcomes  
**Format**: CSV  
**Rows**: 378 (one per sample)  
**Columns**: 3

**Structure**:

| acquisition_id | survival_days | survival_status | rsf_risk_score |
|---|---|---|---|
| UPMC_c001_v001_r001_reg001 | 1001.0 | 1 | 2.456 |
| UPMC_c001_v001_r001_reg002 | 831.0 | 1 | 1.234 |
| UPMC_c001_v001_r001_reg004 | 61.0 | 1 | 3.567 |
| UPMC_c002_v001_r001_reg001 | 1825.0 | 0 | 0.987 |

**Interpretation**:
- `survival_days`: Follow-up time (censoring time)
- `survival_status`: 1 = death event, 0 = censored/alive
- `rsf_risk_score`: Predicted risk (higher = higher risk of death)

---

### File 2: `patient_risk_info.csv`

**Purpose**: Patient-level (not sample-level) risk stratification  
**Format**: CSV  
**Rows**: 7 (one per unique patient)  
**Columns**: 5

**Structure**:
```csv
patient_id,survival_days,survival_status,rsf_risk_score,risk_cohort
2836,1001.0,1,2.456,high
2830,831.0,1,1.234,high
2823,61.0,1,3.567,high
2821,2091.0,0,0.234,low
```

**Interpretation**:
- `risk_cohort`: Categorical assignment
  - **`high`**: Patient has high predicted death risk
  - **`low`**: Patient has low predicted death risk
- Cutoff typically at median or specific percentile

---

### Additional Output Files (if generated)

**Subdirectory: `estimators/`**

Contains saved models:
```
estimators/
├── biomarker_region_cv0_rsf.pkl
├── biomarker_region_cv1_rsf.pkl
├── ...
├── biomarker_region_cv9_rsf.pkl
├── biomarker_cell_cv0_rsf.pkl
├── ...
└── ripley_cv9_rsf.pkl
```

Each `.pkl` file is a trained RandomSurvivalForest model for:
- Feature type (biomarker_region, biomarker_cell, etc.)
- Cross-validation fold (0–9)

---

# DATA FLOW SUMMARY

## Input → Processing → Output

```
dataset_info.tar.gz
├── marker_names.csv ─────┐
├── sample_metadata.csv ──┼─→ Feature Engineering ──→ Feature Matrix
├── cell_locations_and_labels.csv ─┐
└── labeled_arcsinh_norm_data.pkl ─┤

biomarker_expr_summary.tar.gz────→ Feature 1: Biomarker Region (40 features)
                                  Feature 2: Biomarker Cell (40 features)

                                  Feature 3: Cell Type Proportion (16 features)

all_k_neighborhood_mats.tar.gz───→ Feature 4: Neighborhood Matrix (256 features)

k_fns_norm_by_uw_qc_labeled.csv──→ Feature 5: Ripley K Function (80 features)

            ↓
        Feature Matrix
    (307 samples × ~392 features)
            ↓
    Random Survival Forest
       (10-fold CV)
            ↓
rsf_risk_scores.tar.gz
├── sample_risk_info.csv
└── patient_risk_info.csv
```

---

# QUICK REFERENCE

## How to Load Each File Type

### Python Loading Examples

```python
import pandas as pd
import numpy as np
import pickle

# 1. Marker names
markers = pd.read_csv('dataset_info/marker_names.csv', header=None)[0].tolist()

# 2. Metadata
metadata = pd.read_csv('dataset_info/sample_metadata.csv', index_col='acquisition_id')

# 3. Cell data
cells = pd.read_csv('dataset_info/cell_locations_and_labels.csv')

# 4. Cell expression
with open('dataset_info/labeled_arcsinh_norm_data.pkl', 'rb') as f:
    expr_df = pickle.load(f)

# 5. Biomarker summary
biomarker_sum = pd.read_csv('expression_biomarkers/UPMC_c001_v001_r001_reg001_cell_info.csv', index_col=0)

# 6. Neighborhood matrix (k=10)
neigh = np.load('k10/UPMC_c001_v001_r001_reg001.npy')

# 7. K-functions
kfns = pd.read_csv('k_fns_norm_by_uw_qc_labeled.csv', index_col='acquisition_id')

# 8. Results
results = pd.read_csv('rsf_risk_scores/sample_risk_info.csv', index_col='acquisition_id')
```

---

## File Sizes Summary

| Archive | Extracted Size | Purpose |
|---------|---|---------|
| `dataset_info.tar.gz` | ~1 GB | Core data |
| `biomarker_expr_summary.tar.gz` | ~100 MB | Biomarker summaries |
| `all_k_neighborhood_mats.tar.gz` | ~500 MB | Neighborhood matrices |
| `k_fns_norm_by_uw_qc_labeled.csv` | ~50 MB | Ripley functions |
| `patwa_comparisons.tar.gz` | ~500 MB | Comparison features |
| `denvar.tar.gz` | ~200 MB | DenVar analysis |
| `rsf_risk_scores.tar.gz` | ~10 MB | Model outputs |
| **TOTAL** | **~2.5 GB** | |

---

## Dataset Statistics

| Metric | Value |
|--------|-------|
| Patients | 7 |
| Tissue Samples (regions) | 378 |
| Cells | ~570,000+ |
| Biomarkers (IHC markers) | 40 |
| Cell Types | 16 |
| Features Generated | ~392 |
| Cross-Validation Folds | 10 |
| Survival Events (deaths) | ~185 (49%) |
| Follow-up Range | 1–2500+ days |

---

