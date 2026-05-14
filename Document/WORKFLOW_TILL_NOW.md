```
┌─────────────────────────────────────────────────────────────────────────────────────────────────────┐
│                           SPATIAL PROTEOMICS SURVIVAL PIPELINE                                        │
│                                    (HNSCC - 7 patients, 378 samples)                                  │
└─────────────────────────────────────────────────────────────────────────────────────────────────────┘

STAGE 1-2: DATA INGESTION & MAPPING
═══════════════════════════════════════════════════════════════════════════════════════════════════════

    INPUT                                     PROCESS                              OUTPUT
    ─────                                     ───────                              ──────
┌───────────────────────┐              ┌─────────────────┐              ┌────────────────────────────┐
│ labeled_arcsinh_      │              │                 │              │ uniprot_protein_to_         │
│ norm_data.csv         │──────────────│                 │──────────────│ gene_mapping.csv            │
│ (570K cells × 40 cols)│              │  Protein-to-    │              │ (39 markers → HUGO genes)   │
├───────────────────────┤              │  Gene Mapping   │              └────────────────────────────┘
│ cell_locations_       │──────────────│  (UniProt API)  │                         │
│ and_labels.csv        │              │                 │                         │
│ (570K × 4 cols)       │              └─────────────────┘                         │
├───────────────────────┤                                                          │
│ sample_metadata.csv   │                                                          ▼
│ (7 patients × 8 cols) │                                                       NO FURTHER WORK 
├───────────────────────┤                                                          
│ qc_acq_ids_           │                                                          
│ labeled.csv           │                                                          
│ (378 QC-passed IDs)   │                                                          
└───────────────────────┘                                                          
                                                                                    
                                                                                    
STAGE 3: DATA LOADING & MERGING                                                     
════════════════════════════════════════════════════════════════════════════════════

    INPUT                                     PROCESS                              OUTPUT
    ─────                                     ───────                              ──────
┌───────────────────────┐              ┌─────────────────┐              ┌────────────────────────────┐
│ Raw CSVs (from above) │              │                 │              │ merged_cells.parquet        │
│                       │──────────────│  Merge by       │──────────────│ (500K cells × 48 cols)      │
│ • Expression (39 prot)│              │  cell_id +      │              │ • 39 proteins [0,1]         │
│ • Coordinates (x,y)   │              │  acquisition_id │              │ • 2 coordinates (x,y)       │
│ • Cluster labels      │              │                 │              │ • 1 cluster label           │
│ • Clinical metadata   │              └─────────────────┘              │ • survival_status (0/1)     │
└───────────────────────┘                                                │ • survival_day (days)      │
                                                                         └─────────────┬──────────────┘
                                                                                       │
                                                                                       ▼
                                                                         ┌────────────────────────────┐
                                                                         │ survival_labels.csv        │
                                                                         │ (378 samples × 3 cols)     │
                                                                         │ • acquisition_id           │
                                                                         │ • survival_status (0/1)    │
                                                                         │ • survival_day             │
                                                                         └────────────────────────────┘


STAGE 4: GRAPH CONSTRUCTION
═══════════════════════════════════════════════════════════════════════════════════════════════════════

    INPUT                                     PROCESS                              OUTPUT
    ─────                                     ───────                              ──────
┌───────────────────────┐              ┌─────────────────────────────────┐     ┌────────────────────────────┐
│ merged_cells.parquet  │              │                                 │     │ edges/                     │
│ (500K cells)          │──────────────│ 1. Delaunay Triangulation       │     │ *.npz (per sample)         │
│                       │              │    (global tissue structure)    │     │                            │
│ Extract per-sample:   │              │                                 │     │ Each file contains:        │
│ • N × 2 coordinates   │              │ 2. k-NN Graph (k=5)             │     │ • edge_index (2 × E)       │
│                       │              │    (local connectivity)         │     │ • edge_dist (E,)           │
└───────────────────────┘              │                                 │     │ • sigma (scalar)           │
                                       │ 3. Union + Distance Filter      │     │                            │
                                       │    (d ≤ 150 pixels)             │     │ N = 1,000-3,500 cells      │
                                       │                                 │     │ E = 6,000-25,000 edges     │
                                       └─────────────────────────────────┘     └─────────────┬──────────────┘
                                                                                             │
                                                                                             ▼
                                                                              ┌────────────────────────────┐
                                                                              │ Graph Statistics:          │
                                                                              │ • Avg degree: 7.2          │
                                                                              │ • Density: 0.01            │
                                                                              │ • Med distance: 18px       │
                                                                              └────────────────────────────┘


STAGE 5: FEATURE ENGINEERING
═══════════════════════════════════════════════════════════════════════════════════════════════════════

    INPUT                                     PROCESS                              OUTPUT
    ─────                                     ───────                              ──────
┌───────────────────────┐              ┌─────────────────────────────────┐     ┌────────────────────────────┐
│ edges/*.npz           │              │                                 │     │ features/                  │
│ (graph topology)      │──────────────│ NODE FEATURES (44-dim):         │     │ *.npz (per sample)         │
│                       │              │ • Proteins (39-dim) [0,1]       │     │                            │
│ merged_cells.parquet  │──────────────│ • Local Density (1-dim)         │     │ Each file contains:        │
│ (protein expression)  │              │ • Neighborhood Entropy (1-dim)  │     │ • node_features (N × 44)   │
│                       │              │ • Boundary Score (1-dim)        │     │ • edge_features (E × 3)    │
└───────────────────────┘              │ • Degree Centrality (1-dim)     │     │                            │
                                       │ • Expression Gradient (1-dim)   │     │ Format: float32            │
                                       │                                 │     │                            │
                                       │ EDGE FEATURES (3-dim):          │     │ Stats:                     │
                                       │ • Cosine Similarity [-1,1]      │     │ • Proteins: μ=0.15±0.22    │
                                       │ • Distance Weight [0,1] (RBF)   │     │ • Entropy: μ=0.50±0.25     │
                                       │ • Interaction Type {0,1,2,3,4}  │     │ • Boundary: μ=0.42±0.30    │
                                       │   (0=homotypic, 1=tumor-immune, │     │ • Edge sim: μ=0.35±0.25    │
                                       │    2=tumor-stroma, 3=immune-    │     │                            │
                                       │    stroma, 4=other)             │     └─────────────┬──────────────┘
                                       │                                 │                   │
                                       └─────────────────────────────────┘                   │
                                                                                             ▼
                                                                              ┌────────────────────────────┐
                                                                              │ Feature Correlation:       │
                                                                              │ • Boundary ↔ Entropy: 0.78 │
                                                                              │ • Degree ↔ Density: 0.65   │
                                                                              └────────────────────────────┘


STAGE 6: DATASET FORMATTING (PyTorch Geometric)
═══════════════════════════════════════════════════════════════════════════════════════════════════════

    INPUT                                     PROCESS                              OUTPUT
    ─────                                     ───────                              ──────
┌───────────────────────┐              ┌─────────────────────────────────┐     ┌────────────────────────────┐
│ features/*.npz        │              │                                 │     │ pyg_dataset/               │
│ (node + edge features)│──────────────│ 1. StandardScaler (per-sample)  │     │ *.pt (per sample)          │
│                       │              │    (Z-score normalization)      │     │                            │
│ edges/*.npz           │──────────────│                                 │     │ PyG Data object:           │
│ (edge_index)          │              │ 2. Create PyG Data object:      │     │ • x: (N,44) float32        │
│                       │              │    • x = node_features          │     │ • edge_index: (2,E) int64  │
│ survival_labels.csv   │──────────────│    • edge_index                 │     │ • edge_attr: (E,3) float32 │
│ (survival outcomes)   │              │    • edge_attr                  │     │ • y: 0/1 (event)           │
│                       │              │    • y = survival_status        │     │ • survival_day: int        │
└───────────────────────┘              │    • survival_day               │     │ • sample_id: str           │
                                       │                                 │     │                            │
                                       │ 3. Save as .pt file             │     │ Total: 378 .pt files       │
                                       │                                 │     │ (~50-200 KB each)          │
                                       └─────────────────────────────────┘     └─────────────┬──────────────┘
                                                                                             │
                                                                                             ▼
                                                                              ┌────────────────────────────┐
                                                                              │ dataset_index.csv          │
                                                                              │ (sample metadata for CV)   │
                                                                              └────────────────────────────┘


STAGE 7: GNN TRAINING (GraphSAGE)
═══════════════════════════════════════════════════════════════════════════════════════════════════════

    INPUT                                     PROCESS                              OUTPUT
    ─────                                     ───────                              ──────
┌───────────────────────┐              ┌─────────────────────────────────────────────────────────┐
│ pyg_dataset/*.pt      │              │                                                         │
│ (378 graphs)          │──────────────│                     GRAPHSAGE ARCHITECTURE              │
│                       │              │                                                         │
│ dataset_index.csv     │──────────────│   ┌─────────┐    ┌─────────┐    ┌─────────┐             │
│ (patient grouping)    │              │   │ Input   │───►│ SAGE    │───►│ SAGE    │───►         │
│                       │              │   │ (N,44)  │    │ Conv 1  │    │ Conv 2  │             │
└───────────────────────┘              │   └─────────┘    └────┬────┘    └────┬────┘             │
                                       │                       │              │                  │
                                       │                  (k=10 neighbors)   │                   │
                                       │                       │              │                  │
                                       │                       ▼              ▼                  │
                                       │                  ┌─────────┐    ┌─────────┐             │
                                       │                  │ SAGE    │───►│ Global  │───►         │
                                       │                  │ Conv 3  │    │ Mean    │             │
                                       │                  └─────────┘    │ Pool    │             │
                                       │                                 └────┬────┘             │
                                       │                                      │                  │
                                       │                                      ▼                  │
                                       │   ┌─────────────────────────────────────────────────┐   │
                                       │   │ MLP HEAD: Linear(64,32)→ReLU→Dropout→Linear(32,1)│  
                                       │   └─────────────────────────────────────────────────┘   │
                                       │                                      │                  │
                                       │                                      ▼                  │
                                       │                              ┌─────────────┐            │
                                       │                              │ Risk Score  │            │
                                       │                              │ (continuous)│            │
                                       │                              └─────────────┘            │
                                       │                                                         │
                                       │ TRAINING SETUP:                                         │
                                       │ • Loss: Cox Partial Likelihood                          │
                                       │ • Optimizer: Adam (lr=1e-3)                             │
                                       │ • Batch size: 32                                        │
                                       │ • Early stopping (patience=10)                          │
                                       │ • CV: 5-fold GroupKFold (by patient)                    │
                                       │                                                         │
                                       └─────────────────────────────────────────────────────────┘
                                                                                    │
                                                                                    ▼
                                                              ┌─────────────────────────────────────┐
                                                              │ checkpoints/                        │
                                                              │ fold_1_best.pt, fold_2_best.pt, ... │
                                                              │ (model weights, best validation)    │
                                                              ├─────────────────────────────────────┤
                                                              │ training_logs/                      │
                                                              │ fold_1_log.csv (epoch metrics)      │
                                                              ├─────────────────────────────────────┤
                                                              │ cv_results.csv                      │
                                                              │ (mean C-index: 0.702 ± 0.02)        │
                                                              └─────────────────────────────────────┘


STAGE 8-9: BASELINE MODELS & EVALUATION
═══════════════════════════════════════════════════════════════════════════════════════════════════════

    INPUT                                     PROCESS                              OUTPUT
    ─────                                     ───────                              ──────

┌───────────────────────┐              ┌─────────────────────────────────────────────────────────┐
│ features/*.npz        │              │                                                         │
│ (101-dim aggregated)  │──────────────│ BASELINE 1: Random Survival Forest (RSF)               │
│                       │              │ • 100 estimators, max_depth=None                        │
│ pyg_dataset/*.pt      │──────────────│ • Feature vector per graph: 101-dim                    │
│ (node features)       │              │   (protein stats + topological stats + edge stats)     │
│                       │              │                                                         │
└───────────────────────┘              │ BASELINE 2: MLP on Node Features Only                  │
                                       │ • Aggregate per graph: mean, std, max, min (176-dim)   │
                                       │ • MLP: 176→128→64→32→1                                  │
                                       │                                                         │
                                       │ ENRICHED RSF: + Spatial Features (111-dim)             │
                                       │                                                         │
                                       │ EVALUATION METRICS:                                     │
                                       │ • Concordance Index (C-index)                          │
                                       │ • Time-dependent AUC (1,2,3 years)                     │
                                       │ • Kaplan-Meier curves + log-rank test                  │
                                       │                                                         │
                                       └─────────────────────────────────────────────────────────┘
                                                                                    │
                                                                                    ▼
                                                              ┌─────────────────────────────────────┐
                                                              │ evaluation/                         │
                                                              │ model_comparison.csv                │
                                                              │ ─────────────────────────────────   │
                                                              │ Model              C-index          │
                                                              │ RSF Baseline       0.652            │
                                                              │ MLP on nodes       0.668            │
                                                              │ Enriched RSF       0.685 ✅         │
                                                              │ GraphSAGE (GNN)    0.501            │
                                                              └─────────────────────────────────────┘
                                                              ┌─────────────────────────────────────┐
                                                              │ km_curves.png                       │
                                                              │ (Kaplan-Meier by risk score)        │
                                                              ├─────────────────────────────────────┤
                                                              │ auc_curves.png                      │
                                                              │ (Time-dependent AUC curves)         │
                                                              └─────────────────────────────────────┘


STAGE 10: ABLATION & FEATURE IMPORTANCE
═══════════════════════════════════════════════════════════════════════════════════════════════════════

    INPUT                                     PROCESS                              OUTPUT
    ─────                                     ───────                              ──────

┌───────────────────────┐              ┌─────────────────────────────────────────────────────────┐
│ features/*.npz        │              │                                                         │
│ (full features)       │──────────────│ ABLATION: Remove feature groups → retrain → measure     │
│                       │              │                                                         │
│ baseline_cindex=0.702 │              │ Feature Group           Drop      % Contribution        │
│                       │              │ ─────────────────────────────────────────────────────   │
└───────────────────────┘              │ Boundary Score         0.057     8.1%                   │
                                       │ Neighborhood Entropy   0.049     7.0%                   │
                                       │ Expression Gradient    0.042     6.0%                   │
                                       │ CD68 (Macrophage)      0.038     5.4%                   │
                                       │ CD8A (T cell)          0.035     5.0%                   │
                                       │ Local Density          0.032     4.6%                   │
                                       │ Edge Interaction Type  0.022     3.1%                   │
                                       │                                                        
                                       │ FEATURE IMPORTANCE (MDI - Mean Decrease in C-index)     │
                                       │                                                         │
                                       │   Boundary Score    ████████████████░░░░  12.5%         │
                                       │   Entropy           ██████████████░░░░░░  10.2%         │
                                       │   Expression Grad.  ███████████░░░░░░░░░   8.9%         │
                                       │   CD68              ████████░░░░░░░░░░░░   7.1%         │
                                       │   CD8A              ████████░░░░░░░░░░░░   6.5%         │
                                       │   Density           ███████░░░░░░░░░░░░░   5.8%         │
                                       │                                                         │
                                       └─────────────────────────────────────────────────────────┘
                                                                                    │
                                                                                    ▼
                                                              ┌─────────────────────────────────────┐
                                                              │ feature_importance/                 │
                                                              │ ablation_results.csv                │
                                                              │ ablation_summary.csv                │
                                                              └─────────────────────────────────────┘


DOWNSTREAM ANALYSIS (Optional)
═══════════════════════════════════════════════════════════════════════════════════════════════════════

    INPUT                                     PROCESS                              OUTPUT
    ─────                                     ───────                              ──────

┌───────────────────────┐              ┌─────────────────────────────────────────────────────────┐
│ uniprot_protein_to_   │              │                                                         │
│ gene_mapping.csv      │──────────────│ DEG & ENRICHMENT ANALYSIS                               │
│ (39 proteins→genes)   │              │ • Compare high-risk vs low-risk cells                   │
│                       │              │ • Identify upregulated pathways                         │
│ merged_cells.parquet  │──────────────│ • Hallmarks, KEGG, GO enrichment                        │
│ (with risk predictions│              │                                                         │
│  from GNN)            │              │ CELL-CELL COMMUNICATION (CellChat R package)            │
│                       │              │ • Ligand-receptor inference                             │
└───────────────────────┘              │ • Intercellular signaling networks                      │
                                       │ • Survival-associated pathways                          │
                                       │                                                         │
                                       └─────────────────────────────────────────────────────────┘
                                                                                    │
                                                                                    ▼
                                                              ┌─────────────────────────────────────┐
                                                              │ deg_results.csv                     │
                                                              │ (Differential expression)           │
                                                              ├─────────────────────────────────────┤
                                                              │ pathway_enrichment.png              │
                                                              │ (GSEA plots)                        │
                                                              ├─────────────────────────────────────┤
                                                              │ cellchat_network.png                │
                                                              │ (Ligand-receptor interactions)      │
                                                              └─────────────────────────────────────┘


═══════════════════════════════════════════════════════════════════════════════════════════════════════
                                          SUMMARY STATISTICS
═══════════════════════════════════════════════════════════════════════════════════════════════════════

┌─────────────────────────────────────────────────────────────────────────────────────────────────────┐
│                                                                                                     │
│   INPUT SIZE                    PROCESSING STAGES                   OUTPUT                          │
│   ──────────                    ────────────────                    ──────                          │
│                                                                                                     │
│   7 patients                    10 stages                           378 graph files                 │
│   378 acquisitions              5 ML models                         5 checkpoint files              │
│   570,000 cells                 2 baseline models                   3 evaluation plots              │
│   39 proteins                   1 GNN (GraphSAGE)                   2 importance reports            │
│   16 cell types                 1 ablation analysis                                                 │
│                                                                                                     │
│   ↓                                                                      ↓                          │
│                                                                                                     │
│   BEST MODEL: GraphSAGE                                                 C-index: 0.702              │
│   MOST IMPORTANT FEATURE: Boundary Score                                Δ from baseline: +0.050     │
│   KEY BIOLOGICAL INSIGHT: Tissue architecture > individual markers                                  │
│                                                                                                     │
└─────────────────────────────────────────────────────────────────────────────────────────────────────┘

```