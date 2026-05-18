# Cell-Cell Communication Papers: How They Help This Project

 our current pipeline already builds one graph per tissue region, where cells are nodes, spatial neighbors are edges, node features include protein expression plus spatial statistics, and survival models use those graph features. 
 The missing next layer is biological cell-cell communication: not just "which cells are near each other?", but "which nearby cells may be signaling to each other?"

### Paper 1: The good, the bad, and the ugly: opportunities, challenges, and pitfalls in spatial proteomics modeling (2026).

Review paper of spatial proteomics that explains why spatial proteomics is powerful, where it fails, and what can go wrong in computational modeling. 

Most useful points for our project:
- Spatial proteomics can directly measure protein states, which is useful because proteins are closer to function than mRNA.
- Graph models are natural for spatial proteomics: cells can be nodes, protein expression can be node features, and spatial proximity can define edges.
- But proximity alone is not enough. A nearby cell pair may not be biologically communicating.
- Segmentation errors can strongly corrupt cell-type assignment, expression values, neighborhood statistics, and inferred cell-cell interactions.
- Batch effects and sample-preparation artifacts can look like real spatial patterns.
- Any communication feature should be validated or at least sensitivity-tested.

Practical lesson:
Use spatial proximity as the candidate interaction graph, but add biological evidence before calling it communication.

### Paper 2: Spatial proteomics in precision medicine: technologies, bioinformatics, and translational applications (2026).

A review paper that focuses on spatial proteomics technologies, bioinformatics, clinical translation, and precision medicine. It discusses technologies like CODEX/PhenoCycler, MIBI-TOF, IMC, DSP, mass spectrometry imaging, and multi-omics integration.

Most useful points for our project:
- Multiplexed imaging data can support cell-neighborhood analysis, tumor-immune interaction analysis, and prognosis/therapy-response modeling.
- AI/GNN methods are appropriate because spatial proteomics data are graph-like.
- Reproducibility needs per-marker QC, careful segmentation, batch-aware analysis, and interpretable features.
- Spatial multi-omics and ligand-receptor tools can guide interaction modeling, but direct protein panels are usually limited, so missing ligand/receptor genes must be expected.

Practical lesson:
Our pipeline should treat cell-cell communication as an extra graph feature layer, not as a fully certain biological truth.

**We can use these to defend our thesis:**
- graph-based spatial proteomics modeling.
- why cell-cell communication should be spatially constrained.
- Designing QC checks for segmentation, batch effects, and false gradients.
- Choosing primary methods to adapt: CellChat, Squidpy, COMMOT, NATMI, stLearn, MISTy/LIANA+, Giotto, GraphCompass.

## 3. Primary Papers 

### Squidpy

Paper: **Squidpy: a scalable framework for spatial omics analysis**

Link: https://www.nature.com/articles/s41592-021-01358-2

What it does:
Squidpy provides spatial graph construction, neighborhood enrichment, co-occurrence analysis, spatial autocorrelation, and ligand-receptor style analysis through OmniPath/CellPhoneDB-style resources.
It can compute spatial neighborhood enrichment cleanly.

How to adapt:
- Use our AnnData object from `01_protein_to_gene_pipeline.py`.
- Store cell coordinates in `adata.obsm["spatial"]`.
- Use Squidpy-style spatial neighbors to validate our own graph construction.
- Compute cell-type neighborhood enrichment.
- Add enrichment scores as graph-level or edge-level features.

Possible features:
- `neighbor_enrichment_celltype_pair`
- `observed_over_expected_adjacency`
- `celltype_cooccurrence_score`

### COMMOT

Paper: **Screening cell-cell communication in spatial transcriptomics via collective optimal transport**

Link: https://www.nature.com/articles/s41592-022-01728-4

What it does:
COMMOT uses optimal transport to infer spatial communication while considering ligand-receptor expression, competition among signals, and spatial distance.
It is designed for spatial transcriptomics, but the scoring idea is excellent for graph edges.

How to adapt:
- For each candidate edge, compute a ligand-receptor score.
- Weight it by spatial distance.
- Normalize across competing sender/receiver cells so one highly expressed ligand does not unrealistically connect to everyone.


### MISTy / LIANA+

Paper: **LIANA+ provides an all-in-one framework for cell-cell communication inference**

Link: https://www.nature.com/articles/s41556-024-01469-w

What it does:
LIANA+ combines multiple CCC methods and includes spatially weighted metrics. It adapts MISTy-style multi-view modeling to spatial data.

How to adapt:
- Treat each cell as having multiple "views":
  - intrinsic protein expression
  - immediate-neighbor proteins
  - broader-neighborhood cell-type composition
  - ligand-receptor interaction scores
- Use these views as node-level features before graph pooling.


### GraphCompass

Paper: **GraphCompass: spatial metrics for differential analyses of cell organization across conditions**

Link: https://academic.oup.com/bioinformatics/article/40/Supplement_1/i548/7700863

What it does:
GraphCompass compares spatial cell organization between conditions using graph metrics.
We can use similar graph metrics to compare high-risk vs low-risk tissue regions.

Possible features:
- Cell-type-specific subgraph density.
- Tumor-immune edge fraction.
- Immune-stroma edge fraction.
- Centrality of macrophages/CD8 cells/tumor cells.
- Graph distances between high-risk and low-risk samples.

## 4. What Information We Can Use

1. Ligand-receptor databases:
   - CellChatDB
   - OmniPath
   - CellPhoneDB
   - connectomeDB2020

2. Edge-level LR scoring:
   - product: `ligand_i * receptor_j`
   - geometric mean: `sqrt(ligand_i * receptor_j)`
   - min score: `min(ligand_i, receptor_j)`
   - pathway aggregation: sum or mean over LR pairs in pathway

3. Spatial constraint:
   - only score existing graph edges
   - or only score cells within a radius
   - apply distance decay

4. Statistical validation:
   - shuffle cell labels
   - shuffle spatial coordinates within sample
   - compare observed LR edge scores against null scores

5. Features for survival model:
   - per-edge LR score
   - per-edge pathway score
   - sender/receiver cell-type pair score
   - per-node outgoing communication strength
   - per-node incoming communication strength
   - per-sample communication summaries

## 5. Proposed Pipeline for Our Project

### Stage 1: Keep Existing Pipeline


```text
01_data_loading.py
  -> merged cell table

02_graph_construction.py
  -> spatial graph edges

03_feature_engineering.py
  -> protein + spatial node features and simple edge features

04_normalize_and_save.py
  -> PyTorch Geometric graph objects
```

### Stage 2: Add Ligand-Receptor Resource

Create a new file:

```text
spatial_survival/lr_database.py
```

This should contain a small LR database filtered to genes available in our panel.

Start with:

- CellChatDB human
- OmniPath
- NATMI/connectomeDB2020

But after filtering to our panel, many pairs may disappear. That is expected.

Example available genes from our panel:

```text
CD274/PDCD1     # PDL1-PD1 checkpoint axis
CTLA4/CD80-like # CTLA4 is available, but CD80/CD86 may not be in panel
ICOS/ICOSLG     # ICOS available, ICOSLG likely absent
CD47/SIRPA      # CD47 available, SIRPA likely absent
PECAM1/PECAM1   # CD31 endothelial adhesion
PDPN/CLEC2      # PDPN available, receptor may be absent
```

Because our panel is protein-marker focused, we should also include marker-defined contact biology:

- Tumor-CD8 adjacency
- Tumor-macrophage adjacency
- Tumor-stroma adjacency
- CD8-macrophage adjacency
- Vessel-immune adjacency
- PD1+ immune near PDL1+ tumor/macrophage
- Ki67+ tumor near immune exclusion or infiltration zones

### Stage 3: Compute Edge Communication Features

Create a new script:

```text
spatial_survival/12_cell_communication_features.py
```

For each graph edge `i -> j`, compute:

```text
spatial_weight = exp(-distance^2 / (2 * sigma^2))
lr_raw_score = ligand_expr_i * receptor_expr_j
lr_spatial_score = lr_raw_score * spatial_weight
```

Aggregate:

```text
edge_comm_total
edge_comm_max
edge_comm_count_nonzero
edge_checkpoint_score
edge_adhesion_score
edge_immune_activation_score
edge_tumor_immune_score
```

Add these to `edge_attr`.

Current edge features:

```text
[cosine_similarity, distance_weight, interaction_type]
```

Proposed enriched edge features:

```text
[
  cosine_similarity,
  distance_weight,
  interaction_type,
  lr_total_score,
  lr_max_score,
  lr_pair_count,
  checkpoint_score,
  tumor_immune_contact_score
]
```

### Stage 4: Compute Node Communication Features

For each cell/node:

```text
outgoing_comm_strength = sum scores from this cell to neighbors
incoming_comm_strength = sum scores from neighbors to this cell
checkpoint_incoming
checkpoint_outgoing
communication_entropy
top_partner_celltype_fraction
```

Add these to node features.

Current node features:

```text
39 protein features + 5 spatial features = 44
```

Proposed:

```text
39 protein features
+ 5 spatial features
+ 5-10 communication features
```

### Stage 5: Sample-Level Communication Summary

For Random Survival Forest and MLP models, summarize each sample:

```text
mean LR score
max LR score
fraction of high-communication edges
tumor->immune score
immune->tumor score
macrophage->tumor score
T cell checkpoint score
PD1-PDL1 local score
cell-type pair communication matrix flattened
```

This will fit naturally into `10_rsf_feature_importance.py` and `11_mlp_enriched.py`.

