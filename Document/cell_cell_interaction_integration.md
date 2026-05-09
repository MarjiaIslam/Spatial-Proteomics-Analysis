# Cell-Cell Interaction Integration

This note describes how to integrate cell-cell interaction signals into the spatial trajectory analysis workflow. It is intentionally conceptual and maps to the existing pipeline stages.

## Where It Fits in the Pipeline

1) Graph Construction: build the base spatial neighbor graph.
2) Feature Engineering: add interaction-derived edge and node features.
3) Trajectory Inference: use the interaction-weighted graph to learn progression.

## Conceptual Steps

### 1) Build the Base Spatial Graph

- Nodes: single cells (or spots).
- Edges: spatial neighbors (k-NN or radius threshold).
- Spatial weight: a distance decay, e.g. exp(-d^2 / (2 * sigma^2)).

This defines who can interact based on proximity.

### 2) Define Interaction Evidence (CellChat Database)

Use CellChat ligand-receptor (LR) pairs to compute interaction strength between spatial neighbors. The idea is to quantify signaling potential from sender cell i to receiver cell j.

Inputs
- Expression for ligands and receptors (protein or gene-mapped).
- CellChat LR pairs (ligand, receptor, pathway label).

Conceptual scoring
- For each LR pair, score the edge (i, j):
	- s_lr(i, j) = f(L_i, R_j)
	- A simple choice is the product or geometric mean of ligand and receptor expression.
- Aggregate across all LR pairs:
	- w_interaction_ij = sum over LR pairs of s_lr(i, j)
- Optionally apply pathway weights if some pathways are more relevant to the tissue or disease.

This defines who is likely to influence whom biologically based on known LR signaling.

### 3) Combine Spatial and Interaction Weights

Create a combined edge weight for the spatial graph:

- w_ij = alpha * w_spatial_ij + beta * w_interaction_ij

Use alpha and beta to control the balance between proximity and biology. The combined weight is used in downstream inference.

### 4) Add Interaction-Derived Features

Add features that summarize interaction at the edge and node level.

Edge features
- interaction_score(i, j): continuous strength.
- interaction_type(i, j): categorical label for the pair.
- combined_weight(i, j): final edge weight for trajectory.

Node features
- interaction_strength(i): sum or mean of interaction scores to neighbors.
- interaction_entropy(i): diversity of interaction partner types.
- interaction_centrality(i): centrality on the interaction-weighted graph.

These features enrich the existing density, entropy, and centrality features.

### 5) Run Trajectory on the Interaction-Weighted Graph

Convert combined weights to a transition matrix:

- P_ij = w_ij / sum_k w_ik

Trajectory inference (PAGA, Monocle3, diffusion) then follows interaction-supported transitions rather than purely spatial adjacency.

### 6) Validate the Effect of Interactions

Compare trajectories with and without interaction weights:

- Are branch points more biologically plausible?
- Do known micro-niches align with high interaction connectivity?
- Does progression score correlate better with clinical outcomes?

## Minimal Integration Summary

- Use spatial neighbors to define possible interactions.
- Use biological evidence to score how strong those interactions are.
- Inject those scores into edge weights and node summaries.
- Run trajectory inference on the interaction-weighted graph.

## CellChat-Specific Notes

- If you already run CellChat on the expression matrix, you can reuse its LR scores and restrict them to spatial neighbors.
- If you run LR scoring directly, make sure the same gene/protein mapping used in the pipeline is applied to the CellChat pairs.
- Consider filtering to the top LR pairs per cell type pair to reduce noise.

## Related References and How They Connect

The following resources are useful context for integrating cell-cell interactions into spatial trajectory analysis. The descriptions below are intentionally high level and based on the paper titles and public landing pages.

1) Nature Protocols (2024): NicheNet workflow for cell-cell communication
- Link: https://www.nature.com/articles/s41596-024-01121-9
- Focus: infers active ligands from transcriptomics data to explain target gene changes.
- Connection: provides a principled way to score LR signaling and prioritize ligand effects, which can be used to derive interaction weights between neighboring cells.

2) OSTA (Orchestrating Spatial Transcriptomics Analysis with Bioconductor)
- Link: https://lmweber.org/OSTA/
- Focus: a Bioconductor-based, modular spatial analysis framework.
- Connection: offers a structured way to organize spatial data objects, metadata, and analysis steps, which can host the interaction-weighted graph and trajectory stages in a reproducible workflow.

3) Nature Communications (2023): robust mapping of spatiotemporal trajectories and cell-cell interactions
- Link: https://www.nature.com/articles/s41467-023-43120-6
- Focus: integrates trajectories with cell-cell interactions in tissue contexts.
- Connection: supports the idea that interaction-aware graphs better reflect biological progression and can inform where trajectory paths and branch points are inferred.
