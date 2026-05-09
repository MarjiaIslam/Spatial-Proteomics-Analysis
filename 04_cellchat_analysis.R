# CellChat Analysis from MIBI Protein-as-Gene Data
# ================================================
# R script for cell-cell communication inference.
#
# Prerequisites:
#   install.packages("devtools")
#   devtools::install_github("jinworks/CellChat")
#   install.packages(c("reticulate", "data.table"))

library(CellChat)
library(data.table)

# ─────────────────────────────────────────────────
# 1. Load exported data from Python pipeline
# ─────────────────────────────────────────────────
expr_mat  <- fread("output/gene_expression_matrix.csv", header = TRUE)
cell_meta <- fread("output/cell_metadata.csv", header = TRUE)

# Convert to matrix (genes × cells for CellChat)
gene_names <- colnames(expr_mat)[-1]  # drop row index
expr_matrix <- as.matrix(expr_mat[, -1])
rownames(expr_matrix) <- expr_mat[[1]]
expr_matrix <- t(expr_matrix)  # genes × cells

# Metadata
meta <- as.data.frame(cell_meta)
rownames(meta) <- meta$V1
meta$labels <- meta$cell_type

cat("Expression matrix:", nrow(expr_matrix), "genes x", ncol(expr_matrix), "cells\n")
cat("Cell types:", length(unique(meta$labels)), "\n")


# ─────────────────────────────────────────────────
# 2. Create CellChat object
# ─────────────────────────────────────────────────
cellchat <- createCellChat(
  object = expr_matrix,
  meta   = meta,
  group.by = "labels"
)

# Use human CellChatDB
CellChatDB <- CellChatDB.human

# IMPORTANT: With only 40 protein-derived genes, most ligand-receptor
# pairs won't be found. Filter to interactions involving our genes.
our_genes <- rownames(expr_matrix)

# Check which LR pairs overlap with our panel
all_lr <- CellChatDB$interaction
lr_genes <- unique(c(
  all_lr$ligand, all_lr$receptor,
  unlist(strsplit(all_lr$receptor, "_"))
))
overlap <- intersect(our_genes, lr_genes)
cat("\nGenes in CellChatDB LR pairs:", length(overlap), "/", length(our_genes), "\n")
cat("Overlapping genes:", paste(overlap, collapse = ", "), "\n")

# Use full database — CellChat will automatically subset to available genes
cellchat@DB <- CellChatDB

# ─────────────────────────────────────────────────
# 3. Preprocessing
# ─────────────────────────────────────────────────
cellchat <- subsetData(cellchat)

# NOTE: overExpression analysis may find few interactions with 40 genes; this is expected.
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)

# ─────────────────────────────────────────────────
# 4. Communication inference
# ─────────────────────────────────────────────────
# computeCommunProb uses a permutation test
cellchat <- computeCommunProb(cellchat, type = "truncatedMean", trim = 0.1)

# Filter communications with < 10 cells in either group
cellchat <- filterCommunication(cellchat, min.cells = 10)

# Compute pathway-level communication
cellchat <- computeCommunProbPathway(cellchat)

# Aggregate network
cellchat <- aggregateNet(cellchat)

# ─────────────────────────────────────────────────
# 5. Visualization
# ─────────────────────────────────────────────────
dir.create("output/cellchat_results", showWarnings = FALSE, recursive = TRUE)

# Network plot
pdf("output/cellchat_results/interaction_network.pdf", width = 10, height = 8)
netVisual_circle(cellchat@net$count, vertex.weight = table(meta$labels),
                 weight.scale = TRUE, label.edge = FALSE,
                 title.name = "Number of interactions")
dev.off()

# Heatmap
pdf("output/cellchat_results/interaction_heatmap.pdf", width = 10, height = 8)
netVisual_heatmap(cellchat, measure = "weight")
dev.off()

# Signaling roles
pdf("output/cellchat_results/signaling_roles.pdf", width = 12, height = 6)
netAnalysis_signalingRole_heatmap(cellchat, signaling = NULL, pattern = "outgoing")
dev.off()

# LR pairs
lr_df <- subsetCommunication(cellchat)
write.csv(lr_df, "output/cellchat_results/lr_interactions.csv", row.names = FALSE)
cat("\nDetected", nrow(lr_df), "ligand-receptor interactions\n")

# Save object
saveRDS(cellchat, "output/cellchat_results/cellchat_object.rds")
cat("\n✅ CellChat analysis complete. Results in output/cellchat_results/\n")

# ─────────────────────────────────────────────────
# 6. CAVEAT for interpretation
# ─────────────────────────────────────────────────
cat("\n")
cat("╔══════════════════════════════════════════════════════════════╗\n")
cat("║ IMPORTANT LIMITATIONS OF CELLCHAT ON PROTEIN DATA           ║\n")
cat("║                                                              ║\n")
cat("║ 1. Only 40 protein-derived genes — most LR pairs will be    ║\n")
cat("║    missing. Results are INCOMPLETE, not wrong.               ║\n")
cat("║ 2. Protein expression ≠ mRNA expression; some LR pairs may  ║\n")
cat("║    show different strengths than RNA-based analysis.         ║\n")
cat("║ 3. Secreted ligands measured by MIBI reflect surface-bound   ║\n")
cat("║    or intracellular protein, NOT secreted pool.              ║\n")
cat("║ 4. Validate key findings against literature or orthogonal   ║\n")
cat("║    data before biological conclusions.                       ║\n")
cat("╚══════════════════════════════════════════════════════════════╝\n")
