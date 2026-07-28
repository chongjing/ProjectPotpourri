# =============================================================================
# H. glycines gene-centred VST expression clustering (Mfuzz, k=30)
# =============================================================================
# Purpose:
#   Soft-cluster H. glycines genes by VST-normalised expression profiles
#   across 7 developmental stages / conditions (E, pJ2, ppJ2, J3, J4, RF, RM,
#   each with 3 biological replicates = 21 samples total).
#
# Method:
#   Mfuzz (soft clustering via fuzzy c-means) implemented in the ClusterGVis
#   package.  k=30 was selected as the final cluster number based on elbow
#   analysis (k-means WSS, k=1–50).
#
# Input:  03.TPM4cluster.tsv
#   - 14,862 H. glycines genes (rows)
#   - 21 VST-normalised expression values (columns)
#   - Values are mean-centred and variance-stabilised
#
# Outputs:
#   04.getCluster.pdf / .jpeg — Elbow plot (k=1–50, WSS)
#   05.k30.line_heatmap.pdf / .jpeg — Line + heatmap visualisation (k=30)
#   05.cm.30clusters.csv — Cluster assignments (wide format, 25 columns)
#
# Dependencies:
#   R 4.x, Bioconductor, ClusterGVis, Biobase, Mfuzz, factoextra
#   Install: BiocManager::install(c("ClusterGVis", "Biobase", "Mfuzz"))
#            install.packages("factoextra")
# =============================================================================

library(ClusterGVis)   # Mfuzz clustering wrapper + visualisation
library(Biobase)       # ExpressionSet container (required by Mfuzz)
library(Mfuzz)         # Soft clustering (fuzzy c-means)
library(factoextra)    # Elbow plot (fviz_nbclust)

# ---------------------------------------------------------------------------
# 1. Load input data
# ---------------------------------------------------------------------------
# 03.TPM4cluster.tsv: 14,862 genes × 21 samples
# Gene IDs are row names, columns are biological conditions
counts <- read.csv("03.TPM4cluster.tsv",
                    sep      = "\t",
                    row.names = 1)

cat(sprintf("Input: %d genes × %d samples\n", nrow(counts), ncol(counts)))

# ---------------------------------------------------------------------------
# 2. Elbow plot — determine optimal number of clusters
# ---------------------------------------------------------------------------
# k-means WSS (within-cluster sum of squares) for k = 1 to 50
# Optimal k is at the "elbow" — the point where adding more clusters
# yields diminishing returns in variance explained.
# Here we use the second-derivative method as a heuristic.

pdf("04.getCluster.pdf", width = 8, height = 5)
elbow_plot <- fviz_nbclust(counts,
                            kmeans,
                            method  = "wss",
                            k.max   = 50)
print(elbow_plot)
dev.off()

jpeg("04.getCluster.jpeg",
     width = 10, height = 6, units = "in", res = 600)
print(elbow_plot)
dev.off()

cat("Elbow plot saved: 04.getCluster.pdf / .jpeg\n")

# Second-derivative elbow detection (for reference)
wss <- sapply(1:50, function(k) {
  kmeans(counts, centers = k, nstart = 10)$tot.withinss
})
diff1  <- diff(wss)
diff2  <- diff(diff1)
elbow_k <- which.min(diff2) + 1
cat(sprintf("Elbow-detected k: %d\n", elbow_k))

# ---------------------------------------------------------------------------
# 3. Mfuzz soft clustering — k = 30
# ---------------------------------------------------------------------------
# After exploring k = 12, 15, 20, 30, we select k = 30 as the final
# cluster number, providing the best resolution of expression patterns
# across the 7 developmental stages.
#
# Mfuzz (fuzzy c-means) assigns each gene a membership score for every
# cluster (0–1), allowing genes to belong partially to multiple clusters.
# The final hard assignment is the cluster with the highest membership.

cm30 <- clusterData(counts,
                    cluster.method = "mfuzz",
                    cluster.num    = 30)

cat("Mfuzz clustering (k=30) complete\n")

# ---------------------------------------------------------------------------
# 4. Visualisation — line plots + heatmap
# ---------------------------------------------------------------------------
# The "both" plot.type shows:
#   - Top panel: mean expression profile per cluster (line plot, ±SD)
#   - Bottom panel: heatmap of individual gene expression within each cluster
# Colour gradient: green (low) → orange (mid) → red (high)

pdf("05.k30.line_heatmap.pdf", width = 9, height = 15)
visCluster(object        = cm30,
           plot.type     = "both",
           ms.col        = c("green", "orange", "red"),
           column_names_rot = 45)
dev.off()

jpeg("05.k30.line_heatmap.jpeg",
     width = 9, height = 15, units = "in", res = 600)
visCluster(object        = cm30,
           plot.type     = "both",
           ms.col        = c("green", "orange", "red"),
           column_names_rot = 45)
dev.off()

cat("Line+heatmap saved: 05.k30.line_heatmap.pdf / .jpeg\n")

# ---------------------------------------------------------------------------
# 5. Export cluster assignments
# ---------------------------------------------------------------------------
# Wide-format CSV:
#   Columns: gene, E_1–RM_3 (21 samples), cluster, membership
#   cluster   — hard assignment (1–30)
#   membership — membership score for the assigned cluster (0–1)

write.csv(cm30$wide.res,
          "05.cm.30clusters.csv",
          row.names = TRUE, quote = FALSE)

cat(sprintf("Cluster assignments saved: 05.cm.30clusters.csv (%d genes)\n",
            nrow(cm30$wide.res)))
cat("Analysis complete.\n")