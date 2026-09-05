# Methods: Pooled-Control-Centered Heatmaps (v3)

## 1. Data Overview

Phenotype measurements were obtained from two independent effector-tissue assays using the same BM effector library:

- **Root assay**: 72 putative effectors + 3 control genotypes (Col-0 wild-type, EmptyVector, GFP), assessed for *Heterodera schachtii* infection phenotypes
- **Shoot assay**: 74 putative effectors + 3 control genotypes (Col-0 wild-type, EmptyVector, GFP), assessed for *Beth*-tissue infection phenotypes

Control replicate counts: Col-0 (root n=11, shoot n=12), EmptyVector (root n=7, shoot n=7), GFP (root n=11, shoot n=11).

Eight quantitative phenotypes were measured:
- **area_end_mm2**: Final lesion area at endpoint
- **AUC**: Area under the disease progress curve
- **AGR_mean**: Mean absolute growth rate
- **AGR_max**: Maximum absolute growth rate
- **RGR_mean**: Mean relative growth rate
- **RGR_peak**: Peak relative growth rate
- **RGR_final**: Final relative growth rate
- **tt10**: Time to reach 10 mm² (days post-inoculation or hours)

For shoot data, right-censored tt10 values (infections that never reached 10 mm²) were imputed as max(tt10) + 1.0. Rows with zero area (area_end_mm2 = 0) were excluded from analysis.

## 2. Control Homogeneity Verification

Prior to normalization, we verified that the three control genotypes are statistically indistinguishable across all phenotypes. For each phenotype, the Kruskal-Wallis test was applied to compare the three control groups:

| Phenotype | H statistic | p-value |
|-----------|-------------|---------|
| area_end_mm2 | 0.28 | 0.868 |
| AUC | 0.49 | 0.784 |
| AGR_mean | 0.08 | 0.960 |
| AGR_max | 0.40 | 0.820 |
| RGR_mean | 1.03 | 0.357 |
| RGR_peak | 0.49 | 0.780 |
| RGR_final | 1.38 | 0.251 |

All p > 0.05 — controls are homogeneous and can be pooled as a single baseline. EmptyVector showed higher coefficient of variation (2-3x) than Col-0 or GFP but this does not invalidate pooling given the non-significant Kruskal-Wallis result.

## 3. Pooled-Control-Centered Normalization

For each phenotype *p*, the pooled baseline was computed as the **median of all pooled control replicates**:

$$\text{pooled\_median}_p = \text{median}\left(\{x_{i,p} \mid i \in \text{Col-0} \cup \text{EmptyVector} \cup \text{GFP}\}\right)$$

The centered z-score for each effector *g* was then:

$$z_{g,p} = \frac{\text{median}(x_{g,p}) - \text{pooled\_median}_p}{\text{SD}\left(\{\text{median}(x_{j,p}) \mid j \in \text{all genes}\}\right)}$$

where the denominator is the population standard deviation of gene-level medians, computed with Bessel's correction (ddof=0). Using global SD (rather than control SD) preserves the distance structure of the original clustering — the pooled baseline is a constant shift per phenotype, leaving pairwise Euclidean distances between genes unchanged.

**No direction inversion** was applied to tt10. Higher tt10 values indicate slower growth (less virulent); this is correctly reflected in the heatmap (blue = faster, red = slower).

## 4. Cluster Analysis

### 4.1 Dendrogram Computation

Hierarchical clustering was performed using Ward's minimum variance method on Euclidean distances between gene median z-scores:

$$d(a, b) = \sqrt{\sum_{p=1}^{8} (z_{a,p} - z_{b,p})^2}$$

Ward linkage iteratively merges clusters to minimize within-cluster variance increase. The linkage matrix Z was computed with optimal leaf ordering (default scipy implementation).

### 4.2 Cluster Assignment

For the **recommended** partition, the dendrogram was cut at k=3 (root) or k=2 (shoot, combined) using `fcluster` with `criterion="maxclust"`. For the **supplementary** partition, k=8 was used for all three analyses. Cut heights were selected based on the dendrogram gap statistic and silhouette analysis.

### 4.3 Dendrogram Preservation

The current analysis **reuses the original dendrogram structure** from the primary clustering (global-centered z-scores). Because the pooled-control baseline is a constant shift per phenotype, pairwise distances between effectors remain identical:

$$d_{\text{pooled}}(a, b) = d_{\text{global}}(a, b)$$

Therefore, re-running Ward linkage on pooled-centered data produces the exact same dendrogram topology. We verified this computationally and preserved the original leaf ordering for consistency.

### 4.4 Combined Analysis

For the 16-dimensional combined heatmap, root and shoot z-score matrices were concatenated column-wise after intersecting on shared genes (75 genes). The original combined dendrogram (computed on 16-D global-centered data) was preserved.

## 5. Heatmap Rendering

### 5.1 Colormap

A custom diverging colormap with a white noise band was designed to highlight biologically meaningful deviations:

- |z| < 0.25: rendered as **white** (noise level)
- z = -1.0: light blue
- z = -2.5: deep blue
- z = +1.0: light red
- z = +2.5: deep red

A `TwoSlopeNorm` mapping [-2.5, -0.25, 0.25, 2.5] → [0, 0.45, 0.55, 1.0] ensures the white band occupies the central 10% of the color range, collapsing noise to background.

### 5.2 Figure Specifications

- **Resolution**: 1200 DPI raster (TIFF LZW compressed, JPEG quality 95), vector PDF/SVG
- **Font**: DejaVu Sans (or Liberation Sans fallback), 7pt axis labels, 5pt gene labels
- **Size**: Labelled heatmaps 8.5 × 11.2 inches; compact heatmaps 6.5 × 8.2 inches
- **Colorblind-safe**: Palette chosen for deuteranopia/protanopia accessibility

### 5.3 Annotation Tracks

- **Row colors**: Cluster assignment color-coded per CPAL palette
- **Column colors**: Tissue annotation (root = blue, shoot = red)
- **Reference tracks**: Per-control z-profiles plotted as line graphs showing deviation from pooled baseline (for diagnostic validation)

## 6. Software Implementation

- **Python**: 3.13
- **NumPy**: 1.x
- **Pandas**: 2.x
- **SciPy**: 1.x (scipy.cluster.hierarchy, scipy.spatial.distance)
- **Seaborn**: 0.x (clustermap)
- **Matplotlib**: 3.x
- **Code**: `02.v3_pooled_control_heatmaps/scripts/run_ctrl_heatmaps_v3.py`

## 7. Output Summary

| Category | Count | Format |
|----------|-------|--------|
| Main heatmaps | 24 (12 × 2 formats) | PDF, JPEG |
| Diagnostic plots | 10 (5 × 2 formats) | PDF, JPEG |
| Data tables | 8 | CSV |
| Run summary | 1 | Markdown |

The 12 main heatmaps cover: Root (k=3, k=8), Shoot (k=2, k=8), Combined (k=2, k=8), each in labelled and compact variants.
