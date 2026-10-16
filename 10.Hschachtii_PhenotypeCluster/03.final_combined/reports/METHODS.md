# Methods: final-sheet phenotype clustering

## Input and exclusion rules

The finalised root and shoot phenotype summary CSVs were analysed independently of previous clustering outputs. Records marked included=False were excluded. A further exclusion applied when both final sizes equalled zero for the same effector and replicate identifiers. A zero final size in one tissue alone did not trigger this exclusion. Matching identifiers was an operational filtering rule; biological pairing between tissues was not established. A tissue-specific exclusion did not remove an otherwise eligible record from the other tissue. No additional statistical outlier removal was performed. Input files remained unchanged.

## Traits and incomplete measurements

Eight traits were analysed per tissue: final size, area under the curve (AUC), capped time to 10 mm², mean and maximum absolute growth rate (AGR), and mean, peak and final relative growth rate (RGR). Values supplied in the final sheets were used without recomputing growth curves. Negative growth rates were retained.

The data owner confirmed that blank time_to_10mm2 values on included records indicated no threshold crossing during observation. With explicit approval, threshold time was encoded as a day-29-capped timing feature. Recorded times below 29 were unchanged. Recorded times at or above 29, and right-censored records observed through at least day 29, received a value of 29. Thus, 29 means "day 29 or later", not a crossing on day 29. This cap is on the supplied dpg time axis, not on duration_days. It was selected using the shortest available follow-up among retained right-censored records. Later crossing-time resolution is intentionally lost. The supplied left-censored root value for BM00092 replicate 1 was retained at its upper bound of day 8 and flagged. It is not an exact crossing time. Original times, end-of-observation times and censoring classifications were exported alongside the derived timing values. No max+1 replacement was used.

Other missing trait values, including undefined RGR measurements, were not imputed. Each genotype was represented by its median available value for each trait. Contributing replicate counts were exported for every genotype–trait combination. Consequently, RGR summaries describe replicates for which RGR was estimable and may have different contributing replicate sets from size or AUC summaries. Low-count summaries require caution. This analysis does not model censoring beyond the explicitly approved bounded timing representation.

## Scaling and colour references

For each tissue and trait, the reference mean (mu) and population standard deviation (sigma; ddof=0) were calculated from all available retained individual replicate measurements, including controls. The derived capped timing values entered the timing reference calculation. Reference statistics were not calculated from genotype medians. Therefore, genotypes with more retained measurements contribute more observations to these reference statistics.

For genotype median m, the all-replicate-centred value was (m-mu)/sigma. The alternative control-centred value was (m-c)/sigma, where c was the median of pooled retained Col-0, EmptyVector and GFP replicate measurements for that tissue and trait. Controls were pooled at the replicate level without equal-genotype weighting or a claim of statistical equivalence. The two views share the same denominator. Actual control rows were not forced to zero, and their average need not equal zero. A pooled-reference strip was drawn outside the dendrogram. It is zero in the control-centred view and (c-mu)/sigma in the all-replicate-centred view.

## Clustering and PCA

A new Euclidean-distance Ward hierarchy was calculated for each analysis using the scaled genotype profiles. Optimal leaf ordering was applied. Root and shoot were analysed separately; the combined matrix concatenated their trait blocks after matching genotype identifiers. Tissue-specific replicate reference statistics were retained in the combined matrix. Controls participated as genotype rows. No previous z-score matrices, cluster assignments or linkage matrices were used.

Each tree was cut to exactly the requested number of clusters: 7, 8, 9, 10, 11, 12. Cluster labels were assigned in dendrogram order and are local to each cut. The two colour views used the same tree because changing the origin with a fixed denominator preserves Euclidean distances. Matrix-distance equality was checked numerically. This check is an algebraic validation, not evidence of biological cluster stability.

Principal component analysis (PCA) was fitted to the complete scaled genotype matrix after column centring, using full singular-value decomposition. PC1–PC2 plots were coloured by the Ward assignments. All component scores, loadings and explained-variance fractions were exported. PCA was not used to select the clustering space. Silhouette, Calinski–Harabasz and Davies–Bouldin scores, cluster sizes and singleton counts were calculated descriptively for the requested cuts. No optimal k was asserted. No stability, trait-correlation or sensitivity analyses were performed, and no inferential p-values were calculated.

## Figure rendering and verification

Heatmaps used a continuous blue–white–red scale centred on zero with limits -2.5 and +2.5. Values beyond those limits saturated only in the rendering; numerical tables and clustering used untruncated scaled values. No white band or significance threshold was applied. Higher numerical values were red; for timing this means later threshold attainment. A single row permutation was applied and checked against the row labels. Both labelled and compact heatmaps were exported. All figures were saved as vector PDF/SVG and as JPEG/TIFF at 900 dpi. The rendering manifest records row/column order and plotting coordinates for independent PDF verification.

## Software and provenance

numpy 2.1.3, pandas 2.2.3, SciPy 1.15.3, scikit-learn 1.6.1, Matplotlib 3.10.0, Pillow 12.2.0, Python 3.13.5. Source hashes, analysis settings, exclusion logs, trait-level replicate counts, numerical matrices and figure manifests accompany the results. The source CSVs are the only research-data inputs.
