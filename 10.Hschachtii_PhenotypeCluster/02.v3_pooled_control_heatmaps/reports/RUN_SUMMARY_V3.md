# Pooled-control-centered heatmaps v3 — run summary

## Normalization

`z = (gene_median − pooled_control_median) / global_sd`

- Dendrograms and cluster labels are reused from the original Ward clustering (global-centered z-scores).
- Combined 16-D uses the 75-gene intersection SD (BM00068 is shoot-only and is excluded from that SD).
- White band: `|z| < 0.25`.
- `tt10` is not sign-flipped (higher = slower).
- Heatmaps are passed unpermuted; seaborn applies `row_linkage` once. Root k=3 landmarks (BM00079 top/red/C3; BM00084/88/53 bottom/blue/C1) are checked at runtime.

## Outputs in this snapshot

- 12 heatmaps × jpeg+pdf
- 5 diagnostic plots × jpeg+pdf
- 8 z-score / baseline tables
