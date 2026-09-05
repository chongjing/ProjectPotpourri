#!/usr/bin/env python3
"""Control-centered heatmaps v3 — pooled baseline, original dendrograms.

Key features:
- Pooled baseline: z = (gene_median - pooled_control_median) / global_SD
- All 3 controls near z=0 (white) in heatmap
- Original dendrograms reused (cluster assignments preserved)
- White-noise-band colormap (|z|<0.25 rendered white)
- No direction inversion (tt10 raw)
"""

import matplotlib
matplotlib.use("Agg")
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm
from pathlib import Path
from scipy.cluster.hierarchy import linkage
from scipy.spatial.distance import pdist

HERE = Path(__file__).resolve().parent.parent
HSCH = HERE.parent
OUT = HERE
INP = HERE / "inputs"
ROOT_PROJ = HSCH
for d in [OUT/"tables", OUT/"figures/main", OUT/"figures/diagnostics", OUT/"reports", INP]:
    d.mkdir(parents=True, exist_ok=True)

CONTROLS = ["Col-0", "EmptyVector", "GFP"]
NPHE = ["area_end_mm2","AUC","AGR_mean","AGR_max","RGR_mean","RGR_peak","RGR_final","tt10"]
CPAL = {1:"#4C78A8",2:"#F58518",3:"#54A24B",4:"#B279A2",5:"#E45756",6:"#72B7B2",7:"#EECA3B",8:"#9D755D"}

sns.set_theme(context="paper",style="white",font="DejaVu Sans",rc={
    "axes.edgecolor":"#222222","axes.linewidth":0.7,"axes.labelcolor":"#111111",
    "axes.titlesize":9,"axes.labelsize":8,"xtick.labelsize":7,"ytick.labelsize":7,
    "legend.fontsize":7,"legend.title_fontsize":7,"figure.dpi":120,"savefig.dpi":1200,
    "pdf.fonttype":42,"ps.fonttype":42,"svg.fonttype":"none",
})
plt.rcParams["axes.spines.top"]=False
plt.rcParams["axes.spines.right"]=False

def sf(fig, prefix):
    for ext, kw in [("pdf", {}), ("jpeg", {"pil_kwargs": {"quality": 95}})]:
        fig.savefig(prefix.with_suffix(f".{ext}"), bbox_inches="tight", dpi=1200, facecolor="white", **kw)
    plt.close(fig)

def scg(cg, prefix):
    sf(cg.fig, prefix)

# White-band colormap
def make_cmap():
    colors = [
        (0.0, "#2166AC"),   # z = -2.5 (deep blue)
        (0.4, "#4393C3"),   # z = -0.5 (light blue)
        (0.45, "#F7F7F7"),  # z = -0.25 (white)
        (0.55, "#F7F7F7"),  # z = +0.25 (white)
        (0.6, "#D6604D"),   # z = +0.5 (light red)
        (1.0, "#B2182B"),   # z = +2.5 (deep red)
    ]
    return LinearSegmentedColormap.from_list("ctrl_cmap", colors)

CTRL_CMAP = make_cmap()

def load_root():
    df = pd.read_csv(ROOT_PROJ / "phenotype_summary_for_clustering_with_AGRmax_and_days_to_10mm2.csv")
    df = df[df["outliers"] != True].copy()
    df = df.rename(columns={"dpi_to_10mm2":"tt10"})
    return df

def load_shoot():
    df = pd.read_csv(INP / "phenotype_summary_for_clustering_shoot.csv")
    df = df[df["area_end_mm2"].astype(float) > 0].copy()
    max_tt = df["time_to_10mm2"].max()
    mask = (df["area_end_mm2"].astype(float) > 0) & df["time_to_10mm2"].isna()
    df.loc[mask, "time_to_10mm2"] = max_tt + 1.0
    df = df.rename(columns={"time_to_10mm2":"tt10"})
    return df

def pooled_centered_zscore(df, tissue_label):
    """z = (gene_median - pooled_control_median) / global_SD. No direction inversion."""
    med = df.groupby("effector", sort=True)[NPHE].median()
    exclude = [g for g in med.index if med.loc[g].isna().any()]
    keep = [g for g in med.index if g not in exclude]
    
    # Pooled control median (all 3 controls combined)
    ctrl_df = df[df["effector"].isin(CONTROLS)]
    pooled_median = ctrl_df[NPHE].median()
    
    # Global SD from all gene medians
    global_sd = med.loc[keep, NPHE].std(ddof=0)
    
    # Pooled-centered z-score
    z = (med.loc[keep, NPHE] - pooled_median) / global_sd
    z.columns = [f"{tissue_label}_{c}" for c in NPHE]
    
    # Per-control z-profiles (for reference tracks)
    ctrl_z = {}
    for ctrl in CONTROLS:
        ctrl_df = df[df["effector"] == ctrl]
        ctrl_med = ctrl_df[NPHE].median()
        ctrl_z[ctrl] = (ctrl_med - pooled_median) / global_sd
    
    return z, keep, exclude, pooled_median, global_sd, ctrl_z

def load_cluster_labels(csv_path):
    df = pd.read_csv(csv_path, index_col=0)
    return df["cluster"]

def assert_index_aligned(z_matrix, orig_matrix, name):
    if list(z_matrix.index) != list(orig_matrix.index):
        raise RuntimeError(f"{name}: pooled z-score row order != original clustering matrix")

def make_heatmap(z_matrix, Z_linkage, labels, paths, prefix, labelled=True, K=None):
    """Pass the unpermuted matrix. Seaborn applies row_linkage once.

    Pre-ordering rows and then passing the original linkage double-permutes
    the heatmap cells relative to y-tick labels (v3 bug).
    """
    suffix = f"_k{K}" if K else ""
    labels = labels.reindex(z_matrix.index)
    if labels.isna().any():
        missing = labels[labels.isna()].index.tolist()
        raise RuntimeError(f"{prefix}: missing cluster labels for {missing}")
    row_colors = labels.map(CPAL)
    n = z_matrix.shape[1]
    if n == 8:
        col_colors = ["#1B5299"] * 8
    elif n == 16:
        col_colors = ["#1B5299"] * 8 + ["#E45756"] * 8
    else:
        col_colors = None
    figsize = (8.5, 11.2) if labelled else (6.5, 8.2)
    norm = TwoSlopeNorm(vmin=-2.5, vcenter=0, vmax=2.5)
    cg = sns.clustermap(
        z_matrix, row_linkage=Z_linkage, col_cluster=False,
        cmap=CTRL_CMAP, norm=norm,
        row_colors=row_colors, col_colors=col_colors,
        linewidths=0.0, figsize=figsize,
        dendrogram_ratio=(0.12, 0.10), colors_ratio=0.025,
        cbar_pos=(0.02, 0.82, 0.02, 0.12),
        cbar_kws={"label":"Pooled-ctrl z-score"},
        yticklabels=False,
    )
    dendro_order = z_matrix.index[cg.dendrogram_row.reordered_ind].tolist()
    if list(cg.data2d.index) != dendro_order:
        raise RuntimeError(f"{prefix}: heatmap row order != dendrogram leaves")
    cg.ax_heatmap.set_xlabel("")
    cg.ax_heatmap.set_ylabel("")
    cg.ax_heatmap.tick_params(axis="x", rotation=35, labelsize=7)
    if labelled:
        n_genes = len(dendro_order)
        cg.ax_heatmap.set_yticks(np.arange(n_genes) + 0.5)
        cg.ax_heatmap.set_yticklabels(dendro_order, fontsize=5)
    n_cl = int(pd.Series(labels).nunique())
    title_k = f"k={n_cl}" if not K else f"k={K}"
    handles = [mpatches.Patch(color=CPAL[c], label=f"C{c} ({(labels==c).sum()})") for c in sorted(set(labels))]
    cg.ax_heatmap.legend(handles=handles, title=title_k, bbox_to_anchor=(1.02,1), loc="upper left", frameon=False)
    scg(cg, paths / (prefix + suffix))
    return dendro_order, cg.data2d

def check_root_k3(dendro_order, data2d, labels):
    """The three visual landmarks for the recommended root heatmap."""
    errors = []
    if dendro_order[0] != "BM00079":
        errors.append(f"top leaf is {dendro_order[0]}, expected BM00079")
    if dendro_order[-3:] != ["BM00084", "BM00088", "BM00053"]:
        errors.append(f"bottom leaves are {dendro_order[-3:]}, expected BM00084/BM00088/BM00053")
    z79 = data2d.loc["BM00079"]
    if z79["R_area_end_mm2"] <= 1.0 or z79["R_AUC"] <= 1.0:
        errors.append(f"BM00079 growth z is {z79[['R_area_end_mm2','R_AUC']].to_dict()}, expected red (>1)")
    if int(labels.loc["BM00079"]) != 3:
        errors.append(f"BM00079 cluster is {labels.loc['BM00079']}, expected C3")
    if int(labels.loc[dendro_order[0]]) != 3:
        errors.append(f"top-row cluster colour is C{int(labels.loc[dendro_order[0]])}, expected C3")
    for g in dendro_order[-3:]:
        if data2d.loc[g, "R_area_end_mm2"] >= -3.0:
            errors.append(f"{g} area z is {data2d.loc[g,'R_area_end_mm2']:.2f}, expected saturated blue (<-3)")
        if int(labels.loc[g]) != 1:
            errors.append(f"{g} cluster is C{int(labels.loc[g])}, expected C1")
    for g in CONTROLS:
        if int(labels.loc[g]) != 3:
            errors.append(f"{g} cluster is C{int(labels.loc[g])}, expected C3")
    if errors:
        raise RuntimeError("root k=3 landmark checks failed:\n  - " + "\n  - ".join(errors))
    print("root k=3 landmark checks: PASS")
    print(f"  top={dendro_order[0]} area_z={z79['R_area_end_mm2']:.3f} cluster=C{int(labels.loc['BM00079'])}")
    print(f"  bottom={dendro_order[-3:]}")

def make_reference_track(ctrl_z_dict, paths, prefix, tissue_label):
    fig, ax = plt.subplots(figsize=(8.5, 1.5))
    ctrl_colors = {"Col-0":"#888888","EmptyVector":"#333333","GFP":"#BBBBBB"}
    markers = {"Col-0":"o","EmptyVector":"s","GFP":"^"}
    x = np.arange(len(NPHE))
    for ctrl in CONTROLS:
        vals = [ctrl_z_dict[ctrl][phe] for phe in NPHE]
        ax.plot(x, vals, marker=markers[ctrl], color=ctrl_colors[ctrl],
                label=ctrl, linewidth=1.5, markersize=6, zorder=3)
    ax.axhline(0, color="#333333", linewidth=0.8, linestyle="--", alpha=0.5)
    ax.set_xticks(x)
    ax.set_xticklabels(NPHE, rotation=35, ha="right", fontsize=7)
    ax.set_ylabel("Z-score", fontsize=8)
    ax.set_title(f"{tissue_label}: Control z-profiles (pooled baseline)", fontsize=9)
    ax.legend(fontsize=7, frameon=False, loc="upper right")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.set_ylim(-1.5, 1.5)
    sf(fig, paths / prefix)

def make_zdist_plot(z_matrix, paths, prefix, tissue_label):
    ncols = z_matrix.shape[1]
    nrows = (ncols + 3) // 4
    fig, axes = plt.subplots(nrows, 4, figsize=(10, 2.5*nrows))
    axes = axes.flatten()
    for i, phe in enumerate(z_matrix.columns):
        ax = axes[i]
        vals = z_matrix[phe].drop(CONTROLS, errors="ignore")
        ax.hist(vals, bins=20, color="#4C78A8", alpha=0.7, edgecolor="white")
        ax.axvline(-0.25, color="#E45756", linewidth=0.8, linestyle="--", alpha=0.7)
        ax.axvline(0.25, color="#E45756", linewidth=0.8, linestyle="--", alpha=0.7)
        ax.axvline(0, color="#333333", linewidth=0.8, linestyle="-", alpha=0.5)
        ax.set_title(phe, fontsize=8)
        ax.set_xlabel("z", fontsize=7)
        ax.set_ylabel("N", fontsize=7)
        ax.tick_params(labelsize=6)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
    fig.suptitle(f"{tissue_label}: Effector z-distributions (white band = noise |z|<0.25)", fontsize=9)
    fig.tight_layout()
    sf(fig, paths / prefix)

print("Pooled-control-centered heatmap generation v3")
print("=" * 60)

# Load data
root_df = load_root()
shoot_df = load_shoot()
print(f"Root: {len(root_df)} rows, {root_df.effector.nunique()} genes")
print(f"Shoot: {len(shoot_df)} rows, {shoot_df.effector.nunique()} genes")

# Compute pooled-centered z-scores
r_z, r_keep, r_excl, r_pooled, r_sd, r_ctrl_z = pooled_centered_zscore(root_df, "R")
s_z, s_keep, s_excl, s_pooled, s_sd, s_ctrl_z = pooled_centered_zscore(shoot_df, "S")
print(f"\nRoot: {len(r_keep)} genes (excluded: {r_excl})")
print(f"Shoot: {len(s_keep)} genes (excluded: {s_excl})")

# Load original dendrograms (linkage leaf i = row i of these matrices)
r_orig = pd.read_csv(ROOT_PROJ / "01.clustering_analysis/tables/gene_level_median_zscore.csv", index_col=0)
s_orig = pd.read_csv(INP / "shoot_gene_level_median_zscore.csv", index_col=0)
c_orig = pd.read_csv(INP / "combined_gene_level_median_zscore_16d.csv", index_col=0)
assert_index_aligned(r_z, r_orig, "root")
assert_index_aligned(s_z, s_orig, "shoot")
r_Z = linkage(pdist(r_orig.values, metric="euclidean"), method="ward", optimal_ordering=True)
s_Z = linkage(pdist(s_orig.values, metric="euclidean"), method="ward", optimal_ordering=True)
c_Z = linkage(pdist(c_orig.values, metric="euclidean"), method="ward", optimal_ordering=True)

# Load cluster assignments
r_labels_k3 = load_cluster_labels(ROOT_PROJ / "01.clustering_analysis/tables/gene_cluster_assignments_recommended_relabelled.csv")
r_labels_k8 = load_cluster_labels(ROOT_PROJ / "01.clustering_analysis/tables/gene_cluster_assignments_k8_relabelled.csv")
s_labels_k2 = load_cluster_labels(INP / "shoot_gene_cluster_assignments_recommended_relabelled.csv")
s_labels_k8 = load_cluster_labels(INP / "shoot_gene_cluster_assignments_k8_relabelled.csv")
c_labels_k2 = load_cluster_labels(INP / "combined_gene_cluster_assignments_recommended_relabelled.csv")
c_labels_k8 = load_cluster_labels(INP / "combined_gene_cluster_assignments_k8_relabelled.csv")

r_labels_k3 = r_labels_k3.reindex(r_z.index)
r_labels_k8 = r_labels_k8.reindex(r_z.index)
s_labels_k2 = s_labels_k2.reindex(s_z.index)
s_labels_k8 = s_labels_k8.reindex(s_z.index)

# ROOT
print("\n--- ROOT ---")
root_order, root_data2d = make_heatmap(r_z, r_Z, r_labels_k3, OUT/"figures/main", "root_ctrl_heatmap_labelled", labelled=True)
check_root_k3(root_order, root_data2d, r_labels_k3)
make_heatmap(r_z, r_Z, r_labels_k3, OUT/"figures/main", "root_ctrl_heatmap_compact", labelled=False)
make_heatmap(r_z, r_Z, r_labels_k8, OUT/"figures/main", "root_ctrl_heatmap_labelled", labelled=True, K=8)
make_heatmap(r_z, r_Z, r_labels_k8, OUT/"figures/main", "root_ctrl_heatmap_compact", labelled=False, K=8)
make_reference_track(r_ctrl_z, OUT/"figures/diagnostics", "root_control_reference", "Root")
make_zdist_plot(r_z, OUT/"figures/diagnostics", "root_zdist", "Root")

# SHOOT
print("\n--- SHOOT ---")
make_heatmap(s_z, s_Z, s_labels_k2, OUT/"figures/main", "shoot_ctrl_heatmap_labelled", labelled=True)
make_heatmap(s_z, s_Z, s_labels_k2, OUT/"figures/main", "shoot_ctrl_heatmap_compact", labelled=False)
make_heatmap(s_z, s_Z, s_labels_k8, OUT/"figures/main", "shoot_ctrl_heatmap_labelled", labelled=True, K=8)
make_heatmap(s_z, s_Z, s_labels_k8, OUT/"figures/main", "shoot_ctrl_heatmap_compact", labelled=False, K=8)
make_reference_track(s_ctrl_z, OUT/"figures/diagnostics", "shoot_control_reference", "Shoot")
make_zdist_plot(s_z, OUT/"figures/diagnostics", "shoot_zdist", "Shoot")

# COMBINED 16-D
# Recentre shoot on the 75-gene intersection so the plotted matrix is a
# constant shift of the original 16-D clustering matrix (BM00068 is
# shoot-only and must not enter the combined SD).
print("\n--- COMBINED ---")
common = list(c_orig.index)
if sorted(set(r_z.index) & set(s_z.index)) != sorted(common):
    raise RuntimeError("combined gene set != original 16-D clustering matrix")
r_med = root_df.groupby("effector", sort=True)[NPHE].median().loc[common]
s_med = shoot_df.groupby("effector", sort=True)[NPHE].median().loc[common]
r_z_c = (r_med - r_pooled) / r_med.std(ddof=0)
r_z_c.columns = [f"R_{c}" for c in NPHE]
s_z_c = (s_med - s_pooled) / s_med.std(ddof=0)
s_z_c.columns = [f"S_{c}" for c in NPHE]
combined_16d = pd.concat([r_z_c, s_z_c], axis=1)
assert_index_aligned(combined_16d, c_orig, "combined")
shift_col_std = (combined_16d.values - c_orig.values).std(axis=0)
if np.max(shift_col_std) > 1e-10:
    raise RuntimeError(f"combined is not a constant shift of original 16-D; col std={shift_col_std}")
print(f"  Combined: {combined_16d.shape[0]} genes; max shift-std vs original={shift_col_std.max():.2e}")
c_labels_k2_c = c_labels_k2.reindex(combined_16d.index)
c_labels_k8_c = c_labels_k8.reindex(combined_16d.index)
make_heatmap(combined_16d, c_Z, c_labels_k2_c, OUT/"figures/main", "combined_ctrl_heatmap_labelled", labelled=True)
make_heatmap(combined_16d, c_Z, c_labels_k2_c, OUT/"figures/main", "combined_ctrl_heatmap_compact", labelled=False)
make_heatmap(combined_16d, c_Z, c_labels_k8_c, OUT/"figures/main", "combined_ctrl_heatmap_labelled", labelled=True, K=8)
make_heatmap(combined_16d, c_Z, c_labels_k8_c, OUT/"figures/main", "combined_ctrl_heatmap_compact", labelled=False, K=8)
make_zdist_plot(combined_16d, OUT/"figures/diagnostics", "combined_zdist", "Combined")

# Save tables
r_z.to_csv(OUT/"tables"/"root_pooled_centered_zscore.csv", index_label="effector")
s_z.to_csv(OUT/"tables"/"shoot_pooled_centered_zscore.csv", index_label="effector")
combined_16d.to_csv(OUT/"tables"/"combined_pooled_centered_zscore_16d.csv", index_label="effector")
r_pooled.to_frame("value").to_csv(OUT/"tables"/"root_pooled_median.csv", index_label="phenotype")
s_pooled.to_frame("value").to_csv(OUT/"tables"/"shoot_pooled_median.csv", index_label="phenotype")
r_sd.to_frame("value").to_csv(OUT/"tables"/"root_global_sd.csv", index_label="phenotype")
s_sd.to_frame("value").to_csv(OUT/"tables"/"shoot_global_sd.csv", index_label="phenotype")
ctrl_ref = pd.DataFrame({"Root_Col-0": r_ctrl_z["Col-0"], "Root_EV": r_ctrl_z["EmptyVector"],
                          "Root_GFP": r_ctrl_z["GFP"], "Shoot_Col-0": s_ctrl_z["Col-0"],
                          "Shoot_EV": s_ctrl_z["EmptyVector"], "Shoot_GFP": s_ctrl_z["GFP"]})
ctrl_ref.to_csv(OUT/"tables"/"control_reference_profiles.csv", index_label="phenotype")

summary = """# Pooled-Control-Centered Heatmaps v3 - Run Summary

## Normalization: z = (gene_median - pooled_control_median) / global_sd
- Original dendrograms reused (cluster assignments preserved)
- White-noise-band colormap (|z|<0.25 rendered white)
- No direction inversion (tt10 raw)
- Combined 16-D uses the 75-gene intersection SD (not the 76-gene shoot SD)

## Plot fix
- Heatmaps are passed unpermuted; seaborn applies `row_linkage` once.
- Y-tick labels are taken from `dendrogram_row.reordered_ind` (same order as cells).
- Root k=3 landmark checks (BM00079 top/red/C3; low trio bottom/blue/C1) must pass or the run aborts.

## Outputs:
- 12 heatmaps (jpeg + pdf)
- 5 diagnostic plots (jpeg + pdf)
- 8 tables (z-scores, pooled medians, global SDs, control profiles)
"""
(OUT/"reports"/"RUN_SUMMARY_V3.md").write_text(summary)

print("\n" + "=" * 60)
print(f"COMPLETE. Output: {OUT}")
print("=" * 60)
