#!/usr/bin/env python3
"""Final-sheet-only phenotype clustering. Source CSVs are never modified.

Biological root/shoot pairing is unconfirmed. Identifier matching is used only
for the explicitly approved joint-zero exclusion, not for paired inference.
"""
from pathlib import Path
import numpy as np
import pandas as pd

KEYS = ['effector', 'replicate']
CONTROLS = ['Col-0', 'EmptyVector', 'GFP']
TRAITS = ['area_end_mm2', 'AUC', 'time_to_10mm2_capped', 'AGR_mean',
          'AGR_max', 'RGR_mean', 'RGR_peak', 'RGR_final']


def filter_final_inputs(frames):
    """Honour final-sheet exclusions and remove identifier-matched joint zeros.

    A missing counterpart or missing final size never satisfies the zero test.
    Nonzero AUC does NOT rescue a pair whose final sizes are both zero: the
    user explicitly selected final size as the criterion for this analysis.
    """
    for tissue, frame in frames.items():
        if frame[KEYS].isna().any().any() or frame.duplicated(KEYS).any():
            raise ValueError(f'{tissue}: missing or duplicate identifiers')
        if not frame['included'].isin([True, False]).all():
            raise ValueError(f'{tissue}: included must contain explicit booleans')
    zeros = {t: set(map(tuple, d.loc[d.area_end_mm2.eq(0), KEYS].to_numpy()))
             for t, d in frames.items()}
    joint = zeros['root'] & zeros['shoot']
    logs, cleaned = [], {}
    for tissue, frame in frames.items():
        other_name = 'shoot' if tissue == 'root' else 'root'
        other = frames[other_name].set_index(KEYS)['area_end_mm2']
        matches = pd.Series([tuple(x) in joint for x in frame[KEYS].to_numpy()],
                            index=frame.index)
        retained = frame.included.eq(True) & ~matches
        out = frame[KEYS + ['included', 'exclusion_reason', 'area_end_mm2']].copy()
        out.insert(0, 'tissue', tissue)
        out.insert(1, 'source_csv_line', np.arange(len(frame)) + 2)
        out['counterpart_identifier_found'] = [tuple(x) in other.index
                                               for x in frame[KEYS].to_numpy()]
        out['counterpart_final_size'] = [other.get(tuple(x), np.nan)
                                         for x in frame[KEYS].to_numpy()]
        out['both_final_sizes_zero'] = matches
        out['retained'] = retained
        out['disposition'] = np.where(~frame.included.eq(True), 'excluded_by_final_sheet',
                                     np.where(matches, 'excluded_joint_zero', 'retained'))
        logs.append(out)
        cleaned[tissue] = frame.loc[retained].copy()
    return cleaned, pd.concat(logs, ignore_index=True)


def derive_timing(frame, cap=29.):
    """Encode the approved bounded timing feature; retain all raw fields.

    Blank T10 on an included row was explicitly confirmed by Dr. Xia to mean
    no crossing observed. Right-censored follow-up must reach the chosen cap.
    A supplied left-censored bound is retained and flagged, never presented as
    an exact crossing. RGR blanks are not touched by this transformation.
    """
    out = frame.copy()
    time = out['time_to_10mm2']
    right = time.isna()
    left = out['time_to_10mm2_left_censored'].eq(True)
    if (right & (out.end_dpg.isna() | out.end_dpg.lt(cap))).any():
        raise ValueError('Right-censored follow-up ends before the timing cap')
    if (right & left).any():
        raise ValueError('A left-censored record lacks its supplied upper bound')
    if np.isinf(time.to_numpy(float)).any():
        raise ValueError('Non-finite observed threshold time')
    out['time_to_10mm2_capped'] = time.clip(upper=cap).fillna(cap)
    out['time_to_10mm2_status'] = np.select(
        [left, right, time.gt(cap), time.eq(cap)],
        ['left_censored_upper_bound', 'right_censored', 'observed_after_cap', 'observed_at_cap'],
        default='observed')
    out['time_to_10mm2_cap_dpg'] = cap
    return out


def summarise_tissue(frame, prefix):
    """Genotype medians; reference moments use retained individual values.

    Missing RGR values remain missing. Each reference trait is fitted on all
    available retained replicate measurements (including controls), ddof=0.
    Controls define only the alternative origin, not an alternative scale.
    """
    values = frame[TRAITS].astype(float)
    if np.isinf(values.to_numpy()).any():
        raise ValueError('Infinite phenotype measurement')
    means, sd = values.mean(), values.std(ddof=0)
    if sd.isna().any() or sd.le(0).any():
        raise ValueError('Trait has no usable positive reference SD')
    med = frame.groupby('effector', sort=True)[TRAITS].median()
    counts = frame.groupby('effector', sort=True)[TRAITS].count()
    if med.isna().any().any():
        missing = med.index[med.isna().any(axis=1)].tolist()
        raise ValueError(f'Incomplete genotype summaries require a decision: {missing}')
    if not set(CONTROLS).issubset(med.index):
        raise ValueError('A required control genotype is absent')
    ctrl_values = frame.loc[frame.effector.isin(CONTROLS), TRAITS]
    ctrl = ctrl_values.median()
    if ctrl.isna().any():
        raise ValueError('A pooled-control reference is missing')
    global_z, control_z = (med - means) / sd, (med - ctrl) / sd
    ref = pd.DataFrame({'all_replicate_mean': means, 'all_replicate_sd': sd,
                        'pooled_control_median': ctrl,
                        'n_retained_values': values.count(),
                        'n_pooled_control_values': ctrl_values.count()})
    names = [f'{prefix}_{t}' for t in TRAITS]
    for obj in (med, counts, global_z, control_z):
        obj.columns = names
        obj.index.name = 'effector'
    ref.index = names
    ref.index.name = 'trait'
    return {'medians': med, 'counts': counts, 'global_z': global_z,
            'control_z': control_z, 'reference': ref}


def analyse_matrix(global_z, control_z, ks):
    """A fresh Euclidean/Ward tree, exact-k cuts and descriptive PCA.

    Both colour matrices have the same fixed trait scaling. Recentring therefore
    cannot change distances. This algebraic check is QA, not a stability study.
    No bootstrap, trait-correlation, sensitivity or paired-inference tests run.
    """
    from scipy.cluster.hierarchy import linkage, leaves_list, cut_tree
    from scipy.spatial.distance import pdist
    from sklearn.decomposition import PCA
    from sklearn.metrics import silhouette_score, calinski_harabasz_score, davies_bouldin_score
    if not global_z.index.equals(control_z.index) or not global_z.columns.equals(control_z.columns):
        raise ValueError('Colour matrix identifiers are not aligned')
    if not np.isfinite(global_z.to_numpy()).all() or not np.isfinite(control_z.to_numpy()).all():
        raise ValueError('Clustering requires complete finite genotype profiles')
    if any(k < 2 or k >= len(global_z) for k in ks):
        raise ValueError('Requested k must be between 2 and n_genotypes - 1')
    dist = pdist(global_z.to_numpy(), 'euclidean')
    dist2 = pdist(control_z.to_numpy(), 'euclidean')
    np.testing.assert_allclose(dist, dist2, rtol=1e-12, atol=1e-12)
    tree = linkage(dist, method='ward', optimal_ordering=True)
    positions = leaves_list(tree)
    order = global_z.index[positions].tolist()
    assignment = pd.DataFrame(index=global_z.index)
    rows = []
    for k in ks:
        raw = cut_tree(tree, n_clusters=[k]).ravel()
        remap = {old: new+1 for new, old in enumerate(dict.fromkeys(raw[positions]))}
        labels = np.array([remap[x] for x in raw])
        assert len(np.unique(labels)) == k
        assignment[f'k{k}'] = labels
        sizes = pd.Series(labels).value_counts()
        rows.append({'k': k, 'silhouette': silhouette_score(global_z, labels),
                     'calinski_harabasz': calinski_harabasz_score(global_z, labels),
                     'davies_bouldin': davies_bouldin_score(global_z, labels),
                     'smallest_cluster': int(sizes.min()), 'largest_cluster': int(sizes.max()),
                     'singleton_clusters': int(sizes.eq(1).sum())})
    pca = PCA(svd_solver='full')
    coords = pca.fit_transform(global_z.to_numpy())
    pc = [f'PC{i+1}' for i in range(coords.shape[1])]
    scores = pd.DataFrame(coords, index=global_z.index, columns=pc)
    loadings = pd.DataFrame(pca.components_.T, index=global_z.columns, columns=pc)
    # Compare projection through the same fitted basis after each matrix's own
    # column centring: avoids arbitrary independent-PC sign conventions.
    projected = (control_z.to_numpy()-control_z.to_numpy().mean(axis=0)) @ pca.components_.T
    np.testing.assert_allclose(coords, projected, atol=1e-10, rtol=1e-10)
    variance = pd.DataFrame({'PC':pc, 'explained_variance_ratio':pca.explained_variance_ratio_,
                             'eigenvalue':pca.explained_variance_})
    for obj in (assignment, scores):
        obj.index.name = 'effector'
    loadings.index.name = 'trait'
    return {'tree':tree, 'leaf_order':order, 'assignments':assignment,
            'metrics':pd.DataFrame(rows), 'scores':scores, 'loadings':loadings,
            'variance':variance,
            'validation':{'max_distance_difference':float(np.max(np.abs(dist-dist2))),
                          'max_same_basis_pca_difference':float(np.max(np.abs(coords-projected)))}}


PALETTE = ['#4477AA','#EE6677','#228833','#CCBB44','#66CCEE','#AA3377',
           '#BBBBBB','#332288','#44AA99','#999933','#882255','#AA4499']
DISPLAY = {'area_end_mm2':'Final size', 'AUC':'AUC',
           'time_to_10mm2_capped':'T10 (cap 29 dpg)', 'AGR_mean':'Mean AGR',
           'AGR_max':'Max AGR', 'RGR_mean':'Mean RGR', 'RGR_peak':'Peak RGR',
           'RGR_final':'Final RGR'}


def plotting_setup():
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams.update({'font.family':'sans-serif',
                         'font.sans-serif':['Liberation Sans','DejaVu Sans'],
                         'font.size':8, 'pdf.fonttype':42, 'svg.fonttype':'none',
                         'axes.spines.top':False,'axes.spines.right':False})
    return plt


def heatmap_colours():
    from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm
    cmap = LinearSegmentedColormap.from_list('blue_white_red',
                                             ['#2166AC','#FFFFFF','#B2182B'],N=257)
    cmap.set_bad('#BDBDBD')
    norm = TwoSlopeNorm(vmin=-2.5,vcenter=0.,vmax=2.5)
    return cmap,norm


def make_heatmap(matrix, result, k, view, baseline, scope, labelled=True):
    """Apply the leaf permutation exactly once; independent baseline strip.

    Pcolormesh retains vector cells in PDF/SVG. Metadata stores the actual
    plotting rectangle for independent PDF pixel and row-label verification.
    """
    plt = plotting_setup()
    from matplotlib.colors import ListedColormap, BoundaryNorm
    from matplotlib.patches import Patch
    from scipy.cluster.hierarchy import dendrogram
    order = result['leaf_order']
    data = matrix.loc[order]
    groups = result['assignments'].loc[order,f'k{k}']
    n,p = data.shape
    fig = plt.figure(figsize=(10.4 if p==16 else 8.6,11.4) if labelled
                     else (8.7 if p==16 else 7.0,8.8),dpi=100)
    right = .83 if labelled else .91
    x0,y0,w,h = .16,.19,right-.16,.635
    ax = fig.add_axes([x0,y0,w,h],label='heatmap')
    dend = fig.add_axes([.025,y0,.105,h],label='dendrogram')
    strip = fig.add_axes([.137,y0,.012,h],label='cluster_strip')
    base_ax = fig.add_axes([x0,.85,w,.017],label='pooled_reference')
    cmap,norm = heatmap_colours()
    mesh = ax.pcolormesh(data.values,cmap=cmap,norm=norm,shading='flat',rasterized=False)
    ax.set_xlim(0,p); ax.set_ylim(n,0)
    ax.set_xticks(np.arange(p)+.5)
    names = [f'{c[:1]}: {DISPLAY[c[2:]]}' for c in matrix.columns]
    ax.set_xticklabels(names,rotation=45,ha='right',fontsize=7)
    ax.tick_params(axis='x',length=0,pad=3)
    ax.set_yticks(np.arange(n)+.5 if labelled else [])
    if labelled:
        ax.set_yticklabels(order,fontsize=5.7)
        ax.yaxis.tick_right()
        ax.tick_params(axis='y',length=0,pad=3)
        for text in ax.get_yticklabels():
            if text.get_text() in CONTROLS:
                text.set_fontweight('bold')
    for spine in ax.spines.values(): spine.set_visible(False)
    dendrogram(result['tree'],orientation='left',no_labels=True,
               color_threshold=0,above_threshold_color='#555555',ax=dend)
    dend.set_ylim(10*n,0); dend.set_axis_off()
    cnorm = BoundaryNorm(np.arange(.5,13.5),12)
    strip.pcolormesh(groups.to_numpy()[:,None],cmap=ListedColormap(PALETTE),norm=cnorm)
    strip.set_ylim(n,0); strip.set_axis_off()
    base_ax.pcolormesh(baseline.reindex(matrix.columns).values[None,:],cmap=cmap,norm=norm,
                       edgecolors='#BBBBBB',linewidth=.35)
    base_ax.set_xlim(0,p); base_ax.set_ylim(1,0); base_ax.set_axis_off()
    fig.text(x0,.875,'Pooled-control reference (outside clustering)',fontsize=7,ha='left')
    if p==16:
        ax.axvline(8,color='black',lw=.7)
        base_ax.axvline(8,color='black',lw=.7)
    for boundary in np.flatnonzero(np.diff(groups.to_numpy()))+1:
        ax.axhline(boundary,color='black',lw=.45,alpha=.7)
    sizes=groups.value_counts()
    handles=[Patch(facecolor=PALETTE[c-1],label=f'C{c} (n={sizes[c]})') for c in range(1,k+1)]
    fig.legend(handles=handles,loc='upper center',bbox_to_anchor=(.53,.948),
               ncol=6,fontsize=7,frameon=False,handlelength=1.0,columnspacing=1.0)
    centring = 'all-replicate mean' if view=='global' else 'pooled-control median'
    fig.suptitle(f'{scope.capitalize()} | k={k} | centred on {centring}',y=.983,fontsize=11)
    cb_ax=fig.add_axes([.34,.064,.38,.013])
    cb=fig.colorbar(mesh,cax=cb_ax,orientation='horizontal',extend='both')
    cb.set_ticks([-2.5,-1.25,0,1.25,2.5]); cb.ax.tick_params(labelsize=7)
    cb.set_label('Deviation / all-retained-replicate SD',fontsize=8)
    note = 'Controls are bold.' if labelled else 'See labelled version for genotype IDs.'
    fig.text(.04,.017,f'T10 is capped at 29 dpg; higher T10 = later growth. {note}\n'
             'Continuous colours; no significance threshold. Cluster numbers are local to this k.',fontsize=7)
    fig.canvas.draw()
    np.testing.assert_allclose(np.asarray(mesh.get_array()),data.values)
    if labelled:
        assert [t.get_text() for t in ax.get_yticklabels()] == order
    meta={'kind':'heatmap','scope':scope,'k':int(k),'view':view,'labelled':bool(labelled),
          'row_order':order,'trait_order':matrix.columns.tolist(),
          'heatmap_bounds':list(ax.get_position().bounds),
          'reference_bounds':list(base_ax.get_position().bounds),
          'baseline_values':baseline.reindex(matrix.columns).tolist(),
          'figure_inches':fig.get_size_inches().tolist(),
          'vmin':-2.5,'vcenter':0.,'vmax':2.5,
          'clipped_cells':int(np.abs(data.to_numpy()).__gt__(2.5).sum())}
    return fig,meta


def save_figure(fig, stem, dpi=900):
    """Vector PDF/SVG plus one high-resolution raster render for JPEG/TIFF."""
    import gc
    from PIL import Image
    import matplotlib.pyplot as plt
    stem=Path(stem); stem.parent.mkdir(parents=True,exist_ok=True)
    files={ext:str(stem.with_suffix('.'+ext)) for ext in ('pdf','jpeg','svg','tiff')}
    fig.savefig(files['pdf'],facecolor='white')
    fig.savefig(files['svg'],facecolor='white')
    original_dpi=fig.dpi
    fig.set_dpi(dpi); fig.canvas.draw()
    image=Image.fromarray(np.asarray(fig.canvas.buffer_rgba())).convert('RGB')
    image.save(files['jpeg'],format='JPEG',quality=95,subsampling=0,dpi=(dpi,dpi))
    image.save(files['tiff'],format='TIFF',compression='tiff_lzw',dpi=(dpi,dpi))
    image.close(); fig.set_dpi(original_dpi); plt.close(fig); gc.collect()
    for path in files.values():
        if Path(path).stat().st_size==0: raise RuntimeError(f'Empty figure: {path}')
    return files


def make_pca_plot(result,k,scope):
    plt=plotting_setup()
    from matplotlib.lines import Line2D
    fig=plt.figure(figsize=(7.4,5.6),dpi=100)
    ax=fig.add_axes([.12,.16,.62,.70])
    scores=result['scores']; labels=result['assignments'][f'k{k}']
    for cluster in range(1,k+1):
        ids=labels.index[labels.eq(cluster)]
        ordinary=[g for g in ids if g not in CONTROLS]
        controls=[g for g in ids if g in CONTROLS]
        if ordinary:
            ax.scatter(scores.loc[ordinary,'PC1'],scores.loc[ordinary,'PC2'],
                       s=24,c=PALETTE[cluster-1],edgecolor='white',lw=.35)
        if controls:
            ax.scatter(scores.loc[controls,'PC1'],scores.loc[controls,'PC2'],
                       s=62,c=PALETTE[cluster-1],marker='D',edgecolor='black',lw=.7,zorder=4)
    offsets={'Col-0':(10,14),'EmptyVector':(-65,16),'GFP':(10,-16)}
    for gene in CONTROLS:
        if gene in scores.index:
            ax.annotate(gene,(scores.loc[gene,'PC1'],scores.loc[gene,'PC2']),
                        xytext=offsets[gene],textcoords='offset points',fontsize=7,
                        bbox={'facecolor':'white','edgecolor':'none','alpha':.85,'pad':.8},
                        arrowprops={'arrowstyle':'-','color':'#666666','lw':.5})
    var=result['variance'].explained_variance_ratio.to_numpy()*100
    ax.set_xlabel(f'PC1 ({var[0]:.1f}% variance)')
    ax.set_ylabel(f'PC2 ({var[1]:.1f}% variance)')
    ax.margins(.17)
    sizes=labels.value_counts()
    handles=[Line2D([],[],marker='o',ls='',color=PALETTE[c-1],
                     label=f'C{c} (n={sizes[c]})',markersize=5) for c in range(1,k+1)]
    ax.legend(handles=handles,bbox_to_anchor=(1.025,1),loc='upper left',frameon=False,fontsize=8)
    fig.suptitle(f'{scope.capitalize()} | PCA coloured by Ward clusters | k={k}',y=.96,fontsize=11)
    fig.text(.08,.045,'Diamonds: controls. PCA is a projection, not the clustering method.\n'
             'T10 capped at 29 dpg. Both heatmap centre choices share these PCA coordinates.',fontsize=7)
    return fig


def make_comparison(result,scope):
    plt=plotting_setup()
    fig,axes=plt.subplots(1,3,figsize=(10,3.5),dpi=100)
    df=result['metrics']
    axes[0].plot(df.k,df.silhouette,'o-',color='#4477AA')
    axes[0].set_ylabel('Mean silhouette'); axes[0].set_title('Separation (descriptive)',fontsize=9)
    axes[1].plot(df.k,df.singleton_clusters,'o-',color='#AA3377')
    axes[1].set_ylabel('Singleton clusters'); axes[1].set_title('Single-genotype groups',fontsize=9)
    for k in df.k:
        sizes=result['assignments'][f'k{k}'].value_counts().sort_index()
        axes[2].scatter(np.repeat(k,len(sizes)),sizes,s=18,c=PALETTE[:len(sizes)],alpha=.9)
    axes[2].set_ylabel('Genotypes per cluster'); axes[2].set_title('Cluster sizes',fontsize=9)
    for ax in axes:
        ax.set_xticks(df.k); ax.set_xlabel('Number of clusters (k)')
    fig.suptitle(f'{scope.capitalize()} | requested k comparison',fontsize=11)
    fig.subplots_adjust(left=.07,right=.97,bottom=.23,top=.80,wspace=.46)
    fig.text(.07,.035,'Exploratory comparison only. No stability analysis or statistical significance tests.',fontsize=7)
    return fig


def read_final_csv(path):
    frame=pd.read_csv(path,dtype={'effector':str,'replicate':str})
    source_traits=['time_to_10mm2' if t=='time_to_10mm2_capped' else t for t in TRAITS]
    required=KEYS+source_traits+['included','exclusion_reason','end_dpg',
                                 'time_to_10mm2_left_censored','is_control']
    absent=set(required)-set(frame.columns)
    if absent: raise ValueError(f'{path}: missing columns {sorted(absent)}')
    for col in source_traits+['end_dpg']:
        frame[col]=pd.to_numeric(frame[col],errors='raise')
    for col in ['included','is_control','time_to_10mm2_left_censored']:
        raw=frame[col].astype('string').str.lower()
        if col=='time_to_10mm2_left_censored': raw=raw.fillna('false')
        if not raw.isin(['true','false']).all(): raise ValueError(f'Invalid boolean column {col}')
        frame[col]=raw.eq('true').astype(bool)
    if not frame.is_control.equals(frame.effector.isin(CONTROLS)):
        raise ValueError('is_control flags disagree with the declared control names')
    frame['source_csv_line']=np.arange(len(frame))+2
    return frame


def write_json(path,data):
    import json
    path=Path(path); path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(data,indent=2,allow_nan=False)+'\n')


def portable_path(path,base):
    """Serialise filesystem paths relative to the analysis output root."""
    import os
    return Path(os.path.relpath(Path(path).resolve(),Path(base).resolve())).as_posix()


def default_paths():
    root=Path(__file__).resolve().parents[1]
    return root/'data',root


def write_methods(output,summary):
    versions=', '.join(f'{k} {v}' for k,v in summary['software_versions'].items())
    text=f'''# Methods: final-sheet phenotype clustering

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

Each tree was cut to exactly the requested number of clusters: {', '.join(map(str,summary['k_values']))}. Cluster labels were assigned in dendrogram order and are local to each cut. The two colour views used the same tree because changing the origin with a fixed denominator preserves Euclidean distances. Matrix-distance equality was checked numerically. This check is an algebraic validation, not evidence of biological cluster stability.

Principal component analysis (PCA) was fitted to the complete scaled genotype matrix after column centring, using full singular-value decomposition. PC1–PC2 plots were coloured by the Ward assignments. All component scores, loadings and explained-variance fractions were exported. PCA was not used to select the clustering space. Silhouette, Calinski–Harabasz and Davies–Bouldin scores, cluster sizes and singleton counts were calculated descriptively for the requested cuts. No optimal k was asserted. No stability, trait-correlation or sensitivity analyses were performed, and no inferential p-values were calculated.

## Figure rendering and verification

Heatmaps used a continuous blue–white–red scale centred on zero with limits -2.5 and +2.5. Values beyond those limits saturated only in the rendering; numerical tables and clustering used untruncated scaled values. No white band or significance threshold was applied. Higher numerical values were red; for timing this means later threshold attainment. A single row permutation was applied and checked against the row labels. Both labelled and compact heatmaps were exported. All figures were saved as vector PDF/SVG and as JPEG/TIFF at {summary['raster_dpi']} dpi. The rendering manifest records row/column order and plotting coordinates for independent PDF verification.

## Software and provenance

{versions}. Source hashes, analysis settings, exclusion logs, trait-level replicate counts, numerical matrices and figure manifests accompany the results. The source CSVs are the only research-data inputs.
'''
    (output/'reports/METHODS.md').write_text(text)


def run_analysis(input_dir,output_dir,scopes=None,ks=None,dpi=900,tables_only=False,overwrite=False):
    import datetime,hashlib,importlib.metadata,logging,json
    input_dir,output=Path(input_dir).resolve(),Path(output_dir).resolve()
    scopes=list(scopes or ['root','shoot','combined']); ks=list(ks or range(7,13))
    for directory in ['tables','reports','figures','logs']:
        (output/directory).mkdir(parents=True,exist_ok=True)
    if (output/'reports/RUN_SUMMARY.json').exists() and not overwrite:
        raise FileExistsError('Existing run summary: use --overwrite explicitly to regenerate this run')
    sources={t:input_dir/f'phenotype_summary_{t}.csv' for t in ['root','shoot']}
    hashes={t:hashlib.sha256(path.read_bytes()).hexdigest() for t,path in sources.items()}
    frames={t:read_final_csv(path) for t,path in sources.items()}
    cleaned,log=filter_final_inputs(frames)
    log.to_csv(output/'tables/replicate_disposition.csv',index=False)
    prepared={t:derive_timing(d,29.) for t,d in cleaned.items()}
    bundles={t:summarise_tissue(d,'R' if t=='root' else 'S') for t,d in prepared.items()}
    common=sorted(set(bundles['root']['medians'].index)&set(bundles['shoot']['medians'].index))
    bundles['combined']={key:pd.concat([bundles[t][key].loc[common] for t in ['root','shoot']],axis=1)
                         for key in ['medians','counts','global_z','control_z']}
    bundles['combined']['reference']=pd.concat([bundles[t]['reference'] for t in ['root','shoot']])
    qc={}
    for t,d in prepared.items():
        d.to_csv(output/'tables'/f'{t}_retained_replicates_with_timing.csv',index=False)
        raw=frames[t]
        qc[t]={'input_records':len(raw),'source_exclusions':int((~raw.included).sum()),
               'additional_joint_zero_exclusions':int(log.query('tissue == @t').disposition.eq('excluded_joint_zero').sum()),
               'retained_records':len(d),'genotypes':d.effector.nunique(),
               'remaining_zero_final_size':int(d.area_end_mm2.eq(0).sum()),
               'excluded_entire_genotypes':sorted(set(raw.effector)-set(d.effector)),
               'control_replicates':{c:int(d.effector.eq(c).sum()) for c in CONTROLS},
               'timing_status_counts':{str(a):int(b) for a,b in d.time_to_10mm2_status.value_counts().items()},
               'missing_analysis_values':{a:int(b) for a,b in d[TRAITS].isna().sum().items()},
               'low_replicate_genotypes':{str(a):int(b) for a,b in d.groupby('effector').size().items() if b<3}}
    write_json(output/'reports/QC_SUMMARY.json',qc)
    manifests=[]; analyses_summary={}
    for scope in scopes:
        logging.info('Starting %s',scope)
        bundle=bundles[scope]
        result=analyse_matrix(bundle['global_z'],bundle['control_z'],ks)
        td=output/'tables'/scope; td.mkdir(parents=True,exist_ok=True)
        for name,df in bundle.items(): df.to_csv(td/f'{name}.csv')
        for name in ['assignments','metrics','scores','loadings','variance']:
            fname={'assignments':'cluster_assignments','metrics':'k_comparison',
                   'scores':'pca_scores','loadings':'pca_loadings','variance':'pca_variance'}[name]
            result[name].to_csv(td/f'{fname}.csv',index=name not in ['metrics','variance'])
        pd.DataFrame(result['tree'],columns=['left_child','right_child','ward_distance','n_leaves']).to_csv(td/'linkage.csv',index=False)
        pd.DataFrame({'row_position':range(len(result['leaf_order'])),'effector':result['leaf_order']}).to_csv(td/'dendrogram_order.csv',index=False)
        cluster_sizes=[]
        for k in ks:
            counts=result['assignments'][f'k{k}'].value_counts().sort_index()
            cluster_sizes.extend({'k':k,'cluster':int(c),'n_genotypes':int(n)} for c,n in counts.items())
        pd.DataFrame(cluster_sizes).to_csv(td/'cluster_sizes.csv',index=False)
        # Per-cell source data, including actual control values and replicate support.
        for view in ['global','control']:
            mat=bundle[f'{view}_z']
            cells=mat.rename_axis(index='effector',columns='trait').stack().rename('z').reset_index()
            counts=bundle['counts'].rename_axis(index='effector',columns='trait').stack().rename('n_replicates').reset_index()
            cells=cells.merge(counts,on=['effector','trait'],validate='one_to_one')
            cells['is_control']=cells.effector.isin(CONTROLS)
            cells.to_csv(td/f'{view}_heatmap_source.csv',index=False)
        analyses_summary[scope]={'n_genotypes':len(bundle['medians']),
                                 'n_traits':len(bundle['medians'].columns),
                                 'pca_pc1_fraction':float(result['variance'].explained_variance_ratio.iloc[0]),
                                 'pca_pc2_fraction':float(result['variance'].explained_variance_ratio.iloc[1]),
                                 'validation':result['validation']}
        manifest=[]
        if not tables_only:
            for k in ks:
                for view in ['global','control']:
                    ref=bundle['reference']
                    baseline=(ref.pooled_control_median-ref.all_replicate_mean)/ref.all_replicate_sd if view=='global' else pd.Series(0.,index=ref.index)
                    for labelled in [True,False]:
                        label='labelled' if labelled else 'compact'
                        fig,meta=make_heatmap(bundle[f'{view}_z'],result,k,view,baseline,scope,labelled)
                        written=save_figure(fig,output/'figures'/scope/f'k{k}'/f'{scope}_k{k}_{view}_{label}',dpi)
                        meta['files']={ext:portable_path(path,output) for ext,path in written.items()}
                        meta['matrix_csv']=portable_path(td/f'{view}_z.csv',output);meta['raster_dpi']=dpi
                        manifest.append(meta)
                        write_json(output/'reports'/f'{scope}_figure_manifest.json',manifest)
                fig=make_pca_plot(result,k,scope)
                written=save_figure(fig,output/'figures'/scope/f'k{k}'/f'{scope}_k{k}_pca',dpi)
                manifest.append({'kind':'pca','scope':scope,'k':k,'raster_dpi':dpi,
                                 'files':{ext:portable_path(path,output) for ext,path in written.items()}})
                write_json(output/'reports'/f'{scope}_figure_manifest.json',manifest)
                logging.info('%s k=%s complete: paired heatmap views and PCA',scope,k)
            fig=make_comparison(result,scope)
            written=save_figure(fig,output/'figures'/scope/f'{scope}_k_comparison',dpi)
            manifest.append({'kind':'k_comparison','scope':scope,'raster_dpi':dpi,
                             'files':{ext:portable_path(path,output) for ext,path in written.items()}})
        write_json(output/'reports'/f'{scope}_figure_manifest.json',manifest)
        manifests.extend(manifest)
        logging.info('Completed %s',scope)
    versions={name:importlib.metadata.version(pkg) for name,pkg in
              [('numpy','numpy'),('pandas','pandas'),('SciPy','scipy'),
               ('scikit-learn','scikit-learn'),('Matplotlib','matplotlib'),('Pillow','Pillow')]}
    import platform
    versions['Python']=platform.python_version()
    summary={'status':'tables_complete' if tables_only else 'complete',
             'completed_local_time':datetime.datetime.now().astimezone().isoformat(),
             'input_files':{t:portable_path(path,output) for t,path in sources.items()},'input_sha256':hashes,
             'k_values':ks,'raster_dpi':dpi,'timing_cap_dpg':29.,'scopes':scopes,
             'pairing_confirmed':False,'sd_ddof':0,'controls_in_clustering':True,
             'stability_analysis':False,'trait_correlation_analysis':False,'sensitivity_analysis':False,
             'qc':qc,'analyses':analyses_summary,'software_versions':versions,
             'figure_sets':len(manifests),'figure_files':sum(len(m['files']) for m in manifests)}
    assert hashes=={t:hashlib.sha256(path.read_bytes()).hexdigest() for t,path in sources.items()}
    expected=len(scopes)*(len(ks)*5+1) if not tables_only else 0
    if len(manifests)!=expected: raise AssertionError('Figure-set count mismatch')
    for entry in manifests:
        for path in entry['files'].values():
            resolved=(output/path).resolve()
            assert resolved.is_file() and resolved.stat().st_size>0
    write_json(output/'reports/figure_manifest.json',manifests)
    write_methods(output,summary)
    write_json(output/'reports/RUN_SUMMARY.json',summary)
    return summary


def main():
    import argparse,logging,json
    parser=argparse.ArgumentParser(description=__doc__)
    default_input,default_output=default_paths()
    parser.add_argument('--input-dir',type=Path,default=default_input)
    parser.add_argument('--output-dir',type=Path,default=default_output)
    parser.add_argument('--scopes',nargs='+',choices=['root','shoot','combined'],default=['root','shoot','combined'])
    parser.add_argument('--ks',nargs='+',type=int,default=list(range(7,13)))
    parser.add_argument('--dpi',type=int,default=900)
    parser.add_argument('--tables-only',action='store_true')
    parser.add_argument('--overwrite',action='store_true')
    args=parser.parse_args()
    logdir=args.output_dir/'logs';logdir.mkdir(parents=True,exist_ok=True)
    logging.basicConfig(level=logging.INFO,format='%(asctime)s | %(message)s',
                        handlers=[logging.StreamHandler(),logging.FileHandler(logdir/'analysis.log',mode='a')])
    logging.getLogger('fontTools').setLevel(logging.WARNING)
    try:
        summary=run_analysis(args.input_dir,args.output_dir,args.scopes,args.ks,args.dpi,args.tables_only,args.overwrite)
    except Exception:
        logging.exception('Analysis failed; partial outputs must not be treated as complete')
        raise
    print(json.dumps(summary,indent=2))


if __name__=='__main__':
    main()
