#!/usr/bin/env python3
"""Independent verification of final-sheet clustering and exported figures.

Does not import the production pipeline. Recomputes filters, timing encoding,
reference moments and genotype summaries directly from source CSVs. Reads PDF
cell colours and label positions, and checks four-format coverage/raster DPI.
"""
from pathlib import Path
import argparse
import collections
import hashlib
import json
import xml.etree.ElementTree as ET
import numpy as np
import pandas as pd
import fitz
from PIL import Image
from matplotlib.colors import LinearSegmentedColormap, Normalize
from scipy.cluster.hierarchy import leaves_list, cut_tree, linkage
from scipy.spatial.distance import pdist

CONTROLS=['Col-0','EmptyVector','GFP']
KEYS=['effector','replicate']


def resolve_path(base,path):
    path=Path(path)
    return path if path.is_absolute() else (Path(base)/path).resolve()
TRAITS=['area_end_mm2','AUC','time_to_10mm2_capped','AGR_mean','AGR_max',
        'RGR_mean','RGR_peak','RGR_final']


def verify(base):
    base=Path(base).resolve()
    summary=json.loads((base/'reports/RUN_SUMMARY.json').read_text())
    manifest=json.loads((base/'reports/figure_manifest.json').read_text())
    assert summary['status']=='complete'
    frames={t:pd.read_csv(resolve_path(base,p),dtype={'effector':str,'replicate':str})
            for t,p in summary['input_files'].items()}
    for t,path in summary['input_files'].items():
        assert hashlib.sha256(resolve_path(base,path).read_bytes()).hexdigest()==summary['input_sha256'][t]
    zero_keys={t:set(map(tuple,d.loc[d.area_end_mm2.eq(0),KEYS].values)) for t,d in frames.items()}
    joint=zero_keys['root']&zero_keys['shoot']
    expected={}
    for t,d in frames.items():
        joint_mask=np.array([tuple(x) in joint for x in d[KEYS].values])
        keep=d.loc[d.included.eq(True)&~joint_mask].copy()
        right=keep.time_to_10mm2.isna()
        assert keep.loc[right,'end_dpg'].ge(29).all()
        keep['time_to_10mm2_capped']=np.where(right,29.,np.minimum(keep.time_to_10mm2,29.))
        retained_csv=pd.read_csv(base/'tables'/f'{t}_retained_replicates_with_timing.csv',dtype={'effector':str,'replicate':str})
        assert list(map(tuple,retained_csv[KEYS].values))==list(map(tuple,keep[KEYS].values))
        np.testing.assert_allclose(retained_csv[TRAITS].values,keep[TRAITS].values,equal_nan=True)
        np.testing.assert_allclose(retained_csv.time_to_10mm2,keep.time_to_10mm2,equal_nan=True)
        mean,sd=keep[TRAITS].mean(),keep[TRAITS].std(ddof=0)
        ctrl=keep.loc[keep.effector.isin(CONTROLS),TRAITS].median()
        med=keep.groupby('effector',sort=True)[TRAITS].median()
        counts=keep.groupby('effector',sort=True)[TRAITS].count()
        global_z=(med-mean)/sd; control_z=(med-ctrl)/sd
        prefix='R' if t=='root' else 'S'
        for obj in (med,counts,global_z,control_z): obj.columns=[prefix+'_'+c for c in obj.columns]
        expected[t]={'medians':med,'counts':counts,'global_z':global_z,'control_z':control_z}
        assert len(keep)==summary['qc'][t]['retained_records']
        assert int((d.included.eq(True)&joint_mask).sum())==summary['qc'][t]['additional_joint_zero_exclusions']
    common=sorted(set(expected['root']['medians'].index)&set(expected['shoot']['medians'].index))
    expected['combined']={key:pd.concat([expected[t][key].loc[common] for t in ['root','shoot']],axis=1)
                          for key in expected['root']}
    all_k=summary['k_values']; matrix_cache={}
    for scope in summary['scopes']:
        for name,df in expected[scope].items():
            saved=pd.read_csv(base/'tables'/scope/f'{name}.csv',index_col=0)
            assert saved.index.tolist()==df.index.tolist()
            assert saved.columns.tolist()==df.columns.tolist()
            np.testing.assert_allclose(saved.to_numpy(),df.to_numpy(),atol=1e-10,rtol=1e-10)
            if name.endswith('_z'): matrix_cache[(scope,name)]=saved
        X=expected[scope]['global_z']
        Z=linkage(pdist(X.to_numpy()),method='ward',optimal_ordering=True)
        stored_Z=pd.read_csv(base/'tables'/scope/'linkage.csv').values
        np.testing.assert_allclose(Z,stored_Z,atol=1e-10,rtol=1e-10)
        leaves=X.index[leaves_list(Z)].tolist()
        order=pd.read_csv(base/'tables'/scope/'dendrogram_order.csv').effector.tolist()
        assert leaves==order
        labels=pd.read_csv(base/'tables'/scope/'cluster_assignments.csv',index_col=0)
        assert labels.index.tolist()==X.index.tolist()
        for k in all_k:
            raw=cut_tree(Z,n_clusters=[k]).ravel()
            order_pos=leaves_list(Z)
            mapping={old:i+1 for i,old in enumerate(dict.fromkeys(raw[order_pos]))}
            wanted=np.array([mapping[v] for v in raw])
            np.testing.assert_array_equal(labels[f'k{k}'],wanted)
            assert labels[f'k{k}'].nunique()==k
        np.testing.assert_allclose(pdist(expected[scope]['global_z']),pdist(expected[scope]['control_z']),atol=1e-10)
    expected_sets=len(summary['scopes'])*(len(all_k)*5+1)
    assert len(manifest)==expected_sets==summary['figure_sets']
    paths=[resolve_path(base,p) for entry in manifest for p in entry['files'].values()]
    assert len(paths)==len(set(paths))==summary['figure_files']
    actual=set(p.resolve() for p in (base/'figures').rglob('*') if p.suffix in ['.pdf','.jpeg','.svg','.tiff'])
    assert actual==set(p.resolve() for p in paths),'Unlisted or missing figure files'
    cmap=LinearSegmentedColormap.from_list('expected',['#2166AC','#FFFFFF','#B2182B'],N=257)
    norm=Normalize(-2.5,2.5)
    checked_cells=checked_labels=checked_refs=0
    max_pdf_error=max_raster_error=0
    checks=[]
    Image.MAX_IMAGE_PIXELS=None  # Only our manifest-listed, locally generated images.
    for entry in manifest:
        assert set(entry['files'])=={'pdf','jpeg','svg','tiff'}
        for ext,path in entry['files'].items():
            path=resolve_path(base,path)
            assert path.is_file() and path.stat().st_size>0
            if ext=='svg': ET.parse(path)
            if ext in ['jpeg','tiff']:
                with Image.open(path) as image:
                    assert image.format=={'jpeg':'JPEG','tiff':'TIFF'}[ext]
                    assert all(abs(v-summary['raster_dpi'])<.1 for v in image.info['dpi'])
                    if entry['kind']=='heatmap':
                        np.testing.assert_allclose(image.size,np.array(entry['figure_inches'])*summary['raster_dpi'],atol=1)
                    image.verify()
        with fitz.open(resolve_path(base,entry['files']['pdf'])) as doc:
            assert doc.page_count==1
            page=doc[0]
            for word in page.get_text('words'):
                assert word[0]>=-1 and word[1]>=-1 and word[2]<=page.rect.width+1 and word[3]<=page.rect.height+1,('Clipped text',entry['files']['pdf'],word)
            if entry['kind']!='heatmap':
                checks.append({'pdf':entry['files']['pdf'],'kind':entry['kind'],'status':'PASS'})
                continue
            df=matrix_cache[(entry['scope'],entry['view']+'_z')].loc[entry['row_order'],entry['trait_order']]
            expected_rgb=np.rint(cmap(norm(df.values))[...,:3]*255).astype(int)
            n,p=df.shape
            x,y,w,h=entry['heatmap_bounds']
            pix=page.get_pixmap(matrix=fitz.Matrix(3,3),alpha=False)
            img=np.frombuffer(pix.samples,dtype=np.uint8).reshape(pix.height,pix.width,pix.n)[...,:3]
            for i in range(n):
                for j in range(p):
                    ix=int((x+(j+.5)/p*w)*pix.width)
                    iy=int((1-y-h+(i+.5)/n*h)*pix.height)
                    error=int(np.abs(img[iy,ix].astype(int)-expected_rgb[i,j]).max())
                    max_pdf_error=max(max_pdf_error,error)
                    assert error<=3,(entry['files']['pdf'],df.index[i],df.columns[j],df.iloc[i,j],img[iy,ix].tolist(),expected_rgb[i,j].tolist())
                    checked_cells+=1
            if entry['labelled']:
                words=[v for v in page.get_text('words') if v[4] in set(df.index) and v[0]>(x+w)*page.rect.width-1]
                assert len(words)==n,(entry['files']['pdf'],'label count',len(words),n)
                words=sorted(words,key=lambda v:(v[1]+v[3])/2)
                assert [v[4] for v in words]==df.index.tolist()
                for i,v in enumerate(words):
                    expected_y=(1-y-h+(i+.5)/n*h)*page.rect.height
                    assert abs((v[1]+v[3])/2-expected_y)<3.5
                    checked_labels+=1
            bx,by,bw,bh=entry['reference_bounds']
            ref_rgb=np.rint(cmap(norm(entry['baseline_values']))[...,:3]*255).astype(int)
            if entry['view']=='control':
                assert np.all(np.array(entry['baseline_values'])==0)
                assert np.all(ref_rgb==255)
            for j in range(p):
                ix=int((bx+(j+.5)/p*bw)*pix.width); iy=int((1-by-bh/2)*pix.height)
                assert np.max(np.abs(img[iy,ix].astype(int)-ref_rgb[j]))<=3
                checked_refs+=1
            # Sample actual control-cell colours in BOTH requested raster formats.
            for ext in ['jpeg','tiff']:
                with Image.open(resolve_path(base,entry['files'][ext])) as raster:
                    for gene in CONTROLS:
                        if gene not in df.index: continue
                        i=df.index.get_loc(gene)
                        for j in range(p):
                            ix=int((x+(j+.5)/p*w)*raster.width)
                            iy=int((1-y-h+(i+.5)/n*h)*raster.height)
                            value=np.array(raster.getpixel((ix,iy))[:3],dtype=int)
                            error=int(np.abs(value-expected_rgb[i,j]).max())
                            max_raster_error=max(max_raster_error,error)
                            assert error<=(12 if ext=='jpeg' else 2),(ext,gene,df.columns[j],error)
            checks.append({'pdf':entry['files']['pdf'],'kind':'heatmap','cells_checked':n*p,
                           'row_labels_checked':n if entry['labelled'] else 0,'status':'PASS'})
    report={'status':'PASS','verified_figure_sets':len(manifest),'verified_figure_files':len(paths),
            'figures_by_kind':dict(collections.Counter(e['kind'] for e in manifest)),
            'files_by_extension':dict(collections.Counter(p.suffix for p in paths)),
            'pdf_heatmap_cells_checked':checked_cells,'pdf_row_labels_checked':checked_labels,
            'pooled_reference_cells_checked':checked_refs,'max_pdf_rgb_channel_error':max_pdf_error,
            'max_raster_control_rgb_channel_error':max_raster_error,
            'raster_dpi':summary['raster_dpi'],'input_hashes_unchanged':True,
            'independent_numerical_reconstruction_passed':True,'exact_k_cuts_passed':True,
            'checks':checks}
    (base/'reports/VERIFICATION.json').write_text(json.dumps(report,indent=2)+'\n')
    concise={k:v for k,v in report.items() if k!='checks'}
    (base/'reports/VERIFICATION.txt').write_text(json.dumps(concise,indent=2)+'\n')
    print(json.dumps(concise,indent=2))
    return report


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--base',type=Path,default=Path(__file__).resolve().parents[1])
    verify(parser.parse_args().base)
