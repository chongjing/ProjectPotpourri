"""Regression tests; small synthetic fixtures are test data, not research results."""
from pathlib import Path
import importlib.util
import numpy as np
import pandas as pd

SCRIPT = Path(__file__).resolve().parents[1] / 'scripts' / 'run_final_clustering.py'


def pipeline():
    assert SCRIPT.exists(), 'The new standalone pipeline has not been implemented'
    spec = importlib.util.spec_from_file_location('final_clustering', SCRIPT)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_filter_only_joint_zero_and_explicit_exclusion():
    m = pipeline()
    root = pd.DataFrame({'effector': ['A','B','C','D','E'], 'replicate': ['1']*5,
                         'included': [True,True,True,False,True],
                         'area_end_mm2': [0.,0.,4.,np.nan,0.],
                         'exclusion_reason': ['', '', '', 'manual', '']})
    shoot = pd.DataFrame({'effector': ['A','B','C','D'], 'replicate': ['1']*4,
                          'included': [True,True,False,True],
                          'area_end_mm2': [0.,2.,np.nan,5.],
                          'exclusion_reason': ['', '', 'manual', '']})
    original = root.copy(deep=True)
    cleaned, log = m.filter_final_inputs({'root': root, 'shoot': shoot})
    assert cleaned['root'].effector.tolist() == ['B','C','E']
    assert cleaned['shoot'].effector.tolist() == ['B','D']
    assert log.query("disposition == 'excluded_joint_zero'").shape[0] == 2
    assert log.query("disposition == 'excluded_by_final_sheet'").shape[0] == 2
    assert len(log) == len(root) + len(shoot)
    pd.testing.assert_frame_equal(root, original)


def test_cap_times_preserves_originals_and_refuses_early_censoring():
    m = pipeline()
    assert hasattr(m, 'derive_timing'), 'Censor-aware timing encoding is not implemented'
    frame = pd.DataFrame({'time_to_10mm2': [14., 32., np.nan, 8.],
                          'end_dpg': [34., 34., 29., 34.],
                          'time_to_10mm2_left_censored': [False, False, False, True]})
    out = m.derive_timing(frame, 29.)
    np.testing.assert_allclose(out.time_to_10mm2_capped, [14.,29.,29.,8.])
    assert np.isnan(out.loc[2,'time_to_10mm2'])
    assert out.time_to_10mm2_status.tolist() == ['observed', 'observed_after_cap',
                                               'right_censored', 'left_censored_upper_bound']
    assert not frame.columns.str.contains('capped').any()
    import pytest
    with pytest.raises(ValueError, match='before'):
        m.derive_timing(frame.assign(end_dpg=[34.,34.,28.,34.]),29.)


def test_scaling_uses_individual_replicates_and_one_shared_denominator():
    m = pipeline()
    assert hasattr(m, 'summarise_tissue'), 'Replicate-based scaling is not implemented'
    frame = pd.DataFrame({'effector': ['Col-0','EmptyVector','GFP','A','A','B'],
                          'replicate': list('123456')})
    for trait in m.TRAITS:
        frame[trait] = [1.,2.,3.,10.,20.,30.]
    frame.loc[4,'RGR_mean'] = np.nan
    out = m.summarise_tissue(frame, 'R')
    ref = out['reference']
    x = frame.area_end_mm2
    assert ref.loc['R_area_end_mm2','all_replicate_mean'] == x.mean()
    assert ref.loc['R_area_end_mm2','all_replicate_sd'] == x.std(ddof=0)
    assert ref.loc['R_RGR_mean','n_retained_values'] == 5
    assert out['counts'].loc['A','R_RGR_mean'] == 1
    expected_ctrl = frame.loc[frame.effector.isin(m.CONTROLS),'area_end_mm2'].median()
    assert ref.loc['R_area_end_mm2','pooled_control_median'] == expected_ctrl
    np.testing.assert_allclose(out['global_z'].loc['A','R_area_end_mm2'],
                               (15.-x.mean())/x.std(ddof=0))
    np.testing.assert_allclose(out['control_z'].loc['A','R_area_end_mm2'],
                               (15.-expected_ctrl)/x.std(ddof=0))
    from scipy.spatial.distance import pdist
    np.testing.assert_allclose(pdist(out['global_z']),pdist(out['control_z']),atol=1e-12)


def test_fresh_tree_has_exact_requested_cuts_and_pca():
    m = pipeline()
    assert hasattr(m, 'analyse_matrix'), 'Fresh clustering/PCA is not implemented'
    genes = m.CONTROLS + [f'TEST{i:02}' for i in range(12)]
    rng = np.random.default_rng(2026)
    z = pd.DataFrame(rng.normal(size=(15,8)), index=genes,
                     columns=[f'R_{t}' for t in m.TRAITS])
    ctrl = z + np.linspace(-1,1,8)
    result = m.analyse_matrix(z, ctrl, [7,8,9,10,11,12])
    assert set(result['leaf_order']) == set(genes)
    for k in range(7,13):
        assert result['assignments'][f'k{k}'].nunique() == k
        ordered = result['assignments'].loc[result['leaf_order'],f'k{k}'].to_numpy()
        assert np.all(np.diff(ordered) >= 0)
    assert result['scores'].shape == (15,8)
    assert result['loadings'].shape == (8,8)
    assert result['validation']['max_distance_difference'] < 1e-10
    assert result['metrics']['k'].tolist() == [7,8,9,10,11,12]


def test_heatmap_alignment_and_four_format_export(tmp_path):
    m = pipeline()
    assert hasattr(m, 'make_heatmap'), 'Verified heatmap rendering is not implemented'
    assert hasattr(m, 'save_figure'), 'Four-format export is not implemented'
    rng = np.random.default_rng(42)
    genes = m.CONTROLS + [f'TEST{i:02}' for i in range(12)]
    z = pd.DataFrame(rng.normal(size=(15,8)),index=genes,
                     columns=[f'R_{t}' for t in m.TRAITS])
    result = m.analyse_matrix(z,z+1.,[7])
    baseline = pd.Series(0.,index=z.columns)
    fig, meta = m.make_heatmap(z,result,7,'control',baseline,'root',True)
    ax = next(a for a in fig.axes if a.get_label() == 'heatmap')
    np.testing.assert_allclose(ax.collections[0].get_array(),z.loc[result['leaf_order']].values)
    assert [x.get_text() for x in ax.get_yticklabels()] == result['leaf_order']
    assert meta['row_order'] == result['leaf_order']
    files = m.save_figure(fig,tmp_path/'test_heatmap',dpi=90)
    assert set(files) == {'pdf','jpeg','svg','tiff'}
    from PIL import Image
    import fitz
    import xml.etree.ElementTree as ET
    for ext in ['jpeg','tiff']:
        with Image.open(files[ext]) as im:
            assert im.width > 100
            assert abs(im.info['dpi'][0]-90.) < .1
    with fitz.open(files['pdf']) as doc:
        assert doc.page_count == 1
        assert 'Col-0' in doc[0].get_text()
    ET.parse(files['svg'])


def test_end_to_end_pipeline_uses_only_final_csvs(tmp_path):
    m = pipeline()
    assert hasattr(m, 'run_analysis'), 'Standalone input-to-output pipeline is not implemented'
    rng=np.random.default_rng(87)
    genes=m.CONTROLS+[f'TEST{i:02}' for i in range(12)]
    for tissue in ('root','shoot'):
        rows=[]
        for gene in genes:
            for rep in ('1','2'):
                row={'effector':gene,'replicate':rep,'included':True,
                     'is_control':gene in m.CONTROLS,'exclusion_reason':np.nan,
                     'end_dpg':34.,'start_dpg':6.,'time_to_10mm2_left_censored':False}
                for trait in m.TRAITS:
                    row['time_to_10mm2' if trait=='time_to_10mm2_capped' else trait]=rng.uniform(10.,24.)
                rows.append(row)
        pd.DataFrame(rows).to_csv(tmp_path/f'phenotype_summary_{tissue}.csv',index=False)
    out=tmp_path/'outputs'
    summary=m.run_analysis(tmp_path,out,scopes=['root','shoot','combined'],ks=[7],dpi=60)
    assert summary['status']=='complete'
    assert summary['figure_sets']==3*(2*2+1+1)
    assert summary['figure_files']==summary['figure_sets']*4
    for scope in ('root','shoot','combined'):
        table=pd.read_csv(out/'tables'/scope/'cluster_assignments.csv')
        assert table.k7.nunique()==7
        assert len(table)==len(genes)
        assert (out/'reports'/f'{scope}_figure_manifest.json').exists()
    assert (out/'reports'/'METHODS.md').exists()
    verifier=SCRIPT.with_name('verify_final_outputs.py')
    assert verifier.exists(), 'Independent output/PDF verifier is not implemented'
    import subprocess,sys,json
    check=subprocess.run([sys.executable,str(verifier),'--base',str(out)],
                         text=True,capture_output=True)
    assert check.returncode==0,check.stdout+'\n'+check.stderr
    verification=json.loads((out/'reports'/'VERIFICATION.json').read_text())
    assert verification['status']=='PASS'
    assert verification['verified_figure_files']==summary['figure_files']
    assert all(not Path(p).is_absolute() for p in summary['input_files'].values())
    manifest=json.loads((out/'reports'/'figure_manifest.json').read_text())
    assert all(not Path(p).is_absolute() for entry in manifest for p in entry['files'].values())
    # Relative metadata must verify regardless of caller working directory.
    check=subprocess.run([sys.executable,str(verifier),'--base',str(out)],cwd=tmp_path.parent,
                         text=True,capture_output=True)
    assert check.returncode==0,check.stdout+'\n'+check.stderr


def test_repository_defaults_use_data_subdirectory():
    m=pipeline()
    assert hasattr(m,'default_paths'),'Portable project defaults are not implemented'
    input_dir,output_dir=m.default_paths()
    assert input_dir==SCRIPT.resolve().parents[1]/'data'
    assert output_dir==SCRIPT.resolve().parents[1]


def test_repository_snapshot_verifier_checks_manifest_and_forbidden_paths(tmp_path):
    verifier=SCRIPT.with_name('verify_repository_snapshot.py')
    assert verifier.exists(),'Repository snapshot verifier is not implemented'
    spec=importlib.util.spec_from_file_location('snapshot_verifier',verifier)
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    import hashlib,json,pytest
    payload=tmp_path/'payload.txt';payload.write_text('portable data\n')
    digest=hashlib.sha256(payload.read_bytes()).hexdigest()
    manifest=tmp_path/'manifest.json'
    manifest.write_text(json.dumps({'files':[{'path':'payload.txt','sha256':digest,'size':payload.stat().st_size}]})+'\n')
    result=module.verify_repository_snapshot(tmp_path,manifest)
    assert result['status']=='PASS'
    payload.write_text('/home/'+'cx264/private\n')
    with pytest.raises(AssertionError):
        module.verify_repository_snapshot(tmp_path,manifest)
