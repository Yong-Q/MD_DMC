import sys
from pathlib import Path
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from transport_framework import block_current, summarize_values, find_plateau

def test_block_current_units():
 p=np.array([[0,0,0],[1,2,3],[2,4,6]],float)
 j=block_current(p,1.,charge_moment=True)
 np.testing.assert_allclose(j[:,0],np.array([1,2,3])*1.602176634e-17)

def test_replica_error_and_no_censoring():
 r=summarize_values([1.,2.,3.,4.],expected=4)
 assert r['mean']==2.5 and r['n']==4 and r['ci95'][0]<2.5<r['ci95'][1]
 assert summarize_values([1.,2.],expected=4)['complete'] is False
 assert summarize_values([1.,None,3.,4.],expected=4)['mean'] is None

def test_plateau_rejects_subdiffusion():
 t=np.arange(15000.);c={'plateau_max_lag_fraction':.25,'plateau_alpha':[.9,1.1],'plateau_min_width_ps':500,'plateau_min_lag_ratio':2}
 assert find_plateau(t,2*t,c)['status']=='plateau_candidate'
 assert find_plateau(t,np.sqrt(t+1),c)['status']=='no_plateau'

def test_group_uses_one_window_and_all_replicas(tmp_path):
 import csv,json,pandas as pd
 from transport_framework import summarize
 manifest=tmp_path/'manifest.csv';rows=[dict(array_index=str(i),base_name='SYN',dmc_per_lifsi='2',replica=str(i+1),dmc_molecules='100',lifsi_pairs='50') for i in range(4)]
 with manifest.open('w',newline='') as f:
  w=csv.DictWriter(f,fieldnames=rows[0]);w.writeheader();w.writerows(rows)
 root=tmp_path/'analysis';t=np.arange(15000.)
 config={'plateau_max_lag_fraction':.25,'plateau_alpha':[.9,1.1],'plateau_min_width_ps':500,'plateau_min_lag_ratio':2,'expected_replicas':4,'msd_fallback_lag_fraction':[.1,.5],'config_sha256':'test'}
 for i in range(4):
  d=root/'cases'/f'{i:03d}';d.mkdir(parents=True)
  j={'config':config,'metadata':{'base_name':'SYN','dmc_per_lifsi':2,'replica_index':i+1,'velocity_seed':100+i,'box_volume_A3':100000.,'temperature_K':298.,'actual_pairs':50},'provenance':{'li_charge_e':.7,'fsi_charge_e':-.7},'conductivity':{'value_mS_cm':.1+i*.01,'status':'needs_review'}}
  (d/'results.json').write_text(json.dumps(j));pd.DataFrame({'lag_ps':t,'Li':2*t,'FSI_COM':3*t,'collective_charge':4*t}).to_csv(d/'curves.csv',index=False)
 out=tmp_path/'summary';summarize(manifest,root,config,out)
 states=pd.read_csv(out/'by_state.csv');cases=pd.read_csv(out/'per_case.csv');assert len(cases)==4
 for name in ('Li','FSI_COM','collective_charge'):
  assert cases[name+'_lag_start_ps'].nunique()==1 and cases[name+'_lag_end_ps'].nunique()==1
 assert states.loc[states.observable=='collective_STACIE','mean'].iloc[0]==.115
 assert states.loc[states.observable=='collective_STACIE','status'].iloc[0]=='review_required'

def test_plateau_invalid_input_rejected():
 t=np.arange(15000.);y=t.copy();y[500]=np.nan;c={'plateau_max_lag_fraction':.25,'plateau_alpha':[.9,1.1],'plateau_min_width_ps':500,'plateau_min_lag_ratio':2}
 assert find_plateau(t,y,c)['status']=='invalid_curve'

def test_config_rejects_wrong_lag_fraction(tmp_path):
 import json,pytest
 from transport_framework import read_config
 c=json.loads((Path(__file__).resolve().parents[1]/'config/transport.json').read_text());c['msd_fallback_lag_fraction']=[.5,.1];p=tmp_path/'bad.json';p.write_text(json.dumps(c))
 with pytest.raises(ValueError):read_config(p)

def test_cache_reuse_and_stale_detection(tmp_path,monkeypatch):
 import json,pytest,transport_v2
 from transport_framework import prepare_cache
 case=tmp_path/'case';(case/'traj').mkdir(parents=True);(case/'metadata.json').write_text(json.dumps({'expected_chunks':15}));(case/'system.data').write_text('topology')
 for i in range(1,16):(case/'traj'/f'chunk_{i:04d}.dcd').write_bytes(b'fake')
 monkeypatch.setattr(transport_v2,'read_positions',lambda *a:(np.zeros((4,1,3)),np.zeros((4,1,3)),np.zeros((4,3)),np.ones((4,6)),{}));monkeypatch.setattr(transport_v2,'temperature_diagnostic',lambda *a:{'status':'test'})
 cache=tmp_path/'cache.npz';prepare_cache(case,cache)
 def forbidden(*a):raise AssertionError('must not reread trajectory')
 monkeypatch.setattr(transport_v2,'read_positions',forbidden);prepare_cache(case,cache)
 (case/'system.data').write_text('changed topology')
 with pytest.raises(ValueError,match='stale cache'):prepare_cache(case,cache)

def test_merge_rejects_duplicate_and_wrong_config_before_copy(tmp_path):
 import json,pytest
 from transport_framework import merge_outputs
 roots=[]
 for name in ('one','two'):
  r=tmp_path/name;d=r/'cases'/'000';d.mkdir(parents=True);(d/'results.json').write_text(json.dumps({'schema_version':1,'config':{'config_sha256':'same'}}));(d/'curves.csv').write_text('lag_ps,Li\n0,0\n');roots.append(r)
 out=tmp_path/'merged'
 with pytest.raises(ValueError,match='duplicate'):merge_outputs(roots,out,tmp_path/'manifest',{'config_sha256':'same'})
 assert not out.exists()
 with pytest.raises(ValueError,match='mismatch'):merge_outputs(roots[:1],out,tmp_path/'manifest',{'config_sha256':'different'})
 assert not out.exists()

def test_cached_shapes_and_nan_rejected():
 import pytest
 from transport_framework import validate_cached_arrays
 x=np.zeros((4,2,3));m=np.zeros((4,3));meta={'actual_pairs':2,'expected_chunks':1,'chunk_ps':4};prov={'dt_ps':1}
 validate_cached_arrays(x,x,m,meta,prov)
 with pytest.raises(ValueError):validate_cached_arrays(x[:-1],x,m,meta,prov)
 x[0,0,0]=np.nan
 with pytest.raises(ValueError):validate_cached_arrays(x,x,m,meta,prov)
