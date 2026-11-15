import sys
from pathlib import Path
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from transport_v2 import choose_and_validate, summarize_group

def test_linear_windows_survive_independent_holdout():
 t=np.arange(7500.);y=2*t+3;w=choose_and_validate(t,y,y)
 assert w['status']=='candidate' and .9<=w['discovery']['alpha']<=1.1 and w['window_ps']==[750.,1500.]

def test_subdiffusive_does_not_pass():
 t=np.arange(7500.);y=np.sqrt(t+1)
 assert choose_and_validate(t,y,y)['status']=='no_candidate'

def test_changed_holdout_does_not_pass():
 t=np.arange(7500.);d=2*t+3;v=10*t+3
 assert choose_and_validate(t,d,v)['status']=='holdout_rejected'

def test_holdout_subdiffusive_rejected():
 t=np.arange(7500.);y=2*t+3;bad=(t+1)**.5
 assert choose_and_validate(t,y,bad)['status']=='holdout_rejected'

def test_no_mean_for_incomplete_replicas():
 assert summarize_group([{'replica':1,'seed':10,'value':2},{'replica':2,'seed':20,'value':4}],expected=4)['status']=='incomplete'
 r=summarize_group([{'replica':i,'seed':100+i,'value':float(i)} for i in range(1,5)],expected=4)
 assert r['status']=='complete' and r['mean']==2.5 and np.isclose(r['sample_std'],np.std([1,2,3,4],ddof=1))
 assert summarize_group([{'replica':1,'seed':10,'value':2},{'replica':1,'seed':20,'value':4}],expected=2)['status']=='incomplete'


def test_duplicate_seeds_rejected():
 assert summarize_group([{'replica':i,'seed':5,'value':float(i)} for i in range(1,5)],expected=4)['status']=='incomplete'

def test_reports_failed_candidate_metrics():
 t=np.arange(7500.);r=choose_and_validate(t,np.sqrt(t+1),np.sqrt(t+1))
 assert len(r['discovery_assessments'])==4
 assert all('alpha' in x['discovery'] for x in r['discovery_assessments'])

def test_end_to_end_figures_and_full_fit(tmp_path,monkeypatch):
 import json
 import transport_v2 as m
 case=tmp_path/'case';case.mkdir();out=tmp_path/'new_output'
 meta=dict(expected_chunks=15,frame_dt_ps=1.,box_volume_A3=100000.,temperature_K=298.,actual_pairs=4,case_name='synthetic',base_name='SYN',dmc_per_lifsi=2,replica_index=1,velocity_seed=123)
 (case/'metadata.json').write_text(json.dumps(meta));n=6000
 pos=np.zeros((n,4,3));moment=np.zeros((n,3));dims=np.ones((n,6))
 provenance=dict(dt_ps=1.,li_charge_e=.7,fsi_charge_e=-.7)
 monkeypatch.setattr(m,'read_positions',lambda case,meta:(pos,pos,moment,dims,provenance))
 monkeypatch.setattr(m,'curves',lambda x,kind:np.arange(len(x),dtype=float)+1)
 monkeypatch.setattr(sys,'argv',['transport_v2.py','--case',str(case),'--output-dir',str(out)])
 m.main();result=json.loads((out/'results.json').read_text())
 assert all(result[k]['status']=='first_half_provisional' for k in ('Li','FSI_COM','collective_charge'))
 assert result['nernst_einstein_mS_cm']>0 and result['collective_charge']['conductivity_mS_cm']>0
 assert (out/'diagnostics.png').stat().st_size>1000 and (out/'diagnostics.pdf').stat().st_size>1000

def test_first_half_mode_does_not_average_rejected_holdout(tmp_path,monkeypatch):
 import json
 import transport_v2 as m
 case=tmp_path/'case';case.mkdir()
 meta=dict(frame_dt_ps=1.,box_volume_A3=100000.,temperature_K=298.,actual_pairs=4,case_name='synthetic',base_name='SYN',dmc_per_lifsi=2,replica_index=1,velocity_seed=123)
 (case/'metadata.json').write_text(json.dumps(meta));pos=np.zeros((6000,4,3));mom=np.zeros((6000,3));dims=np.ones((6000,6));prov=dict(dt_ps=1.,li_charge_e=.7,fsi_charge_e=-.7)
 monkeypatch.setattr(m,'read_positions',lambda case,meta:(pos,pos,mom,dims,prov));monkeypatch.setattr(m,'curves',lambda x,kind:np.arange(len(x),dtype=float)+1)
 choice={'status':'holdout_rejected','window_ps':[100.,500.],'discovery':{'status':'ok','alpha':1.,'slope':2.,'intercept':1.,'r2':1.},'holdout':{'status':'ok','alpha':.5,'slope':20.,'intercept':1.,'r2':1.},'holdout_slope_ratio':10.}
 monkeypatch.setattr(m,'choose_and_validate',lambda *args:choice.copy())
 out=tmp_path/'first';monkeypatch.setattr(sys,'argv',['transport_v2.py','--case',str(case),'--output-dir',str(out)]);m.main();j=json.loads((out/'results.json').read_text())
 assert j['Li']['status']=='first_half_provisional' and j['Li']['holdout_status']=='holdout_rejected' and np.isclose(j['Li']['diffusion_cm2_s'],2e-4/6)
 out=tmp_path/'verified';monkeypatch.setattr(sys,'argv',['transport_v2.py','--case',str(case),'--output-dir',str(out),'--mode','verified']);m.main();j=json.loads((out/'results.json').read_text())
 assert j['Li']['status']=='holdout_rejected' and 'diffusion_cm2_s' not in j['Li']

def test_kinisi_first_half_never_uses_late_frames(tmp_path,monkeypatch):
 import json
 import transport_v2 as m
 case=tmp_path/'case';case.mkdir();out=tmp_path/'out'
 meta=dict(frame_dt_ps=1.,box_volume_A3=100000.,temperature_K=298.,actual_pairs=4,case_name='synthetic',base_name='SYN',dmc_per_lifsi=2,replica_index=1,velocity_seed=123)
 (case/'metadata.json').write_text(json.dumps(meta));pos=np.zeros((6000,4,3));mom=np.zeros((6000,3));dims=np.ones((6000,6));seen=[]
 monkeypatch.setattr(m,'read_positions',lambda case,meta:(pos,pos,mom,dims,dict(dt_ps=1.,li_charge_e=.7,fsi_charge_e=-.7)))
 monkeypatch.setattr(m,'curves',lambda x,kind:np.arange(len(x),dtype=float)+1)
 def check(li,box,meta,window):seen.append((len(li),len(box)));return {'status':'test_only'}
 monkeypatch.setattr(m,'kinisi_check',check)
 monkeypatch.setattr(sys,'argv',['transport_v2.py','--case',str(case),'--output-dir',str(out),'--kinisi']);m.main()
 assert seen==[(3000,3000)]

def test_temperature_diagnostic(tmp_path):
 from transport_v2 import temperature_diagnostic
 logs=tmp_path/'logs';logs.mkdir()
 for i in (1,2):
  (logs/f'thermo_chunk_{i:04d}.csv').write_text('step,temp_K\n'+'\n'.join(f'{k},298' for k in range(4))+'\n')
 meta={'expected_chunks':2,'chunk_ps':4,'frame_dt_ps':1,'temperature_K':298}
 r=temperature_diagnostic(tmp_path,meta)
 assert r['first_status']=='pass' and r['second_status']=='pass'
 (logs/'thermo_chunk_0002.csv').unlink()
 assert temperature_diagnostic(tmp_path,meta)['status']=='unavailable'

def test_batch_summary_uses_first_half_and_temperature(tmp_path):
 import csv,json,subprocess
 root=Path(__file__).resolve().parents[1];manifest=tmp_path/'m.csv';rows=[dict(array_index=str(i),base_name='SYN',dmc_per_lifsi='2',replica=str(i+1)) for i in range(4)]
 with manifest.open('w',newline='') as f:
  w=csv.DictWriter(f,fieldnames=rows[0]);w.writeheader();w.writerows(rows)
 results=tmp_path/'r';results.mkdir()
 for i in range(4):
  d=results/f'{i:03d}';d.mkdir();v=dict(status='first_half_provisional',diffusion_cm2_s=1.+i,holdout_status='holdout_rejected')
  j=dict(base_name='SYN',dmc_per_lifsi=2,replica_index=i+1,velocity_seed=i+100,method_version='v2 test',analysis_mode='first-half',temperature_K=298,volume_A3=1000,temperature_diagnostic={'first_status':'pass','second_status':'pass'},Li=v,FSI_COM={'status':'no_candidate'},collective_charge={'status':'no_candidate'})
  (d/'results.json').write_text(json.dumps(j))
 output=tmp_path/'summary.csv';subprocess.run([sys.executable,str(root/'scripts/summarize_transport_v2.py'),'--manifest',str(manifest),'--results-root',str(results),'--output',str(output)],check=True)
 records=list(csv.DictReader(output.open()));li=next(r for r in records if r['observable']=='Li')
 assert li['status']=='complete' and float(li['mean'])==2.5 and 'holdout_rejected' in li['holdout_statuses']
 file=results/'003/results.json';j=json.loads(file.read_text());j['temperature_diagnostic']={};file.write_text(json.dumps(j));output=tmp_path/'missing.csv'
 subprocess.run([sys.executable,str(root/'scripts/summarize_transport_v2.py'),'--manifest',str(manifest),'--results-root',str(results),'--output',str(output)],check=True)
 li=next(r for r in csv.DictReader(output.open()) if r['observable']=='Li');assert li['status']=='incomplete' and li['mean']==''
