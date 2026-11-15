import sys
from pathlib import Path
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from transport_v5 import select_stable

def test_linear_diffusion_correct():
 t=np.arange(15000.);r=select_stable(t,2*t+3)
 assert r['status']=='candidate' and abs(r['fit']['slope']-2)<1e-9 and r['window_ps'][1]<=.5*t[-1]

def test_reject_anomalous_diffusion():
 t=np.arange(15000.)
 for y in (np.sqrt(t+1),t*t+1,np.ones_like(t)):
  assert select_stable(t,y)['status']=='no_candidate'

def test_short_curve_and_invalid_data():
 t=np.arange(1000.);assert select_stable(t,t+1)['status']=='no_candidate'
 t=np.arange(15000.);y=t+1;y[300]=np.nan;assert select_stable(t,y)['status']=='invalid_curve'

def test_end_to_end_uses_full_trajectory(tmp_path,monkeypatch):
 import json
 import transport_v5 as m
 case=tmp_path/'case';case.mkdir();out=tmp_path/'out';n=6000
 meta=dict(frame_dt_ps=1.,box_volume_A3=100000.,temperature_K=298.,actual_pairs=4,case_name='synthetic',base_name='SYN',dmc_per_lifsi=2,replica_index=1,velocity_seed=123)
 (case/'metadata.json').write_text(json.dumps(meta));x=np.zeros((n,4,3));moment=np.zeros((n,3));dims=np.ones((n,6));seen=[]
 monkeypatch.setattr(m,'read_positions',lambda case,meta:(x,x,moment,dims,dict(dt_ps=1.,li_charge_e=.7,fsi_charge_e=-.7)))
 def fake_curve(a,kind):seen.append(len(a));return np.arange(len(a),dtype=float)+1
 monkeypatch.setattr(m,'curves',fake_curve);monkeypatch.setattr(sys,'argv',['transport_v5.py','--case',str(case),'--output-dir',str(out)]);m.main();j=json.loads((out/'results.json').read_text())
 assert seen==[n,n,n] and j['Li']['status']=='candidate' and j['collective_charge']['conductivity_mS_cm']>0
 assert (out/'diagnostics.png').stat().st_size>1000

def test_export_reports_missing_and_rejects_wrong_case(tmp_path):
 import csv,json,subprocess
 root=Path(__file__).resolve().parents[1];manifest=tmp_path/'manifest.csv';r=dict(array_index='2',base_name='SYN',qmof_id='q',replica='1',dmc_per_lifsi='2',lifsi_pairs='50',dmc_molecules='100',effective_pore_A3='100000',salt_concentration_per_pore_M='1')
 with manifest.open('w',newline='') as f:w=csv.DictWriter(f,fieldnames=r);w.writeheader();w.writerow(r)
 outputs=tmp_path/'cases';outputs.mkdir();dest=tmp_path/'out.csv';cmd=[sys.executable,str(root/'scripts/export_transport_v5.py'),'--manifest',str(manifest),'--results-root',str(outputs),'--output',str(dest)];subprocess.run(cmd,check=True);row=next(csv.DictReader(dest.open()));assert row['analysis_status']=='not_analyzed' and row['collective_sigma_mS_cm']==''
 d=outputs/'002';d.mkdir();(d/'results.json').write_text(json.dumps({'base_name':'WRONG','method_version':'v5 x'}));cmd[-1]=str(tmp_path/'wrong.csv');p=subprocess.run(cmd,capture_output=True);assert p.returncode!=0 and b'mismatch' in p.stderr

def test_crossover_can_mimic_diffusion_and_requires_review():
 # A sum of sub- and super-linear power laws is not globally linear,
 # yet can pass any local slope screen around a crossover.
 t=np.arange(15000.);z=(t+1)/1500.;r=select_stable(t,z**.2+z**1.8)
 assert r['status']=='candidate' and abs(r['fit']['alpha']-1)<.1

def test_batch_only_lists_complete_and_unprocessed(tmp_path):
 import csv,shutil,subprocess
 root=Path(__file__).resolve().parents[1];(tmp_path/'scripts').mkdir();shutil.copy(root/'scripts/list_completed_for_v5.py',tmp_path/'scripts')
 with (tmp_path/'manifest.csv').open('w',newline='') as f:
  w=csv.DictWriter(f,fieldnames=['array_index']);w.writeheader();w.writerows([{'array_index':0},{'array_index':1}])
 for i in (0,1):
  d=tmp_path/'runs'/f'{i:03d}';(d/'traj').mkdir(parents=True);(d/'state').mkdir()
  for n in range(1,16):
   (d/'traj'/f'chunk_{n:04d}.dcd').write_bytes(b'x');(d/'state'/f'restart_{n:04d}.bin').write_bytes(b'x')
 (tmp_path/'analysis_v5'/'001').mkdir(parents=True)
 p=subprocess.run([sys.executable,str(tmp_path/'scripts/list_completed_for_v5.py')],capture_output=True,text=True,check=True)
 assert p.stdout.strip()=='0' and 'Skipping 1' in p.stderr

def test_conductivity_conversion_matches_si():
 from robust_transport import conductivity_from_slope
 e=1.602176634e-19;k=1.380649e-23
 assert np.isclose(conductivity_from_slope(1,1000,300), e**2*1e-8/(6*1000e-30*k*300))
