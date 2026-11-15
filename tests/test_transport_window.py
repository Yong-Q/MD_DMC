import sys
from pathlib import Path
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from robust_transport import select_window, displacement_msd, conductivity_from_slope

def test_msd_matches_direct():
    rng=np.random.default_rng(123);x=np.cumsum(rng.normal(size=(35,3)),axis=0)
    actual=displacement_msd(x)
    expected=np.array([np.mean(np.sum((x[k:]-x[:len(x)-k])**2,axis=1)) for k in range(len(x))])
    np.testing.assert_allclose(actual,expected,atol=1e-10)

def test_selects_linear_and_rejects_ballistic():
    t=np.arange(0,15001,dtype=float)
    result=select_window(t,3*t,15000)
    assert result['status']=='selected' and abs(result['alpha']-1)<1e-8 and abs(result['slope']-3)<1e-8
    assert select_window(t,t*t,15000)['status']=='no_diffusive_window'
    assert select_window(t,np.zeros_like(t),15000)['status']=='no_diffusive_window'

def test_rejects_short_and_nonfinite():
    t=np.arange(30.);assert select_window(t,t,30)['status']=='no_diffusive_window'
    t=np.arange(1000.);y=t.copy();y[200]=np.nan
    assert select_window(t,y,1000)['status']=='invalid_curve'

def test_si_units():
    e=1.602176634e-19;k=1.380649e-23
    assert np.isclose(conductivity_from_slope(1,1000,300), e**2*1e-8/(6*1000*1e-30*k*300))

def test_summary_checks_independent_replica_identity(tmp_path):
    import csv,json,subprocess
    root=Path(__file__).resolve().parents[1]; manifest=tmp_path/'manifest.csv';results=tmp_path/'results';results.mkdir()
    rows=[dict(array_index='0',base_name='MOF',dmc_per_lifsi='2',replica='1'),dict(array_index='1',base_name='MOF',dmc_per_lifsi='2',replica='2')]
    with manifest.open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=rows[0]);w.writeheader();w.writerows(rows)
    for i,val in enumerate((1.,3.)):
        d=results/f'{i:03d}';d.mkdir();(d/'results.json').write_text(json.dumps(dict(base_name='MOF',dmc_per_lifsi=2,replica_index=i+1,provenance={'source':'raw DCD'},collective_charge={'conductivity_mS_cm':val})))
    out=tmp_path/'summary.csv';subprocess.run([sys.executable,str(root/'scripts/summarize_transport.py'),'--manifest',str(manifest),'--results-root',str(results),'--output',str(out)],check=True)
    r=next(csv.DictReader(out.open()));assert r['n_qualified']=='2' and float(r['conductivity_mean_mS_cm'])==2 and np.isclose(float(r['conductivity_sample_std_mS_cm']),np.sqrt(2))
    bad=results/'001/results.json';j=json.loads(bad.read_text());j['replica_index']=1;bad.write_text(json.dumps(j))
    proc=subprocess.run([sys.executable,str(root/'scripts/summarize_transport.py'),'--manifest',str(manifest),'--results-root',str(results),'--output',str(tmp_path/'bad.csv')],capture_output=True)
    assert proc.returncode!=0 and b'mismatch' in proc.stderr
