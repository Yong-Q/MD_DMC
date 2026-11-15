#!/usr/bin/env python3
"""Predeclared-window exploratory transport; held-out time-block checks."""
import argparse,json,math,sys,traceback
from pathlib import Path
import numpy as np
from robust_transport import displacement_msd,conductivity_from_slope

# In ps. Predeclared before looking at each case; a single discovery-selected
# window is evaluated on holdout, with no retry after failure.
WINDOWS=((750.,1500.),(500.,1000.),(250.,750.),(100.,500.))

def metrics(t,y,window):
 t=np.asarray(t,float);y=np.asarray(y,float);start,end=window
 mask=(t>=start)&(t<=end);x=t[mask];v=y[mask]
 if len(x)<20 or not np.isfinite(v).all() or np.any(v<=0):return {'status':'invalid'}
 sample=np.linspace(np.log(x[0]),np.log(x[-1]),min(60,len(x)))
 alpha=float(np.polyfit(sample,np.interp(sample,np.log(x),np.log(v)),1)[0]);slope,intercept=np.polyfit(x,v,1)
 den=np.sum((v-v.mean())**2);r2=float(1-np.sum((v-(slope*x+intercept))**2)/den) if den>0 else -1.
 if not all(np.isfinite(z) for z in (alpha,slope,intercept,r2)):return {'status':'invalid'}
 return {'status':'ok','alpha':alpha,'slope':float(slope),'intercept':float(intercept),'r2':r2}

def passes(m):return m['status']=='ok' and m['slope']>0 and .9<=m['alpha']<=1.1 and m['r2']>=.98

def choose_and_validate(t,discovery,holdout):
 t=np.asarray(t,float);maxlag=t[-1]*.2; assessments=[]
 for win in WINDOWS:
  if win[1]>maxlag+1:continue
  d=metrics(t,discovery,win);assessments.append({'window_ps':list(win),'discovery':d})
  if not passes(d):continue
  h=metrics(t,holdout,win);valid=passes(h)
  ratio=h['slope']/d['slope'] if h['status']=='ok' and d['slope']>0 else None
  if ratio is None or not .5<=ratio<=2.:valid=False
  return {'status':'candidate' if valid else 'holdout_rejected','window_ps':list(win),'discovery':d,'holdout':h,'holdout_slope_ratio':ratio,'discovery_assessments':assessments}
 return {'status':'no_candidate','discovery_assessments':assessments}

def summarize_group(rows,expected=4):
 valid=[r for r in rows if r.get('value') is not None and np.isfinite(r['value'])]
 ids=[r['replica'] for r in valid];seeds=[r.get('seed') for r in valid]
 if len(ids)!=expected or len(set(ids))!=expected or None in seeds or len(set(seeds))!=expected:return {'status':'incomplete','n_valid':len(ids),'n_expected':expected}
 vals=[r['value'] for r in valid]
 return {'status':'complete','n_valid':len(vals),'n_expected':expected,'mean':float(np.mean(vals)),'sample_std':float(np.std(vals,ddof=1)) if expected>1 else None}

def temperature_diagnostic(case,meta):
 import csv
 try:
  n=int(meta['expected_chunks']);expected=round(n*float(meta['chunk_ps'])/float(meta['frame_dt_ps']));vals=[]
  for i in range(1,n+1):
   path=case/'logs'/f'thermo_chunk_{i:04d}.csv'
   with path.open(newline='') as f:
    for row in csv.DictReader(f):
     if row.get('temp_K'):vals.append(float(row['temp_K']))
  if len(vals)<.8*expected or not np.isfinite(vals).all():return {'status':'incomplete_or_invalid','samples':len(vals),'expected':expected}
  half=len(vals)//2;target=float(meta['temperature_K']);first=float(np.mean(vals[:half]));second=float(np.mean(vals[half:]));return {'status':'checked','samples':len(vals),'first_mean_K':first,'second_mean_K':second,'target_K':target,'first_status':'pass' if abs(first-target)<=10 else 'outside_10K','second_status':'pass' if abs(second-target)<=10 else 'outside_10K','note':'mean-temperature check only; no energy, structure, equilibration or drift proof'}
 except (OSError,ValueError,KeyError,ZeroDivisionError) as exc:return {'status':'unavailable','reason':str(exc)}

def read_positions(case,meta):
 import MDAnalysis as mda
 from MDAnalysis.transformations import NoJump
 n=int(meta['expected_chunks']); dt=float(meta['frame_dt_ps']);duration=float(meta['chunk_ps']);trajs=[case/'traj'/f'chunk_{i:04d}.dcd' for i in range(1,n+1)]
 if n!=15 or any(not x.is_file() or x.stat().st_size==0 for x in trajs):raise ValueError('requires 15 nonempty trajectory chunks')
 if any(not (case/'state'/f'restart_{i:04d}.bin').is_file() for i in range(1,n+1)):raise ValueError('missing completed restart')
 u=mda.Universe(str(case/'system.data'),[str(x) for x in trajs],format='LAMMPS',atom_style='id resid type charge x y z');u.trajectory.add_transformations(NoJump())
 types=np.asarray(u.atoms.types,dtype=int);li=u.atoms[types==int(meta['li_type'])];anchors=u.atoms[types==int(meta['fsi_anchor_type'])];count=int(meta['actual_pairs'])
 if len(li)!=count or len(anchors)!=count:raise ValueError('Li/FSI count mismatch')
 indices=[];weights=[];charges=[]
 for a in anchors:
  group=a.residue.atoms; ts=np.asarray(group.types,dtype=int)
  if len(group)!=9 or not np.all(np.isin(ts,meta['guest_types']['fsi_all'])):raise ValueError('cannot identify 9-atom FSI molecules')
  indices.append(group.indices);weights.append(group.masses/group.masses.sum());charges.append(group.charges.sum())
 indices=np.asarray(indices);weights=np.asarray(weights,float)
 if len(np.unique(indices))!=indices.size:raise ValueError('FSI molecules overlap')
 qli=np.asarray(li.charges,float);qfs=np.asarray(charges,float);q=np.asarray(u.atoms[int(meta['framework_atom_count']):].charges,float)
 if abs(q.sum())>1e-3 or np.ptp(qli)>1e-5 or np.ptp(qfs)>1e-5 or abs(qli.mean()+qfs.mean())>1e-3:raise ValueError('ion or mobile charge check failed')
 frames=len(u.trajectory);expected=round(n*duration/dt)
 if frames!=expected or frames%2:raise ValueError(f'frame count {frames}, expected even {expected}')
 apos=np.empty((frames,count,3));bpos=np.empty_like(apos);moment=np.empty((frames,3));dimensions=np.empty((frames,6));mobile=u.atoms[int(meta['framework_atom_count']):]
 for k,ts in enumerate(u.trajectory):
  apos[k]=li.positions;pos=u.atoms.positions[indices];bpos[k]=np.einsum('ijk,ij->ik',pos,weights);moment[k]=np.einsum('i,ij->j',q,mobile.positions);dimensions[k]=ts.dimensions
 if not np.allclose(dimensions,dimensions[0],rtol=1e-6,atol=1e-4):raise ValueError('variable cell not supported')
 return apos,bpos,moment,dimensions,{'n_frames':frames,'dt_ps':dt,'li_charge_e':float(qli.mean()),'fsi_charge_e':float(qfs.mean()),'mobile_charge_e':float(q.sum()),'fsi_representation':'9-atom mass-weighted COM','frame_count_check':'passed','pbc_continuity':'minimum-image assumption, not independently proven','stationarity':'not independently established'}

def curves(x,kind):
 import tidynamics
 if kind=='collective':return displacement_msd(x)
 return np.mean([tidynamics.msd(x[:,i,:]) for i in range(x.shape[1])],axis=0)

def kinisi_check(li,dimensions,meta,window):
 """Optional numerical sensitivity diagnostic, not an automatic correctness proof."""
 import MDAnalysis as mda
 from MDAnalysis.coordinates.memory import MemoryReader
 import scipp as sc
 from kinisi.analyze import DiffusionAnalyzer
 u=mda.Universe.empty(li.shape[1],trajectory=True);u.add_TopologyAttr('types',[str(meta['li_type'])]*li.shape[1]);u.load_new(li.astype(np.float32),format=MemoryReader,dimensions=dimensions.astype(np.float32),dt=float(meta['frame_dt_ps']))
 samples=[];errors=[]
 for n in (16,31):
  try:
   lags=np.unique(np.rint(np.linspace(window[0],window[1],n)/float(meta['frame_dt_ps'])).astype(int))*float(meta['frame_dt_ps'])
   a=DiffusionAnalyzer.from_universe(u,specie=str(meta['li_type']),time_step=float(meta['frame_dt_ps'])*sc.Unit('ps'),step_skip=1*sc.Unit('dimensionless'),dt=sc.array(dims=['time interval'],values=lags,unit='ps'),progress=False)
   for seed in (101,202):
    a.diffusion(window[0]*sc.Unit('ps'),recondition=True,n_samples=600,n_burn=400,n_thin=5,n_walkers=32,progress=False,random_state=np.random.RandomState(seed));v=np.asarray(a.D.values,float).ravel();condition=float(np.linalg.cond(a.diff.covariance_matrix.values))
    if len(v)<100 or not np.isfinite(v).all():raise ValueError('invalid posterior')
    samples.append({'lag_points':n,'seed':seed,'mean_cm2_s':float(v.mean()),'ci95_cm2_s':[float(vv) for vv in np.percentile(v,[2.5,97.5])],'covariance_condition':condition})
  except Exception as exc:errors.append(f'{type(exc).__name__}: {exc}')
 means=[r['mean_cm2_s'] for r in samples];ratio=max(means)/min(means) if means and min(means)>0 else None
 return {'status':'numerically_stable' if len(samples)==4 and ratio is not None and ratio<=1.5 and max(x['covariance_condition'] for x in samples)<1e10 else 'unstable_or_failed','runs':samples,'max_mean_ratio':ratio,'errors':errors,'limitations':'posterior conditional on model/window; no MCMC convergence or coverage certificate'}

def main():
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--case',type=Path,required=True);p.add_argument('--output-dir',type=Path,required=True);p.add_argument('--kinisi',action='store_true');p.add_argument('--mode',choices=('first-half','verified'),default='first-half');a=p.parse_args();case=a.case.resolve();out=a.output_dir.resolve();out.mkdir(parents=True,exist_ok=False)
 try:
  meta=json.loads((case/'metadata.json').read_text());li,fsi,moment,dimensions,provenance=read_positions(case,meta);dt=provenance['dt_ps'];split=len(li)//2;t=np.arange(split)*dt;tf=np.arange(len(li))*dt
  arrays={'Li':(li,'self'),'FSI_COM':(fsi,'self'),'collective_charge':(moment,'collective')};results={};allcurves={}
  for name,(x,kind) in arrays.items():
   first=curves(x[:split],kind);second=curves(x[split:],kind);full=curves(x,kind);choice=choose_and_validate(t,first,second);allcurves[name]=(first,second,full)
   choice['holdout_status']=choice['status']
   if 'window_ps' in choice and (a.mode=='first-half' or choice['status']=='candidate'):
    m=choice['discovery'] if a.mode=='first-half' else metrics(tf,full,choice['window_ps']);choice['fit']=m
    if m['status']!='ok' or m['slope']<=0:choice['status']='fit_invalid'
    else:
     if a.mode=='first-half':choice['status']='first_half_provisional'
     if name=='collective_charge':choice['conductivity_mS_cm']=10*conductivity_from_slope(m['slope'],float(meta['box_volume_A3']),float(meta['temperature_K']))
     else:choice['diffusion_cm2_s']=m['slope']*1e-4/6
   results[name]=choice
  if results['Li']['status'] in ('candidate','first_half_provisional') and results['FSI_COM']['status'] in ('candidate','first_half_provisional'):
   ql=provenance['li_charge_e'];qf=provenance['fsi_charge_e'];factor=10*int(meta['actual_pairs'])*1.602176634e-19**2*1e-8/(6*float(meta['box_volume_A3'])*1e-30*1.380649e-23*float(meta['temperature_K']))
   results['nernst_einstein_mS_cm']=factor*(ql**2*results['Li']['fit']['slope']+qf**2*results['FSI_COM']['fit']['slope'])
  if a.kinisi and results['Li']['status'] in ('candidate','first_half_provisional'):
   kl=li[:split] if a.mode=='first-half' else li;kd=dimensions[:split] if a.mode=='first-half' else dimensions
   results['Li']['kinisi_crosscheck']=kinisi_check(kl,kd,meta,results['Li']['window_ps'])
  elif a.kinisi:results['Li']['kinisi_crosscheck']={'status':'not_run_no_candidate'}
  results.update(temperature_diagnostic=temperature_diagnostic(case,meta),case_name=meta['case_name'],base_name=meta['base_name'],dmc_per_lifsi=meta['dmc_per_lifsi'],replica_index=meta['replica_index'],velocity_seed=meta['velocity_seed'],temperature_K=meta['temperature_K'],volume_A3=meta['box_volume_A3'],provenance=provenance,analysis_mode=a.mode,method_version='v2 predeclared discovery/holdout exploratory',thresholds={'windows_ps':WINDOWS,'log_slope':[.9,1.1],'r2_min':.98,'holdout_slope_ratio':[.5,2.],'max_lag_fraction_each_half':.2},warning='first-half estimates are provisional; block differences may indicate state changes, sampling noise or analysis artifacts; stationarity and model accuracy are not established')
  import pandas as pd
  df={'lag_ps':tf}
  for name,(first,second,full) in allcurves.items():df[name+'_full']=full;df[name+'_discovery']=np.pad(first,(0,len(full)-len(first)),constant_values=np.nan);df[name+'_holdout']=np.pad(second,(0,len(full)-len(second)),constant_values=np.nan)
  pd.DataFrame(df).to_csv(out/'curves.csv',index=False)
  import matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt
  fig,axes=plt.subplots(3,2,figsize=(12,11),constrained_layout=True)
  for i,(name,(first,second,full)) in enumerate(allcurves.items()):
   result=results[name];mask=(t>0)&(t<=1500)&(first>0)&(second>0)
   for ax in axes[i]:
    ax.plot(t[mask],first[mask],label='discovery');ax.plot(t[mask],second[mask],label='holdout');ax.set_xlabel('Lag (ps)');ax.set_ylabel(name+' MSD');ax.grid(alpha=.2)
    if 'window_ps' in result:ax.axvspan(*result['window_ps'],color='orange',alpha=.25)
   axes[i,0].set_xscale('log');axes[i,0].set_yscale('log')
   title=name+' | '+result['status']
   if result.get('discovery',{}).get('status')=='ok' and result.get('holdout',{}).get('status')=='ok':title+=f" | α: {result['discovery']['alpha']:.2f} / {result['holdout']['alpha']:.2f}"
   axes[i,0].set_title(title);axes[i,0].legend(fontsize=8)
   if result['status'] in ('candidate','first_half_provisional'):
    m=result['fit'];win=result['window_ps'];maskfit=(t>=win[0])&(t<=win[1]);axes[i,1].plot(t[maskfit],m['slope']*t[maskfit]+m['intercept'],color='crimson',label='first-half fit' if a.mode=='first-half' else 'full fit')
    value=(f"σ={result['conductivity_mS_cm']:.3g} mS/cm" if name=='collective_charge' else f"D={result['diffusion_cm2_s']:.3g} cm²/s");axes[i,1].set_title(('FIRST-HALF provisional | ' if a.mode=='first-half' else 'Full | ')+value);axes[i,1].legend(fontsize=8)
   else:axes[i,1].set_title('No screened estimate')
  fig.savefig(out/'diagnostics.png',dpi=160);fig.savefig(out/'diagnostics.pdf');plt.close(fig)
  tmp=out/'results.json.tmp';tmp.write_text(json.dumps(results,indent=2,allow_nan=False)+'\n');tmp.replace(out/'results.json')
  print(json.dumps({k:results[k]['status'] for k in arrays}),flush=True)
 except Exception:
  (out/'failure.json').write_text(json.dumps({'status':'error','traceback':traceback.format_exc()},indent=2)+'\n');raise
if __name__=='__main__':main()
