#!/usr/bin/env python3
"""Unified cached trajectory analysis: STACIE conductivity and replica-level MSD QC."""
import argparse,csv,hashlib,json,sys,traceback,warnings
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
E=1.602176634e-19;KB=1.380649e-23

def block_current(positions,dt_ps,charge_moment=False):
 x=np.asarray(positions,float)
 if x.ndim not in (2,3) or not np.isfinite(x).all() or dt_ps<=0:raise ValueError('invalid position/time array')
 diff=np.diff(x,axis=0)*1e-10/(dt_ps*1e-12)
 if charge_moment:diff*=E
 return diff.reshape(len(diff),-1).T

def summarize_values(values,expected=4):
 from scipy.stats import t
 n=len(values)
 if n==0 or any(x is None or not np.isfinite(x) for x in values):return {'n':n,'complete':False,'mean':None,'sd':None,'ci95':None}
 mean=float(np.mean(values));sd=float(np.std(values,ddof=1)) if n>1 else None
 error=float(t.ppf(.975,n-1)*sd/np.sqrt(n)) if n>1 else None
 return {'n':n,'complete':n==expected,'mean':mean,'sd':sd,'ci95':[mean-error,mean+error] if error is not None else None,'uncertainty':'Student-t across independent replicas; no filtering of successful values'}

def find_plateau(t,y,config):
 t=np.asarray(t,float);y=np.asarray(y,float)
 if len(t)!=len(y) or len(t)<2 or not np.isfinite(t).all() or not np.isfinite(y).all() or np.any(np.diff(t)<=0):return {'status':'invalid_curve'}
 maxlag=t[-1]*config['plateau_max_lag_fraction'];mask=(t>0)&(t<=maxlag)&(y>0)
 if mask.sum()<40 or maxlag<=50:return {'status':'no_plateau'}
 grid=np.geomspace(50,maxlag,80);values=np.interp(grid,t[mask],y[mask]);a=np.full(80,np.nan)
 for i in range(3,77):a[i]=np.polyfit(np.log(grid[i-3:i+4]),np.log(values[i-3:i+4]),1)[0]
 good=(a>=config['plateau_alpha'][0])&(a<=config['plateau_alpha'][1]);runs=[];start=None
 for i,v in enumerate(np.r_[good,False]):
  if v and start is None:start=i
  elif not v and start is not None:
   end=i-1
   if grid[end]-grid[start]>=config['plateau_min_width_ps'] and grid[end]/grid[start]>=config['plateau_min_lag_ratio']:runs.append((start,end))
   start=None
 info={'local_lag_ps':grid.tolist(),'local_alpha':[float(x) if np.isfinite(x) else None for x in a]}
 if not runs:return dict(status='no_plateau',**info)
 begin,end=max(runs,key=lambda ij:grid[ij[1]]-grid[ij[0]])
 return dict(status='plateau_candidate',window_ps=[float(grid[begin]),float(grid[end])],**info)

def read_config(path):
 c=json.loads(path.read_text())
 if c.get('schema_version')!=1 or c.get('expected_replicas',0)<2:raise ValueError('invalid schema/replica count')
 lo,hi=c['msd_fallback_lag_fraction']
 if not 0<lo<hi<=.5 or not 0<c['plateau_max_lag_fraction']<=.5:raise ValueError('invalid lag fractions')
 if not 0<c['plateau_alpha'][0]<c['plateau_alpha'][1] or c['plateau_min_width_ps']<=0 or c['plateau_min_lag_ratio']<=1:raise ValueError('invalid plateau criteria')
 c['config_sha256']=hashlib.sha256(path.read_bytes()).hexdigest();return c

def source_signature(case,meta):
 files=[case/'metadata.json',case/'system.data']+[case/'traj'/f'chunk_{i:04d}.dcd' for i in range(1,int(meta['expected_chunks'])+1)]
 return [{'name':str(f.relative_to(case)),'size':f.stat().st_size,'mtime_ns':f.stat().st_mtime_ns,**({'sha256':hashlib.sha256(f.read_bytes()).hexdigest()} if f.name in ('metadata.json','system.data') else {})} for f in files]

def prepare_cache(case,cache):
 from transport_v2 import read_positions,temperature_diagnostic
 case=case.resolve();meta=json.loads((case/'metadata.json').read_text());sig=source_signature(case,meta)
 if cache.exists():
  with np.load(cache,allow_pickle=False) as z:
   if json.loads(str(z['source_signature']))!=sig:raise ValueError(f'stale cache {cache}; preserve it and select a new output root')
  return
 li,fsi,moment,dims,prov=read_positions(case,meta);prov['temperature_diagnostic']=temperature_diagnostic(case,meta);cache.parent.mkdir(parents=True,exist_ok=True);tmp=cache.with_suffix('.tmp.npz')
 np.savez_compressed(tmp,li=li,fsi=fsi,moment=moment,metadata=json.dumps(meta),provenance=json.dumps(prov),source_signature=json.dumps(sig));tmp.replace(cache)

def spectral_axes(moment,meta,dt,config,out):
 from importlib.metadata import version
 from stacie import compute_spectrum,estimate_acint,ExpPolyModel,plot_results,UnitConfig
 if version('stacie')!=config['stacie_version']:raise ValueError('STACIE version differs from locked configuration')
 currents=block_current(moment,dt,True);axes=[]
 uc=UnitConfig(acint_unit_str='S/m',time_unit=1e-12,time_unit_str='ps',freq_unit=1e12,freq_unit_str='THz')
 prefactor=1/(float(meta['box_volume_A3'])*1e-30*KB*float(meta['temperature_K']))
 for axis in range(3):
  try:
   with warnings.catch_warnings(record=True) as ws:
    warnings.simplefilter('always');s=compute_spectrum(currents[axis:axis+1],prefactors=prefactor,timestep=dt*1e-12,include_zero_freq=False)
    r=estimate_acint(s,ExpPolyModel(config['spectral_powers']),rng=np.random.default_rng(config['random_seed']))
   npars=len(r.props['pars']);qc=[]
   if r.neff<20*npars:qc.append('too_few_effective_frequencies')
   if abs(r.props['cost_zscore'])>2 or abs(r.props['criterion_zscore'])>2:qc.append('spectral_model_warning')
   if r.acint<=0 or r.acint_std/r.acint>config['relative_spectral_error_warning']:qc.append('large_relative_uncertainty')
   axes.append(dict(axis='xyz'[axis],value_mS_cm=float(r.acint)*10,std_mS_cm=float(r.acint_std)*10,neff=float(r.neff),cost_z=float(r.props['cost_zscore']),criterion_z=float(r.props['criterion_zscore']),qc=qc,warnings=[str(w.message) for w in ws]));plot_results(str(out/f'spectrum_{axis}.pdf'),r,uc=uc)
  except Exception as exc:axes.append(dict(axis='xyz'[axis],error=f'{type(exc).__name__}: {exc}'))
 if any('error' in r for r in axes):return {'status':'failed','axes':axes,'value_mS_cm':None}
 # Sum of marginal standard errors /3 bounds the trace-average SE under unknown axis covariance.
 return {'status':'spectral_candidate' if all(not r['qc'] for r in axes) else 'needs_review','value_mS_cm':float(np.mean([r['value_mS_cm'] for r in axes])),'conditional_se_upper_bound_mS_cm':float(sum(r['std_mS_cm'] for r in axes)/3),'axes':axes,'note':'single-axis spectra; independent-replica uncertainty is reported separately'}

def validate_cached_arrays(li,fsi,moment,meta,prov):
 n=len(li);pairs=int(meta['actual_pairs']);dt=float(prov['dt_ps'])
 expected=round(int(meta['expected_chunks'])*float(meta['chunk_ps'])/dt) if dt>0 else -1
 if li.shape!=(n,pairs,3) or fsi.shape!=li.shape or moment.shape!=(n,3) or n!=expected or n<2:raise ValueError('cached array shapes or frame count disagree with metadata')
 if not all(np.isfinite(x).all() for x in (li,fsi,moment)):raise ValueError('non-finite cached coordinates')

def analyze(cache,out,config):
 from transport_v2 import curves,metrics
 out.mkdir(parents=True,exist_ok=False)
 try:
  with np.load(cache,allow_pickle=False) as z:
   meta=json.loads(str(z['metadata']));prov=json.loads(str(z['provenance']));li=z['li'];fsi=z['fsi'];moment=z['moment']
  validate_cached_arrays(li,fsi,moment,meta,prov)
  dt=float(prov['dt_ps']);t=np.arange(len(li))*dt;window=[float(x)*t[-1] for x in config['msd_fallback_lag_fraction']];data={};result={'schema_version':1,'config':config,'metadata':meta,'provenance':prov,'cache':str(cache),'msd_window_role':'prespecified diagnostic; final common plateau fitted at state level'}
  for name,x,kind in [('Li',li,'self'),('FSI_COM',fsi,'self'),('collective_charge',moment,'collective')]:
   y=curves(x,kind);data[name]=y;result[name]={'diagnostic_window_ps':window,'diagnostic_fit':metrics(t,y,window)}
  result['conductivity']=spectral_axes(moment,meta,dt,config,out)
  temp=prov.get('temperature_diagnostic',{});result['quality_alerts']=[]
  if temp.get('first_status')!='pass' or temp.get('second_status')!='pass':result['quality_alerts'].append('temperature_not_verified')
  result['quality_alerts'].append('energy_structure_and_equilibration_require_review')
  import pandas as pd
  pd.DataFrame({'lag_ps':t,**data}).to_csv(out/'curves.csv',index=False)
  import matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt
  fig,axs=plt.subplots(3,2,figsize=(11,10),constrained_layout=True)
  for i,(name,y) in enumerate(data.items()):
   mask=(t>0)&(y>0)&(t<=window[1]);m=result[name]['diagnostic_fit']
   for ax in axs[i]:ax.plot(t[mask],y[mask]);ax.axvspan(*window,color='orange',alpha=.2);ax.set_xlabel('Lag (ps)');ax.set_ylabel(name+' MSD')
   axs[i,0].set_xscale('log');axs[i,0].set_yscale('log');axs[i,0].set_title(f"Diagnostic α={m.get('alpha',float('nan')):.3f}")
   if m['status']=='ok':
    select=(t>=window[0])&(t<=window[1]);axs[i,1].plot(t[select],m['slope']*t[select]+m['intercept'],color='red')
  fig.savefig(out/'msd_diagnostics.png',dpi=160);fig.savefig(out/'msd_diagnostics.pdf');plt.close(fig)
  tmp=out/'results.tmp';tmp.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n');tmp.replace(out/'results.json');return result
 except Exception:
  (out/'failure.json').write_text(json.dumps({'error':traceback.format_exc()},indent=2));raise

def summarize(manifest,output_root,config,dest):
 import pandas as pd
 from transport_v2 import metrics
 rows=list(csv.DictReader(manifest.open()));groups={};case_table=[];state_table=[];dest.mkdir(parents=True,exist_ok=False)
 for row in rows:groups.setdefault((row['base_name'],int(row['dmc_per_lifsi'])),[]).append(row)
 for (base,ratio),members in groups.items():
  available=[];errors=[]
  for row in members:
   folder=output_root/'cases'/f"{int(row['array_index']):03d}";f=folder/'results.json'
   if not f.exists():errors.append(f"{row['array_index']}:missing");continue
   j=json.loads(f.read_text());m=j['metadata']
   if j['config']['config_sha256']!=config['config_sha256'] or m['base_name']!=base or int(m['dmc_per_lifsi'])!=ratio or int(m['replica_index'])!=int(row['replica']):raise ValueError(f'case or configuration mismatch {f}')
   available.append((row,j,pd.read_csv(folder/'curves.csv')))
  if not available:continue
  ids=[a[1]['metadata']['velocity_seed'] for a in available]
  if len(ids)!=len(set(ids)):raise ValueError('duplicate velocity seeds')
  volumes=[float(j['metadata']['box_volume_A3']) for _,j,_ in available];temps=[float(j['metadata']['temperature_K']) for _,j,_ in available]
  if np.ptp(volumes)>max(volumes)*1e-5 or np.ptp(temps)>1e-3:raise ValueError('incompatible states')
  time=available[0][2]['lag_ps'].to_numpy();
  if any(not np.array_equal(time,d['lag_ps']) for _,_,d in available):raise ValueError('incompatible time grids')
  common={};statuses={};fit_results={};expected=config['expected_replicas'];complete=len(available)==expected and len(members)==expected
  temp_ok=all(j['provenance'].get('temperature_diagnostic',{}).get('first_status')=='pass' and j['provenance'].get('temperature_diagnostic',{}).get('second_status')=='pass' for _,j,_ in available)
  for name in ('Li','FSI_COM','collective_charge'):
   mean=np.mean([d[name].to_numpy() for _,_,d in available],axis=0);p=find_plateau(time,mean,config) if complete else {'status':'incomplete_replicas'};statuses[name]=p['status'];win=p.get('window_ps',[time[-1]*f for f in config['msd_fallback_lag_fraction']]);common[name]=win
   fit_results[name]=[metrics(time,d[name].to_numpy(),win) for _,_,d in available]
   if p['status']=='plateau_candidate':
    if any(m['status']!='ok' or not .8<=m['alpha']<=1.2 or m['slope']<=0 or m['r2']<.98 for m in fit_results[name]):statuses[name]='replica_disagreement'
   if statuses[name]=='plateau_candidate' and not temp_ok:statuses[name]='temperature_review_required'
   values=[m.get('slope') for m in fit_results[name]]
   factor=1e-4/6 if name!='collective_charge' else 10*E**2*1e-8/(6*volumes[0]*1e-30*KB*temps[0]);stats=summarize_values([v*factor if v is not None else None for v in values],expected)
   state_table.append(dict(mof=base,dmc_per_lifsi=ratio,observable=name+'_Einstein',status=statuses[name],lag_start_ps=win[0],lag_end_ps=win[1],**flat_stats(stats),unit='cm2/s' if name!='collective_charge' else 'mS/cm'))
  sigma=[j['conductivity']['value_mS_cm'] for _,j,_ in available]
  factor=10*E**2*1e-8/(6*volumes[0]*1e-30*KB*temps[0]);eh=[m.get('slope')*factor if m['status']=='ok' else None for m in fit_results['collective_charge']]
  paired=summarize_values([e-s if e is not None and s is not None else None for e,s in zip(eh,sigma)],expected);ci=paired['ci95'];consistent=bool(ci and ci[0]<=0<=ci[1])
  state_table.append(dict(mof=base,dmc_per_lifsi=ratio,observable='collective_STACIE',status='review_required' if not complete or not temp_ok or not consistent or statuses['collective_charge']!='plateau_candidate' or any(j['conductivity']['status']!='spectral_candidate' for _,j,_ in available) else 'spectral_candidate',**flat_stats(summarize_values(sigma,expected)),unit='mS/cm',paired_EH_minus_STACIE_mean=paired['mean'],paired_difference_ci95_low=ci[0] if ci else None,paired_difference_ci95_high=ci[1] if ci else None))
  ne_values=[]
  for i,(row,j,d) in enumerate(available):
   rec=dict(array_index=int(row['array_index']),mof=base,dmc_per_lifsi=ratio,replica=int(row['replica']),DMC_molecules=int(row['dmc_molecules']),LiFSI_pairs=int(row['lifsi_pairs']),sigma_STACIE_mS_cm=j['conductivity']['value_mS_cm'],spectral_status=j['conductivity']['status'],missing_replicas=';'.join(errors),config_sha256=config['config_sha256'])
   for name in ('Li','FSI_COM','collective_charge'):
    m=fit_results[name][i];win=common[name];factor=1e-4/6 if name!='collective_charge' else 10*E**2*1e-8/(6*volumes[0]*1e-30*KB*temps[0]);rec[name+'_diagnostic_value']=m['slope']*factor if m['status']=='ok' else None;rec[name+'_status']=statuses[name];rec[name+'_alpha']=m.get('alpha');rec[name+'_lag_start_ps']=win[0];rec[name+'_lag_end_ps']=win[1]
   if all(fit_results[k][i]['status']=='ok' for k in ('Li','FSI_COM')):
    ql=j['provenance']['li_charge_e'];qf=j['provenance']['fsi_charge_e'];N=j['metadata']['actual_pairs'];rec['sigma_NE_diagnostic_mS_cm']=10*N*E**2/(KB*temps[0]*volumes[0]*1e-30)*1e-4*(ql**2*rec['Li_diagnostic_value']+qf**2*rec['FSI_COM_diagnostic_value'])
   ne_values.append(rec.get('sigma_NE_diagnostic_mS_cm'));case_table.append(rec)
  state_table.append(dict(mof=base,dmc_per_lifsi=ratio,observable='NE_diagnostic',status='diffusion_review_required' if statuses['Li']!='plateau_candidate' or statuses['FSI_COM']!='plateau_candidate' else 'diffusion_candidate',**flat_stats(summarize_values(ne_values,expected)),unit='mS/cm'))
  plot_state(dest,base,ratio,time,available,common,statuses)
 for r in state_table:
  r['screened_mean']=r['mean'] if r['complete'] and r['status'] in ('plateau_candidate','spectral_candidate','diffusion_candidate') and r['ci95_low'] is not None and r['ci95_low']>0 else None
  r['quality_note']='descriptive mean includes every available replica; screened_mean is a statistical screen, not physical validation'
 pd.DataFrame(case_table).to_csv(dest/'per_case.csv',index=False);pd.DataFrame(state_table).to_csv(dest/'by_state.csv',index=False);(dest/'config.json').write_text(json.dumps(config,indent=2))

def flat_stats(s):return dict(n=s['n'],complete=s['complete'],mean=s['mean'],sample_sd=s['sd'],ci95_low=s['ci95'][0] if s['ci95'] else None,ci95_high=s['ci95'][1] if s['ci95'] else None)

def plot_state(dest,base,ratio,t,available,common,statuses):
 import matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt
 from transport_v2 import metrics
 fig,axs=plt.subplots(3,2,figsize=(12,11),constrained_layout=True)
 for i,name in enumerate(('Li','FSI_COM','collective_charge')):
  w=common[name];mask=(t>0)&(t<=t[-1]*.5)
  for _,j,d in available:
   y=d[name].to_numpy();good=mask&(y>0)
   for ax in axs[i]:ax.plot(t[good],y[good],alpha=.6);ax.set_xlabel('Lag (ps)');ax.set_ylabel(name+' MSD')
   m=metrics(t,y,w)
   if m['status']=='ok':fit=(t>=w[0])&(t<=w[1]);axs[i,1].plot(t[fit],m['slope']*t[fit]+m['intercept'],lw=1,color='black')
  for ax in axs[i]:ax.axvspan(*w,color='orange',alpha=.2)
  axs[i,0].set_xscale('log');axs[i,0].set_yscale('log');axs[i,0].set_title(name+' | '+statuses[name]);axs[i,1].set_title('One common lag window for all replicas')
 fig.savefig(dest/f'{base}_ratio{ratio:02d}.png',dpi=160);fig.savefig(dest/f'{base}_ratio{ratio:02d}.pdf');plt.close(fig)

def run_case_task(args):
 case,cache,out,config=args
 try:prepare_cache(case,cache);analyze(cache,out,config);return {'case':str(case),'status':'completed'}
 except Exception:return {'case':str(case),'status':'failed','error':traceback.format_exc()}

def merge_outputs(inputs,output,manifest,config):
 import shutil
 seen=set();planned=[]
 for source in inputs:
  for f in sorted((source/'cases').glob('*/results.json')):
   index=f.parent.name
   if index in seen:raise ValueError(f'duplicate case across machines: {index}')
   j=json.loads(f.read_text())
   if j['config']['config_sha256']!=config['config_sha256'] or j.get('schema_version')!=1:raise ValueError('method/config mismatch while merging')
   if not (f.parent/'curves.csv').is_file():raise ValueError('missing curves while merging')
   seen.add(index);planned.append((index,f.parent))
 if not planned:raise ValueError('no completed analyses to merge')
 output.mkdir(parents=True,exist_ok=False)
 for index,source in planned:
  dest=output/'cases'/index;dest.mkdir(parents=True)
  for name in ('results.json','curves.csv'):shutil.copy2(source/name,dest/name)
 summarize(manifest,output,config,output/'summary')

def main():
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('command',choices=('cache','analyze','batch','summarize','merge'));p.add_argument('--config',type=Path,default=ROOT/'config/transport.json');p.add_argument('--inputs',nargs='+',type=Path);p.add_argument('--case',type=Path);p.add_argument('--cache',type=Path);p.add_argument('--cache-root',type=Path);p.add_argument('--output',type=Path,required=True);p.add_argument('--manifest',type=Path,default=ROOT/'manifest.csv');p.add_argument('--runs-root',type=Path,default=ROOT/'runs');p.add_argument('--workers',type=int,default=1);p.add_argument('--resume',action='store_true');a=p.parse_args();config=read_config(a.config)
 if a.command=='cache' and a.case is None:p.error('--case required for cache')
 if a.command=='analyze' and a.cache is None:p.error('--cache required for analyze')
 if a.command=='cache':prepare_cache(a.case,a.output)
 elif a.command=='analyze':analyze(a.cache,a.output,config)
 elif a.command=='summarize':summarize(a.manifest,a.runs_root,config,a.output)
 elif a.command=='merge':
  if not a.inputs:p.error('--inputs required for merge')
  merge_outputs(a.inputs,a.output,a.manifest,config)
 else:
  import os,concurrent.futures
  if a.workers<1 or a.workers>int(os.environ.get('SLURM_CPUS_PER_TASK',a.workers)):raise ValueError('worker count exceeds allocation')
  if a.output.exists():
   if not a.resume:raise FileExistsError('output exists; use --resume to validate and continue')
   previous=json.loads((a.output/'config.json').read_text())
   if previous['config_sha256']!=config['config_sha256']:raise ValueError('configuration changed; use a new output root')
  else:
   a.output.mkdir(parents=True);(a.output/'config.json').write_text(json.dumps(config,indent=2))
  tasks=[];skipped=[]
  for r in csv.DictReader(a.manifest.open()):
   i=int(r['array_index']);case=a.runs_root/f'{i:03d}';files=[case/'traj'/f'chunk_{j:04d}.dcd' for j in range(1,16)]+[case/'state'/f'restart_{j:04d}.bin' for j in range(1,16)]
   if not all(x.is_file() and x.stat().st_size>0 for x in files):continue
   cache=(a.cache_root or a.output/'cache')/f'{i:03d}.npz';folder=a.output/'cases'/f'{i:03d}';result=folder/'results.json'
   if result.exists():
    old=json.loads(result.read_text());meta=old['metadata']
    if old['config']['config_sha256']!=config['config_sha256'] or meta['base_name']!=r['base_name'] or int(meta['replica_index'])!=int(r['replica']):raise ValueError('existing result identity/config mismatch')
    if not (folder/'curves.csv').is_file():raise ValueError('existing result lacks curves')
    prepare_cache(case,cache);skipped.append({'case':str(case),'status':'already_complete'});continue
   if folder.exists():
    import time
    archive=a.output/'failed_attempts'/f'{i:03d}_{time.time_ns()}';archive.parent.mkdir(parents=True,exist_ok=True);folder.rename(archive)
   tasks.append((case,cache,folder,config))
  if not tasks and not skipped:raise ValueError('no complete cases found')
  with concurrent.futures.ProcessPoolExecutor(max_workers=a.workers) as pool:
   results=skipped
   for r in pool.map(run_case_task,tasks):
    results.append(r);print(r['case'],r['status'],flush=True)
    tmp=a.output/'batch_status.tmp';tmp.write_text(json.dumps(results,indent=2));tmp.replace(a.output/'batch_status.json')
  (a.output/'batch_status.json').write_text(json.dumps(results,indent=2))
  summary=a.output/'summary'
  if summary.exists():
   import time
   summary=a.output/f'summary_{time.time_ns()}'
  summarize(a.manifest,a.output,config,summary);(a.output/'latest_summary.txt').write_text(summary.name+'\n');print('SUMMARY',summary,flush=True)
  if any(r['status']=='failed' for r in results):raise SystemExit(1)
if __name__=='__main__':main()
