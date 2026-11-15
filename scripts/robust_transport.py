#!/usr/bin/env python3
"""Exploratory transport estimates with explicit diffusive-window rejection."""
import argparse, json, sys, hashlib
from pathlib import Path
import numpy as np

E=1.602176634e-19
KB=1.380649e-23

def displacement_msd(x):
    x=np.asarray(x,dtype=float);n=len(x)
    if x.ndim!=2 or n<2 or not np.isfinite(x).all(): raise ValueError('invalid displacement series')
    # Translation reduces FFT cancellation without altering displacement.
    x=x-x[0]; nfft=1<<(2*n-1).bit_length()
    corr=np.fft.irfft(np.sum(abs(np.fft.rfft(x,n=nfft,axis=0))**2,axis=1),n=nfft)[:n]
    sq=np.sum(x*x,axis=1);prefix=np.r_[0.,np.cumsum(sq)];lag=np.arange(n)
    out=(prefix[n-lag]+prefix[n]-prefix[lag]-2*corr)/(n-lag);out[0]=0
    return np.maximum(out,0.)

def conductivity_from_slope(slope, volume_A3, temperature_K):
    return slope*E**2*1e-8/(6*volume_A3*1e-30*KB*temperature_K)

def select_window(t,y,total_ps):
    """Deterministic exploratory gate; widths/origins and failures remain visible."""
    t=np.asarray(t,float);y=np.asarray(y,float)
    if len(t)!=len(y) or len(t)<2 or not np.isfinite(t).all() or not np.isfinite(y).all() or not np.all(np.diff(t)>0):return {'status':'invalid_curve'}
    # At most 20% of trajectory, avoid overselling long-lag sample counts.
    maxlag=min(3000.,float(total_ps)*0.2,t[-1]); minstart=50.; minwidth=200.
    candidates=[]
    for start in np.arange(minstart,maxlag-minwidth+1,50.):
        for end in np.arange(start+minwidth,maxlag+1,50.):
            idx=np.flatnonzero((t>=start)&(t<=end));
            if len(idx)<20 or np.any(y[idx]<=0):continue
            # Uniform log-time sample: keep the whole interval, not just the best local point.
            n=min(60,len(idx));logt=np.linspace(np.log(t[idx[0]]),np.log(t[idx[-1]]),n)
            yl=np.interp(logt,np.log(t[idx]),np.log(y[idx])); alpha=float(np.polyfit(logt,yl,1)[0])
            if not .8<=alpha<=1.2:continue
            a,b=np.polyfit(t[idx],y[idx],1); pred=a*t[idx]+b;den=np.sum((y[idx]-np.mean(y[idx]))**2)
            r2=float(1-np.sum((y[idx]-pred)**2)/den) if den>0 else -1
            half=len(idx)//2; a1=np.polyfit(t[idx[:half]],y[idx[:half]],1)[0];a2=np.polyfit(t[idx[half:]],y[idx[half:]],1)[0]
            sensitivity=float(abs(a1-a2)/abs(a)) if a>0 else float('inf')
            # Half-window positive slope prevents a single averaged crossing from qualifying.
            if a<=0 or a1<=0 or a2<=0 or r2<.98 or sensitivity>.5:continue
            candidates.append({'start_ps':float(start),'end_ps':float(end),'alpha':alpha,'slope':float(a),'intercept':float(b),'r2':r2,'half_slope_relative_difference':sensitivity})
    if not candidates:return {'status':'no_diffusive_window','max_lag_ps':float(maxlag),'n_candidates':0}
    # Longest interval, then latest start; no conductivity optimisation.
    best=max(candidates,key=lambda c:(c['end_ps']-c['start_ps'],c['start_ps']))
    best.update(status='selected',max_lag_ps=float(maxlag),n_candidates=len(candidates),early_window_warning=best['start_ps']<200.)
    return best

def read_curves(case, cached):
    import pandas as pd
    meta=json.loads((case/'metadata.json').read_text())
    if cached:
        dfs=[pd.read_csv(cached/name) for name in ('li_msd.csv','fsi_msd.csv','charge_moment_msd.csv')]
        t=dfs[0].iloc[:,0].to_numpy();ys=[d.iloc[:,1].to_numpy() for d in dfs]
        if any(len(d)!=len(t) or not np.allclose(d.iloc[:,0],t) for d in dfs):raise ValueError('cached time grids differ')
        return meta,t,ys,{'source':'existing analysis curves','warning':'cached curve provenance/charge convention must be verified by user'}
    import MDAnalysis as mda
    from MDAnalysis.transformations import NoJump
    import tidynamics
    count=int(meta['expected_chunks']);trajs=[case/'traj'/f'chunk_{i:04d}.dcd' for i in range(1,count+1)]
    if count!=15 or any(not p.is_file() or p.stat().st_size==0 for p in trajs):raise ValueError('need 15 non-empty ordered production chunks')
    if any(not (case/'state'/f'restart_{i:04d}.bin').is_file() for i in range(1,count+1)):raise ValueError('missing completed restart: do not analyze incomplete chunks')
    u=mda.Universe(str(case/'system.data'),[str(p) for p in trajs],format='LAMMPS'); u.trajectory.add_transformations(NoJump())
    alltypes=np.asarray(u.atoms.types,dtype=int);liq=np.asarray(u.atoms.charges,float);li=u.atoms[alltypes==int(meta['li_type'])];fsi=u.atoms[alltypes==int(meta['fsi_anchor_type'])]
    mobile=u.atoms[int(meta['framework_atom_count']):]; q=np.asarray(mobile.charges,float)
    if abs(q.sum())>1e-3:raise ValueError('mobile net charge exceeds 1e-3 e')
    pairs=int(meta['actual_pairs']); fsitypes=meta['guest_types']['fsi_all'];fsp=np.isin(alltypes,fsitypes)
    if len(li)!=pairs or len(fsi)!=pairs:raise ValueError('Li or FSI count mismatch')
    liq_unique=np.unique(np.round(np.asarray(li.charges),8));qfsi=float(liq[fsp].sum()/pairs)
    if len(liq_unique)!=1 or abs(float(liq_unique[0])+qfsi)>1e-3:raise ValueError('inconsistent scaled ion charges')
    n=len(u.trajectory);expected=round(count*float(meta['chunk_ps'])/float(meta['frame_dt_ps']))
    if n!=expected:raise ValueError(f'frames {n} != expected {expected}; inspect DCD boundaries')
    a=np.empty((n,len(li),3));b=np.empty((n,len(fsi),3));moment=np.empty((n,3))
    for i,_ in enumerate(u.trajectory):
        a[i]=li.positions;b[i]=fsi.positions;moment[i]=np.einsum('i,ij->j',q,mobile.positions)
    def ion_msd(coords):return np.mean([tidynamics.msd(coords[:,i,:]) for i in range(coords.shape[1])],axis=0)
    t=np.arange(n)*float(meta['frame_dt_ps']);ys=[ion_msd(a),ion_msd(b),displacement_msd(moment)]
    return meta,t,ys,{'source':'raw DCD','li_charge_e':float(liq_unique[0]),'fsi_net_charge_e':qfsi,'mobile_charge_e':float(q.sum()),'n_frames':n,'fsi_diffusion':'FSI nitrogen anchor proxy; not center of mass'}

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--case',type=Path,required=True);p.add_argument('--output-dir',type=Path,required=True);p.add_argument('--cached-curves',type=Path,help='exploratory only; three legacy CSVs')
    a=p.parse_args();case=a.case.resolve();out=a.output_dir.resolve();out.mkdir(parents=True,exist_ok=False)
    meta,t,ys,provenance=read_curves(case,a.cached_curves); total_ps=len(t)*float(meta['frame_dt_ps']); names=['Li','FSI_anchor','collective_charge'];results={}
    for name,y in zip(names,ys):results[name]=select_window(t,y,total_ps)
    li,fsi,coll=(results[x] for x in names);volume=float(meta['box_volume_A3']);temp=float(meta['temperature_K'])
    if coll['status']=='selected':coll['conductivity_mS_cm']=10*conductivity_from_slope(coll['slope'],volume,temp)
    if li['status']=='selected':li['diffusion_cm2_s']=li['slope']*1e-4/6
    if fsi['status']=='selected':fsi['diffusion_cm2_s']=fsi['slope']*1e-4/6
    if li['status']==fsi['status']=='selected' and 'li_charge_e' in provenance:
        ql=provenance['li_charge_e'];qf=provenance['fsi_net_charge_e'];n=meta['actual_pairs'];
        results['nernst_einstein_mS_cm']=10*n*E**2*(ql**2*li['slope']+qf**2*fsi['slope'])*1e-8/(6*volume*1e-30*KB*temp)
    results.update(provenance=provenance,case_name=meta['case_name'],base_name=meta['base_name'],dmc_per_lifsi=meta['dmc_per_lifsi'],replica_index=meta['replica_index'],velocity_seed=meta['velocity_seed'],temperature_K=temp,volume_A3=volume,frame_dt_ps=meta['frame_dt_ps'],algorithm='exploratory window screening v1; not a convergence certificate',thresholds={'alpha':[0.8,1.2],'linear_r2_min':.98,'half_slope_relative_difference_max':.5,'min_start_ps':50,'min_width_ps':200,'max_lag_fraction':.2,'max_lag_ps':3000})
    (out/'results.json').write_text(json.dumps(results,indent=2,allow_nan=False)+'\n')
    import pandas as pd
    pd.DataFrame({'lag_ps':t,'Li_msd_A2':ys[0],'FSI_anchor_msd_A2':ys[1],'collective_charge_msd_e2_A2':ys[2]}).to_csv(out/'curves.csv',index=False)
    import matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt
    fig,axes=plt.subplots(3,2,figsize=(11,11),constrained_layout=True)
    for i,(name,y) in enumerate(zip(names,ys)):
        res=results[name]; mask=(t>0)&(t<=min(3000,total_ps*.2))&(y>0)
        axes[i,0].loglog(t[mask],y[mask],lw=1);axes[i,1].plot(t[mask],y[mask],lw=1)
        for ax in axes[i]:
            ax.set_xlabel('Lag time (ps)');ax.set_ylabel(name+' MSD');ax.grid(alpha=.2)
            if res['status']=='selected':
                ax.axvspan(res['start_ps'],res['end_ps'],color='orange',alpha=.25)
                sub=(t>=res['start_ps'])&(t<=res['end_ps']);ax.plot(t[sub],res['slope']*t[sub]+res['intercept'],color='crimson',lw=1.5)
        axes[i,0].set_title(name+' | '+res['status']+(f" | α={res['alpha']:.3f}" if res['status']=='selected' else ''))
        if res['status']=='selected':
            val=(f"σ = {res['conductivity_mS_cm']:.3g} mS/cm" if name=='collective_charge' else f"D = {res['diffusion_cm2_s']:.3g} cm²/s")
            axes[i,1].set_title(f"{res['start_ps']:.0f}–{res['end_ps']:.0f} ps | {val}")
        else:
            axes[i,1].set_title('No eligible interval — no value reported')
    fig.savefig(out/'diffusion_windows.png',dpi=180);fig.savefig(out/'diffusion_windows.pdf');plt.close(fig)
    print(json.dumps({k:results[k]['status'] for k in names}),flush=True)
if __name__=='__main__':main()
