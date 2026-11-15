#!/usr/bin/env python3
"""Optional seeded Brownian/confinement sanity check; execute in a Slurm allocation."""
import json,sys
from pathlib import Path
import numpy as np
import tidynamics
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from transport_v5 import select_stable

def main():
 records=[];n=6000;particles=50;target=1e-6;step=(2*target/1e-4)**.5
 for typ in ('brownian','confined'):
  for seed in range(8):
   rng=np.random.default_rng(9081+seed);noise=rng.normal(0,step,size=(n,particles,3))
   if typ=='brownian':x=np.cumsum(noise,axis=0)
   else:
    x=np.zeros_like(noise);a=np.exp(-1/150)
    for i in range(1,n):x[i]=a*x[i-1]+noise[i]
   y=np.mean([tidynamics.msd(x[:,i,:]) for i in range(particles)],axis=0);r=select_stable(np.arange(n,dtype=float),y)
   records.append(dict(type=typ,seed=seed,status=r['status'],D_est_cm2_s=r['fit']['slope']*1e-4/6 if r['status']=='candidate' else None))
 print(json.dumps({'target_D_cm2_s':target,'records':records},indent=2))
 brown=[r for r in records if r['type']=='brownian'];confined=[r for r in records if r['type']=='confined']
 assert all(r['status']=='candidate' and .5<r['D_est_cm2_s']/target<1.5 for r in brown)
 assert all(r['status']=='no_candidate' for r in confined)
if __name__=='__main__':main()
