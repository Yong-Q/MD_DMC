#!/usr/bin/env python3
"""Actual STACIE end-to-end unit/replica test. Execute in a Slurm allocation."""
import argparse,json,sys
from pathlib import Path
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from transport_framework import spectral_axes,read_config,summarize_values

def main():
 p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);a=p.parse_args();a.output.mkdir(parents=True,exist_ok=False)
 cfg=read_config(Path(__file__).resolve().parents[1]/'config/transport.json');meta={'box_volume_A3':100000.,'temperature_K':298.};scale=np.array([1.,2.,3.]);e=1.602176634e-19;k=1.380649e-23
 # Analytic zero-frequency conductivity for independent white charge increments.
 expected=10*e**2*1e-20*np.mean(scale**2)/(2*1e-12*100000e-30*k*298);records=[]
 for seed in range(4):
  rng=np.random.default_rng(500+seed);moment=np.cumsum(rng.normal(size=(12000,3))*scale,axis=0);out=a.output/f'seed_{seed}';out.mkdir();r=spectral_axes(moment,meta,1.,cfg,out);records.append(r);print(seed,r['value_mS_cm'],r['status'],flush=True)
 stats=summarize_values([r['value_mS_cm'] for r in records]);(a.output/'results.json').write_text(json.dumps({'analytic_sigma_mS_cm':expected,'records':records,'replica_statistics':stats},indent=2))
 assert all(r['value_mS_cm'] is not None and .65<r['value_mS_cm']/expected<1.35 for r in records)
 assert stats['ci95'][0]<expected<stats['ci95'][1]
 print('PASS analytic value and independent-replica interval')
if __name__=='__main__':main()
