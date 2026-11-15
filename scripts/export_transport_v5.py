#!/usr/bin/env python3
"""Make an auditable per-case CSV from completed v5 results and the campaign manifest."""
import argparse,csv,json
from pathlib import Path
p=argparse.ArgumentParser(description=__doc__);p.add_argument('--manifest',type=Path,required=True);p.add_argument('--results-root',type=Path,required=True);p.add_argument('--output',type=Path,required=True);args=p.parse_args()
fields=['array_index','mof','qmof_id','replica','dmc_per_lifsi','lifsi_pairs','dmc_molecules','pore_volume_A3_proxy','lifsi_concentration_pore_M_proxy','dmc_concentration_pore_M_proxy','analysis_status','li_status','li_D_cm2_s','li_lag_start_ps','li_lag_end_ps','li_alpha','li_alpha_near_boundary','li_window_shift_sensitivity','fsi_com_status','fsi_com_D_cm2_s','fsi_lag_start_ps','fsi_lag_end_ps','fsi_alpha','fsi_alpha_near_boundary','fsi_window_shift_sensitivity','collective_status','collective_sigma_mS_cm','collective_lag_start_ps','collective_lag_end_ps','collective_alpha','collective_alpha_near_boundary','collective_window_shift_sensitivity','nernst_einstein_sigma_mS_cm','temperature_first_mean_K','temperature_second_mean_K','temperature_status','note']
args.output.parent.mkdir(parents=True,exist_ok=True)
with args.output.open('x',newline='') as f:
 w=csv.DictWriter(f,fieldnames=fields,lineterminator='\n');w.writeheader()
 for r in csv.DictReader(args.manifest.open()):
  i=int(r['array_index']);p=args.results_root/f'{i:03d}'/'results.json';j=json.loads(p.read_text()) if p.exists() else {};ratio=int(r['dmc_per_lifsi']);sal=float(r['salt_concentration_per_pore_M']);
  if j and (j.get('base_name')!=r['base_name'] or int(j.get('dmc_per_lifsi',-1))!=ratio or int(j.get('replica_index',-1))!=int(r['replica']) or not str(j.get('method_version','')).startswith('v5 ')):raise ValueError(f'result/manifest mismatch: {p}')
  q=dict(array_index=i,mof=r['base_name'],qmof_id=r['qmof_id'],replica=r['replica'],dmc_per_lifsi=ratio,lifsi_pairs=r['lifsi_pairs'],dmc_molecules=r['dmc_molecules'],pore_volume_A3_proxy=r['effective_pore_A3'],lifsi_concentration_pore_M_proxy=sal,dmc_concentration_pore_M_proxy=sal*ratio,analysis_status='completed' if j else ('failed' if (args.results_root/f'{i:03d}'/'failure.json').exists() else 'not_analyzed'),note='screening candidate; not a convergence certificate' if j else 'no v5 estimate')
  for prefix,key,valuekey in [('li','Li','diffusion_cm2_s'),('fsi','FSI_COM','diffusion_cm2_s'),('collective','collective_charge','conductivity_mS_cm')]:
   x=j.get(key,{});alpha=x.get('fit',{}).get('alpha');window=x.get('window_ps',['','']);status_key={'li':'li_status','fsi':'fsi_com_status','collective':'collective_status'}[prefix];value_key={'li':'li_D_cm2_s','fsi':'fsi_com_D_cm2_s','collective':'collective_sigma_mS_cm'}[prefix]
   q.update({status_key:x.get('status','not_analyzed'),value_key:x.get(valuekey,''),prefix+'_lag_start_ps':window[0],prefix+'_lag_end_ps':window[1],prefix+'_alpha':alpha if alpha is not None else '',prefix+'_alpha_near_boundary':(alpha<.92 or alpha>1.08) if alpha is not None else '',prefix+'_window_shift_sensitivity':x.get('shift_slope_max_relative_difference','')})
  temp=j.get('temperature_diagnostic',{});q.update(nernst_einstein_sigma_mS_cm=j.get('nernst_einstein_mS_cm',''),temperature_first_mean_K=temp.get('first_mean_K',''),temperature_second_mean_K=temp.get('second_mean_K',''),temperature_status=temp.get('status','not_analyzed'));w.writerow(q)
print(args.output)
