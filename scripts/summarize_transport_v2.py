#!/usr/bin/env python3
"""Report means only when all expected independent replica screenings pass."""
import argparse,csv,json,sys
from collections import defaultdict
from pathlib import Path
from transport_v2 import summarize_group
p=argparse.ArgumentParser(description=__doc__);p.add_argument('--manifest',type=Path,required=True);p.add_argument('--results-root',type=Path,required=True);p.add_argument('--output',type=Path,required=True);p.add_argument('--mode',choices=('first-half','verified'),default='first-half');a=p.parse_args()
groups=defaultdict(list)
for row in csv.DictReader(a.manifest.open()):
 idx=int(row['array_index']);f=a.results_root/f'{idx:03d}'/'results.json';j=json.loads(f.read_text()) if f.is_file() else None
 if j and (j.get('base_name')!=row['base_name'] or int(j.get('dmc_per_lifsi',-1))!=int(row['dmc_per_lifsi']) or int(j.get('replica_index',-1))!=int(row['replica']) or not str(j.get('method_version','')).startswith('v2 ') or j.get('analysis_mode')!=a.mode):raise ValueError(f'wrong case or method: {f}')
 groups[(row['base_name'],row['dmc_per_lifsi'])].append((row,j))
a.output.parent.mkdir(parents=True,exist_ok=True)
with a.output.open('x',newline='') as f:
 fields=['base_name','dmc_per_lifsi','mode','observable','n_expected','n_candidate','status','mean','sample_std','units','case_statuses','holdout_statuses','kinisi_statuses']
 w=csv.DictWriter(f,fieldnames=fields);w.writeheader()
 for (base,ratio),cases in sorted(groups.items(),key=lambda x:(x[0][0],int(x[0][1]))):
  for name,key,units in [('Li','diffusion_cm2_s','cm2/s'),('FSI_COM','diffusion_cm2_s','cm2/s'),('collective_charge','conductivity_mS_cm','mS/cm')]:
   vals=[];statuses=[];holdouts=[];kinisi_statuses=[];target='first_half_provisional' if a.mode=='first-half' else 'candidate'
   for row,j in cases:
    obs=j.get(name,{}) if j else {};status=obs.get('status','missing_or_failed');temp=j.get('temperature_diagnostic',{}) if j else {};temp_ok=temp.get('first_status')=='pass' and (a.mode!='verified' or temp.get('second_status')=='pass');
    if status==target and not temp_ok:status='temperature_unchecked_or_outside'
    statuses.append(f"{row['array_index']}:{status}");holdouts.append(f"{row['array_index']}:{obs.get('holdout_status','missing')}");kinisi_statuses.append(f"{row['array_index']}:{obs.get('kinisi_crosscheck',{}).get('status','not_requested')}" if name=='Li' else '');vals.append({'replica':int(row['replica']),'seed':j.get('velocity_seed') if j else None,'value':obs.get(key) if status==target else None})
   r=summarize_group(vals,expected=4); present=[j for _,j in cases if j]; vols=[float(j['volume_A3']) for j in present];temps=[float(j['temperature_K']) for j in present];
   if r['status']=='complete' and (max(vols)-min(vols)>max(vols)*1e-5 or max(temps)-min(temps)>1e-3):r={'status':'conditions_mismatch','n_valid':r['n_valid'],'n_expected':r['n_expected']}
   w.writerow({'base_name':base,'dmc_per_lifsi':ratio,'mode':a.mode,'observable':name,'n_expected':r['n_expected'],'n_candidate':r['n_valid'],'status':r['status'],'mean':r.get('mean',''),'sample_std':r.get('sample_std',''),'units':units,'case_statuses':';'.join(statuses),'holdout_statuses':';'.join(holdouts),'kinisi_statuses':';'.join(x for x in kinisi_statuses if x)})
print(a.output)
