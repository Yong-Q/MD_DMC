#!/usr/bin/env python3
"""Summarize independently run case results without treating fit windows as replicas."""
import argparse,csv,json,statistics
from collections import defaultdict
from pathlib import Path
p=argparse.ArgumentParser(description=__doc__)
p.add_argument('--manifest',type=Path,required=True);p.add_argument('--results-root',type=Path,required=True,help='directory containing INDEX/results.json, e.g. 000/results.json')
p.add_argument('--output',type=Path,required=True);args=p.parse_args()
rows=list(csv.DictReader(args.manifest.open()));groups=defaultdict(list)
for row in rows:
 idx=int(row['array_index']);result=args.results_root/f'{idx:03d}'/'results.json';v=None
 if result.is_file():
  j=json.loads(result.read_text())
  if (j.get('base_name')!=row['base_name'] or int(j.get('dmc_per_lifsi',-1))!=int(row['dmc_per_lifsi']) or int(j.get('replica_index',-1))!=int(row['replica'])):raise ValueError(f'result/manifest mismatch: {result}')
  v=j.get('collective_charge',{}).get('conductivity_mS_cm')
  # Cached legacy curves lack independent charge/trajectory verification.
  if j.get('provenance',{}).get('source')!='raw DCD':v=None
 groups[(row['base_name'],row['dmc_per_lifsi'])].append((idx,v))
args.output.parent.mkdir(parents=True,exist_ok=True)
with args.output.open('x',newline='') as f:
 fields=['base_name','dmc_per_lifsi','n_expected','n_qualified','qualified_indices','conductivity_mean_mS_cm','conductivity_sample_std_mS_cm','status']
 w=csv.DictWriter(f,fieldnames=fields);w.writeheader()
 for (base,ratio),values in sorted(groups.items(),key=lambda z:(z[0][0],int(z[0][1]))):
  good=[(i,x) for i,x in values if x is not None];xs=[x for i,x in good]
  w.writerow({'base_name':base,'dmc_per_lifsi':ratio,'n_expected':len(values),'n_qualified':len(good),'qualified_indices':';'.join(str(i) for i,_ in good),'conductivity_mean_mS_cm':statistics.mean(xs) if xs else '', 'conductivity_sample_std_mS_cm':statistics.stdev(xs) if len(xs)>1 else '', 'status':'all_replicas_qualified' if len(good)==len(values) else 'incomplete_or_rejected'})
print(args.output)
