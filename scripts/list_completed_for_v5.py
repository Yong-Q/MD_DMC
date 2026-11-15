#!/usr/bin/env python3
"""Print runnable, not-yet-analyzed complete case indices as an sbatch array."""
import csv,sys
from pathlib import Path
root=Path(__file__).resolve().parents[1];ready=[];excluded=[]
for r in csv.DictReader((root/'manifest.csv').open()):
 i=int(r['array_index']);case=root/'runs'/f'{i:03d}'
 if not case.is_dir():continue
 files=[case/'traj'/f'chunk_{n:04d}.dcd' for n in range(1,16)]+[case/'state'/f'restart_{n:04d}.bin' for n in range(1,16)]
 if not all(x.is_file() and x.stat().st_size>0 for x in files):continue
 if (root/'analysis_v5'/f'{i:03d}').exists():excluded.append(i);continue
 ready.append(str(i))
if excluded:print(f'Skipping {len(excluded)} existing output directories; inspect failures separately.',file=sys.stderr)
print(','.join(ready))
