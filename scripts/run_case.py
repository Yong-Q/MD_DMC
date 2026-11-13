#!/usr/bin/env python3
"""Run one built system in 15 restartable production chunks."""
import argparse, csv, json, os, shutil, subprocess, sys, time
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
p=argparse.ArgumentParser()
p.add_argument('--index',type=int,required=True);p.add_argument('--lammps',default=os.environ.get('LAMMPS_EXE','lmp'))
p.add_argument('--ranks',type=int,default=int(os.environ.get('SLURM_CPUS_PER_TASK','1')))
p.add_argument('--launcher',choices=['srun','mpirun','none'],default='srun')
p.add_argument('--resume',action='store_true');args=p.parse_args()
rows=list(csv.DictReader((ROOT/'manifest.csv').open())); matches=[r for r in rows if int(r['array_index'])==args.index]
if len(matches)!=1: sys.exit('invalid case index')
source=ROOT/matches[0]['case_dir']; run=ROOT/'runs'/f'{args.index:03d}'
if any(c.isspace() for c in str(run)): sys.exit('repository path must not contain whitespace')
if args.resume:
    if not run.is_dir(): sys.exit('run does not exist; omit --resume for a fresh run')
else:
    run.mkdir(parents=True,exist_ok=False)
    for d in ('inputs','analysis','logs','state','traj'): (run/d).mkdir()
    for name in ('system.data','metadata.json'): shutil.copy2(source/name,run/name)
    shutil.copy2(source/'inputs/forcefield_coeffs.in',run/'inputs/forcefield_coeffs.in')
    (run/'inputs/init_nvt.in').write_text((source/'inputs/init_nvt.in').read_text().replace('__CASE_DIR__',str(run)))
meta=json.loads((run/'metadata.json').read_text()); completed=0
for i in range(1,16):
    log=run/'logs'/('init.log' if i==1 else f'chunk_{i:04d}.log')
    files=[run/'state'/f'restart_{i:04d}.bin',run/'traj'/f'chunk_{i:04d}.dcd',log]
    if all(x.is_file() and x.stat().st_size for x in files) and f'COMPLETED_CHUNK {i}' in log.read_text(errors='replace'):
        completed=i
    else: break
if args.resume and completed==0: sys.exit('No complete first chunk. Preserve or rename this run directory before starting fresh.')
if completed==15:
    status=run/'analysis/status.json'; status.write_text(json.dumps({'status':'completed','completed_chunks':15,'max_chunks':15},indent=2)+'\n')
    sys.exit('All 15 chunks already complete.')
if args.resume and completed<15:
    # Preserve incomplete output, including autosaves, rather than mixing it with a new attempt.
    i=completed+1; candidates=[run/'traj'/f'chunk_{i:04d}.dcd',run/'state'/f'restart_{i:04d}.bin',run/'state'/f'data_{i:04d}.data',run/'logs'/f'chunk_{i:04d}.log',run/'logs'/f'thermo_chunk_{i:04d}.csv']+list((run/'state').glob('autosave.*.bin'))
    candidates=[x for x in candidates if x.exists()]
    if candidates:
        archive=run/'interrupted'/f'{time.time_ns()}';archive.mkdir(parents=True)
        for x in candidates: shutil.move(str(x),str(archive/(x.parent.name+'_'+x.name)))
launcher={'srun':['srun','--mpi='+os.environ.get('MPI_TYPE','pmi2'),'-n',str(args.ranks),'-c','1'],'mpirun':['mpirun','-np',str(args.ranks)],'none':[]}[args.launcher]
for i in range(completed+1,16):
    inp=run/'inputs'/('init_nvt.in' if i==1 else f'chunk_{i:04d}.in')
    if i>1:
        prev=i-1; f=meta['framework_atom_count']
        inp.write_text(f'''log {run}/logs/chunk_{i:04d}.log
units real
boundary p p p
atom_style full
box tilt large
neighbor 2.0 bin
neigh_modify every 1 delay 0 check yes
read_restart {run}/state/restart_{prev:04d}.bin
pair_style lj/cut/coul/long 12.0
pair_modify mix geometric
special_bonds lj/coul 0.0 0.0 0.5
bond_style harmonic
angle_style harmonic
dihedral_style opls
improper_style cvff
kspace_style pppm 1.0e-4
include {run}/inputs/forcefield_coeffs.in
group framework id 1:{f}
group mobile subtract all framework
compute tmobile mobile temp
thermo_style custom step c_tmobile pe ke etotal
thermo_modify temp tmobile flush yes
variable stepv equal step
variable tempv equal c_tmobile
variable pev equal pe
variable kev equal ke
variable etv equal etotal
fix freeze framework setforce 0.0 0.0 0.0
velocity framework set 0.0 0.0 0.0
fix recenter mobile momentum 100 linear 1 1 1
dump traj all dcd 1000 {run}/traj/chunk_{i:04d}.dcd
restart 50000 {run}/state/autosave.*.bin
fix nvt_prod mobile nvt temp 298.0 298.0 100.0
fix thermo_prod all print 1000 "${{stepv}},${{tempv}},${{pev}},${{kev}},${{etv}}" file {run}/logs/thermo_chunk_{i:04d}.csv screen no title "step,temp_K,pe_kcal_mol,ke_kcal_mol,etotal_kcal_mol"
run 1000000
unfix thermo_prod
unfix nvt_prod
undump traj
write_restart {run}/state/restart_{i:04d}.bin
write_data {run}/state/data_{i:04d}.data
print "COMPLETED_CHUNK {i}"
''')
    subprocess.run(launcher+[args.lammps,'-in',str(inp)],check=True)
    log=run/'logs'/('init.log' if i==1 else f'chunk_{i:04d}.log')
    for x in (log,run/'state'/f'restart_{i:04d}.bin',run/'traj'/f'chunk_{i:04d}.dcd'):
        if not x.is_file() or not x.stat().st_size: sys.exit(f'missing chunk output: {x}')
    if f'COMPLETED_CHUNK {i}' not in log.read_text(errors='replace'): sys.exit(f'chunk {i} lacks completion marker')
    status=run/'analysis/status.json';tmp=status.with_suffix('.tmp');tmp.write_text(json.dumps({'status':'completed' if i==15 else 'continue','completed_chunks':i,'max_chunks':15},indent=2)+'\n');tmp.replace(status)
    print(f'Finished {i}/15',flush=True)
