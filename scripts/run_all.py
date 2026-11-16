"""Offline frozen-data replication. No figure redraw, manuscript edit or data retrieval."""
from pathlib import Path
import subprocess,sys,os,argparse
R=Path(__file__).resolve().parents[1];a=argparse.ArgumentParser();a.add_argument('--with-inference',action='store_true',help='Replay the original4000-replicate bootstrap and original model specifications for verification only');args=a.parse_args()
env=os.environ.copy();env['PYTHONUTF8']='1';env['OPENBLAS_NUM_THREADS']='1';env['OMP_NUM_THREADS']='1';env['MKL_NUM_THREADS']='1'
stages=['scripts/data_preparation/prepare_frozen.py','scripts/decoupling/recompute_targeted.py']
if args.with_inference:stages+=['scripts/statistical_analysis/analyze.py','scripts/statistical_analysis/additional_checks.py']
stages+=['scripts/validate_archive.py']
logs=R/'.reproduction/logs';logs.mkdir(parents=True,exist_ok=True)
for stage in stages:
    print('Running '+stage,flush=True)
    with (logs/(Path(stage).stem+'.log')).open('w',encoding='utf-8') as f:
        p=subprocess.run([sys.executable,'-I',str(R/stage)],cwd=R,env=env,stdout=f,stderr=subprocess.STDOUT)
    if p.returncode:raise SystemExit('Failed '+stage+';see .reproduction/logs/'+Path(stage).stem+'.log')
print('Replication checks completed;frozen results/figures never overwritten.')
