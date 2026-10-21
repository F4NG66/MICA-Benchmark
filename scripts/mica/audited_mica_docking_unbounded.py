"""Run pinned MICA docking with explicit per-domain outcome checks.

Search implementation, masking, ordering and merging stay in pinned MICA.
CPU count follows the pinned native MICA rule. Outcome checks and user-authorized
time/memory guards remain; installed Phenix version is explicitly recorded.
"""
from pathlib import Path
import argparse, importlib.util, json, os, signal, subprocess, time, multiprocessing

ROOT=Path(__import__('os').environ['MICA_WORK_ROOT'])

def classify(returncode, has_model, stop_reason=None):
    if stop_reason or returncode!=0:return 'operational_failure'
    return 'docked' if has_model else 'completed_no_model'

def main():
    ap=argparse.ArgumentParser()
    for a in ['density_map','fasta_file','AF3_results_path','phenix_act']:ap.add_argument('--'+a,required=True)
    ap.add_argument('--resolution',type=float,required=True)
    ap.add_argument('--contour_level',type=float,required=True)
    ap.add_argument('--per_domain_timeout',type=int,default=3600)
    args=ap.parse_args()
    requested_timeout=args.per_domain_timeout
    budget_file=ROOT/'config/mica_unbounded_budget.json'
    budget=json.loads(budget_file.read_text()) if budget_file.exists() else {}
    args.per_domain_timeout=int(budget.get('per_domain_timeout_seconds',args.per_domain_timeout))
    assert args.per_domain_timeout==0, 'User requested no time limit'
    spec=importlib.util.spec_from_file_location('native_mica_docking',ROOT/'repos/MICA_merizo_fixed/utils/dock_in_map.py')
    native=importlib.util.module_from_spec(spec);spec.loader.exec_module(native)
    native_nproc=min(max(int(multiprocessing.cpu_count()*0.6),1),48)
    base=Path(args.AF3_results_path).parent
    records=base/'domain_outcomes.jsonl'
    if records.exists():raise RuntimeError('Existing attempt; refusing overwrite')
    outcomes=[]
    class AuditedProcessor(native.PhenixDockingProcessor):
        def run_phenix_dock_in_map(self,density_map,atomic_model,resolution,save_file,fasta_file):
            name=Path(atomic_model).name
            out=Path(save_file)
            if out.exists():raise RuntimeError('Unexpected existing domain output')
            command=['bash','-c','source "$1"; shift; exec "$@"','--',self.phenix_command,
                str(Path(os.environ['MICA_PHENIX_PREFIX'])/'bin/python'),'-u',str(ROOT/'scripts/observe_phaser_docking.py'),
                str(Path(os.environ['MICA_PHENIX_PREFIX'])/'lib/python3.11/site-packages/phenix/command_line/dock_in_map.py'),
                'search_model='+atomic_model,'map_file='+density_map,'nproc='+str(native_nproc),'quick=False',
                'resolution='+str(resolution),'pdb_out='+save_file,'dock_chains_individually=True','sequence='+fasta_file]
            rec=dict(domain=name,command=command,started_unix=time.time(),status='running',requested_timeout_seconds=requested_timeout,effective_timeout_seconds=args.per_domain_timeout,budget_file=str(budget_file))
            log=Path(self.log_directory)/('audited_'+name+'.log')
            stop=None
            with log.open('x') as handle:
                proc=subprocess.Popen(command,stdout=handle,stderr=subprocess.STDOUT,start_new_session=True)
                rec['pid']=proc.pid
                def terminate_owned_child(signum,frame):
                    if proc.poll() is None:
                        os.killpg(proc.pid,signal.SIGTERM)
                        try:proc.wait(timeout=10)
                        except subprocess.TimeoutExpired:os.killpg(proc.pid,signal.SIGKILL);proc.wait()
                    raise SystemExit(128+signum)
                old_term=signal.signal(signal.SIGTERM,terminate_owned_child)
                with records.open('a') as f:f.write(json.dumps(rec)+'\n')
                try:
                    while proc.poll() is None:
                        available=int(next(x.split()[1] for x in Path('/proc/meminfo').read_text().splitlines() if x.startswith('MemAvailable:')))*1024
                        if args.per_domain_timeout and time.time()-rec['started_unix']>args.per_domain_timeout:stop='per_domain_timeout';break
                        if available<20*2**30:stop='host_memory_reserve_20GiB';break
                        time.sleep(5)
                finally:
                    if proc.poll() is None:
                        os.killpg(proc.pid,signal.SIGTERM)
                        try:proc.wait(timeout=10)
                        except subprocess.TimeoutExpired:os.killpg(proc.pid,signal.SIGKILL);proc.wait()
                signal.signal(signal.SIGTERM,old_term)
            has_model=False
            if out.exists() and out.stat().st_size:
                # Any coordinate residue is sufficient; exact identity/coverage is scored later.
                model=self.parser.get_structure(name,str(out))
                has_model=any(a.id=='CA' for a in model.get_atoms())
            rec.update(returncode=proc.returncode,stop_reason=stop,output_has_ca=has_model,
                wall_seconds=time.time()-rec['started_unix'],status=classify(proc.returncode,has_model,stop))
            outcomes.append(rec)
            with records.open('a') as f:f.write(json.dumps(rec)+'\n')
            if rec['status']=='operational_failure':raise RuntimeError('Phenix operational failure: '+name+'; '+str(rec))
            return has_model
    processor=AuditedProcessor(args.phenix_act,args.AF3_results_path,log_directory=str(base/'docking_logs'))
    expected={f.name for f in (base/'AF3_domains').glob('*.pdb')}
    result=processor.process_docking(args.density_map,args.contour_level,args.resolution,args.fasta_file,'domainwise')
    assert {r['domain'] for r in outcomes}==expected,'Not all searches were attempted'
    assert all(r['status'] in ['docked','completed_no_model'] for r in outcomes)
    if not result['successful']:raise RuntimeError('All searches completed without a docked model')
    if not processor.combine_af3_docked_results():raise RuntimeError('Native model merge failed')
    summary=dict(status='success',total_domains=len(expected),docked=sum(r['status']=='docked' for r in outcomes),
        no_model=sum(r['status']=='completed_no_model' for r in outcomes),nproc=native_nproc,
        note='Completed no-hit searches are allowed; killed/timeout/error searches are not counted as normal partial docking.',
        outcomes=outcomes)
    (base/'docking_outcome_summary.json').write_text(json.dumps(summary,indent=2))

if __name__=='__main__':main()
