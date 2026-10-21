"""Run the six benchmark workflows in separate, preinstalled environments."""
from pathlib import Path
import argparse,gzip,hashlib,json,os,shutil,subprocess,sys,time
ROOT=Path(__file__).resolve().parents[1]
METHODS=['modelangelo','cryoatom','cryoatom2','emprot','emodelx','mica']
def plan(method,t,cfg,data,job):
 pid=t['pdb_id'];m=cfg['methods'][method];raw=job/pid;mp=job.parent.parent/'inputs'/pid/'input.map';inp=data/'data/targets'/pid;prep=job/'prepared';resolution=str(t['resolution']);cmds=[]
 mode='mica' if method=='mica' else 'monomers' if method in ['emprot','emodelx'] else None
 if method=='modelangelo':cmds=[[m['executable'],'build','-v',str(mp),'-pf',str(inp/'unique_sequences.fasta'),'-o',str(raw),'--device','0']]
 elif method in ['cryoatom','cryoatom2']:cmds=[[m['executable'],'build','-v',str(mp),'-s' if method=='cryoatom' else '-ps',str(inp/'unique_sequences.fasta'),'-o',str(raw),'-d','cuda:0','-k']]
 elif method=='emprot':
  chains=[str(prep/'templates'/(pid+'_'+e['entity_id'])/'ranked_0.pdb') for e in t['entities']]
  cmds=[[m['executable'],'build','--map',str(mp),'--seq',str(prep/'full_copy_sequences.fasta'),'--resolution',resolution,'--output',str(raw),'--device','0','--nt',str(cfg.get('threads',4)),'--chain']+chains]
 elif method=='emodelx':cmds=[[m['python'],str(Path(m['repo'])/'run.py'),'--protocol=temp_flex','--EM_map='+str(mp),'--fasta='+str(prep/(pid+'.fasta')),'--template_dir='+str(prep/'templates'),'--output_dir='+str(raw),'--run_pulchra','--pulchra_path='+m['pulchra'],'--mul_proc_num='+str(cfg.get('threads',4))]]
 elif method=='mica':
  work=Path(m['work_root']);repo=Path(m['repo'])
  if repo!=work/'repos/MICA_merizo_fixed':raise ValueError('MICA repo must be WORK/repos/MICA_merizo_fixed; see docs/modeling.md')
  contour=next(x['level'] for x in t['map_contours'] if x.get('primary'))
  cmds=[[m['python'],str(work/'scripts/audited_mica_af3_processing_fixed_unbounded.py'),'-f',str(prep/(pid+'.fasta')),'-a',str(prep/'AF3_results'),'-d','cpu'],[m['python'],str(work/'scripts/audited_mica_docking_unbounded.py'),'--density_map',str(mp),'--contour_level',str(contour),'--resolution',resolution,'--fasta_file',str(prep/(pid+'.fasta')),'--AF3_results_path',str(prep/'AF3_results'),'--phenix_act',cfg['phenix_activation'],'--per_domain_timeout','0'],[m['python'],str(repo/'run.py'),'-m',str(mp),'-f',str(prep/(pid+'.fasta')),'-a',str(prep/'AF3_results'),'-o',str(raw),'--run_pulchra','--pulchra_path='+m['pulchra'],'--resolution='+resolution,'--device=cuda']]
 else:raise ValueError(method)
 return {'pdb_id':pid,'method':method,'af3_mode':mode,'commands':cmds,'map':str(mp),'raw_output':str(raw),'prepared':str(prep),'requires_refinement':method in ['mica','emodelx'],'native_final':str(raw/(pid+'.cif')) if method in ['modelangelo','cryoatom','cryoatom2'] else str(raw/'output.cif') if method=='emprot' else None,'cwd':m.get('repo') or str(job)}
def sha(p):
 h=hashlib.sha256()
 with Path(p).open('rb') as f:
  for b in iter(lambda:f.read(2**20),b''):h.update(b)
 return h.hexdigest()
def run_one(p,t,cfg,data,job):
 import gemmi
 from prepare_inputs import prepare
 from restore_formal_prediction_frame import restore
 from select_refined_export import select
 if job.exists():raise FileExistsError('Use a fresh output directory: '+str(job))
 job.mkdir(parents=True);sf=job/'status.json';state={'status':'running','started_unix':time.time(),'plan':p,'stages':[]}
 def save():sf.write_text(json.dumps(state,indent=2))
 save();env={**os.environ,'CUDA_VISIBLE_DEVICES':str(cfg.get('gpu',0)),'OMP_NUM_THREADS':str(cfg.get('threads',4)),'OPENBLAS_NUM_THREADS':str(cfg.get('threads',4)),'MKL_NUM_THREADS':str(cfg.get('threads',4)),'PYTHONDONTWRITEBYTECODE':'1'}
 method=p['method'];m=cfg['methods'][method]
 if m.get('python'):env['PATH']=str(Path(m['python']).parent)+os.pathsep+env.get('PATH','')
 if method=='mica':env.update(MICA_WORK_ROOT=m['work_root'],MICA_PHENIX_PREFIX=cfg['phenix_prefix'],PYTHONPATH=m['repo']+os.pathsep+env.get('PYTHONPATH',''))
 def stage(command,cwd):
  index=len(state['stages']);rec={'command':command,'started_unix':time.time()};state['stages'].append(rec);save()
  with (job/f'stage_{index:02}.log').open('x') as log:
   r=subprocess.run(command,cwd=cwd,env=env,stdout=log,stderr=subprocess.STDOUT)
  rec.update(returncode=r.returncode,wall_seconds=time.time()-rec['started_unix']);save();r.check_returncode()
 try:
  mp=Path(p['map']);mp.parent.mkdir(parents=True,exist_ok=True);source=data/'data/targets'/t['pdb_id']/'canonical.map.gz'
  expected=json.loads((data/'data_manifest.json').read_text());row=next(x for x in expected if x['path']==str(source.relative_to(data)))
  if sha(source)!=row['sha256']:raise ValueError('Canonical map archive hash mismatch')
  if not mp.exists():
   with gzip.open(source,'rb') as src,mp.open('xb') as dst:shutil.copyfileobj(src,dst)
  if row.get('uncompressed_sha256') and sha(mp)!=row['uncompressed_sha256']:raise ValueError('Canonical map hash mismatch')
  if p['af3_mode']:prepare(data,t['pdb_id'],p['af3_mode'],Path(p['prepared']))
  if method=='mica':
   for src in (ROOT/'scripts/mica').glob('*.py'):
    dst=Path(m['work_root'])/'scripts'/src.name;dst.parent.mkdir(parents=True,exist_ok=True)
    if dst.exists() and dst.read_bytes()!=src.read_bytes():raise ValueError('Different MICA wrapper exists: '+str(dst))
    if not dst.exists():shutil.copyfile(src,dst)
  for command in p['commands']:stage(command,p['cwd'])
  if p['requires_refinement']:
   files=list(Path(p['raw_output']).glob('*_all_atom_model.pdb'))
   if len(files)!=1:raise ValueError('Expected one all-atom model')
   source=files[0];refine=job/'refinement';refine.mkdir()
   cmd=['bash','-c','source "$1"; shift; exec phenix.real_space_refine "$@"','--',cfg['phenix_activation'],str(source),p['map'],str(ROOT/'historical/config/mica_phenix_2_2_compat.eff'),'resolution='+str(t['resolution']),'dry_run=False','output.overwrite=False']
   if method=='emodelx':cmd+=['nproc='+str(cfg.get('threads',4))]
   stage(cmd,str(refine));chosen,proof=select(refine,source);state['refinement_export']=proof
  else:chosen=Path(p['native_final'])
  st=gemmi.read_structure(str(chosen))
  if not any(a.name=='CA' for model in st for c in model for r in c for a in r):raise ValueError('No final protein CA atoms')
  state['frame_restore']=restore(chosen,data/'data/targets'/t['pdb_id']/'map_transform.json',job/'final.cif');state.update(status='success',final_model=str(job/'final.cif'))
 except BaseException as ex:
  state.update(status='failed',error=repr(ex));raise
 finally:state.update(finished_unix=time.time(),wall_seconds=time.time()-state['started_unix']);save()
def main():
 ap=argparse.ArgumentParser();ap.add_argument('--config',type=Path,required=True);ap.add_argument('--data',type=Path,required=True);ap.add_argument('--output',type=Path,required=True);ap.add_argument('--methods',nargs='+',default=['all']);ap.add_argument('--targets',nargs='+',default=['all']);ap.add_argument('--execute',action='store_true');a=ap.parse_args()
 cfg=json.loads(a.config.read_text());targets=json.loads((ROOT/'manifest/targets.json').read_text())['targets'];methods=METHODS if a.methods==['all'] else a.methods
 if any(m not in METHODS for m in methods):ap.error('Unknown method')
 if a.targets!=['all']:
  wanted=set(a.targets);targets=[t for t in targets if t['pdb_id'] in wanted]
  if {t['pdb_id'] for t in targets}!=wanted:ap.error('Target outside the 16 common cases')
 data=a.data.resolve();output=a.output.resolve();jobs=[(plan(m,t,cfg,data,output/m/t['pdb_id']),t,output/m/t['pdb_id']) for t in targets for m in methods]
 if not a.execute:print(json.dumps([p for p,t,j in jobs],indent=2));return
 if any(j.exists() for p,t,j in jobs):raise FileExistsError('Some attempts already exist; use another --output')
 failures=[]
 for p,t,j in jobs:
  try:run_one(p,t,cfg,data,j)
  except Exception as ex:failures.append({'method':p['method'],'pdb_id':p['pdb_id'],'error':repr(ex)})
 print(json.dumps({'jobs':len(jobs),'failed':failures},indent=2))
 if failures:raise SystemExit(1)
if __name__=='__main__':main()
