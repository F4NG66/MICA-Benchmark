"""Repeat scoring of the exact supplied final models; no new inference."""
from pathlib import Path
import argparse,json
from evaluate_structures_v2 import evaluate
NAME={'MICA + AF3':'mica_af3','ModelAngelo':'modelangelo','EModelX + AF3':'emodelx_af3','CryoAtom':'cryoatom','EMProt + AF3':'emprot_af3','CryoAtom2 v2.0':'cryoatom2_v2_0'}
if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--data',type=Path,required=True);ap.add_argument('--output',type=Path,required=True);ap.add_argument('--usalign',required=True);ap.add_argument('--phenix',required=True);ap.add_argument('--pdb');a=ap.parse_args()
 root=Path(__file__).resolve().parents[1];rows=json.loads((root/'results/per_target_metrics.json').read_text());results=[]
 for r in rows:
  pid=r['pdb_id'];method=NAME[r['method']]
  if a.pdb and pid!=a.pdb:continue
  pred=list((a.data/'data/predictions'/method).glob(pid+'.*'))
  if len(pred)!=1:raise ValueError('Expected one final prediction: '+method+'/'+pid)
  ref=a.data/'data/targets'/pid/'reference.cif'
  result=evaluate(pred[0].resolve(),ref.resolve(),(a.output/method/pid).resolve(),Path(a.usalign).resolve(),Path(a.phenix).resolve())
  if result['prediction_sha256']!=r['prediction_sha256'] or result['reference_sha256']!=r['reference_sha256']:raise ValueError('Wrong structure version')
  u,p=result['usalign'],result['phenix'];v={'tm_score':u['tm_reference'],'aligned_ca_length':u['aligned_ca'],'sequence_identity':u['aligned_sequence_identity'],'ca_match':p['ca_match'],'ca_quality':p['ca_quality'],'sequence_match':p['sequence_match']}
  delta={k:v[k]-r[k] for k in v};results.append({'pdb_id':pid,'method':r['method'],'observed':v,'difference_from_archived':delta});print(pid,method,delta,flush=True)
 a.output.mkdir(parents=True,exist_ok=True);p=a.output/'comparison.json'
 if p.exists():raise FileExistsError(p)
 p.write_text(json.dumps(results,indent=2))
