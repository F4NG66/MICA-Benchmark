"""Prepare entity AF3 priors and copy-aware FASTA; never read reference coordinates."""
from pathlib import Path
import argparse,json,hashlib,shutil,string

def prepare(data,pdb,mode,out):
 import gemmi
 root=Path(__file__).resolve().parents[1];targets=json.loads((root/'manifest/targets.json').read_text())['targets'];t=next(t for t in targets if t['pdb_id']==pdb)
 if out.exists():raise FileExistsError(out)
 out.mkdir(parents=True);headers=[];proof=[];letters=iter(string.ascii_uppercase+string.ascii_lowercase+string.digits)
 for e in t['entities']:
  seq=e['sequence'];copies=len(e['asym_ids']);n=copies if mode=='mica' else 1;label='seq_'+hashlib.sha256(seq.encode()).hexdigest()[:16]+'_n'+str(n)
  src=data/'data/af3/selected'/(label+'.cif');st=gemmi.read_structure(str(src));ss=[''.join(gemmi.find_tabulated_residue(r.name).one_letter_code for r in c if r.find_atom('CA','*') is not None) for c in st[0]]
  if ss!=[seq]*n:raise ValueError('Wrong AF3 sequence/copies')
  name=pdb+'_'+e['entity_id'];chains=[next(letters) for _ in range(copies)];headers.append('>'+name+'|Chains '+', '.join(chains)+'|Deposited full construct\n'+seq+'\n')
  dst=out/'AF3_results'/('entity_'+e['entity_id'])/('fold_'+name+'_model_0.cif') if mode=='mica' else out/'templates'/name/'ranked_0.pdb';dst.parent.mkdir(parents=True)
  if mode=='mica':shutil.copyfile(src,dst)
  else:st.write_pdb(str(dst))
  proof.append({'entity_id':e['entity_id'],'copies':n,'source_job':label,'source_sha256':hashlib.sha256(src.read_bytes()).hexdigest(),'prepared_sha256':hashlib.sha256(dst.read_bytes()).hexdigest()})
 (out/(pdb+'.fasta')).write_text(''.join(headers));shutil.copyfile(data/'data/targets'/pdb/'full_copy_sequences.fasta',out/'full_copy_sequences.fasta');(out/'adapter_manifest.json').write_text(json.dumps({'reference_coordinates_used':False,'pdb_id':pdb,'mode':mode,'entities':proof},indent=2))
if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--data',type=Path,required=True);ap.add_argument('--pdb',required=True);ap.add_argument('--mode',choices=['mica','monomers'],required=True);ap.add_argument('--output',type=Path,required=True);a=ap.parse_args();prepare(a.data,a.pdb,a.mode,a.output)
