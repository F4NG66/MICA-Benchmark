"""Validate all companion data files by bytes and SHA-256."""
from pathlib import Path
import json,hashlib,argparse

def digest(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for chunk in iter(lambda:f.read(2**20),b''):h.update(chunk)
 return h.hexdigest()
def verify(data):
 files=json.loads((data/'data_manifest.json').read_text());seen=set()
 for row in files:
  rel=Path(row['path'])
  if rel.is_absolute() or '..' in rel.parts:raise ValueError('Unsafe manifest path')
  if str(rel) in seen:raise ValueError('Duplicate manifest path')
  seen.add(str(rel));p=data/rel
  if not p.is_file() or p.stat().st_size!=row['bytes'] or digest(p)!=row['sha256']:raise ValueError('Missing/corrupt file: '+str(p))
 models=[x for x in seen if x.startswith('data/predictions/')];refs=[x for x in seen if x.endswith('/reference.cif')]
 if len(models)!=96 or len(refs)!=16:raise ValueError('Expected 96 models and 16 references')
 return len(files)
if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--data',type=Path,required=True);a=ap.parse_args();print('Verified',verify(a.data),'files')
