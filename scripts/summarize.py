"""Regenerate the six-metric table; Python standard library only."""
from pathlib import Path
import json,statistics
ROOT=Path(__file__).resolve().parents[1]
KEYS=['tm_score','ca_match','ca_quality','aligned_ca_length','sequence_identity','sequence_match']
def summarize(rows):
 methods=list(dict.fromkeys(r['method'] for r in rows)); result=[]
 for m in methods:
  rr=[r for r in rows if r['method']==m];ids=[r['pdb_id'] for r in rr]
  if len(rr)!=16 or len(set(ids))!=16:raise ValueError('Require 16 distinct common targets per method')
  if any(r.get('status')!='valid' or any(r.get(k) is None for k in KEYS) for r in rr):raise ValueError('Missing/invalid measurement')
  result.append(dict(method=m,n=16,**{k:{'mean':statistics.mean(r[k] for r in rr),'sd':statistics.stdev(r[k] for r in rr)} for k in KEYS}))
 if len({tuple(sorted(r['pdb_id'] for r in rows if r['method']==m)) for m in methods})!=1:raise ValueError('Different target sets')
 return result
if __name__=='__main__':
 rows=json.loads((ROOT/'results/per_target_metrics.json').read_text())
 print('| Method | n | TM-score | CA match (%) | CA quality (%) | Aligned CA | Sequence identity (%) | Sequence match (%) |')
 print('|---|---|---|---|---|---|---|---|')
 for s in summarize(rows):
  out=[]
  for k in KEYS:
   scale=100 if k in ['ca_match','ca_quality','sequence_identity','sequence_match'] else 1
   digits=4 if k=='tm_score' else (1 if k=='aligned_ca_length' else 2)
   out.append(f"{s[k]['mean']*scale:.{digits}f} ± {s[k]['sd']*scale:.{digits}f}")
  print('| '+' | '.join([s['method'],str(s['n'])]+out)+' |')
