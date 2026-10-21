"""Select a deterministic export after verifying equivalent PDB/mmCIF outputs."""
from pathlib import Path
import gemmi
def composition(path):
    st=gemmi.read_structure(str(path))
    return sorted((c.name,r.seqid.num,str(r.seqid.icode),r.name) for c in st[0] for r in c if any(a.name=='CA' for a in r))

def atoms(p):
 st=gemmi.read_structure(str(p));rows=[]
 for c in st[0]:
  for r in c:
   for a in r:rows.append(((c.name,r.seqid.num,str(r.seqid.icode),r.name,a.name,str(a.altloc)),(a.pos.x,a.pos.y,a.pos.z)))
 return sorted(rows)
def select(folder,source):
 candidates=sorted(list(Path(folder).glob('*real_space_refined*.pdb'))+list(Path(folder).glob('*real_space_refined*.cif')))
 valid=[p for p in candidates if p.stat().st_size and composition(p)==composition(source)]
 assert valid,'Missing valid refinement export'
 assert len({p.stem for p in valid})==1,'Multiple distinct refined models'
 chosen=next((p for p in valid if p.suffix=='.cif'),valid[0]);base=atoms(chosen);errors={}
 for p in valid:
  rows=atoms(p);assert len(rows)==len(base) and all(a[0]==b[0] for a,b in zip(base,rows)),'Atom identities differ between exports'
  error=max(abs(a[1][i]-b[1][i]) for a,b in zip(base,rows) for i in range(3));assert error<=.002,'Coordinates differ beyond format rounding'
  errors[p.name]=error
 return chosen,{'selection_rule':'Prefer mmCIF after same model stem, source residue identity, atom identity and coordinate-equivalence verification; no scores inspected','candidate_exports':[str(p) for p in valid],'max_coordinate_difference_A':errors}
