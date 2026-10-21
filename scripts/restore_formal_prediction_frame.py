"""Restore only the translation encoded in the original map header."""
from pathlib import Path
import json,hashlib,math
import gemmi

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def restore(source,transform,output):
 source,transform,output=map(Path,[source,transform,output]);d=json.loads(transform.read_text());shift=d['original_map_coordinate_offset_xyz_A'];st=gemmi.read_structure(str(source))
 if output.exists():raise FileExistsError(output)
 before=[(c.name,r.name,r.seqid.num,str(r.seqid.icode),a.name,a.pos.x,a.pos.y,a.pos.z) for c in st[0] for r in c for a in r]
 for model in st:
  for c in model:
   for r in c:
    for a in r:a.pos=gemmi.Position(a.pos.x+shift[0],a.pos.y+shift[1],a.pos.z+shift[2])
 st.make_mmcif_document().write_file(str(output));back=gemmi.read_structure(str(output));after=[(c.name,r.name,r.seqid.num,str(r.seqid.icode),a.name,a.pos.x,a.pos.y,a.pos.z) for c in back[0] for r in c for a in r]
 assert len(before)==len(after) and all(a[:5]==b[:5] for a,b in zip(before,after)),'Coordinate export changed atom identities'
 err=max(abs(b[i]-a[i]-shift[i-5]) for a,b in zip(before,after) for i in [5,6,7]);assert err<1e-3,'Unexpected shift/precision'
 return {'source_sha256':sha(source),'restored_sha256':sha(output),'transform_sha256':sha(transform),'offset_xyz_A':shift,'max_coordinate_roundtrip_error_A':err,'atom_count':len(before),'identity_order_preserved':True,'reference_used_for_transformation':False}
