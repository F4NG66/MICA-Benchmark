"""Metadata-only reversible full-grid map normalization; no reference input."""
from pathlib import Path
import hashlib,json
import mrcfile,numpy as np

def sha(p):
 h=hashlib.sha256()
 with open(p,'rb') as f:
  for b in iter(lambda:f.read(2**20),b''):h.update(b)
 return h.hexdigest()

def normalize(source,destination):
 source,destination=Path(source),Path(destination)
 if destination.exists():raise FileExistsError(destination)
 with mrcfile.open(source,permissive=False) as m:
  h=m.header;axes=[int(h.mapc)-1,int(h.mapr)-1,int(h.maps)-1]
  if sorted(axes)!=[0,1,2]:raise ValueError('Invalid map axis labels')
  voxel=np.array([float(h.cella.x),float(h.cella.y),float(h.cella.z)])/np.array([int(h.mx),int(h.my),int(h.mz)])
  origin=np.array([float(h.origin.x),float(h.origin.y),float(h.origin.z)])
  starts=np.zeros(3);starts[axes]=[int(h.nxstart),int(h.nystart),int(h.nzstart)]
  if np.any(origin) and np.any(starts):raise ValueError('Ambiguous simultaneous origin/starts')
  if not np.allclose([h.cellb.alpha,h.cellb.beta,h.cellb.gamma],[90,90,90]):raise ValueError('Nonorthogonal grid')
  offset=origin+starts*voxel
  labels=axes[::-1];order=[labels.index(i) for i in (2,1,0)]
  arr=np.ascontiguousarray(m.data.transpose(order))
  if not np.isfinite(arr).all():raise ValueError('Nonfinite map values')
  raw_hash=hashlib.sha256(arr.tobytes()).hexdigest()
 with mrcfile.new(destination,overwrite=False) as n:
  n.set_data(arr);n.voxel_size=tuple(voxel);n.header.mapc=1;n.header.mapr=2;n.header.maps=3;n.header.origin=(0.,0.,0.);n.header.nxstart=0;n.header.nystart=0;n.header.nzstart=0
 with mrcfile.open(destination,permissive=False) as n:
  if not np.array_equal(n.data,arr):raise ValueError('Density array changed during write')
  actual=np.array([n.voxel_size.x,n.voxel_size.y,n.voxel_size.z],dtype=float)
  if np.max(np.abs(actual-voxel))*max(arr.shape)>1e-4:raise ValueError('Voxel metadata precision loss')
 return {'original_map_sha256':sha(source),'canonical_map_sha256':sha(destination),'canonical_density_array_sha256':raw_hash,'original_axis_crs':[i+1 for i in axes],'array_transpose_order':order,'voxel_xyz_A':voxel.tolist(),'original_map_coordinate_offset_xyz_A':offset.tolist(),'shape_zyx':list(arr.shape),'value_preservation_exact':True,'rule':'Original coordinates = canonical predicted coordinates + offset. Array axis permutation and header-only origin translation; no interpolation or reference input.'}

if __name__=='__main__':
 import argparse
 p=argparse.ArgumentParser();p.add_argument('source');p.add_argument('destination');p.add_argument('record');a=p.parse_args()
 r=normalize(a.source,a.destination);Path(a.record).write_text(json.dumps(r,indent=2))
