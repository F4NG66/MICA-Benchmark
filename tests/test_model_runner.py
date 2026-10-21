from pathlib import Path
import unittest,json,sys,tempfile,gzip,hashlib,importlib.util
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'scripts'))
from run_models import plan,METHODS
class PlanTest(unittest.TestCase):
 def setUp(self):self.cfg=json.loads((ROOT/'config.example.json').read_text());self.targets=json.loads((ROOT/'manifest/targets.json').read_text())['targets']
 def test_all_96_plans(self):
  count=0
  for t in self.targets:
   for m in METHODS:
    p=plan(m,t,self.cfg,Path('/data'),Path('/runs')/m/t['pdb_id']);self.assertTrue(p['commands']);self.assertEqual(p['requires_refinement'],m in ['mica','emodelx']);count+=1
  self.assertEqual(count,96)
 def test_distinct_cryoatom_flags(self):
  t=self.targets[0]
  for m,flag in [('cryoatom','-s'),('cryoatom2','-ps')]:self.assertIn(flag,plan(m,t,self.cfg,Path('/data'),Path('/runs')/m/t['pdb_id'])['commands'][0])
 def test_full_copy_emprot_and_native_mica(self):
  t=next(t for t in self.targets if t['pdb_id']=='9OKP')
  p=plan('emprot',t,self.cfg,Path('/data'),Path('/runs/emprot/9OKP'));cmd=p['commands'][0];self.assertTrue(cmd[cmd.index('--seq')+1].endswith('full_copy_sequences.fasta'));self.assertEqual(len(cmd[cmd.index('--chain')+1:]),2)
  p=plan('mica',t,self.cfg,Path('/data'),Path('/runs/mica/9OKP'));self.assertEqual(p['af3_mode'],'mica');self.assertEqual(len(p['commands']),3);self.assertEqual(p['commands'][1][-2:],['--per_domain_timeout','0'])
@unittest.skipUnless(importlib.util.find_spec('gemmi'),'Optional Gemmi execution tests')
class ExecutionTest(unittest.TestCase):
 def test_export_translation_and_existing_attempt(self):
  import gemmi
  from run_models import run_one
  with tempfile.TemporaryDirectory() as tmp:
   b=Path(tmp);data=b/'dataset';inp=data/'data/targets/TEST';inp.mkdir(parents=True);raw=b'canonical test map';compressed=inp/'canonical.map.gz';compressed.write_bytes(gzip.compress(raw));(data/'data_manifest.json').write_text(json.dumps([{'path':'data/targets/TEST/canonical.map.gz','sha256':hashlib.sha256(compressed.read_bytes()).hexdigest(),'uncompressed_sha256':hashlib.sha256(raw).hexdigest()}]));(inp/'map_transform.json').write_text(json.dumps({'original_map_coordinate_offset_xyz_A':[2,3,4]}));(inp/'unique_sequences.fasta').write_text('>TEST\nAAA\n')
   st=gemmi.Structure();model=gemmi.Model('1');chain=gemmi.Chain('A');res=gemmi.Residue();res.name='ALA';res.seqid=gemmi.SeqId(1,' ');a=gemmi.Atom();a.name='CA';a.element=gemmi.Element('C');a.pos=gemmi.Position(1,2,3);res.add_atom(a);chain.add_residue(res);model.add_chain(chain);st.add_model(model);source=b/'source.cif';st.make_mmcif_document().write_file(str(source))
   fake=b/'fake_modelangelo';fake.write_text('#!'+sys.executable+'\nfrom pathlib import Path\nimport sys,shutil\no=Path(sys.argv[sys.argv.index("-o")+1]);o.mkdir();shutil.copyfile('+repr(str(source))+',o/(o.name+".cif"))\n');fake.chmod(0o755)
   cfg={'methods':{'modelangelo':{'executable':str(fake),'repo':str(b)}}};t={'pdb_id':'TEST','resolution':3.0};job=b/'runs/modelangelo/TEST';p=plan('modelangelo',t,cfg,data,job);run_one(p,t,cfg,data,job);state=json.loads((job/'status.json').read_text());self.assertEqual(state['status'],'success');pos=gemmi.read_structure(str(job/'final.cif'))[0][0][0][0].pos;self.assertEqual([pos.x,pos.y,pos.z],[3,5,7])
   with self.assertRaises(FileExistsError):run_one(p,t,cfg,data,job)
 def test_refined_export_helper_imports_without_private_module(self):
  from select_refined_export import composition,select
  self.assertTrue(callable(composition));self.assertTrue(callable(select))
if __name__=='__main__':unittest.main()
