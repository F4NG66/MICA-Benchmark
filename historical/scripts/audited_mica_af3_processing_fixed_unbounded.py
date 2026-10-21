"""Observe original MICA conversions/Merizo calls and propagate operational failure.

No inference, segmentation, filtering or coordinate algorithms are changed.
A successful native search that emits no domain is recorded, not invented as an error.
"""
from pathlib import Path
import importlib.util,json,hashlib,time,sys,argparse,shlex

def classify(record):
 return 'completed' if record.get('success') is True and record.get('returncode')==0 else 'operational_failure'

def main():
 root=Path('/ORIGINAL_BENCHMARK_ROOT');source=root/'repos/MICA_merizo_fixed/utils/process_AF3_results.py'
 parser=argparse.ArgumentParser(add_help=False);parser.add_argument('-a','--AF3_results_path',required=True);args,_=parser.parse_known_args()
 af=Path(args.AF3_results_path);parent=af.parent;audit=parent/'domain_processing_audit';audit.mkdir(exist_ok=False)
 spec=importlib.util.spec_from_file_location('mica_original_af3_processing',source);module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
 from Bio.PDB import MMCIFParser
 cif_inputs=sorted(af.glob('*/*_model_0.cif'));reader=MMCIFParser(QUIET=True)
 expected_chains=sum(len(list(reader.get_structure('prediction',str(f))[0].get_chains())) for f in cif_inputs)
 cls=module.AF3ResultsProcessor;original_run=cls.run_merizo_command;original_convert=cls.convert_cif_to_pdb;calls=[];conversions=[]
 def observed_run(self,command,capture_output=True,timeout=None):
  start=time.time();result=original_run(self,command,capture_output=capture_output,timeout=None)
  i=len(calls);tokens=shlex.split(command);chain=tokens[tokens.index('-i')+1] if '-i' in tokens else None
  rec={'command':command,'chain_input':chain,'started_unix':start,'wall_seconds':time.time()-start,'returncode':result.get('returncode'),'status':classify(result),'timeout_seconds':None}
  (audit/f'merizo_{i:03d}.stdout.log').write_text(result.get('stdout') or '');(audit/f'merizo_{i:03d}.stderr.log').write_text(result.get('stderr') or '')
  calls.append(rec)
  with (audit/'calls.jsonl').open('a') as f:f.write(json.dumps(rec)+'\n')
  return result
 def observed_convert(self,*a,**kw):
  result=original_convert(self,*a,**kw);conversions.append(result);return result
 cls.run_merizo_command=observed_run;cls.convert_cif_to_pdb=observed_convert
 failure=None
 try:module.main()
 except BaseException as e:failure=repr(e)
 chain_files=sorted((parent/'AF3_chains').glob('*chain*.pdb'));domains=sorted((parent/'AF3_domains').glob('*.pdb'))
 reasons=[]
 if failure:reasons.append('native_processing_exception:'+failure)
 if not conversions or any(r.get('failed') for r in conversions):reasons.append('CIF_conversion_failed_or_missing')
 if not chain_files or len(chain_files)!=expected_chains or len(calls)!=len(chain_files):reasons.append('Not_all_extracted_chains_were_processed')
 if any(r['status']!='completed' for r in calls):reasons.append('Merizo_operational_failure')
 rec={'status':'failed' if reasons else 'success','reasons':reasons,'source_sha256':hashlib.sha256(source.read_bytes()).hexdigest(),'expected_predicted_chains':expected_chains,'chain_inputs':len(chain_files),'merizo_calls':len(calls),'domain_files':len(domains),'domain_file_paths':[str(f) for f in domains],'conversions':conversions,'completed_unix':time.time(),'interpretation':'Native zero-domain outcome with rc0 is allowed. All extracted chains must be processed; nonzero/timeout and conversion errors block workflow qualification. Does not prove every predicted residue receives a domain.'}
 (audit/'summary.json').write_text(json.dumps(rec,indent=2));print(json.dumps(rec),flush=True)
 if reasons:raise SystemExit(2)
if __name__=='__main__':main()
