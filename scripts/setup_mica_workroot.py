from pathlib import Path
import argparse,shutil
if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('workroot',type=Path);a=ap.parse_args();g=Path(__file__).resolve().parents[1]
 for src in (g/'scripts/mica').glob('*.py'):
  dst=a.workroot/'scripts'/src.name;dst.parent.mkdir(parents=True,exist_ok=True)
  if dst.exists():raise FileExistsError(dst)
  shutil.copyfile(src,dst)
 print('Clone patched MICA into',a.workroot/'repos/MICA_merizo_fixed','and export MICA_WORK_ROOT and MICA_PHENIX_PREFIX.')
