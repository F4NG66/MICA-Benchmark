from pathlib import Path
import unittest,json,sys,copy
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'scripts'))
from summarize import summarize
class SummaryTest(unittest.TestCase):
 def setUp(self):self.rows=json.loads((ROOT/'results/per_target_metrics.json').read_text())
 def test_recorded_statistics(self):
  expected={s['method']:s for s in json.loads((ROOT/'results/summary.json').read_text())}
  for row in summarize(self.rows):
   self.assertEqual(row['n'],16)
   for k in ['tm_score','ca_match','ca_quality','aligned_ca_length','sequence_identity','sequence_match']:
    for stat in ['mean','sd']:self.assertAlmostEqual(row[k][stat],expected[row['method']][k][stat],places=12)
 def test_missing_does_not_become_zero(self):
  rows=copy.deepcopy(self.rows);rows[0]['tm_score']=None
  with self.assertRaises(ValueError):summarize(rows)
 def test_duplicate_target_rejected(self):
  rows=copy.deepcopy(self.rows);same=[r for r in rows if r['method']==rows[0]['method']];same[1]['pdb_id']=same[0]['pdb_id']
  with self.assertRaises(ValueError):summarize(rows)
if __name__=='__main__':unittest.main()
