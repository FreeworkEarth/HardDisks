"""Permanent checks for the analysis screens (no simulator or provider calls)."""
import math
import csv
import json
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from analyze_sound_structure import METRICS, grouped, screen, parse_metrics, load_sources


def row(hold, end, eta=.7, neighbors=5):
    r = {k: 0.5 for k in METRICS}
    r.update(psi6_global_hold=hold, psi6_global_end=end,
             neighbors_hold=neighbors, neighbors_end=neighbors,
             eta=eta, family="test", N_total=100, M=50)
    return r


class Screens(unittest.TestCase):
    def test_documented_no_neighbor_case(self):
        r=row(.2, float('nan'))
        r['psi6_local_end']=float('nan')
        r['neighbors_end']=0
        parsed, reason=parse_metrics(r)
        self.assertIsNone(parsed['psi6_global_end'])
        self.assertEqual(reason, 'undefined_psi6_end_no_neighbors')

    def test_unexpected_nonfinite_rejected(self):
        with self.assertRaises(ValueError):
            parse_metrics(row(.2, float('nan')))
        with self.assertRaises(ValueError):
            parse_metrics(row(.2, float('inf')))

    def test_out_of_range_rejected(self):
        with self.assertRaises(ValueError):
            parse_metrics(row(.2, 1.1))

    def test_paired_sem(self):
        result = screen([row(.1,.2),row(.8,.9)])
        self.assertAlmostEqual(result['delta_mean'], .1)
        self.assertAlmostEqual(result['delta_sem'], 0)

    def test_large_shift(self):
        self.assertTrue(screen([row(.2,.7)]*10)['endpoint_shift_flag'])

    def test_cancelling_changes_remain_visible(self):
        result=screen([row(.1,.8),row(.8,.1)])
        self.assertAlmostEqual(result['delta_mean'],0)
        self.assertEqual(result['fraction_abs_delta_gt_0p1'],1)

    def test_dilute_neighbors_not_square_flag(self):
        self.assertFalse(screen([row(.1,.1,eta=.1,neighbors=2)])['dense_low_neighbors_flag'])

    def test_dense_neighbor_screen(self):
        self.assertTrue(screen([row(.1,.1,neighbors=4)])['dense_low_neighbors_flag'])

    def test_dual_tails_need_sample_size(self):
        self.assertFalse(screen([row(.2,.1),row(.2,.9)])['dual_tail_screen'])
        self.assertTrue(screen([row(.2,.1)]*5+[row(.2,.9)]*5)['dual_tail_screen'])

    def test_single_sample_error_missing(self):
        self.assertIsNone(screen([row(.2,.3)])['delta_sem'])

    def test_mass_groups_stay_separate(self):
        a,b=row(.2,.3),row(.2,.8)
        b['M']=100
        self.assertEqual(len(grouped([a,b], ('family','N_total','eta','M'))),2)


class SourceLoading(unittest.TestCase):
    def setUp(self):
        self.tmp=tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root=Path(self.tmp.name)
        self.parent=self.root/'campaign'/'eta_0p700'
        self.parent.mkdir(parents=True)
        self.scope=patch('analyze_sound_structure.CAMPAIGNS', {'campaign':'Test'})
        self.scope.start()
        self.addCleanup(self.scope.stop)
        self.raw=row(.2,.3)
        self.raw.update(L0=10,wall_mass_factor=50,repeat=0,seed=123)
        self.write_rows([self.raw])
        (self.parent/'00_COMMAND.md').write_text('--particles=100 --edmd-acc=0')
        (self.parent/'run.log').write_text('finished')
        self.status={'complete':True,'requested_runs':1,'valid_runs':1,'invalid_runs':0}
        self.write_status()

    def write_rows(self, rows):
        with (self.parent/'speed_of_sound_psi6.csv').open('w',newline='') as f:
            w=csv.DictWriter(f,fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)

    def write_status(self):
        (self.parent/'speed_of_sound_batch_status.json').write_text(json.dumps(self.status))

    def test_leaf_provenance_and_merged_copy(self):
        merged=self.parent/'merged'; merged.mkdir()
        (merged/'speed_of_sound_psi6.csv').write_text('not a source')
        records,ledgers,manifest=load_sources(self.root)
        self.assertEqual(len(records),1)
        self.assertEqual(records[0]['source_line'],2)
        self.assertTrue(records[0]['eligible_structure'])
        self.assertEqual(len(manifest),4)

    def test_health_overrides_valid_ledger(self):
        (self.parent/'run.log').write_text('[EDMD-HEALTH] L0=10 M=50 run=0 seed=123: forced_advance=0 wall_clamp_repairs=1')
        records,_,_=load_sources(self.root)
        self.assertFalse(records[0]['eligible_structure'])
        self.assertIn('wall_clamp_repairs=1',records[0]['health_reason'])

    def test_unmatched_health_fails(self):
        (self.parent/'run.log').write_text('[EDMD-HEALTH] L0=10 M=50 run=0 seed=456: wall_overdue=1')
        with self.assertRaisesRegex(ValueError,'Unmatched health'):
            load_sources(self.root)

    def test_incomplete_ledger_not_eligible(self):
        self.status['complete']=False; self.write_status()
        self.assertFalse(load_sources(self.root)[0][0]['eligible_structure'])

    def test_count_mismatch_fails(self):
        self.status['valid_runs']=2; self.write_status()
        with self.assertRaisesRegex(ValueError,'count mismatch'):
            load_sources(self.root)

    def test_duplicate_trajectory_fails(self):
        self.write_rows([self.raw,self.raw])
        self.status.update(valid_runs=2,requested_runs=2); self.write_status()
        with self.assertRaisesRegex(ValueError,'Duplicate'):
            load_sources(self.root)

    def test_source_change_fails(self):
        with patch('analyze_sound_structure.digest',return_value='changed'):
            with self.assertRaisesRegex(ValueError,'Source changed'):
                load_sources(self.root)

    def test_accelerated_core_fails(self):
        (self.parent/'00_COMMAND.md').write_text('--particles=100 --edmd-acc=1')
        with self.assertRaisesRegex(ValueError,'Accelerated'):
            load_sources(self.root)


if __name__ == '__main__':
    unittest.main()
