#!/usr/bin/env python3
"""##CHRIS 2026-10-09 (261012 sec. 4.7.14, M3 amendment e): unit tests of the gen3 run-record guard (edmd3_health_guard.py) and
of its call from edmd_acc_guard.guard(): clean and dirty records, missing records (builds without a record, no record at all,
a per-run file without its own record), the engine named by command records (nearest wins), the experiments root, gen2 data.
usage (from hspist3/): python3 -m unittest validation/test_edmd3_health_guard.py -v"""
import os, sys, tempfile, unittest
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import edmd3_health_guard as G3
import edmd_acc_guard as G

BUILT = "[EDMD3] built #1 at t = 0 (internal units): N = 100, box 240 x 240 px, cell 32 px, dividers 1, pistons 0/0; v_ref 1\n"


def rec(rid, clean="1"):
    return (f"[EDMD3-HEALTH] {rid}: clean={clean} overlap_repair=0 wall_overdue=0 past_event=0 validator_every=60 "
            f"hash=0123456789abcdef t_end=1 engine=gen3 build=\"abc1234\" target=mac-O3 run_s=1 engine_s=0.9 driver_share=0.1\n")


SOS = "L0=10.0 M=300 run={r} seed=4242{r}"


class Gen3GuardTest(unittest.TestCase):
    def setUp(self):
        G3.clear_cache(); G.clear_cache()
        self.tmp = tempfile.TemporaryDirectory()
        self.root = os.path.join(self.tmp.name, "experiments_test")
        self.cell = os.path.join(self.root, "testG", "gen3", "m_300")
        os.makedirs(self.cell)
        self.trace = self.data("wall_x_positions_L0_10_wallmassfactor_300_run1.csv")
        self.summary = self.data("summary.csv")

    def tearDown(self):
        self.tmp.cleanup()

    def data(self, name, rel=None):
        p = os.path.join(self.root, rel, name) if rel else os.path.join(self.cell, name)
        os.makedirs(os.path.dirname(p), exist_ok=True)
        with open(p, "w") as fh:
            fh.write("Time,Displacement(σ)\n0,0\n")
        return p

    def write(self, rel, text):
        p = os.path.join(self.root, rel)
        os.makedirs(os.path.dirname(p), exist_ok=True)
        with open(p, "w") as fh:
            fh.write(text)
        G3.clear_cache(); G.clear_cache()
        return p

    def log(self, text, rel="testG/gen3/m_300/run.log"):
        return self.write(rel, text)

    def refused(self, path):
        with self.assertRaises(G3.Gen3RunError):
            G3.check(path)
        with self.assertRaises(G3.Gen3RunError):    # and through the wired guard
            G.guard(path)

    def passed(self, path):
        self.assertEqual(G3.check(path), path)
        self.assertEqual(G.guard(path), path)

    # ---- gen2 data is untouched
    def test_gen2_no_markers(self):
        self.log("seed 1\n⚠️ [EDMD-HEALTH] L0=10.0 M=300 run=1 seed=1: forced_advance=0 wall_clamp_repairs=1 overlap_repairs=0 "
                 "wall_overdue=0 past_events=0\n")
        self.write("testG/gen3/m_300/00_COMMAND.md", "./00ALLINONE --mode=edmd --engine=gen2\n")
        self.passed(self.trace); self.passed(self.summary)

    def test_no_records_at_all(self):
        self.passed(self.trace)

    # ---- clean gen3
    def test_clean_records(self):
        self.log("".join(BUILT + rec(SOS.format(r=r)) for r in range(3)))
        self.passed(self.trace); self.passed(self.summary)

    def test_clean_energy_transfer_and_replay_ids(self):
        self.log(BUILT + rec("energy-transfer seed=7") + BUILT + rec("replay m2_held"))
        self.passed(self.summary)

    # ---- rule 1: dirty
    def test_dirty_own_log(self):
        self.log(BUILT + rec(SOS.format(r=0)) + BUILT + rec(SOS.format(r=1), clean="0"))
        self.refused(self.trace); self.refused(self.summary)

    def test_dirty_other_run_refuses_whole_directory(self):
        self.log(BUILT + rec(SOS.format(r=1)) + BUILT + rec(SOS.format(r=2), clean="0"))
        self.refused(self.trace)

    def test_dirty_ancestor_log(self):
        self.log(BUILT + rec(SOS.format(r=1)))
        self.log(BUILT + rec("energy-transfer seed=3", clean="0"), rel="testG/gen3/campaign.log")
        self.refused(self.trace)

    def test_never_advanced(self):
        self.log(BUILT + "[EDMD3-HEALTH] L0=10.0 M=300 run=1 seed=42421: clean=0 no gen3 state (the run never advanced)\n")
        self.refused(self.trace)

    def test_unreadable_clean_value(self):
        for v in ("x", "1.0", "01", "true", ""):
            with self.subTest(clean=v):
                self.log(BUILT + rec(SOS.format(r=1), clean=v))
                self.refused(self.trace)

    # ---- rule 2: missing records
    def test_build_without_record(self):
        self.log(BUILT + rec(SOS.format(r=0)) + BUILT + "the run crashed here\n")
        self.refused(self.trace); self.refused(self.summary)

    def test_command_names_gen3_no_log(self):
        self.write("testG/gen3/m_300/00_COMMAND.md", "./00ALLINONE --mode=edmd --engine=gen3 --particles=100\n")
        self.refused(self.trace); self.refused(self.summary)

    def test_command_names_gen3_log_without_records(self):
        self.write("testG/gen3/m_300/command.txt", "./00ALLINONE --engine gen3\n")
        self.log("no gen3 lines in this log\n")
        self.refused(self.summary)

    # ---- rule 3: a per-run file needs its own record
    def test_per_run_file_without_its_record(self):
        self.log(BUILT + rec(SOS.format(r=0)) + BUILT + rec(SOS.format(r=2)))
        self.refused(self.trace)                     # run 1 has no record
        self.passed(self.summary)                    # the directory's records are complete and clean
        self.passed(self.data("wall_x_positions_L0_10_wallmassfactor_300_run2.csv"))
        self.passed(self.data("psi6_t_L0_10_wallmassfactor_300_run0.csv"))

    def test_per_run_file_other_mass_or_L0(self):
        self.log(BUILT + rec(SOS.format(r=1)))
        self.refused(self.data("wall_x_positions_L0_10_wallmassfactor_1500_run1.csv"))
        self.refused(self.data("wall_x_positions_L0_20_wallmassfactor_300_run1.csv"))
        self.passed(self.trace)

    # ---- the engine named by command records: the nearest directory decides
    def test_nearest_engine_record_decides(self):
        self.write("testG/00_COMMAND.md", "tasks: --engine=gen2 ... --engine=gen3\n")
        self.write("testG/gen3/m_300/00_COMMAND.md", "./00ALLINONE --engine=gen2\n")
        self.passed(self.trace)                      # own record says gen2, no gen3 log
        self.write("testG/gen3/m_300/00_COMMAND.md", "./00ALLINONE --engine=gen3\n")
        self.refused(self.trace)                     # own record says gen3, no run record

    def test_parent_engine_record_applies(self):
        self.write("testG/00_COMMAND.md", "./00ALLINONE --engine=gen3\n")
        self.refused(self.trace)
        self.log(BUILT + rec(SOS.format(r=1)))
        self.passed(self.trace)

    # ---- the experiments root is never scanned
    def test_experiments_root_not_scanned(self):
        self.write("run.log", BUILT + rec("energy-transfer seed=1", clean="0"))
        self.write("00_COMMAND.md", "--engine=gen3\n")
        self.passed(self.trace)

    def test_sibling_campaign_does_not_decide(self):
        self.write("testG/gen2/m_300/run.log", BUILT + rec(SOS.format(r=1), clean="0"))
        self.passed(self.trace)


if __name__ == "__main__":
    unittest.main()
