#!/usr/bin/env python3
"""##CHRIS 2026-10-08 (261012 sec. 4.7.4, decision 2): unit tests of the loader provenance guard (edmd_acc_guard.py), one per
spelling of the flag, plus the scope rules (own directory, ancestor run records, experiments root, sibling campaigns, log
lines, command columns). usage (from hspist3/): python3 -m unittest validation/test_edmd_acc_guard.py -v"""
import os, sys, tempfile, unittest
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import edmd_acc_guard as G

REFUSED = ["--edmd-acc=1", "--edmd-acc=2", "--edmd-acc=9", "--edmd-acc=10", "--edmd-acc=yes", "--edmd-acc=YES",
           "--edmd-acc=true", "--edmd-acc=True", "--edmd-acc=on", "--edmd-acc", "--edmd-acc --particles=100",
           "--edmd-acc 1", "--edmd-acc=1,", "'--edmd-acc=1'", "--edmd-acc=maybe"]
PASSED = ["--edmd-acc=0", "--edmd-acc=false", "--edmd-acc=no", "--edmd-acc=off", "--edmd-acc=OFF", "--edmd-acc 0",
          "--edmd-acc=0,", "'--edmd-acc=0'", "--particles=100 (no flag)"]


class GuardTest(unittest.TestCase):
    def setUp(self):
        G.clear_cache()
        self.tmp = tempfile.TemporaryDirectory()
        self.root = os.path.join(self.tmp.name, "experiments_test")              # the experiments root (never scanned)
        self.cell = os.path.join(self.root, "campaign_a", "eta_0p39", "m_300")    # own directory of the data
        os.makedirs(self.cell)
        self.trace = os.path.join(self.cell, "wall_x_positions_L0_100_wallmassfactor_300_run0.csv")
        with open(self.trace, "w") as fh:
            fh.write("Time,Displacement(σ)\n0,0\n")

    def tearDown(self):
        self.tmp.cleanup()

    def write(self, rel, text):
        p = os.path.join(self.root, rel)
        os.makedirs(os.path.dirname(p), exist_ok=True)
        with open(p, "w") as fh:
            fh.write(text)
        G.clear_cache()
        return p

    def test_refused_spellings_in_own_record(self):
        for s in REFUSED:
            with self.subTest(spelling=s):
                self.write("campaign_a/eta_0p39/m_300/00_COMMAND.md", f"./00ALLINONE --mode=edmd {s}\n")
                with self.assertRaises(G.AcceleratedRunError):
                    G.guard(self.trace)

    def test_passed_spellings_in_own_record(self):
        for s in PASSED:
            with self.subTest(spelling=s):
                self.write("campaign_a/eta_0p39/m_300/00_COMMAND.md", f"./00ALLINONE --mode=edmd {s}\n")
                self.assertEqual(G.guard(self.trace), self.trace)

    def test_refused_spellings_in_command_txt_and_run_params(self):
        for name in ("command.txt", "x.command.txt", "run_params.json", "01_PLOT_COMMANDS.md", "00_COMMAND_v2.md"):
            with self.subTest(record=name):
                self.write(f"campaign_a/eta_0p39/m_300/{name}", '{"argv": "--edmd-acc=1"}\n')
                with self.assertRaises(G.AcceleratedRunError):
                    G.guard(self.trace)
                os.remove(os.path.join(self.cell, name)); G.clear_cache()

    def test_ancestor_run_record_refuses(self):
        self.write("campaign_a/00_COMMAND.md", "--edmd-acc=1\n")
        with self.assertRaises(G.AcceleratedRunError):
            G.guard(self.trace)

    def test_ancestor_log_refuses(self):
        self.write("campaign_a/driver.log", "./00ALLINONE --edmd-acc=1\n")
        with self.assertRaises(G.AcceleratedRunError):
            G.guard(self.trace)

    def test_ancestor_summary_with_command_column_refuses(self):
        self.write("campaign_a/summary.csv", "timestamp,command\nx,./00ALLINONE --edmd-acc=1\n")
        with self.assertRaises(G.AcceleratedRunError):
            G.guard(self.trace)

    def test_experiments_root_is_never_scanned(self):
        self.write("00_COMMAND.md", "--edmd-acc=1\n")
        self.write("energy_transfer_runs_2walls.csv", "timestamp,command\nx,--edmd-acc=1\n")
        self.assertEqual(G.guard(self.trace), self.trace)

    def test_sibling_campaign_does_not_decide(self):
        self.write("campaign_b/m_300/00_COMMAND.md", "--edmd-acc=1\n")
        self.assertEqual(G.guard(self.trace), self.trace)

    def test_own_log_backend_line(self):
        self.write("campaign_a/eta_0p39/m_300/run.log", "EDMD backend: default\n")
        self.assertEqual(G.guard(self.trace), self.trace)
        self.write("campaign_a/eta_0p39/m_300/run.log", "EDMD backend: default\nEDMD backend: accelerated\n")
        with self.assertRaises(G.AcceleratedRunError):
            G.guard(self.trace)

    def test_own_summary_csv_with_command_column(self):
        self.write("campaign_a/eta_0p39/m_300/summary_9700.csv", "timestamp,command\nx,./00ALLINONE --edmd-acc=1\n")
        with self.assertRaises(G.AcceleratedRunError):
            G.guard(self.trace)

    def test_csv_without_command_column_is_not_a_record(self):
        self.write("campaign_a/eta_0p39/m_300/red_9700.csv", "nu,note\n1.0,--edmd-acc=1 mentioned in a note\n")
        self.assertEqual(G.guard(self.trace), self.trace)

    def test_notes_are_not_records(self):
        self.write("campaign_a/eta_0p39/m_300/README.md", "this campaign never used --edmd-acc=1\n")
        self.assertEqual(G.guard(self.trace), self.trace)

    def test_directory_argument(self):
        self.write("campaign_a/eta_0p39/m_300/00_COMMAND.md", "--edmd-acc=on\n")
        with self.assertRaises(G.AcceleratedRunError):
            G.guard(self.cell)

    def test_returns_the_path_unchanged(self):
        self.assertIs(G.guard(self.trace), self.trace)


if __name__ == "__main__":
    unittest.main()
