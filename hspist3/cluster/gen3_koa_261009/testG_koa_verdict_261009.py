#!/usr/bin/env python3
"""##CHRIS 2026-10-09 (261012 sec. 4.7.20; stage G): PREPARED, NOT RUN, NOT REGISTERED. The verdict of TEST G ON KOA: the REGISTERED Mac
verdict script validation/testG_verdict_261009.py, unchanged and imported, with three of its constants set for the KOA data:
  BUILD  the KOA binary's version line (given here at the registration: "00ALLINONE  git <commit>  target koa");
  LOC    the fetched KOA data (experiments_gen3_gate_koa_261009/, the merged tree of merge_testG_koa.py);
  TASKS  the KOA task list (cluster/gen3_koa_261009/tasks_testG_koa_261009.txt), and the KOA extension lists.
Everything else -- the rule, the inventory, the false-alarm rates, the information rows -- is the registered script's own code.
usage (from hspist3/, on the Mac after the fetch): python3 cluster/gen3_koa_261009/testG_koa_verdict_261009.py
"""
import glob, os, sys
HERE = os.path.dirname(os.path.abspath(__file__)); HS = os.path.dirname(os.path.dirname(HERE))
sys.path.insert(0, os.path.join(HS, "validation"))
BUILD_KOA = "00ALLINONE  git <TO BE SET AT THE REGISTRATION>  target koa"
import testG_verdict_261009 as V
V.BUILD = BUILD_KOA
V.LOC = os.path.join(HS, "experiments_gen3_gate_koa_261009")
V.TASKS = os.path.join(HERE, "tasks_testG_koa_261009.txt")
V.GG = HERE                                             # the inventory reads the extension lists tasks_testG_ext_M*_261009.txt from here
if __name__ == "__main__":
    if "<TO BE SET" in BUILD_KOA: sys.exit("STOP: not registered -- BUILD_KOA is set at the registration commit")
    V.main()
