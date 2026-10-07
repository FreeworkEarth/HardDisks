#!/usr/bin/env python3
"""##CHRIS 2026-10-07 (261012 sec. 4.4.10 E3, gate version 2): the task list of TEST T.
Cell epi8_H_H10_L10, method B, the campaign's protocol (geometry, stride and command from the campaign's own task list
tasks_B_epi8_H_H10_L10.txt); 100 FRESH seeds per mass, run_seed(20261007, 0, mass index, r), r = 0..99 (checked against every
campaign seed and the smoke/pilot seeds by validation/resched_testT_design_261007.py, which prints the list's SHA-256); the SAME
seed on both policies; each (mass, r) as two consecutive lines, minimal then legacy, so the two policies share nodes.
Output: cluster/resched_gate_261005/tasks_T_epi8_H_H10_L10.txt, 1800 lines in the worker's B format plus the policy field:
  B <rel>/<policy>/epi8_H_H10_L10/m_<M> <M> <r> <seed> <L0> <H> <Ns> <stride> 20261007 <policy>
with rel = experiments_speed_of_sound/EDMD/mode1_normalized_units/00_eta_sweep_ROMAN/resched_testT_261007.
usage (from hspist3/): python3 cluster/resched_gate_261005/gen_testT_261007.py
"""
import hashlib, os, sys
HERE = os.path.dirname(os.path.abspath(__file__)); CL = os.path.dirname(HERE); HS = os.path.dirname(CL)
sys.path.insert(0, os.path.join(HS, "validation")); sys.path.insert(0, HS)
import tests_20260913 as T

CID, BASE_T, NSEED = "epi8_H_H10_L10", 20261007, 100
REL = "experiments_speed_of_sound/EDMD/mode1_normalized_units/00_eta_sweep_ROMAN/resched_testT_261007"
geo = {}
for l in open(os.path.join(CL, "confinement_20261013", f"tasks_B_{CID}.txt")):
    f = l.split(); geo.setdefault(int(f[2]), (f[5], f[6], f[7], f[8]))      # L0 H Ns stride, per mass
lines = []
for mi, M in enumerate(T.A1_MASSES):
    L0, H, Ns, stride = geo[M]
    for r in range(NSEED):
        seed = T.run_seed(BASE_T, 0, mi, r)
        for pol in ("minimal", "legacy"):
            lines.append(f"B {REL}/{pol}/{CID}/m_{M} {M} {r} {seed} {L0} {H} {Ns} {stride} {BASE_T} {pol}\n")
out = os.path.join(HERE, f"tasks_T_{CID}.txt")
open(out, "w").write("".join(lines))
seedtxt = "".join(f"{M} {r} {T.run_seed(BASE_T, 0, mi, r)}\n" for mi, M in enumerate(T.A1_MASSES) for r in range(NSEED))
print(f"{out}: {len(lines)} lines; seed list SHA-256 {hashlib.sha256(seedtxt.encode()).hexdigest()}; "
      f"task list SHA-256 {hashlib.sha256(''.join(lines).encode()).hexdigest()}")
