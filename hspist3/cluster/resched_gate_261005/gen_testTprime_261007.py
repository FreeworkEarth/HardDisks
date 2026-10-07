#!/usr/bin/env python3
"""##CHRIS 2026-10-07 (261012 sec. 4.4.13, gate v3, third plan-author decision of 2026-10-07, item 2): the task list of TEST T-PRIME.
Cell epi8_H_H10_L10, method B, the campaign's protocol (geometry, stride and command from the campaign's own task list
tasks_B_epi8_H_H10_L10.txt, as gen_testT_261007.py); masses M = 300 and M = 1500 only; 400 FRESH seeds per policy per mass:
run_seed(20261007, 1, mass index, r), r = 0..399 -- Test T's base with stream index l = 1 (Test T used l = 0), so no seed can
coincide by construction of the stream, and validation/resched_testTprime_261007.py --design checks it against every seed of
the campaign task lists, the smoke/pilot seeds, the A-fixed seeds and Test T; the SAME seed on both policies; each (mass, r) as
two consecutive lines, minimal then legacy, so the two policies share nodes.
Output: cluster/resched_gate_261005/tasks_Tprime_epi8_H_H10_L10.txt, 1600 lines in the worker's B format plus the policy field:
  B <rel>/<policy>/epi8_H_H10_L10/m_<M> <M> <r> <seed> <L0> <H> <Ns> <stride> 20261007 <policy>
with rel = experiments_speed_of_sound/EDMD/mode1_normalized_units/00_eta_sweep_ROMAN/resched_testTprime_261007.
usage (from hspist3/): python3 cluster/resched_gate_261005/gen_testTprime_261007.py
"""
import hashlib, os, sys
HERE = os.path.dirname(os.path.abspath(__file__)); CL = os.path.dirname(HERE); HS = os.path.dirname(CL)
sys.path.insert(0, os.path.join(HS, "validation")); sys.path.insert(0, HS)
import tests_20260913 as T

CID, BASE, STREAM, NSEED, MASSES = "epi8_H_H10_L10", 20261007, 1, 400, (300, 1500)
REL = "experiments_speed_of_sound/EDMD/mode1_normalized_units/00_eta_sweep_ROMAN/resched_testTprime_261007"


def seeds():
    return {(M, r): T.run_seed(BASE, STREAM, T.A1_MASSES.index(M), r) for M in MASSES for r in range(NSEED)}


def seed_text(S):
    return "".join(f"{M} {r} {S[(M, r)]}\n" for M in MASSES for r in range(NSEED))


def main():
    geo = {}
    for l in open(os.path.join(CL, "confinement_20261013", f"tasks_B_{CID}.txt")):
        f = l.split(); geo.setdefault(int(f[2]), (f[5], f[6], f[7], f[8]))      # L0 H Ns stride, per mass
    S = seeds(); lines = []
    for M in MASSES:
        L0, H, Ns, stride = geo[M]
        for r in range(NSEED):
            for pol in ("minimal", "legacy"):
                lines.append(f"B {REL}/{pol}/{CID}/m_{M} {M} {r} {S[(M, r)]} {L0} {H} {Ns} {stride} {BASE} {pol}\n")
    out = os.path.join(HERE, f"tasks_Tprime_{CID}.txt")
    if os.path.exists(out) and open(out).read() != "".join(lines): sys.exit(f"STOP: {out} exists with other content -- not overwriting")
    open(out, "w").write("".join(lines))
    print(f"{out}: {len(lines)} lines; seed list SHA-256 {hashlib.sha256(seed_text(S).encode()).hexdigest()}; "
          f"task list SHA-256 {hashlib.sha256(''.join(lines).encode()).hexdigest()}")


if __name__ == "__main__":
    main()
