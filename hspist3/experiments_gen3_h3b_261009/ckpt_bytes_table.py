#!/usr/bin/env python3
"""##CHRIS 2026-10-09 (261012 sec. 4.7.25; plan-author decision 12, part 3b): the acceptance table of the checkpoint-file fix.
Checkpoint files must not contain memory addresses: two checkpoint files of the same state are byte-identical. Over the evidence
folder of stage H's runner (experiments_gen3_h_261009/stageH_evidence.py, unchanged, run with the binaries of the fix):
  1. every pair of checkpoints of the same state written by two processes: C (wrote and went on) and R1 (wrote and stopped) at each
     checkpoint of the three trajectories; and the chain: C's checkpoint at record:25000 against R2b's, written by a process that
     had RESTARTED at hold:300000 -- byte for byte, with their SHA-256;
  2. the engine test's output (H2, both builds) against stage H's output of ab80304, byte for byte (the fix changes no state and
     no file size, so the test must print the same bytes);
  3. the stage H tables of this run (stageH_tables.py, unchanged) are printed separately.
usage (from hspist3/ of the engine-gen3 worktree):
  python3 experiments_gen3_h3b_261009/ckpt_bytes_table.py --out <evidence dir> --stageH <stage H's evidence dir> > ckpt_bytes_table_output.txt
"""
import argparse, hashlib, os


def sha(p): return hashlib.sha256(open(p, "rb").read()).hexdigest()


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--out", required=True); ap.add_argument("--stageH", required=True)
    a = ap.parse_args(); O = os.path.abspath(a.out); H = os.path.abspath(a.stageH)
    print("# Decision 12, part 3b: checkpoint files of the same state, byte for byte (261012 sec. 4.7.25), printed by "
          "experiments_gen3_h3b_261009/ckpt_bytes_table.py\n")
    print("## 1. Pairs of checkpoints of the same state, written by different processes\n")
    print("| trajectory | state | file 1 (writer) | file 2 (writer) | bytes | SHA-256 of file 1 | identical |\n|---|---|---|---|---|---|---|")
    pairs = []
    for cell, cks in (("A", ("hold1", "hold300000", "record0", "record25000")), ("B", ("hold60000", "record6000")), ("C", ("hold15000", "record3000"))):
        for t in cks:
            pairs.append((cell, t, f"ck_C_{t}.bin", "C: wrote and went on", f"ck_R1_{t}.bin", "R1: wrote and stopped"))
    pairs.append(("A", "record25000", "ck_C_record25000.bin", "C: the uninterrupted run", "ck_R2b_chain.bin", "R2b: restarted at hold:300000"))
    ok = 0
    for cell, t, f1, w1, f2, w2 in pairs:
        p1, p2 = os.path.join(O, "h3", cell, f1), os.path.join(O, "h3", cell, f2)
        b1, b2 = open(p1, "rb").read(), open(p2, "rb").read()
        same = b1 == b2; ok += same
        print(f"| {cell} | {t} | {f1} ({w1}) | {f2} ({w2}) | {len(b1)} / {len(b2)} | {sha(p1)} | {'IDENTICAL' if same else '**DIFFERENT**'} |")
    print(f"\npairs byte-identical: {ok} of {len(pairs)}")
    print("\n## 2. The engine test (H2) against stage H's output (ab80304), byte for byte\n\n| build | this run | stage H | identical |\n|---|---|---|---|")
    ok2 = 0
    for b in ("default", "ld"):
        p, q = os.path.join(O, "h2", b, "checkpoint_test_output.txt"), os.path.join(H, "h2", b, "checkpoint_test_output.txt")
        same = open(p, "rb").read() == open(q, "rb").read(); ok2 += same
        print(f"| {b} | {sha(p)[:16]} | {sha(q)[:16]} | {'IDENTICAL' if same else '**DIFFERENT**'} |")
    print(f"\nACCEPTANCE (part 3b, first item): {'PASS' if ok == len(pairs) else '**NOT MET**'} -- two checkpoint files of the same state are "
          f"byte-identical ({ok} of {len(pairs)} pairs); the engine test's output unchanged in {ok2} of 2 builds")


if __name__ == "__main__":
    main()
