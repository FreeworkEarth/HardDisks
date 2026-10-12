#!/usr/bin/env python3
"""##CHRIS 2026-10-11 (plan-author decision 29, task P; 261012 sec. 4.7.59): cluster/gen3_testR_261011/test_stopped_chunks.py (b5378f4)
with the script, the log prefix and the sbatch named in the expected lines changed to the preview's, and case 6 also holding a Test R
chunk log. The unit test of preview_stopped_chunks.sh (read only: it
prints the resubmission line of decision 25's rule), on synthetic chunk logs written as testR_koa.sbatch writes them (a STOP line names
its node: "STOP: [<node>] ..."), in a new temporary folder (under TMPDIR; kept).
  1 complete         12 chunks, one done log each: "nothing to resubmit"
  2 stopped once     chunks 5, 6, 7 and 11 stopped in their start-up checks on two nodes: the line --array=5-7,11 and
                     --exclude=cn-05-08-07,<the two nodes> (each node once)
  3 resubmitted      chunk 5 stopped, then complete (the resubmission ran): not in the line; chunk 6 stopped twice: not in the line,
                     listed as needing the plan author
  4 started          chunk 8 started and did not finish: not in the line, listed
  5 no log           chunk 12 has no log at all: counted under "no log" with its number, not in the line
  6 other prefixes   the merge job's log and Test G's and Test R's chunk logs in the same folder are ignored
usage: python3 hspist3/cluster/gen3_preview_261011/test_preview_stopped_chunks.py      (a table; the last line ALL PASS or FAIL)
"""
import os, subprocess, sys, tempfile
HERE = os.path.dirname(os.path.abspath(__file__)); SCRIPT = os.path.join(HERE, "preview_stopped_chunks.sh"); PFX = "gen3-preview"
ROWS = []


def done(j): return f"env: gcc (GCC) 14.3.0 | ...\nchunk {j} on cn-03-38-07 (E5): lines {(j - 1) * 100 + 1}-{j * 100}; 00ALLINONE  git a643bc7  target koa\nchunk {j} done\n"
def stopped(node): return f"STOP: [{node}] the preview inputs are not the delivered ones (cluster/gen3_preview_261011/preview_inputs.sha256)\n"
def started(j): return f"env: ...\nchunk {j} on cn-03-38-07 (E5): lines {(j - 1) * 100 + 1}-{j * 100}; 00ALLINONE  git a643bc7  target koa\n"


def run(base, name, logs):
    g = os.path.join(base, name); os.makedirs(g)
    for n, t in logs.items(): open(os.path.join(g, n), "w").write(t)
    p = subprocess.run(["bash", SCRIPT, g, PFX, "12"], capture_output=True, text=True)
    return p.returncode, p.stdout + p.stderr


def case(n, e, ok, ev): ROWS.append((n, e, "PASS" if ok else "**FAIL**", ev))


def main():
    base = tempfile.mkdtemp(prefix="stopped_chunks_test_")
    full = {f"{PFX}_100_{j}.out": done(j) for j in range(1, 13)}
    rc, out = run(base, "c1", full)
    case("1 complete", "nothing to resubmit", rc == 0 and "nothing to resubmit" in out and "sbatch" not in out, f"exit {rc}")
    L = {k: v for k, v in full.items() if k not in (f"{PFX}_100_5.out", f"{PFX}_100_6.out", f"{PFX}_100_7.out", f"{PFX}_100_11.out")}
    L |= {f"{PFX}_100_5.out": stopped("cn-05-08-07"), f"{PFX}_100_6.out": stopped("cn-07-01-01"), f"{PFX}_100_7.out": stopped("cn-05-08-07"),
          f"{PFX}_100_11.out": stopped("cn-07-01-01")}
    rc, out = run(base, "c2", L)
    want = "sbatch --array=5-7,11 --exclude=cn-05-08-07,cn-07-01-01 cluster/gen3_preview_261011/preview_koa.sbatch"
    case("2 chunks 5, 6, 7 and 11 stopped once on two nodes", want, rc == 0 and want in out and out.count("stopped in the start-up checks (ran nothing)") == 4,
         f"exit {rc}; {[l for l in out.splitlines() if l.startswith('sbatch')]}")
    L3 = dict(L); L3[f"{PFX}_200_5.out"] = done(5); L3[f"{PFX}_200_6.out"] = stopped("cn-07-01-02")
    rc, out = run(base, "c3", L3)
    want3 = "sbatch --array=7,11 --exclude=cn-05-08-07,cn-07-01-01 cluster/gen3_preview_261011/preview_koa.sbatch"
    case("3 chunk 5 resubmitted and done; chunk 6 stopped twice", want3 + "; chunk 6 listed for the plan author",
         rc == 0 and want3 in out and "chunk 6: stopped in the start-up checks 2 times" in out, f"exit {rc}")
    L4 = dict(full); L4[f"{PFX}_100_8.out"] = started(8)
    rc, out = run(base, "c4", L4)
    case("4 chunk 8 started and did not finish", "not in the line; listed", rc == 0 and "chunk 8: started and did not finish" in out
         and "nothing to resubmit" in out, f"exit {rc}")
    L5 = {k: v for k, v in full.items() if k != f"{PFX}_100_12.out"}
    rc, out = run(base, "c5", L5)
    case("5 chunk 12 without any log", "counted under no log (12); nothing to resubmit", rc == 0 and "no log 1 (12)" in out and "nothing to resubmit" in out, f"exit {rc}")
    L6 = dict(full); L6["gen3-preview-merge_300.out"] = stopped("cn-x"); L6["gen3-testG_15043056_3.out"] = stopped("cn-y"); L6["gen3-testR_15050001_3.out"] = stopped("cn-z")
    rc, out = run(base, "c6", L6)
    case("6 the merge job's log, a Test G log and a Test R log in the folder", "ignored: nothing to resubmit", rc == 0 and "nothing to resubmit" in out and "cn-x" not in out and "cn-z" not in out, f"exit {rc}")
    print("| case | expected | result | evidence |\n|---|---|---|---|")
    for r in ROWS: print("| " + " | ".join(r) + " |")
    print(f"\n(the synthetic logs are kept in {base})")
    nb = sum(r[2] != "PASS" for r in ROWS)
    print(f"\n{'ALL PASS' if nb == 0 else f'FAIL ({nb} of {len(ROWS)})'} ({len(ROWS)} cases)")
    sys.exit(0 if nb == 0 else 1)


if __name__ == "__main__":
    main()
