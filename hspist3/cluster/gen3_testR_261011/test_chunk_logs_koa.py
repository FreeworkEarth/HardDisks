#!/usr/bin/env python3
"""##CHRIS 2026-10-11 (plan-author decision 28, task R0; 261012 sec. 4.7.57): fd8d8a2's unit test test_testG_chunk_logs_261010.py for
the generalised chunk check chunk_logs_koa.sh: the same cases with Test R's log names (gen3-testR_<id>_<j>.out, the prefix passed as the
new argument) and its STOP line, and one new case:
  8 other prefixes the merge job's own log (gen3-testR-merge_<id>.out) and Test G's chunk logs in the same folder are not counted
fd8d8a2's description follows.
##CHRIS 2026-10-10 (plan-author decision 24, point 2b; 261012 sec. 4.7.52): the unit test of testG_chunk_logs_koa.sh, the merge job's
chunk check, on synthetic log and chunk folders in a new temporary folder (under TMPDIR; kept). The logs are written as
testG_koa.sbatch writes them: the environment line, then either a STOP line (a start-up check) or "chunk <j> on <node> (<cpu>): lines
...", the workers' lines and "chunk <j> done".
  1 complete        12 chunks, one done log each: exit 0, no earlier attempt (and chunk 1 is not confused with 11 or 12)
  2 decision 24     chunks 2 and 11: a first attempt that stopped in a start-up check, then a complete one: exit 0; both earlier
                    attempts listed by name with their STOP lines
  3 stopped only    chunk 3 has only the stopped attempt (which makes no folder): exit 2, "3(no folder)"; 3b the same with a folder
                    present: exit 2, "3(no log with its done line)"
  4 two done logs   chunk 4 has two logs with its done line: exit 2, "4(2 logs with its done line)"
  5 interrupted     chunk 5: an attempt that started and did not finish, then a complete one: exit 0; listed as started, not finished
  6 no folder       chunk 6 has its done log but no folder: exit 2
  7 no log          chunk 1 has no log at all (chunks 11 and 12 have theirs): exit 2, "1(no log with its done line)"
usage: python3 hspist3/cluster/gen3_testR_261011/test_chunk_logs_koa.py      (a table; the last line ALL PASS or FAIL)
"""
import os, subprocess, sys, tempfile

HERE = os.path.dirname(os.path.abspath(__file__)); CHECK = os.path.join(HERE, "chunk_logs_koa.sh"); PFX = "gen3-testR"
ENV = "env: gcc (GCC) 14.3.0 | Python 3.12.14 at /home/charing/envs/hd/bin/python3 | git version 2.52.0\n"
ROWS = []


def done_log(j):
    return ENV + f"chunk {j} on cn-03-38-07 (Intel(R) Xeon(R) CPU E5-2670 v3 @ 2.30GHz): lines {(j - 1) * 100 + 1}-{j * 100}; 00ALLINONE  git a643bc7  target koa\nchunk {j} done\n"


def stopped_log():
    return "STOP: the Test R inputs are not the delivered ones (cluster/gen3_testR_261011/testR_inputs.sha256)\n"


def started_log(j):
    return ENV + f"chunk {j} on cn-03-38-07 (Intel(R) Xeon(R) CPU E5-2670 v3 @ 2.30GHz): lines {(j - 1) * 100 + 1}-{j * 100}; 00ALLINONE  git a643bc7  target koa\n"


def tree(base, name, nch, logs, folders):
    """logs: {file name: text}; folders: the chunk numbers that have a folder"""
    r = os.path.join(base, name); g, c = os.path.join(r, "logs"), os.path.join(r, "chunks")
    os.makedirs(g); os.makedirs(c)
    for n, t in logs.items(): open(os.path.join(g, n), "w").write(t)
    for j in folders: os.makedirs(os.path.join(c, f"chunk_{j}"))
    p = subprocess.run(["bash", CHECK, g, PFX, c, str(nch)], capture_output=True, text=True)
    return p.returncode, p.stdout + p.stderr


def case(name, expect, ok, ev): ROWS.append((name, expect, "PASS" if ok else "**FAIL**", ev))


def main():
    base = tempfile.mkdtemp(prefix="chunk_logs_test_")
    full = {f"gen3-testR_15043056_{j}.out": done_log(j) for j in range(2, 13)} | {"gen3-testR_15043031_1.out": done_log(1)}
    rc, out = tree(base, "c1", 12, full, range(1, 13))
    case("1 complete, 12 chunks", "exit 0; 12 of 12 done; no earlier attempt", rc == 0 and "chunks done: 12 of 12" in out and "earlier attempt of" not in out, f"exit {rc}")
    logs = dict(full); logs["gen3-testR_15043056_2.out"] = stopped_log(); logs["gen3-testR_15043056_11.out"] = stopped_log()
    logs["gen3-testR_15050000_2.out"] = done_log(2); logs["gen3-testR_15050000_11.out"] = done_log(11)
    rc, out = tree(base, "c2", 12, logs, range(1, 13))
    case("2 decision 24: chunks 2 and 11 stopped in a start-up check, then complete", "exit 0; 12 of 12 done; both earlier attempts listed by name with their STOP line",
         rc == 0 and "chunks done: 12 of 12" in out and "chunk logs: 14;" in out and "2 stopped in the start-up checks" in out
         and "earlier attempt of chunk 2, stopped in the start-up checks (ran nothing): gen3-testR_15043056_2.out: STOP: the Test R inputs are not the delivered ones" in out
         and "earlier attempt of chunk 11, stopped in the start-up checks (ran nothing): gen3-testR_15043056_11.out" in out, f"exit {rc}")
    logs = dict(full); logs["gen3-testR_15043056_3.out"] = stopped_log()
    rc, out = tree(base, "c3", 12, logs, [j for j in range(1, 13) if j != 3])
    case("3 chunk 3 only stopped", "exit 2; '3(no folder)' (a stopped attempt makes no folder)", rc == 2 and "NOT DONE: 3(no folder)" in out, f"exit {rc}; {out.splitlines()[-1]}")
    rc, out = tree(base, "c3b", 12, logs, range(1, 13))
    case("3b chunk 3 only stopped, a folder present", "exit 2; '3(no log with its done line)'", rc == 2 and "NOT DONE: 3(no log with its done line)" in out, f"exit {rc}; {out.splitlines()[-1]}")
    logs = dict(full); logs["gen3-testR_15050000_4.out"] = done_log(4)
    rc, out = tree(base, "c4", 12, logs, range(1, 13))
    case("4 chunk 4 with two done logs", "exit 2; '4(2 logs with its done line)'", rc == 2 and "NOT DONE: 4(2 logs with its done line)" in out, f"exit {rc}; {out.splitlines()[-1]}")
    logs = dict(full); logs["gen3-testR_15043056_5.out"] = started_log(5); logs["gen3-testR_15050000_5.out"] = done_log(5)
    rc, out = tree(base, "c5", 12, logs, range(1, 13))
    case("5 chunk 5 interrupted, then complete", "exit 0; listed as started and not finished", rc == 0 and "chunks done: 12 of 12" in out
         and "earlier attempt of chunk 5, started and not finished: gen3-testR_15043056_5.out" in out and "1 started and not finished" in out, f"exit {rc}")
    rc, out = tree(base, "c6", 12, full, [j for j in range(1, 13) if j != 6])
    case("6 chunk 6 without its folder", "exit 2; '6(no folder)'", rc == 2 and "NOT DONE: 6(no folder)" in out, f"exit {rc}; {out.splitlines()[-1]}")
    logs = {k: v for k, v in full.items() if k != "gen3-testR_15043031_1.out"}
    rc, out = tree(base, "c7", 12, logs, range(1, 13))
    case("7 chunk 1 without a log (chunks 11 and 12 have theirs)", "exit 2; '1(no log with its done line)'", rc == 2 and "NOT DONE: 1(no log with its done line)" in out,
         f"exit {rc}; {out.splitlines()[-1]}")
    logs = dict(full); logs["gen3-testR-merge_15060000.out"] = "chunk 1 done\n"; logs["gen3-testG_15043031_1.out"] = done_log(1)
    rc, out = tree(base, "c8", 12, logs, range(1, 13))
    case("8 the merge job's log and a Test G chunk log in the same folder", "exit 0; 12 of 12 done; 12 logs counted", rc == 0
         and "chunk logs: 12;" in out and "chunks done: 12 of 12" in out, f"exit {rc}")
    print("| case | expected | result | evidence |\n|---|---|---|---|")
    for r in ROWS: print("| " + " | ".join(r) + " |")
    print(f"\n(the synthetic folders are kept in {base})")
    n_bad = sum(r[2] != "PASS" for r in ROWS)
    print(f"\n{'ALL PASS' if n_bad == 0 else f'FAIL ({n_bad} of {len(ROWS)})'} ({len(ROWS)} cases)")
    sys.exit(0 if n_bad == 0 else 1)


if __name__ == "__main__":
    main()
