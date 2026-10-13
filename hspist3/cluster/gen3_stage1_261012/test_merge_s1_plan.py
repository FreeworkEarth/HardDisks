#!/usr/bin/env python3
"""##CHRIS 2026-10-12 (plan-author decision 39 item 5 f; gate G1): the unit test of merge_s1_koa.py's changes against Test R's
merge_koa.py (whose own twelve cases, test_merge_koa_plan.py, cover the unchanged rules), on synthetic chunk trees in a new temporary
folder (under TMPDIR; kept, nothing else is written). Each case runs the merge as KOA runs it (a separate python3 process,
--no-reduce). "Nothing moved" is checked byte for byte (every path, size and SHA-256 under chunks/; no merged/).
The synthetic task list holds two groups: G1 = chunks 1-2 (cell X, two masses), G2 = chunks 3-5 (cell Y, one mass); every run leaves
the four stage 1 files (trace, psi6, virial blocks; a snapshot file for run 0).
  1 clean G2      the group G2 (chunks 3, 4, 5 -- not starting at 1): exit 0; the trace, psi6, virial and snapshot files of every run moved
                  to merged/stage1/gen3/Y/m_<M>/ (S2), run.log appended in chunk order, the build line once, .merge_complete
  2 clean G1      the group G1 in its own root: exit 0
  3 missing       G2 without chunk_4: exit 2, the missing chunk printed, nothing moved
  4 extra         G2 with a chunk_1 folder (a chunk of G1) beside it: exit 2, nothing moved
  5 virial dup    the same virial-block name (run index) from two chunks: exit 2, nothing moved
  6 unknown name  a file snapshots_x.dat at the cell level: exit 2 (no rule), nothing moved
  7 bad group     --group G9 (not in the table): exit 3 at the start, nothing changes
  8 plan only     --plan-only on a clean G2 tree: exit 0, nothing moved
usage: python3 hspist3/cluster/gen3_stage1_261012/test_merge_s1_plan.py     (prints a table; the last line ALL PASS or FAIL)
"""
import hashlib, os, subprocess, sys, tempfile
HERE = os.path.dirname(os.path.abspath(__file__)); MERGE = os.path.join(HERE, "merge_s1_koa.py")
BUILD = "00ALLINONE  git f72d764  target koa"
ROWS = []


def w(p, text):
    os.makedirs(os.path.dirname(p), exist_ok=True)
    with open(p, "w") as fh: fh.write(text)


def run_files(cell, r, snap=False):
    w(os.path.join(cell, ".build_git"), BUILD + "\n")
    w(os.path.join(cell, f"wall_x_positions_L0_100_wallmassfactor_50_run{r}.csv"), f"t,x\n0,{r}\n")
    for n in (f"psi6_t_L0_100_wallmassfactor_50_run{r}.csv.gz", f"virial_blocks_L0_100_wallmassfactor_50_run{r}.csv.gz"):
        w(os.path.join(cell, n), f"gz {n}\n")
    if snap: w(os.path.join(cell, f"snapshots_L0_100_wallmassfactor_50_run{r}.bin.gz"), "gz snap\n")
    with open(os.path.join(cell, "run.log"), "a") as fh: fh.write(f"\n##RUN 2026-10-12 12:00:0{r} run {r} seed {1000 + r} node cn-x (2 s)\nrun = {r},\n")
    with open(os.path.join(cell, ".sha256_uncompressed"), "a") as fh: fh.write(f"{r:064x}  psi6_t_L0_100_wallmassfactor_50_run{r}.csv\n")
    w(os.path.join(cell, ".done_runs", f"run{r}", "stdout.log"), f"run {r}\n")


def tree(root, chunks):
    """chunks: {j: [(cell, M, r), ...]}"""
    for j, runs in chunks.items():
        for cell, M, r in runs:
            run_files(os.path.join(root, "chunks", f"chunk_{j}", "stage1", "gen3", cell, f"m_{M}"), r, snap=(r == 0))


def lists(d):
    tl = os.path.join(d, "tasks.txt"); ct = os.path.join(d, "chunks.tsv")
    L = [f"B stage1/gen3/X/m_50 50 {r} {9000 + r} 10 10 50 26 1 60000 1 0 default B\n" for r in range(2)]
    L += [f"B stage1/gen3/X/m_100 100 {r} {9100 + r} 10 10 50 26 1 60000 1 0 default B\n" for r in range(2)]
    L += [f"B stage1/gen3/Y/m_50 50 {r} {9200 + r} 10 10 50 26 1 60000 1 0 default B\n" for r in range(3)]
    w(tl, "".join(L))
    w(ct, "chunk\tfirst\tlast\tlines\tgroup\tcell\tM\tpredicted_s_E5\tlimit_min\tclass_min\n"
          "1\t1\t2\t2\tG1\tX\t50\t1\t30\t30\n2\t3\t4\t2\tG1\tX\t100\t1\t30\t30\n"
          "3\t5\t5\t1\tG2\tY\t50\t1\t30\t30\n4\t6\t6\t1\tG2\tY\t50\t1\t30\t30\n5\t7\t7\t1\tG2\tY\t50\t1\t30\t30\n")
    return tl, ct


def snap(root):
    out = {}
    for dp, dn, fn in os.walk(os.path.join(root, "chunks")):
        for f in fn:
            p = os.path.join(dp, f); out[os.path.relpath(p, root)] = (os.path.getsize(p), hashlib.sha256(open(p, "rb").read()).hexdigest())
    return out


def merge(root, tl, ct, grp, *extra):
    r = subprocess.run([sys.executable, MERGE, root, "--tasks", tl, "--chunks", ct, "--group", grp, "--no-reduce", *extra], capture_output=True, text=True)
    return r.returncode, r.stdout + r.stderr


def case(n, what, ok, detail): ROWS.append((n, what, ok, detail))


def g2(root): tree(root, {3: [("Y", 50, 0)], 4: [("Y", 50, 1)], 5: [("Y", 50, 2)]})


def main():
    base = tempfile.mkdtemp(prefix="test_merge_s1_"); tl, ct = lists(base)
    # 1
    R = os.path.join(base, "c1"); g2(R); rc, out = merge(R, tl, ct, "G2")
    M = os.path.join(R, "merged", "stage1", "gen3", "Y", "m_50")
    moved = all(os.path.exists(os.path.join(M, f"{k}_L0_100_wallmassfactor_50_run{r}.{e}")) for r in range(3)
                for k, e in (("wall_x_positions", "csv"), ("psi6_t", "csv.gz"), ("virial_blocks", "csv.gz")))
    sn = os.path.exists(os.path.join(M, "snapshots_L0_100_wallmassfactor_50_run0.bin.gz"))
    order = [l for l in open(os.path.join(M, "run.log")).read().splitlines() if l.startswith("##RUN")] if os.path.exists(os.path.join(M, "run.log")) else []
    ok = rc == 0 and moved and sn and [o.split()[4] for o in order] == ["0", "1", "2"] and os.path.exists(os.path.join(R, "merged", ".merge_complete")) \
        and open(os.path.join(M, ".build_git")).read().strip() == BUILD
    case(1, "clean G2 (chunks 3-5)", ok, f"exit {rc}; traces, psi6, virial moved: {moved}; snapshot moved: {sn}; run.log order {[o.split()[4] for o in order]}")
    # 2
    R = os.path.join(base, "c2"); tree(R, {1: [("X", 50, 0), ("X", 50, 1)], 2: [("X", 100, 0), ("X", 100, 1)]}); rc, out = merge(R, tl, ct, "G1")
    case(2, "clean G1 (chunks 1-2)", rc == 0 and os.path.exists(os.path.join(R, "merged", "stage1", "gen3", "X", "m_100", "virial_blocks_L0_100_wallmassfactor_50_run1.csv.gz")), f"exit {rc}")
    # 3
    R = os.path.join(base, "c3"); tree(R, {3: [("Y", 50, 0)], 5: [("Y", 50, 2)]}); b = snap(R); rc, out = merge(R, tl, ct, "G2")
    case(3, "missing chunk_4", rc == 2 and "missing [4]" in out and snap(R) == b and not os.path.exists(os.path.join(R, "merged")), f"exit {rc}")
    # 4
    R = os.path.join(base, "c4"); g2(R); tree(R, {1: [("X", 50, 0)]}); b = snap(R); rc, out = merge(R, tl, ct, "G2")
    case(4, "extra chunk_1 of another group", rc == 2 and "extra [1]" in out and snap(R) == b and not os.path.exists(os.path.join(R, "merged")), f"exit {rc}")
    # 5
    R = os.path.join(base, "c5"); g2(R)
    w(os.path.join(R, "chunks", "chunk_5", "stage1", "gen3", "Y", "m_50", "virial_blocks_L0_100_wallmassfactor_50_run1.csv.gz"), "dup\n")
    b = snap(R); rc, out = merge(R, tl, ct, "G2")
    case(5, "virial file of run 1 in two chunks", rc == 2 and "two entries for one destination" in out and snap(R) == b, f"exit {rc}")
    # 6
    R = os.path.join(base, "c6"); g2(R); w(os.path.join(R, "chunks", "chunk_4", "stage1", "gen3", "Y", "m_50", "snapshots_x.dat"), "x\n")
    b = snap(R); rc, out = merge(R, tl, ct, "G2")
    case(6, "unknown name snapshots_x.dat", rc == 2 and "no rule" in out and snap(R) == b, f"exit {rc}")
    # 7
    R = os.path.join(base, "c7"); g2(R); b = snap(R); rc, out = merge(R, tl, ct, "G9")
    case(7, "group G9 not in the table", rc == 3 and "no chunk of group G9" in out and snap(R) == b, f"exit {rc}")
    # 8
    R = os.path.join(base, "c8"); g2(R); b = snap(R); rc, out = merge(R, tl, ct, "G2", "--plan-only")
    case(8, "--plan-only", rc == 0 and snap(R) == b and not os.path.exists(os.path.join(R, "merged")), f"exit {rc}")
    print(f"# merge_s1_koa.py unit test (synthetic trees under {base})\n\n| case | what | result | detail |\n|---|---|---|---|")
    for n, what, ok, det in ROWS: print(f"| {n} | {what} | {'PASS' if ok else 'FAIL'} | {det} |")
    print(f"\n{'ALL PASS' if all(r[2] for r in ROWS) else 'FAIL'} ({sum(r[2] for r in ROWS)} of {len(ROWS)} cases)")


if __name__ == "__main__":
    main()
