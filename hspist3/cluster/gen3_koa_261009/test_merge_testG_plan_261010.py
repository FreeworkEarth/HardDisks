#!/usr/bin/env python3
"""##CHRIS 2026-10-10 (plan-author decision 21, item 2; 261012 sec. 4.7.49): the unit test of merge_testG_koa.py's plan-then-execute
merge, on synthetic chunk trees in a new temporary folder (under TMPDIR; kept, nothing else is written). Each case runs the merge script as KOA runs it
(a separate python3 process, --no-reduce: the synthetic traces are not real). "Nothing moved" is checked byte for byte: every path, type,
size and SHA-256 under chunks/ before and after, and no merged/ folder.
  1 clean           three chunks; m_50 of both engines from chunks 1-2, x_0 of both engines from chunks 2-3: exit 0, every treatment
                    checked in merged/ (append order, per-chunk names, the runs' own folders, the build line, the marker)
  2 second start    after case 1: exit 3 (a completed merge is never repeated); nothing changes
  3 planted dup     an unknown name (01_PLOT_COMMANDS.md at the cell level) in the LAST chunk and in chunk 1, same destination: exit 2,
                    both printed, zero files moved, no merged/
  4 restart         case 3's tree with the planted files removed: exit 0 (a start after a stop in the scan phase is allowed)
  5 build differs   .build_git of the last chunk differs: exit 2, both lines printed, nothing moved
  6 trace dup       the same trace name (run index) from two chunks: exit 2, nothing moved
  7 layout          a stray file in chunks/ and an unknown cell folder: exit 2, nothing moved
  8 plan only       --plan-only on a clean tree: exit 0, nothing moved
  9 half merged     merged/ without .merge_complete: exit 3, nothing changes
usage: python3 hspist3/cluster/gen3_koa_261009/test_merge_testG_plan_261010.py     (prints a table; the last line ALL PASS or FAIL)
"""
import hashlib, os, subprocess, sys, tempfile

HERE = os.path.dirname(os.path.abspath(__file__)); MERGE = os.path.join(HERE, "merge_testG_koa.py")
BUILD = "00ALLINONE  git a643bc7  target koa"
ROWS = []


def w(p, text):
    os.makedirs(os.path.dirname(p), exist_ok=True)
    with open(p, "w") as fh: fh.write(text)


def b_run(cell, eng, r):
    """one finished speed-of-sound run as testG_worker_koa.sh leaves it in its cell folder"""
    w(os.path.join(cell, ".build_git"), BUILD + "\n")
    w(os.path.join(cell, f"wall_x_positions_L0_100_wallmassfactor_50_run{r}.csv"), f"t,x\n0,{r}\n")
    with open(os.path.join(cell, "run.log"), "a") as fh: fh.write(f"\n##RUN 2026-10-10 12:00:0{r} run {r} seed {1000 + r} node cn-x (2 s)\nrun = {r},\n")
    if eng == "gen3":
        w(os.path.join(cell, f"psi6_t_L0_100_wallmassfactor_50_run{r}.csv.gz"), f"gz {r}\n")
        with open(os.path.join(cell, ".sha256_uncompressed"), "a") as fh: fh.write(f"{r:064x}  psi6_t_L0_100_wallmassfactor_50_run{r}.csv\n")
    d = os.path.join(cell, ".done_runs", f"run{r}")
    for n in ("00_COMMAND.md", "01_PLOT_COMMANDS.md", "00_COMMAND_AND_PLOT_COMMAND.md", "speed_of_sound_psi6.csv",
              "speed_of_sound_batch_status.json", "stdout.log"):
        w(os.path.join(d, n), f"{n} of run {r}\n")


def af_run(x0, seed, j):
    """one finished A-fixed run (KEEP_EV=1); 00_COMMAND.md and run_params.json are rewritten by every run of the chunk"""
    w(os.path.join(x0, ".build_git"), BUILD + "\n")
    for n in (f"ev_{seed}.csv.gz", f"red_{seed}.csv", f"run_{seed}.log", f"summary_{seed}.csv", f"tr_{seed}.csv"): w(os.path.join(x0, n), f"{n}\n")
    w(os.path.join(x0, "00_COMMAND.md"), f"command of seed {seed} (chunk {j})\n"); w(os.path.join(x0, "run_params.json"), f'{{"seed": {seed}}}\n')
    with open(os.path.join(x0, ".sha256_uncompressed"), "a") as fh: fh.write(f"{seed:064x}  ev_{seed}.csv\n")


def tree(root):
    """chunk 1: m_50 runs 0-1 (both engines); chunk 2: m_50 run 2, x_0 seed 111; chunk 3: x_0 seed 222, a failed seed, a stale lock"""
    C = lambda j: os.path.join(root, "chunks", f"chunk_{j}", "testG")
    for eng in ("gen2", "gen3"):
        for r in (0, 1): b_run(os.path.join(C(1), eng, "epi8_H_H10_L10", "m_50"), eng, r)
        b_run(os.path.join(C(2), eng, "epi8_H_H10_L10", "m_50"), eng, 2)
        af_run(os.path.join(C(2), eng, "epi8_H_H10_L10", "x_0"), 111, 2)
        af_run(os.path.join(C(3), eng, "epi8_H_H10_L10", "x_0"), 222, 3)
    w(os.path.join(C(3), "gen2", "epi8_H_H10_L10", "x_0", ".failed_333_20261010_120000", "run_333.log"), "failed\n")
    os.makedirs(os.path.join(C(3), "gen3", "epi8_H_H10_L10", "x_0", ".guard.lock"))


def snap(root):
    out = {}
    for dp, dn, fn in os.walk(os.path.join(root, "chunks")):
        for d in dn: out[os.path.relpath(os.path.join(dp, d), root)] = "dir"
        for f in fn:
            p = os.path.join(dp, f); out[os.path.relpath(p, root)] = (os.path.getsize(p), hashlib.sha256(open(p, "rb").read()).hexdigest())
    return out


def merge(root, *extra):
    p = subprocess.run([sys.executable, MERGE, root, "--no-reduce", *extra], capture_output=True, text=True)
    return p.returncode, p.stdout + p.stderr


def case(name, expect, ok, evidence):
    ROWS.append((name, expect, "PASS" if ok else "**FAIL**", evidence))


def unchanged(root, before, rc, want_rc):
    return rc == want_rc and snap(root) == before and not os.path.exists(os.path.join(root, "merged"))


def main():
    base = tempfile.mkdtemp(prefix="merge_plan_test_")
    mk = lambda n: (lambda d: (tree(d), d)[1])(os.path.join(base, n))
    # 1 clean
    r1 = mk("clean"); before = snap(r1); rc, out = merge(r1); M = os.path.join(r1, "merged", "testG")
    g2, g3 = os.path.join(M, "gen2", "epi8_H_H10_L10"), os.path.join(M, "gen3", "epi8_H_H10_L10")
    log_ok = open(os.path.join(g2, "m_50", "run.log")).read() == "".join(
        open(os.path.join(r1, "chunks", f"chunk_{j}", "testG", "gen2", "epi8_H_H10_L10", "m_50", "run.log")).read() if os.path.exists(
            os.path.join(r1, "chunks", f"chunk_{j}", "testG", "gen2", "epi8_H_H10_L10", "m_50", "run.log")) else "" for j in (1, 2))
    # the chunks' run.log files are still there (APPEND copies; nothing is deleted), so the merged log must equal chunk 1 + chunk 2
    checks = {
        "exit 0": rc == 0,
        "marker": os.path.exists(os.path.join(r1, "merged", ".merge_complete")),
        "run.log = chunk 1 + chunk 2": log_ok and open(os.path.join(g2, "m_50", "run.log")).read().count("##RUN") == 3,
        "traces moved": all(os.path.exists(os.path.join(g, "m_50", f"wall_x_positions_L0_100_wallmassfactor_50_run{r}.csv")) for g in (g2, g3) for r in (0, 1, 2))
                        and not os.path.exists(os.path.join(r1, "chunks", "chunk_1", "testG", "gen2", "epi8_H_H10_L10", "m_50", "wall_x_positions_L0_100_wallmassfactor_50_run0.csv")),
        "psi6 moved": all(os.path.exists(os.path.join(g3, "m_50", f"psi6_t_L0_100_wallmassfactor_50_run{r}.csv.gz")) for r in (0, 1, 2)),
        "run folders per chunk": sorted(os.listdir(os.path.join(g2, "m_50", ".done_runs"))) == ["run0_chunk1", "run1_chunk1", "run2_chunk2"]
                                 and os.path.exists(os.path.join(g2, "m_50", ".done_runs", "run2_chunk2", "speed_of_sound_batch_status.json")),
        "AF commands per chunk": open(os.path.join(g2, "x_0", "00_COMMAND_chunk2.md")).read() == "command of seed 111 (chunk 2)\n"
                                 and open(os.path.join(g2, "x_0", "00_COMMAND_chunk3.md")).read() == "command of seed 222 (chunk 3)\n"
                                 and os.path.exists(os.path.join(g3, "x_0", "run_params_chunk3.json")),
        "sha files appended": open(os.path.join(g3, "x_0", ".sha256_uncompressed")).read().count("  ev_") == 2
                              and open(os.path.join(g3, "m_50", ".sha256_uncompressed")).read().count("psi6_t_") == 3,
        "failed folder kept": os.path.isdir(os.path.join(g2, "x_0", ".failed_333_20261010_120000_chunk3")),
        "lock left": os.path.isdir(os.path.join(r1, "chunks", "chunk_3", "testG", "gen3", "epi8_H_H10_L10", "x_0", ".guard.lock")),
        "one build line": all(open(os.path.join(g, c, ".build_git")).read() == BUILD + "\n" for g in (g2, g3) for c in ("m_50", "x_0")),
        "summary: two-chunk folders": "| testG/gen2/epi8_H_H10_L10/m_50 | 1, 2 |" in out and "| testG/gen3/epi8_H_H10_L10/x_0 | 2, 3 |" in out,
        "summary: treatments": "run.log: append in chunk order" in out and ".build_git: keep one after an equality check" in out
                               and ".done_runs/run<r>: rename per chunk" in out and "00_COMMAND.md: rename per chunk" in out,
    }
    bad = [k for k, v in checks.items() if not v]
    case("1 clean (3 chunks)", "exit 0; every treatment as planned", not bad, f"exit {rc}; {len(checks) - len(bad)}/{len(checks)} checks" + (f"; failed: {bad}" if bad else ""))
    # 2 second start
    s_before = snap(r1); rc, out = merge(r1)
    case("2 second start after case 1", "exit 3; nothing changes", rc == 3 and "never repeated" in out and snap(r1) == s_before, f"exit {rc}")
    # 3 planted unknown duplicate in the last chunk (and chunk 1)
    r3 = mk("planted")
    p_last = os.path.join(r3, "chunks", "chunk_3", "testG", "gen2", "epi8_H_H10_L10", "m_50", "01_PLOT_COMMANDS.md")
    b_run(os.path.dirname(p_last), "gen2", 3); w(p_last, "planted, last chunk\n")
    p_first = os.path.join(r3, "chunks", "chunk_1", "testG", "gen2", "epi8_H_H10_L10", "m_50", "01_PLOT_COMMANDS.md"); w(p_first, "planted, chunk 1\n")
    before = snap(r3); rc, out = merge(r3)
    case("3 planted unknown duplicate in the last chunk", "exit 2; both printed; zero files moved; no merged/",
         unchanged(r3, before, rc, 2) and out.count("no rule for a file named '01_PLOT_COMMANDS.md'") == 2 and "NOTHING MOVED" in out,
         f"exit {rc}; problems printed {out.count('no rule for')}; chunks/ byte-identical {snap(r3) == before}; merged/ {os.path.exists(os.path.join(r3, 'merged'))}")
    # 4 restart after the stop
    os.remove(p_last); os.remove(p_first); rc, out = merge(r3)
    case("4 restart after the scan-phase stop", "exit 0 once the cause is removed", rc == 0 and os.path.exists(os.path.join(r3, "merged", ".merge_complete")), f"exit {rc}")
    # 5 build line differs in the last chunk
    r5 = mk("build"); w(os.path.join(r5, "chunks", "chunk_3", "testG", "gen3", "epi8_H_H10_L10", "x_0", ".build_git"), "00ALLINONE  git 0000000  target koa\n")
    before = snap(r5); rc, out = merge(r5)
    case("5 .build_git differs in the last chunk", "exit 2; both lines printed; nothing moved",
         unchanged(r5, before, rc, 2) and out.count(".build_git differs") == 2, f"exit {rc}; {out.count('.build_git differs')} lines")
    # 6 the same trace name from two chunks
    r6 = mk("tracedup"); b_run(os.path.join(r6, "chunks", "chunk_3", "testG", "gen2", "epi8_H_H10_L10", "m_50"), "gen2", 1)
    before = snap(r6); rc, out = merge(r6)
    case("6 one trace name from two chunks", "exit 2; nothing moved", unchanged(r6, before, rc, 2) and "two entries for one destination" in out
         and "wall_x_positions_L0_100_wallmassfactor_50_run1.csv" in out, f"exit {rc}")
    # 7 layout
    r7 = mk("layout"); w(os.path.join(r7, "chunks", "stray.txt"), "x\n")
    os.makedirs(os.path.join(r7, "chunks", "chunk_2", "testG", "gen2", "epi8_H_H10_L10", "m_abc"))
    before = snap(r7); rc, out = merge(r7)
    case("7 stray file and unknown cell folder", "exit 2; both printed; nothing moved", unchanged(r7, before, rc, 2) and "chunks/stray.txt" in out and "m_abc" in out, f"exit {rc}")
    # 8 plan only
    r8 = mk("planonly"); before = snap(r8); rc, out = merge(r8, "--plan-only")
    case("8 --plan-only on a clean tree", "exit 0; nothing moved; no merged/", unchanged(r8, before, rc, 0) and "--plan-only" in out, f"exit {rc}")
    # 9 half merged
    r9 = mk("half"); os.makedirs(os.path.join(r9, "merged")); before = snap(r9); rc, out = merge(r9)
    case("9 merged/ without .merge_complete", "exit 3; nothing changes", rc == 3 and "inspect by hand" in out and snap(r9) == before, f"exit {rc}")
    print("| case | expected | result | evidence |\n|---|---|---|---|")
    for r in ROWS: print("| " + " | ".join(r) + " |")
    print(f"\n(the synthetic trees are kept in {base})")
    n_bad = sum(r[2] != "PASS" for r in ROWS)
    print(f"\n{'ALL PASS' if n_bad == 0 else f'FAIL ({n_bad} of {len(ROWS)})'} ({len(ROWS)} cases)")
    sys.exit(0 if n_bad == 0 else 1)


if __name__ == "__main__":
    main()
