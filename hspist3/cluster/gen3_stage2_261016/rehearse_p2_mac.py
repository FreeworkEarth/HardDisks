#!/usr/bin/env python3
"""##CHRIS 2026-10-15 (plan-author decision 57 item 2; notes 261012 sec. 4.7.135): THE MAC REHEARSAL OF THE P2 PACKAGE: the runsheet's KOA
blocks A to E run as Chris types them (the runsheet's own lines, read from runsheet_p2_261016.txt), with the stage 1 rehearsal's shims and
these REHEARSAL SUBSTITUTIONS, each listed in the output:
  STAND-IN of the package commit: a shared scratch clone of <pkg> with (i) the chunk table cut to the four rehearsal chunks -- the task
    lines 221, 271, 321, 371: eta 0.78, alpha 5, s = 0, one per (G, pin) (decision 57 item 2) -- (ii) the Lustre root in p2_koa.sbatch,
    p2_merge.sbatch and fetch_p2_koa.sh made local and the fetch's DTN source made local (rehearse_k2_mac.rewrite_lustre), (iii)
    p2_inputs.sha256 made again (make_p2_manifests.py), committed there as a local commit; the 400-line task list is unchanged;
  BLOCK A: GitHub is not used: the clone line becomes a shared clone of this repository (--no-checkout, core.abbrev 7) and its fetch of
    8aaf5fef is not needed (the objects are local); the kissfft clone comes from the main checkout's hspist3/kissfft (febd4ca); the
    fetch of the package commit fetches the stand-in from its clone and checks it out; "exit" is left out;
  BUILD SHIMS (block A): SLURM_JOB_ID set; module (nothing); gcc (--version prints 14.3.0; else cc without -march=x86-64-v2 -mtune=generic);
    make (the real make with the Mac CFLAGS_KOA and BUILD_CF, as rehearse_s1_build_mac.py); ldd (otool -L);
  RUN SHIMS (rehearse_k2_mac.shims): sbatch runs the job at once (array tasks one after the other), squeue empty, srun, lscpu, hostname;
  BLOCK B: --array=1-400%64 becomes --array=221,271,321,371, and the guard's package commit the stand-in's (the build record names it);
    the block is typed twice (the second paste must refuse);
  BLOCKS D, D2: the merge job, its log, then the merge submitted again (it must refuse: exit 3);
  BLOCK E: the fetch runs from the stand-in clone (its DST is the stand-in's own hspist3/experiments_gen3_stage2_koa_261016/P2).
Then p2_tables_mac.py on the fetched folder. The record length is the full one (200 predicted periods), not 10: the stage 1 reducers
(unchanged, decision 57) need the 200-period prefix (TD) and a 20-period ACF window. Writes into --root (scratch: clones, data) and --out
(a new folder: every block's transcript, the fetched tables, the P2 tables).
usage (from hspist3/): python3 cluster/gen3_stage2_261016/rehearse_p2_mac.py --pkg <40 hex> --root <scratch folder> --out <new folder>
"""
import argparse, os, re, shlex, shutil, subprocess, sys, time
HERE = os.path.dirname(os.path.abspath(__file__)); HS = os.path.dirname(os.path.dirname(HERE)); REPO = os.path.dirname(HS)
sys.path.insert(0, os.path.join(HS, "cluster", "gen3_koa_261009"))
import rehearse_k2_mac as K2R
P = "cluster/gen3_stage2_261016"; BUILD = "8aaf5fef27c92391fc4fadd5133324740c29ef06"; MINI = (221, 271, 321, 371)
STEPS = []


def step(lab, ok, det): STEPS.append((lab, bool(ok), det)); print(f"[{'PASS' if ok else 'FAIL'}] {lab}: {det}", flush=True)


def kissfft(dest):
    """kissfft at febd4ca, cloned from the main checkout's hspist3/kissfft (the consolidation moved the gen3 worktree's clone there; its
    local change to kiss_fft_log.h is uncommitted, so the clone holds febd4ca's files)"""
    K2R.sh(["git", "clone", "-q", os.path.join(HS, "kissfft"), os.path.join(dest, "hspist3", "kissfft")])
    K2R.sh(["git", "-C", os.path.join(dest, "hspist3", "kissfft"), "-c", "advice.detachedHead=false", "checkout", "-q", K2R.KISS])


def build_shims(R):
    d = os.path.join(R, "build_shims"); os.makedirs(d)
    real_make = shutil.which("make"); mac = "-O2 -ffp-contract=off"
    for n, body in (("module", "exit 0\n"),
                    ("gcc", 'if [ "${1:-}" = --version ]; then echo "gcc (GCC) 14.3.0"; exit 0; fi\n'
                            'a=(); for x in "$@"; do case "$x" in -march=x86-64-v2|-mtune=generic) ;; *) a+=("$x");; esac; done\nexec cc "${a[@]}"\n'),
                    ("make", f'exec {real_make} "$@" \'CFLAGS_KOA={mac} $(INCLUDE)\' "BUILD_CF={mac} [mac rehearsal]"\n'),
                    ("ldd", 'exec otool -L "$@"\n')):
        p = os.path.join(d, n); open(p, "w").write("#!/bin/bash\n# rehearsal build shim\n" + body); os.chmod(p, 0o755)
    return d


def blocks(path):
    """the runsheet's paste blocks: header line 'P2 BLOCK X (...)' -> its command lines (up to the first empty line)"""
    out, cur = {}, None
    for l in open(path).read().split("\n"):
        m = re.match(r"^P2 BLOCK (A0|A|B|C|D|D2|E) \(", l)
        if m: cur = m.group(1); out[cur] = []; continue
        if cur is not None:
            if not l.strip(): cur = None if out[cur] else cur; continue
            if l.startswith("  "): cur = None; continue
            out[cur].append(l)
    return out


def run_bash(lines, env, cwd, log, title):
    """the block's lines in one bash, each echoed as '$ <line>' before it runs (a terminal's transcript; stderr into stdout)"""
    script = "".join(f"printf '%s\\n' {shlex.quote('$ ' + l)}\n{l}\n" for l in lines)
    t0 = time.time(); p = subprocess.run(["/bin/bash"], input=script, env=env, cwd=cwd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
    txt = p.stdout
    open(log, "w").write(f"### {title}\n### exit {p.returncode}, {time.time() - t0:.0f} s\n" + txt)
    print(f"  {title}: exit {p.returncode}, {time.time() - t0:.0f} s, log {os.path.basename(log)}", flush=True)
    return p.returncode, txt


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--pkg", required=True); ap.add_argument("--root", required=True); ap.add_argument("--out", required=True)
    a = ap.parse_args(); R = os.path.abspath(a.root); O = os.path.abspath(a.out)
    for d in (R, O):
        if os.path.exists(d): sys.exit(f"REFUSED: {d} exists (a new folder only)")
    os.makedirs(os.path.join(R, "logs")); os.makedirs(O); home = os.path.join(R, "home"); os.makedirs(home)
    LUS = os.path.join(R, "lustre"); runsheet = os.path.join(HS, P, "runsheet_p2_261016.txt"); B = blocks(runsheet)
    step("the runsheet's blocks", set(B) == {"A0", "A", "B", "C", "D", "D2", "E"}, ", ".join(f"{k}: {len(v)} lines" for k, v in B.items()))
    # the stand-in of the package commit
    S = os.path.join(R, "standin"); sh_ = K2R.clone(S, a.pkg, K2R.SPARSE_KOA); kissfft(S)
    ct = os.path.join(S, "hspist3", P, "chunks_p2_koa_261016.tsv"); L = open(ct).read().split("\n")
    keep = [L[0]] + [l for l in L[1:] if l and int(l.split("\t")[0]) in MINI]; open(ct, "w").write("\n".join(keep) + "\n")
    done = K2R.rewrite_lustre(S, R, local_fetch=True)
    mk = subprocess.run([sys.executable, "-B", os.path.join(S, "hspist3", P, "make_p2_manifests.py")], cwd=os.path.join(S, "hspist3"), capture_output=True, text=True)
    for c in (["config", "user.email", "rehearsal@localhost"], ["config", "user.name", "P2 rehearsal"], ["add", "-A", "hspist3/" + P],
              ["commit", "-q", "-m", "P2 rehearsal stand-in: the mini chunk table, local Lustre and fetch, the input manifest again"], ["branch", "standin"]):
        K2R.sh(["git", "-C", S] + c)
    SC = K2R.sh(["git", "-C", S, "rev-parse", "HEAD"]).stdout.strip()
    step("the stand-in commit", sh_ == a.pkg and len(keep) == 5 and mk.returncode == 0,
         f"{SC[:12]} on {a.pkg[:12]}: chunk table {len(keep) - 1} chunks {MINI}; Lustre/fetch made local in {', '.join(done)}; manifests: {mk.stdout.strip().splitlines()[1] if mk.returncode == 0 else mk.stderr[-200:]}")
    shim = K2R.shims(R); bsh = build_shims(R)
    env = K2R.env_for(R, home, shim, REHEARSAL_ROOT=R); envA = dict(env, PATH=f"{bsh}:{env['PATH']}", SLURM_JOB_ID="930001")
    # block A
    sub = {f"cd ~ && git clone --depth 1 --branch main --filter=blob:none --sparse https://github.com/FreeworkEarth/HardDisks.git harddisks_gen3_p2":
               f"cd ~ && git clone -q --no-checkout --shared {REPO} harddisks_gen3_p2 && git -C harddisks_gen3_p2 config core.abbrev 7 && git -C harddisks_gen3_p2 sparse-checkout init --no-cone",
           f"cd ~/harddisks_gen3_p2 && git fetch --depth 1 origin {BUILD}": "cd ~/harddisks_gen3_p2",
           "git clone https://github.com/mborgerding/kissfft.git hspist3/kissfft": f"git clone -q {os.path.join(HS, 'kissfft')} hspist3/kissfft",
           f"git fetch --depth 1 origin {a.pkg}": f"git fetch -q {S} standin",
           f"git -c advice.detachedHead=false checkout {a.pkg}": f"git -c advice.detachedHead=false checkout {SC}"}
    linesA = [sub.get(l, l).replace("/mnt/lustre/koa/scratch/charing", LUS) for l in B["A"] if l != "exit"]
    nsub = sum(1 for l in B["A"] if l in sub)
    rc, out = run_bash(linesA, envA, home, os.path.join(O, "block_A.txt"), "block A (the sandbox shell)")
    ok = ("BUILD OK" in out and "BUILD RECORD OK:" in out and "the merge and reducer imports: OK" in out and "00ALLINONE  git 8aaf5fe  target koa" in out
          and re.search(r"build_inputs (\d+) of \1 equal", out) and "No such file or directory" in out.split("\nthe merge and reducer imports: OK")[-1])
    step("block A: clone, build at 8aaf5fef, package checkout, record, checks", rc == 0 and bool(ok), f"{nsub} lines substituted; exit {rc}; "
         + "; ".join(x for x in re.findall(r"^(commit_build \S+|version .*|build_inputs .*|\d+$)", out, re.M)[:6]))
    H2 = os.path.join(home, "harddisks_gen3_p2", "hspist3")
    # block B, twice
    linesB = [l.replace("--array=1-400%64", "--array=" + ",".join(map(str, MINI))).replace("/mnt/lustre/koa/scratch/charing", LUS).replace(a.pkg, SC).replace(a.pkg[:8], SC[:8])
              for l in B["B"]]                                      # the guard's commit: the stand-in's (the build record names it)
    rc, out = run_bash(linesB, env, home, os.path.join(O, "block_B_first.txt"), "block B (the array; the shim runs the 4 tasks one after the other)")
    exits = open(os.path.join(R, "logs", "sbatch_shim_exits.txt")).read() if os.path.exists(os.path.join(R, "logs", "sbatch_shim_exits.txt")) else ""
    step("block B: the guarded submit", "CHECKS GOOD" in out and "P2: THE ARRAY IS SUBMITTED" in out and exits.count(" 0\n") >= 4,
         "; ".join(l for l in out.splitlines() if l.startswith(("CHECKS", "Submitted", "P2:"))) + f"; task exits: {[l.split()[-1] for l in exits.splitlines() if '_' in l.split()[0]]}")
    rc, out = run_bash(linesB, env, home, os.path.join(O, "block_B_second.txt"), "block B pasted again")
    step("block B a second time: refused", "REFUSED:" in out and "P2: THIS PASTE DID NOT SUBMIT THE ARRAY" in out and "Submitted" not in out,
         "; ".join(l for l in out.splitlines() if l.startswith(("REFUSED", "P2:"))))
    # block C
    linesC = [l.replace("/mnt/lustre/koa/scratch/charing", LUS) for l in B["C"]]
    rc, out = run_bash(linesC, env, home, os.path.join(O, "block_C.txt"), "block C (the checks)")
    step("block C: 4 logs, 4 done, 0 FAILED/REFUSED/STOP, the chunk check", re.search(r"chunks done: 4 of 4", out) is not None and
         re.findall(r"^\s*(\d+)\s*$", out, re.M)[:3] == ["4", "4", "0"], " | ".join(l.strip() for l in out.splitlines() if l.strip() and not l.startswith(("squeue", "cd ", "ls ", "grep ", "cat ", "bash "))))
    # block D, then the merge again without the guard
    rc, out = run_bash(B["D"], env, home, os.path.join(O, "block_D.txt"), "block D (the merge job; the shim runs it at once)")
    step("block D: the merge submitted", "Submitted batch job" in out, "; ".join(l for l in out.splitlines() if l.startswith("Submitted")))
    rc, out = run_bash(B["D2"], env, home, os.path.join(O, "block_D2.txt"), "block D2 (the merge's log)")
    step("block D2: the merge and p2_reduce", "P2 merge and reduction done" in out and re.search(r"p2_reduce: 4 folders \(0 without the stage 1 tables\), 4 trajectories, 0 with an error", out) is not None,
         "; ".join(l[:120] for l in out.splitlines() if l.startswith(("merged ", "p2_reduce:", "P2 merge"))))
    rc, out = run_bash([f"cd {H2} && sbatch {P}/p2_merge.sbatch", f"tail -2 \"$(ls -t {H2}/logs/gen3-p2-merge_*.out | head -1)\""], env, home,
                       os.path.join(O, "block_D_merge_again.txt"), "the merge submitted again (no guard)")
    step("the merge a second time: refused (exit 3)", "exists -- a merge ran before; not merging twice" in out,
         "; ".join(l[:160] for l in out.splitlines() if "STOP" in l or "exists" in l))
    # block E
    linesE = [l.replace("cd ~/Desktop/CCS_complex_coupled_systems/Repo/HardDisks", f"cd {S}") for l in B["E"]]
    rc, out = run_bash(linesE, env, home, os.path.join(O, "block_E.txt"), "block E (the fetch, local)")
    F = os.path.join(S, "hspist3", "experiments_gen3_stage2_koa_261016", "P2")
    step("block E: the fetch", "fetched into" in out and "p2_reduce: 4 folders" in out, out.strip().splitlines()[-1][:220])
    big = [os.path.join(dp, f) for dp, dn, fn in os.walk(F) for f in fn if f.startswith(("wall_x_positions_", "snapshots_"))]
    dumps = [os.path.join(dp, f) for dp, dn, fn in os.walk(F) for f in fn if f == "initial_state.txt"]
    step("the fetch's filter: no trace, no snapshot; the 4 initial states", not big and len(dumps) == 4, f"traces/snapshots fetched: {len(big)}; initial states: {len(dumps)}")
    shutil.copytree(os.path.join(F, "merged"), os.path.join(O, "fetched_merged"), ignore=shutil.ignore_patterns(".done_runs"))
    t = subprocess.run([sys.executable, "-B", os.path.join(HS, P, "p2_tables_mac.py"), "--data", F, "--out", os.path.join(O, "p2_tables")], capture_output=True, text=True)
    open(os.path.join(O, "p2_tables_stdout.txt"), "w").write(t.stdout + t.stderr)
    step("p2_tables_mac.py on the fetched folder", t.returncode == 0, f"exit {t.returncode}; {len(t.stdout.splitlines())} lines")
    print("\n== the rehearsal\n\n| step | result | detail |\n|---|---|---|")
    rows = [f"| {lab} | {'PASS' if ok else 'FAIL'} | {det} |" for lab, ok, det in STEPS]; print("\n".join(rows))
    tail = f"\nP2 REHEARSAL: {'PASS' if all(s[1] for s in STEPS) else 'FAIL'} ({sum(s[1] for s in STEPS)} of {len(STEPS)})"
    print(tail); open(os.path.join(O, "REHEARSAL.md"), "w").write("| step | result | detail |\n|---|---|---|\n" + "\n".join(rows) + "\n" + tail + "\n")
    print(f"data under {R}: {K2R.sh(['du', '-sh', R]).stdout.split()[0]}; under {O}: {K2R.sh(['du', '-sh', O]).stdout.split()[0]}")


if __name__ == "__main__":
    main()
