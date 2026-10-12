#!/usr/bin/env python3
"""##CHRIS 2026-10-11 (plan-author decisions 29 and 32, task P; 261012 sec. 4.7.59 and 4.7.63): THE MAC REHEARSAL OF THE PREVIEW SWEEP --
PREVIEW, INFORMATION. Every script the KOA jobs run, from the runsheet's paste blocks as Chris types them, before anything is typed on
KOA. NOT RUN before sec. 4.8 (the stage 1 registration with A1 and A2) is committed (decisions 29 and 32): the rehearsal runs preview
trajectories. It STOPS unless the notes at HEAD hold sec. 4.8. Nothing outside --root is written; the clone borrows the repository's
objects (git clone --shared). At most one simulation process at a time.
THIS FILE IS cluster/gen3_testR_261011/rehearse_testR_mac.py (Test R's, run 2 with its fixed block D check) adapted to the preview:
  - the runsheet runsheet_preview_261011.txt, "PREVIEW BLOCK <X> (" to the first empty line;
  - the homes as on KOA after Test R: home/harddisks_gen3 (the Mac build of a643bc7 under the koa target's name, built or --reuse-binary,
    and the build-time hash file) and home/harddisks_gen3_k2 at Test R's delivery b5378f4 with its sparse set, its logs/ holding a
    pretended Test R merge log and a Test R chunk log; block B UPDATES IT IN PLACE to the preview's delivery commit;
  - the MINI LIST: the registered lines with r = 0 and 1 at M = 50 for every cell, plus pi/8 at M = 100 and 200, in their registered
    order: 18 lines; chunk size 9 (two chunks, every cell folder receiving from both); so E0's slope has three masses at pi/8 and the
    recorded-length checks see all seven cells, the three truncated ones included;
  - after block B, in the updated clone: the mini list and its SHA-256 file, the chunk size 9, Test R's build record with the Mac
    binary's SHA-256 (the preview reads that file), the preview's two manifests rewritten by make_preview_manifests.py; the Lustre path
    rewritten in the preview's job scripts; the fetch's DTN source made local;
  - block E's --array=2-25 becomes --array=2; block F's chunk count 25 becomes 2; block D's "masses 9" lines become the mini list's;
  - after the fetch: E0 (e0_preview_261011.py) and the figure (fig_preview_261011.py, labelled PREVIEW REHEARSAL) on the mini data;
  - decision 32, item 3 (information): one trajectory (the mini list's first pi/8 line) run twice with the worker's B command, with and
    without HD_CONTACT_AUDIT=1; its event hash and its trace compared.
THE SHIMS (the residual differences, printed in the table): the Slurm commands (squeue, sacct, srun, sbatch, sinfo), module, gcc
--version, lscpu, hostname, sha256sum (fails with REHEARSAL_STARTUP_FAIL: a pretended faulty node) -- Test R's, unchanged.
usage (from the repository root): python3 hspist3/cluster/gen3_preview_261011/rehearse_preview_mac.py --root <empty dir> [--pkg <commit>]
                                    [--reuse-binary <Mac koa-name binary of a643bc7>]
"""
import argparse, hashlib, os, re, shutil, subprocess, sys, time
HERE = os.path.dirname(os.path.abspath(__file__)); HS = os.path.dirname(os.path.dirname(HERE)); REPO = os.path.dirname(HS)
GEN3_WT = os.path.join(os.path.dirname(REPO), "HardDisks_gen3")
sys.path.insert(0, os.path.join(os.path.dirname(HERE), "gen3_testR_261011"))
import rehearse_testR_mac as RT                                   # Test R's shims, env, build and run_block, imported unchanged
A643, K2_NOW = RT.A643, "b5378f49bbf3caf39df761aa49fe19cee693e3a5"   # KOA's ~/harddisks_gen3_k2 since Test R's block B
PKG_DEFAULT = "c70d314e927580e35267045d58a22fdf10fb9850"
LUSTRE, SPARSE_K2 = RT.LUSTRE, RT.SPARSE_K2
P, RP = "cluster/gen3_preview_261011", "cluster/gen3_testR_261011"
NOTES = "0000_PLAN_OVERALL/ALL_MARKDOWNS/261012_paper1_confinement.md"
MINI_K = 9
CELLS = ("eta_0p052360", "eta_0p112200", "eta_0p196350", "eta_0p261799", "eta_0p392699", "eta_0p523599", "eta_0p590")
RES = []
sh, sha, L = RT.sh, RT.sha, RT.L


def res(script, crit, ok, ev):
    RES.append((script, crit, "PASS" if ok else "**FAIL**", ev)); print(f"    -> {script}: {'PASS' if ok else '**FAIL**'} ({ev})", flush=True)


def mini_lines(reg):
    out = []
    for l in reg:
        p = l.split(); cid, M, r = p[1].split("/")[2], int(p[2]), int(p[3])
        if r in (0, 1) and (M == 50 or (cid == "eta_0p392699" and M in (100, 200))): out.append(l)
    return out


def k2_clone(home):
    """KOA's ~/harddisks_gen3_k2 after Test R: Test R's delivery b5378f4, its sparse set, a pretended Test R merge log and chunk log"""
    d = os.path.join(home, "harddisks_gen3_k2")
    sh(["git", "clone", "-q", "--no-checkout", "--shared", REPO, d]); sh(["git", "-C", d, "config", "core.abbrev", "7"])
    sh(["git", "-C", d, "config", "core.sparseCheckout", "true"])
    open(os.path.join(d, ".git", "info", "sparse-checkout"), "w").write("\n".join(SPARSE_K2) + "\n")
    sh(["git", "-C", d, "-c", "advice.detachedHead=false", "checkout", "-q", K2_NOW])
    g = os.path.join(d, "hspist3", "logs"); os.makedirs(g)
    open(os.path.join(g, "gen3-testR-merge_15060000.out"), "w").write("merge of .../testR_koa by the scripts of ~/harddisks_gen3_k2 ...\nmerge done\n")
    open(os.path.join(g, "gen3-testR_15050001_1.out"), "w").write("env: ...\nchunk 1 on cn-03-38-07 (E5): lines 1-100; 00ALLINONE  git a643bc7  target koa\nchunk 1 done\n")
    return d


def blocks_of(path):
    out, cur, shas, last = {}, None, {}, None
    for l in open(path).read().splitlines():
        m = re.match(r"^PREVIEW BLOCK ([A-Z]) \(", l)
        if m: cur = last = m.group(1); out[cur] = []; shas[cur] = []; continue
        if cur is not None and not l.strip(): cur = None; continue
        if cur is not None: out[cur].append(l); continue
        q = re.match(r"^\s+([0-9a-f]{64})  (\S+)$", l)
        if last and q: shas[last].append(f"{q.group(1)}  {q.group(2)}")
    return out, shas


def substitute(R, hs, mac_sha):
    reg = open(os.path.join(hs, P, "tasks_preview_koa_261011.txt")).read().splitlines(keepends=True)
    mini = "".join(mini_lines(reg))
    open(os.path.join(hs, P, "tasks_preview_koa_261011.txt"), "w").write(mini)
    open(os.path.join(hs, P, "tasks_preview_koa_261011.sha256"), "w").write(hashlib.sha256(mini.encode()).hexdigest() + "\n")
    open(os.path.join(hs, P, "preview_chunk_size.txt"), "w").write(f"{MINI_K}\n")
    rec = os.path.join(hs, RP, "build_record_a643bc7.txt"); s = open(rec).read()
    open(rec, "w").write(re.sub(r"^sha256 [0-9a-f]{64}$", f"sha256 {mac_sha}", s, flags=re.M))
    done = []
    for f in ("preview_koa.sbatch", "preview_merge.sbatch", "fetch_preview_koa.sh"):
        p = os.path.join(hs, P, f); t = open(p).read()
        t2 = t.replace(LUSTRE, os.path.join(R, "lustre")).replace("charing@koa-dtn.its.hawaii.edu:", "")
        if t2 != t: open(p, "w").write(t2); done.append(f)
    out = sh([sys.executable, os.path.join(P, "make_preview_manifests.py")], cwd=hs).stdout
    return mini.splitlines(), done, out


def expected_D(runsheet):
    """block D's expected per-cell lines, with the mini list's counts: 1 mass per cell (pi/8: 3), 1 trajectory each"""
    lines = open(runsheet).read().splitlines(); i = lines.index(next(l for l in lines if l.startswith("PREVIEW BLOCK D")))
    want = {}
    for l in lines[i:]:
        m = re.match(r"^\s+(eta_\S+) masses 9: traces 9 runs 9 gen3-records-clean 9 build-lines-ok 9$", l)
        if m: cur = m.group(1); continue
        q = re.match(r"^\s+9 (box=.*)$", l)
        if q: want[cur] = q.group(1)
        if l.startswith("  STOP RULE"): break
    return want


def check_D(out, want):
    Lo = L(out); okc = 0
    for c in CELLS:
        n = 3 if c == "eta_0p392699" else 1
        okc += (f"{c} masses 9: traces {n} runs {n} gen3-records-clean {n} build-lines-ok {n}" in Lo) and (f"{n} {want[c]}" in Lo)
    return (okc == 7 and "chunk 1 done" in out and "lines 1-9; 00ALLINONE  git a643bc7  target koa" in out
            and not re.search(r"FAILED|REFUSED|^STOP", out, re.M)), okc


def audit_pair(R, binp, line):
    """decision 32 item 3: the worker's B command run twice, HD_KE_TRACE=1 with and without HD_CONTACT_AUDIT=1; event hash and trace"""
    p = line.split(); M, seed, L0, H, NS, stride, base = p[2], p[4], p[5], p[6], p[7], p[8], p[9]
    out = {}
    for tag, extra in (("audit", {"HD_CONTACT_AUDIT": "1"}), ("plain", {})):
        d = os.path.join(R, "audit_pair", tag); os.makedirs(d)
        env = {k: v for k, v in os.environ.items() if k not in ("HD_CONTACT_AUDIT",)}; env.update(HD_KE_TRACE="1", **extra)
        cmd = [binp, "--mode=edmd", "--experiment=speed_of_sound", "--headless", "--kbt1", "--seed-drift-order=drift-first", "--edmd-acc=0",
               f"--particles={2 * int(NS)}", f"--particles-boxes={NS},{NS}", f"--height={H}", "--particle-radius=0.5", "--wall-thickness=0.05",
               "--wall-thickness-vis=0.05", f"--lengths={L0}", f"--wall-masses={M}", "--repeats=1", f"--seed={base}", "--wall-hold-steps=2000",
               "--fixed-dt=0.4", "--target-oscillations=200", "--oscillation-safety=1.0", "--oscillation-min-steps=10000",
               "--oscillation-max-steps=400000000", f"--speed-sound-log-stride={stride}", f"--speed-sound-run-dir={d}",
               f"--speed-sound-exact-seed={seed}", "--engine=gen3"]
        q = subprocess.run(cmd, cwd=d, env=env, capture_output=True, text=True)
        open(os.path.join(d, "stdout.log"), "w").write(q.stdout + q.stderr)
        h = re.search(r"^\[EDMD3-HEALTH\] .*?clean=(\d).*? hash=([0-9a-f]+)", q.stdout, re.M)
        tr = [f for f in os.listdir(d) if f.startswith("wall_x_positions_")]
        out[tag] = dict(rc=q.returncode, clean=h.group(1) if h else None, hash=h.group(2) if h else None,
                        trace=sha(os.path.join(d, tr[0])) if len(tr) == 1 else None,
                        audit_lines=len(re.findall(r"^\[EDMD-CONTACT\]|^\[EDMD3-GAP\]", q.stdout, re.M)))
    return out


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--root", required=True); ap.add_argument("--pkg", default=PKG_DEFAULT)
    ap.add_argument("--runsheet", default=os.path.join(HERE, "runsheet_preview_261011.txt")); ap.add_argument("--reuse-binary", default=None)
    a = ap.parse_args(); R = os.path.abspath(a.root)
    notes = sh(["git", "-C", REPO, "show", f"HEAD:{NOTES}"]).stdout
    if "\n### 4.8 " not in notes: sys.exit("STOP: sec. 4.8 (the stage 1 registration with A1 and A2) is not committed: no preview trajectory "
                                           "before it (decisions 29 and 32)")
    if not os.path.isdir(R) or os.listdir(R): sys.exit(f"STOP: {R} must be an existing empty folder")
    PKG = sh(["git", "-C", REPO, "rev-parse", a.pkg]).stdout.strip()
    for d in ("home", "lustre", "logs"): os.makedirs(os.path.join(R, d))
    home = os.path.join(R, "home"); shim = RT.shims(R)
    blocks, shas = blocks_of(a.runsheet)
    print(f"# PREVIEW REHEARSAL: package {PKG} (engine-gen3), runsheet {a.runsheet} (blocks {''.join(sorted(blocks))}), root {R}\n", flush=True)
    m = re.search(r"git fetch --depth 1 origin ([0-9a-f]{40})", "\n".join(blocks.get("B", [])))
    res("runsheet block B", "fetches the package commit", bool(m) and m.group(1) == PKG, f"block B fetches {m.group(1) if m else None}")
    if a.reuse_binary:
        hs0 = os.path.join(home, "harddisks_gen3", "hspist3"); os.makedirs(os.path.join(hs0, "logs"))
        os.symlink(os.path.abspath(a.reuse_binary), os.path.join(hs0, "00ALLINONE"))
        mac_sha = sha(a.reuse_binary); open(os.path.join(hs0, "logs", "BUILD_KOA_LAST.hash"), "w").write(f"{mac_sha}  00ALLINONE\n")
        ver = sh([os.path.join(hs0, "00ALLINONE"), "--version"]).stdout.splitlines()[0]
        print(f"  the binary: {a.reuse_binary} (symlinked, not rebuilt; SHA-256 {mac_sha[:12]}...)", flush=True)
    else:
        ver, mac_sha = RT.build_binary(R, home)
    binp = os.path.join(home, "harddisks_gen3", "hspist3", "00ALLINONE")
    res("the Mac binary (a643bc7, koa target name)", "version line '00ALLINONE  git a643bc7  target koa'", ver == "00ALLINONE  git a643bc7  target koa", ver)
    k2 = k2_clone(home); k2hs = os.path.join(k2, "hspist3")
    E = RT.env_for(R, home, shim)
    rc, out = RT.run_block(R, "A", blocks["A"], E)
    res("block A", "squeue, srun (shims)", rc == 0 and "rehearsal shim: srun" in out, f"exit {rc}")
    rc, out = RT.run_block(R, "B", blocks["B"], RT.env_for(R, home, shim, SLURM_JOB_ID=950001))
    got = re.findall(r"^[0-9a-f]{64}  cluster/gen3_\S+$", out, re.M)
    ok_b = ["0" in L(out), f"HEAD is now at {PKG[:7]}" in out, got == shas["B"], L(out).count("20") >= 2,
            "00ALLINONE  git a643bc7  target koa" in out, out.count("ALL PASS") == 3, "the E0 and nu_d imports: OK" in out,
            "No such file or directory" in out]
    res("block B (the update in place, the checks, the unit tests)", "as the runsheet's good lines, but KOA's four binaries (not on the Mac)",
        all(ok_b), f"checks {ok_b}; the eight SHA lines {'equal' if got == shas['B'] else 'DIFFER'}; KOA's binaries: "
        f"{len(re.findall(r': FAILED', out))} FAILED here (expected: the Mac build and three absent files)")
    mini, rew, man = substitute(R, k2hs, mac_sha)
    res("the rehearsal substitutions", "mini list (18 lines, chunk size 9), build record, manifests, Lustre path, local fetch",
        len(mini) == 18 and len(rew) == 3, f"rewritten {rew}; manifests: {man.count('files)')} written")
    rc, out = RT.run_block(R, "C", blocks["C"], E)
    res("block C (chunk 1 alone)", "Submitted batch job", rc == 0 and "Submitted batch job" in out, f"exit {rc}")
    rc, out = RT.run_block(R, "D", blocks["D"], E)
    ok_d, n_d = check_D(out, expected_D(a.runsheet))
    res("block D (chunk 1's checks)", "chunk 1 on ... lines 1-9; chunk 1 done; per cell its counts and the runsheet's recorded box and eta; "
        "no FAILED/REFUSED/STOP", ok_d, f"{n_d}/7 cells")
    linesE = [l.replace("--array=2-25", "--array=2") for l in blocks["E"]]
    rc, out = RT.run_block(R, "E", linesE, RT.env_for(R, home, shim, REHEARSAL_STARTUP_FAIL=1, REHEARSAL_HOST="cn-07-07-07"))
    lg = [p for p in os.listdir(os.path.join(k2hs, "logs")) if p.startswith("gen3-preview_") and p.endswith("_2.out")]
    l2 = open(os.path.join(k2hs, "logs", lg[0])).read() if len(lg) == 1 else ""
    res("block E (the rest, --array=2 here; the first attempt on a pretended faulty node)", "STOP: [cn-07-07-07] ... before 'chunk 2 on'; no folder",
        "STOP: [cn-07-07-07] the preview inputs are not the delivered ones" in l2 and "chunk 2 on" not in l2
        and not os.path.exists(os.path.join(R, "lustre", "harddisks_gen3", "preview_koa", "chunks", "chunk_2")), l2.strip().splitlines()[-1] if l2 else "no log")
    linesF = [l.replace("preview_stopped_chunks.sh logs gen3-preview 25", "preview_stopped_chunks.sh logs gen3-preview 2") for l in blocks["F"]]
    rc, out = RT.run_block(R, "F", linesF, E)
    rl = re.findall(r"^sbatch --array=\S+ --exclude=\S+ cluster/gen3_preview_261011/preview_koa\.sbatch$", out, re.M)
    res("block F (checks, the resubmission line)", "2 logs, 1 done, 0 FAILED; the line --array=2 --exclude=cn-05-08-07,cn-07-07-07",
        rl == ["sbatch --array=2 --exclude=cn-05-08-07,cn-07-07-07 cluster/gen3_preview_261011/preview_koa.sbatch"], f"{rl}")
    rc, out = RT.run_block(R, "R", ["squeue -u charing", f"cd ~/harddisks_gen3_k2/hspist3 && {rl[0]}"] if rl else ["true"], E)
    rc, out = RT.run_block(R, "F", linesF, E, "2")
    res("block F again (after the resubmission)", "3 logs, 2 done, 0 FAILED, nothing to resubmit",
        "3" in L(out) and "2" in L(out) and "0" in L(out) and "nothing to resubmit" in out, "")
    q = os.path.join(R, "logs", "pretend_squeue.txt"); open(q, "w").write("930999_7 PENDING\n")
    rc, out = RT.run_block(R, "G", blocks["G"], RT.env_for(R, home, shim, REHEARSAL_SQUEUE=q), "_pending")
    ex = open(os.path.join(R, "logs", "sbatch_shim_exits.txt")).read().splitlines()[-1]
    mlog = sorted(p for p in os.listdir(os.path.join(k2hs, "logs")) if p.startswith("gen3-preview-merge_"))
    t = open(os.path.join(k2hs, "logs", mlog[-1])).read() if mlog else ""
    res("the merge job: a pending gen3-preview job", "exit 2, nothing moved", ex.endswith(" 2") and "still pending or running" in t
        and not os.path.exists(os.path.join(R, "lustre", "harddisks_gen3", "preview_koa", "merged")), ex)
    lg2 = sorted(p for p in os.listdir(os.path.join(k2hs, "logs")) if re.match(r"^gen3-preview_\d+_2\.out$", p))
    done2 = [p for p in lg2 if "chunk 2 done" in open(os.path.join(k2hs, "logs", p)).read()]
    hold = os.path.join(R, "logs", "held_aside_" + done2[0]); shutil.move(os.path.join(k2hs, "logs", done2[0]), hold)
    rc, out = RT.run_block(R, "G", blocks["G"], E, "_heldaside")
    ex = open(os.path.join(R, "logs", "sbatch_shim_exits.txt")).read().splitlines()[-1]
    mlog = sorted(p for p in os.listdir(os.path.join(k2hs, "logs")) if p.startswith("gen3-preview-merge_"))
    t = open(os.path.join(k2hs, "logs", mlog[-1])).read()
    shutil.move(hold, os.path.join(k2hs, "logs", done2[0]))
    res("the merge job: chunk 2's done log held aside", "exit 2, 'NOT DONE: 2(no log with its done line)', nothing moved",
        ex.endswith(" 2") and "NOT DONE: 2(no log with its done line)" in t
        and not os.path.exists(os.path.join(R, "lustre", "harddisks_gen3", "preview_koa", "merged")), ex)
    rc, out = RT.run_block(R, "G", blocks["G"], E)
    ex = open(os.path.join(R, "logs", "sbatch_shim_exits.txt")).read().splitlines()[-1]
    rc, out = RT.run_block(R, "H", blocks["H"], E)
    okH = [ex.endswith(" 0"), "merge done; next: the fetch" in out, "scan: 2 chunks, 108 plan entries, 0 problem(s)" in out,
           "plan entries by treatment: APPEND 36, EQUAL 18, MOVE 36, PER CHUNK 18" in out,
           "merged 2 chunks: 36 files moved, 36 files appended, 18 kept per chunk (0 failed-run folders), 0 incomplete run folders, "
           "0 locks or temp files left in place, 18 build files equal; build '00ALLINONE  git a643bc7  target koa'" in out,
           out.count(": exit 0 2 trajectories: nud_runs.csv (0 fit failures), runs_health.csv (2 ##RUN sections), failed_runs.csv (0 failed runs)") == 9,
           out.count("2 traces compressed") == 9, "earlier attempt of chunk 2, stopped in the start-up checks (ran nothing)" in out]
    res("blocks G and H (the merge, the reduction, the compression; its record)", "exit 0; 108 plan entries; the merged line; 9 reductions; "
        "9 compressions; chunk 2's earlier attempt listed", all(okH), f"checks {okH}")
    rc, out = RT.run_block(R, "G", blocks["G"], E, "_second")
    ex = open(os.path.join(R, "logs", "sbatch_shim_exits.txt")).read().splitlines()[-1]
    res("a second merge submission", "exit 3 (merged/ exists)", ex.endswith(" 3"), ex)
    p = sh(["bash", os.path.join(P, "fetch_preview_koa.sh")], cwd=k2hs, env=E, log=os.path.join(R, "logs", "fetch.out"), check=False)
    F = os.path.join(k2hs, "experiments_gen3_preview_koa_261011")
    fetched = [os.path.relpath(os.path.join(dp, f), F) for dp, dn, fn in os.walk(F) for f in fn]
    raw = [f for f in fetched if "wall_x_positions_" in f or "psi6_t_" in f or ".done_runs" in f or f.startswith("chunks/")]
    res("the fetch (local rsync)", "exit 0; the merge record, the tables, the job logs; no trace, psi6 or run folder",
        p.returncode == 0 and "merged/.merge_complete" in fetched and not raw and sum(f.endswith("red_nu.csv") for f in fetched) == 9,
        f"{len(fetched)} files; raw {len(raw)}")
    # E0 and the figure on the mini data (from the repository, as after KOA's fetch)
    e0csv = os.path.join(R, "e0_mini.csv")
    e = sh([sys.executable, os.path.join(HS, P, "e0_preview_261011.py"), "--preview", os.path.join(F, "merged"), "--tasks",
            os.path.join(k2hs, P, "tasks_preview_koa_261011.txt"), "--out", e0csv], cwd=HS, log=os.path.join(R, "logs", "e0_mini.out"), check=False)
    import pandas as pd
    rec = pd.read_csv(e0csv.replace(".csv", "_recorded.csv")) if os.path.exists(e0csv.replace(".csv", "_recorded.csv")) else None
    rows = pd.read_csv(e0csv) if os.path.exists(e0csv) else None
    chk_cols = [c for c in (rec.columns if rec is not None else []) if c.startswith(("the run's", "box width", "box eta", "exact_box", "N = "))]
    res("E0 on the mini data (e0_preview_261011.py --preview)", "exit 0; the recorded-length checks pass in all 9 (cell, mass) folders "
        "(the three truncated cells included); one E0 row (pi/8: three masses)",
        e.returncode == 0 and rec is not None and len(rec) == 9 and all(rec[c].all() for c in chk_cols) and rows is not None and len(rows) == 1,
        f"exit {e.returncode}; recorded rows {0 if rec is None else len(rec)}; E0 rows {0 if rows is None else len(rows)}")
    f = sh([sys.executable, os.path.join(HS, P, "fig_preview_261011.py"), "--preview", e0csv, "--out", os.path.join(R, "figure"), "--label",
            "PREVIEW REHEARSAL (mini list, 2 seeds; not data)"], cwd=HS, log=os.path.join(R, "logs", "figure_mini.out"), check=False)
    ft = open(os.path.join(R, "logs", "figure_mini.out")).read()
    res("the figure and the table on the mini data", "exit 0; the pi/8 row 'the same cell' as A1 v2's", f.returncode == 0 and "the same cell" in ft, f"exit {f.returncode}")
    # decision 32, item 3: HD_CONTACT_AUDIT on and off
    pi8 = [l for l in mini if l.split()[1].split("/")[2] == "eta_0p392699"][0]
    ap_ = audit_pair(R, binp, pi8)
    same = ap_["audit"]["hash"] is not None and ap_["audit"]["hash"] == ap_["plain"]["hash"] and ap_["audit"]["trace"] == ap_["plain"]["trace"]
    RES.append(("decision 32 item 3 (information): one trajectory with and without HD_CONTACT_AUDIT=1", "event hash and trace compared",
                "INFORMATION", f"{pi8.split()[1]} r {pi8.split()[3]}: hash {ap_['audit']['hash']} / {ap_['plain']['hash']}; trace SHA-256 "
                f"{(ap_['audit']['trace'] or '-')[:12]} / {(ap_['plain']['trace'] or '-')[:12]}; equal: {'yes' if same else 'NO'}; clean "
                f"{ap_['audit']['clean']} / {ap_['plain']['clean']}; audit lines {ap_['audit']['audit_lines']} / {ap_['plain']['audit_lines']}"))
    print(f"    -> audit on/off: equal {'yes' if same else 'NO'}", flush=True)
    print("\n| step (script) | criterion | result | evidence |\n|---|---|---|---|")
    for r in RES: print("| " + " | ".join(str(x) for x in r) + " |")
    print("\nRESIDUAL DIFFERENCES: the shims (Slurm commands, module, gcc --version, lscpu, hostname, sha256sum); the Mac binary (its SHA-256 "
          f"in Test R's build record; KOA's four binaries absent{'; symlinked from an earlier run, not rebuilt' if a.reuse_binary else ''}); "
          "the mini list (18 lines, chunk size 9) with its SHA-256 file and the manifests rewritten; the Lustre path; the local fetch; block "
          "E's --array=2; block F's chunk count 2; block D's counts the mini list's")
    nb = sum(r[2] not in ("PASS", "INFORMATION") for r in RES)
    print(f"\nPREVIEW REHEARSAL: {'PASS' if nb == 0 else 'FAIL'} ({sum(r[2] == 'PASS' for r in RES)} of {len(RES) - 1}; plus one INFORMATION row)")


if __name__ == "__main__":
    main()
