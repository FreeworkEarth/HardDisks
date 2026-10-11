#!/usr/bin/env python3
"""##CHRIS 2026-10-10 (plan-author decision 20, task 3; 261012 sec. 4.7.47): THE FULL MAC REHEARSAL of every script a KOA job executes,
end to end, before Chris types it (decision 20's new rule), in scratch clones of the commits KOA uses. Nothing outside --root is
written; the clones borrow the repository's objects (git clone --shared) and are checked out with KOA's sparse sets.
  home/harddisks_gen3       engine-gen3 a643bc7, KOA's step-1 sparse set and the pinned kissfft -- as ~/harddisks_gen3 on KOA; its
                            binaries built here (below). Runs what K2 runs from that clone: the cross-node task and its compare (as in
                            step 4) and Test G with the merge (K2-3, K2-4: the registered scripts, unchanged).
  home/harddisks_resched2   7b08827, built the same way: the gate's E0 reference binary, as ~/harddisks_resched2.
  home/harddisks_gen3_k2    the FIX commit with the K2-1 sparse set (scripts only): runs the K2-2 job (the gate's (c) and (d) with
                            the binaries of home/harddisks_gen3), the fetch and the verdict, as on KOA and on the Mac.
  home_gate/harddisks_gen3  the FIX commit with KOA's step-1 sparse set and its own builds: runs the gate job (a)-(d) as committed.
MAC BUILDS: the Makefile's koa and koa-ld targets with CFLAGS_KOA="-O2 -ffp-contract=off $(INCLUDE)" (the koa flags without
-march=x86-64-v2 -mtune=generic, which arm64 does not take), so the version line reads "target koa" / "koa-ld" as on KOA; the M1/M2
harnesses with build_gen3_koa.sh's flags without -march/-mtune; the hash files written exactly as build_gen3_koa.sh writes them.
core.abbrev=7 in every clone (KOA's shallow clones abbreviate to 7 digits; the full local repository to 8).
SHIMS (every one is listed in the printed table of residual differences): the Slurm variables (SLURM_SUBMIT_DIR, SLURM_JOB_ID,
SLURM_ARRAY_JOB_ID, SLURM_ARRAY_TASK_ID) and the log file names Slurm would write (logs/%x_%j.out, %x_%A_%a.out); `module` (a no-op),
`gcc --version` (prints 14.3.0; nothing is compiled inside a job), `lscpu` (the Mac's CPU), `hostname` (a pretended node name for the
cross-node tasks); HOME; the Lustre path /mnt/lustre/koa/scratch/charing rewritten to <root>/lustre in the clones' job scripts; the
fetch's DTN source (ssh) made local; in home/harddisks_gen3 the registered Test G task list and its SHA-256 replaced by a MINI list (52
of the registered lines, in their order: main M = 50 r 0-13, M = 300 and 1500 r 0-1, info_m r 0-1, dense r 0-1, n400 r 0-1, afix the
first two seeds; three chunks of 25, every block type); in home/harddisks_gen3_k2 the expected SHA-256 of K2-2 replaced by the Mac
builds' (k2_binaries_a643bc7.sha256 holds KOA's).
ORDER (at most 8 simulation processes at any time): builds; the unit test of the fix; the gate job; K2-2 (and two STOP checks); the
cross-node array (4 tasks, 2 pretended nodes) and the compare; Test G stage 1 (chunk 1 alone), stage 2 (chunks 2-3) and a chunk
beyond the list (STOP check); the merge; the fetch; the verdict (A: as committed; B: its inventory's verdict forced to clean, so the rest
of the registered main() runs on the mini data -- a component check, information).
usage (from the repository root): python3 hspist3/cluster/gen3_koa_261009/rehearse_k2_mac.py --root <empty dir> --fix <engine-gen3 commit>
"""
import argparse, glob, hashlib, os, re, shutil, subprocess, sys, time

HERE = os.path.dirname(os.path.abspath(__file__)); HS = os.path.dirname(os.path.dirname(HERE)); REPO = os.path.dirname(HS)
GEN3_WT = os.path.join(os.path.dirname(REPO), "HardDisks_gen3")          # the engine-gen3 worktree (its kissfft clone is the source)
GATE_C, REF_C, KISS = "a643bc7edb1b41574c0f9af0c5fea199dc797680", "7b08827", "febd4caeed32e33ad8b2e0bb5ea77542c40f18ec"
SPARSE_KOA = ["/hspist3/*.c", "/hspist3/*.h", "/hspist3/*.py", "/hspist3/Makefile", "/hspist3/edmd_core/", "/hspist3/cluster/",
              "/hspist3/validation/", "/hspist3/kissfft"]                  # runsheet step 1
SPARSE_K2 = ["/hspist3/*.py", "/hspist3/cluster/", "/hspist3/validation/"]  # runsheet addendum K2-1
LUSTRE = "/mnt/lustre/koa/scratch/charing"
MACF = "-O2 -ffp-contract=off"
K = "cluster/gen3_koa_261009"
RES = []                                                                  # (script, run, criterion, result, evidence)


def sh(cmd, cwd=None, env=None, log=None, timeout=None, check=True):
    if log:
        with open(log, "w") as out:
            p = subprocess.run(cmd, cwd=cwd, env=env, stdout=out, stderr=subprocess.STDOUT, timeout=timeout)
    else:
        p = subprocess.run(cmd, cwd=cwd, env=env, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, timeout=timeout)
    if check and p.returncode != 0:
        sys.exit(f"STOP: {' '.join(cmd)} (in {cwd}) exit {p.returncode}" + (f"; see {log}" if log else f": {(p.stdout or '')[-800:]}"))
    return p


def sha(p): return hashlib.sha256(open(p, "rb").read()).hexdigest()


def clone(dest, commit, sparse):
    sh(["git", "clone", "-q", "--no-checkout", "--shared", REPO, dest])
    sh(["git", "-C", dest, "config", "core.abbrev", "7"])
    sh(["git", "-C", dest, "config", "core.sparseCheckout", "true"])
    with open(os.path.join(dest, ".git", "info", "sparse-checkout"), "w") as fh: fh.write("\n".join(sparse) + "\n")
    sh(["git", "-C", dest, "-c", "advice.detachedHead=false", "checkout", "-q", commit])
    return sh(["git", "-C", dest, "rev-parse", "HEAD"]).stdout.strip()


def kissfft(dest):
    sh(["git", "clone", "-q", os.path.join(GEN3_WT, "hspist3", "kissfft"), os.path.join(dest, "hspist3", "kissfft")])
    sh(["git", "-C", os.path.join(dest, "hspist3", "kissfft"), "-c", "advice.detachedHead=false", "checkout", "-q", KISS])


def build(hs, R, ld=True, harness=True):
    """the Mac equivalents of build_gen3_koa.sh's builds, the same hash files"""
    os.makedirs(os.path.join(hs, "logs"), exist_ok=True); os.makedirs(os.path.join(hs, "bin"), exist_ok=True)
    tag = os.path.relpath(os.path.dirname(hs), R).replace(os.sep, "_")
    mk = lambda t, cf: sh(["make", "-B", t, f"CFLAGS_KOA={MACF} $(INCLUDE)", f"BUILD_CF={cf}"], cwd=hs, log=os.path.join(R, "logs", f"make_{t}_{tag}.txt"))
    if ld:
        mk("koa-ld", f"{MACF} -DEDMD3_LONG_DOUBLE [mac rehearsal]")
        shutil.move(os.path.join(hs, "00ALLINONE"), os.path.join(hs, "bin", "00ALLINONE_ld"))
        open(os.path.join(hs, "logs", "BUILD_KOA_LD.hash"), "w").write(f"{sha(os.path.join(hs, 'bin', '00ALLINONE_ld'))}  00ALLINONE_ld\n")
    mk("koa", f"{MACF} [mac rehearsal]")
    open(os.path.join(hs, "logs", "BUILD_KOA_LAST.hash"), "w").write(f"{sha(os.path.join(hs, '00ALLINONE'))}  00ALLINONE\n")
    if harness:
        for n, src in (("gen3_m1", "gen3_m1_harness.c"), ("gen3_m2", "gen3_m2_harness.c")):
            sh(["cc", "-std=gnu11", "-O2", "-ffp-contract=off", "-o", f"bin/{n}", f"edmd_core/tests/{src}", "edmd_core/edmd_gen3.c", "edmd_core/edmd.c", "-lm"], cwd=hs)
        open(os.path.join(hs, "logs", "BUILD_KOA_HARNESS.hash"), "w").write("".join(f"{sha(os.path.join(hs, 'bin', n))}  {n}\n" for n in ("gen3_m1", "gen3_m2")))
    return [sh([os.path.join(hs, b), "--version"]).stdout.splitlines()[0] for b in (["00ALLINONE"] + (["bin/00ALLINONE_ld"] if ld else []))]


def shims(R):
    b = os.path.join(R, "shim", "bin"); os.makedirs(b)
    files = {
        "module": "#!/bin/bash\n# rehearsal shim: Lmod is not on the Mac; koa_env.sh's 'module purge' / 'module load compiler/GCC/14.3.0' do nothing\nexit 0\n",
        "gcc": "#!/bin/bash\n# rehearsal shim: koa_env.sh checks 'gcc --version' for 14.3; no job compiles\n"
               "if [ \"${1:-}\" = --version ]; then echo \"gcc (GCC) 14.3.0 [rehearsal shim; here: $(/usr/bin/cc --version | head -1)]\"; exit 0; fi\n"
               "echo 'rehearsal shim: gcc is not used inside a job' >&2; exit 97\n",
        "lscpu": "#!/bin/bash\n# rehearsal shim: Linux's lscpu\necho \"Model name:          $(sysctl -n machdep.cpu.brand_string) (Mac rehearsal)\"\n",
        "hostname": "#!/bin/bash\n# rehearsal shim: a pretended node name for the cross-node tasks (REHEARSAL_HOST), else the Mac's\n"
                    "if [ -n \"${REHEARSAL_HOST:-}\" ]; then echo \"$REHEARSAL_HOST\"; else exec /bin/hostname \"$@\"; fi\n",
        # ##CHRIS 2026-10-10 (sec. 4.7.47): added during the first run on 9d4490c, before its fetch step; written here since
        "rsync": "#!/bin/bash\n# rehearsal shim: the fetch's rsync with --exclude=e0_*/ added, so the gate folder's E0 outputs (A-fixed event logs,\n"
                 "# about 0.2 GB) are not copied a second time; every other file goes through /usr/bin/rsync with the fetch's own arguments\n"
                 "echo \"[rehearsal rsync shim] /usr/bin/rsync --exclude=e0_*/ $*\" >&2\nexec /usr/bin/rsync \"--exclude=e0_*/\" \"$@\"\n",
    }
    for n, t in files.items():
        p = os.path.join(b, n); open(p, "w").write(t); os.chmod(p, 0o755)
    return b


def rewrite_lustre(clone_dir, R, local_fetch=False):
    done = []
    for p in sorted(glob.glob(os.path.join(clone_dir, "hspist3", "cluster", "**", "*"), recursive=True)):
        if not os.path.isfile(p) or not p.endswith((".sbatch", ".sh", ".py")): continue
        s = open(p, errors="surrogateescape").read()
        if LUSTRE not in s: continue
        t = s.replace(LUSTRE, os.path.join(R, "lustre"))
        if local_fetch: t = t.replace("charing@koa-dtn.its.hawaii.edu:", "")
        open(p, "w", errors="surrogateescape").write(t); done.append(os.path.relpath(p, clone_dir))
    return done


def env_for(R, home, shim, **slurm):
    e = {k: v for k, v in os.environ.items() if not k.startswith("SLURM_")}
    e.update(HOME=home, PATH=f"{shim}:{os.environ['PATH']}", **{k: str(v) for k, v in slurm.items()})
    e.pop("REHEARSAL_HOST", None)
    return e


def job(name, script, args, hs, env, log, timeout=7200):
    t0 = time.time()
    p = sh(["/bin/bash", script] + list(args), cwd=hs, env=env, log=log, timeout=timeout, check=False)
    txt = open(log, errors="ignore").read()
    print(f"  {name}: exit {p.returncode}, {time.time() - t0:.0f} s, log {os.path.relpath(log, os.path.dirname(os.path.dirname(hs)))}", flush=True)
    return p.returncode, txt


def res(script, run, crit, ok, ev):
    RES.append((script, run, crit, "PASS" if ok else "**FAIL**", ev)); print(f"    -> {script}: {'PASS' if ok else '**FAIL**'} ({ev})", flush=True)


def table_rows(txt, head):
    i = txt.find(head)
    if i < 0: return []
    rows = []
    for l in txt[i:].splitlines()[2:]:
        if not l.startswith("| "): break
        rows.append([c.strip() for c in l.strip().strip("|").split("|")])
    return rows


def mini_list(src):
    out, af = [], []
    for l in open(src):
        f = l.split()
        if f[0] == "B":
            blk, M, r = f[11], int(f[2]), int(f[3])
            keep = (blk == "main" and M == 50 and r < 14) or (blk == "main" and M in (300, 1500) and r < 2) or (blk in ("info_m", "dense", "n400") and r < 2)
        else:
            if f[3] not in af and len(af) < 2: af.append(f[3])
            keep = f[3] in af
        if keep: out.append(l)
    return out


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--root", required=True); ap.add_argument("--fix", required=True)
    ap.add_argument("--only", default="unit,gate,k2,xnode,testg", help="parts to run (a later run after a fix: the parts whose scripts changed)")
    a = ap.parse_args(); PARTS = set(a.only.split(","))
    R = os.path.abspath(a.root)
    if not os.path.isdir(R) or os.listdir(R): sys.exit(f"STOP: {R} must be an existing empty folder")
    FIX = sh(["git", "-C", REPO, "rev-parse", a.fix]).stdout.strip()
    for d in ("logs", "home", "home_gate", "lustre"): os.makedirs(os.path.join(R, d))
    t00 = time.time()
    print(f"# Mac rehearsal of the KOA job scripts {time.strftime('%Y-%m-%d %H:%M:%S %Z')}: root {R}; fix commit {FIX}; gate commit {GATE_C}; parts {sorted(PARTS)}")
    print(f"tools: {sh(['/bin/bash', '--version']).stdout.splitlines()[0]} | {sh([sys.executable, '--version']).stdout.strip()} | "
          f"{sh(['cc', '--version']).stdout.splitlines()[0]} | {sh(['git', '--version']).stdout.strip()} | {(sh(['/usr/bin/rsync', '--version'], check=False).stdout or '?').splitlines()[0]}")
    shim = shims(R)
    # ---------------------------------------------------------------- clones and builds
    H, HG = os.path.join(R, "home"), os.path.join(R, "home_gate")
    A, REF, K2C, G = (os.path.join(H, "harddisks_gen3"), os.path.join(H, "harddisks_resched2"), os.path.join(H, "harddisks_gen3_k2"),
                      os.path.join(HG, "harddisks_gen3"))
    print("== clones and Mac builds")
    CL = [(A, GATE_C, SPARSE_KOA, True), (K2C, FIX, SPARSE_K2, False)] + ([(REF, REF_C, SPARSE_KOA, True), (G, FIX, SPARSE_KOA, True)] if "gate" in PARTS else [])
    for d, c, sp, kf in CL:
        head = clone(d, c, sp)
        if kf: kissfft(d)
        print(f"  {os.path.relpath(d, R)}: HEAD {head}; files {sh(['git', '-C', d, 'ls-files', '-t']).stdout.count(chr(10))} tracked, "
              f"{len([x for x in sh(['git', '-C', d, 'ls-files', '-v']).stdout.splitlines() if x.startswith('H ')])} checked out")
    vA = build(os.path.join(A, "hspist3"), R); vR = vG = "(not built: no gate part)"
    if "gate" in PARTS:
        os.symlink(REF, os.path.join(HG, "harddisks_resched2"))
        vR = build(os.path.join(REF, "hspist3"), R, ld=False, harness=False); vG = build(os.path.join(G, "hspist3"), R)
    print(f"  versions: gen3 a643bc7 {vA}; reference {vR}; gate clone (fix commit) {vG}")
    os.makedirs(os.path.join(K2C, "hspist3", "logs"))
    for d in [x[0] for x in CL if x[0] != K2C]:
        os.makedirs(os.path.join(d, "hspist3", "logs"), exist_ok=True)
    # ---------------------------------------------------------------- shims in the clones (after the builds: the version lines stay clean)
    print("== shims in the clones")
    for d, lf in [(A, False), (K2C, True)] + ([(G, False)] if "gate" in PARTS else []):
        print(f"  Lustre path rewritten in {os.path.relpath(d, R)}: {', '.join(rewrite_lustre(d, R, lf))}")
    reg = os.path.join(A, "hspist3", K, "tasks_testG_koa_261009.txt"); assert sha(reg) == "cff8e6f9ef3d95440f54b4c621fff776b64e48700083e7258e472351fc6412b9"
    ml = mini_list(reg); txt = "".join(ml)
    open(reg, "w").write(txt); open(os.path.join(A, "hspist3", K, "tasks_testG_koa_261009.sha256"), "w").write(hashlib.sha256(txt.encode()).hexdigest() + "\n")
    from collections import Counter
    print(f"  mini Test G list in home/harddisks_gen3: {len(ml)} of the registered lines ({len(ml) // 25 + (len(ml) % 25 > 0)} chunks of 25); per block "
          + ", ".join(f"{k} {v}" for k, v in Counter(l.split()[-1] for l in ml).items()))
    kfile = os.path.join(K2C, "hspist3", K, "k2_binaries_a643bc7.sha256")
    koa_k2 = open(kfile).read()
    open(kfile, "w").write("".join(f"{sha(os.path.join(A, 'hspist3', p))}  {p}\n" for p in ("00ALLINONE", "bin/00ALLINONE_ld", "bin/gen3_m1", "bin/gen3_m2")))
    print(f"  k2_binaries_a643bc7.sha256 in the K2 clone: the Mac builds' SHA-256 (committed: KOA's, {koa_k2.split()[0][:8]}...)")
    # ---------------------------------------------------------------- 0. the unit test of the fix
    os.makedirs(os.path.join(R, "tmp"))
    if "unit" in PARTS: unit(R, K2C)
    if "gate" in PARTS: gate(R, G, HG, shim)
    if "k2" in PARTS: k2(R, K2C, H, HG, shim)
    if "xnode" in PARTS: xnode(R, A, H, shim)
    if "testg" in PARTS: testg(R, A, K2C, H, shim, PARTS)
    table(R, t00)


def unit(R, K2C):
    print("== the unit test of the fix (task 1)")
    p = sh([sys.executable, os.path.join(K2C, "hspist3", K, "test_g3_rate_ld_record_261010.py")], env=dict(os.environ, TMPDIR=os.path.join(R, "tmp")),
           log=os.path.join(R, "logs", "unit_test.txt"), check=False)
    t = open(os.path.join(R, "logs", "unit_test.txt")).read()
    res("test_g3_rate_ld_record_261010.py", "Mac python3", "exit 0, ALL PASS (old pattern fails as on KOA, record() parses)", p.returncode == 0 and "ALL PASS" in t,
        f"{t.count('| PASS |')} cases PASS")


def gate(R, G, HG, shim):
    print("== 1. the gate job g3_gate_koa.sbatch (a)-(d), clone at the fix commit")
    hs = os.path.join(G, "hspist3"); JG = 910001
    rc, t = job("gate", f"{K}/g3_gate_koa.sbatch", [], hs, env_for(R, HG, shim, SLURM_SUBMIT_DIR=hs, SLURM_JOB_ID=JG), os.path.join(hs, "logs", f"gen3-gate_{JG}.out"))
    e0 = re.findall(r"^ (ctrl_min|ctrl_leg|afix|afix_leg)  audit vs plain:  (\S+)  plain vs ref:  (\S+)", t, re.M)
    e2 = re.findall(r"^E2 \(amended\): (.*)$", t, re.M)
    rc_c = table_rows(t, "| N | eta | events | run [s]"); rc_d = table_rows(t, "| N | M | seed | build | clean")
    ok_a = len(e0) == 8 and all(x[1] == x[2] == "IDENTICAL" for x in e0)
    ok_e2 = len(e2) == 2 and all(x.startswith("NOT TESTED by this job") and x.endswith("PASS") for x in e2)
    ok_b = "(m1 exit 0;m2 quick exit 0;m2 exit 0;)" in t
    ok_c = len(rc_c) == 12 and all(r[-1] == "1" for r in rc_c)
    ok_d = len(rc_d) == 16 and all(r[4] == "1" and r[7] not in ("nan", "") for r in rc_d)
    res("g3_gate_koa.sbatch", f"job {JG}", "exit 0; (a) 8 x IDENTICAL/IDENTICAL; E2 labelled NOT TESTED, E0 PASS; (b) 3 exits 0; (c) 12 rows clean 1; (d) 16 rows clean 1; done line",
        rc == 0 and ok_a and ok_e2 and ok_b and ok_c and ok_d and "gen3 gate (KOA, one node) done" in t,
        f"exit {rc}; E0 {sum(x[1] == x[2] == 'IDENTICAL' for x in e0)}/8; E2 lines {e2[:1]}; (b) {ok_b}; (c) {len(rc_c)} rows; (d) {len(rc_d)} rows")
    res("audit_runs_261007.py (run, report)", f"inside job {JG}", "E0 PASS twice", ok_a and ok_e2, f"{len(e0)} E0 rows")
    res("g3_rate_ld_koa.py", f"inside job {JG}", "(c) and (d) tables complete", ok_c and ok_d, f"(c) {len(rc_c)}, (d) {len(rc_d)}")


def k2(R, K2C, H, HG, shim):
    print("== 2. K2-2: g3_rate_ld_k2.sbatch from the K2 clone, the binaries of home/harddisks_gen3")
    hs = os.path.join(K2C, "hspist3"); JK = 910002
    rc, t = job("K2-2", f"{K}/g3_rate_ld_k2.sbatch", [], hs, env_for(R, H, shim, SLURM_SUBMIT_DIR=hs, SLURM_JOB_ID=JK), os.path.join(hs, "logs", f"gen3-rate-ld-k2_{JK}.out"))
    rc_c = table_rows(t, "| N | eta | events | run [s]"); rc_d = table_rows(t, "| N | M | seed | build | clean")
    res("g3_rate_ld_k2.sbatch", f"job {JK}", "exit 0; binaries' SHA-256 and version lines accepted; (c) 12 rows clean 1; (d) 16 rows clean 1; done line",
        rc == 0 and "STOP" not in t and len(rc_c) == 12 and all(r[-1] == "1" for r in rc_c) and len(rc_d) == 16 and all(r[4] == "1" for r in rc_d)
        and "K2-2 (gate parts (c) and (d)) done" in t, f"exit {rc}; (c) {len(rc_c)} rows; (d) {len(rc_d)} rows")
    os.symlink(K2C, os.path.join(H, "harddisks_gen3_wrongname"))
    hw = os.path.join(H, "harddisks_gen3_wrongname", "hspist3")
    rc1, t1 = job("K2-2 STOP check: another clone name", f"{K}/g3_rate_ld_k2.sbatch", [], hw, env_for(R, H, shim, SLURM_SUBMIT_DIR=hw, SLURM_JOB_ID=910012),
                  os.path.join(R, "logs", "k2_stop_name.out"))
    rc2, t2 = job("K2-2 STOP check: other binaries in ~/harddisks_gen3", f"{K}/g3_rate_ld_k2.sbatch", [], hs, env_for(R, HG, shim, SLURM_SUBMIT_DIR=hs, SLURM_JOB_ID=910013),
                  os.path.join(R, "logs", "k2_stop_sha.out"))
    res("g3_rate_ld_k2.sbatch (STOP rules)", "jobs 910012, 910013", "exit 2 with the STOP line, before any run: wrong clone name; binaries not the recorded ones",
        rc1 == 2 and "STOP: the K2 jobs run from" in t1 and rc2 == 2 and "STOP: the binaries of ~/harddisks_gen3 are not" in t2,
        f"exits {rc1}, {rc2}")


def xnode(R, A, H, shim):
    print("== 3. the cross-node array g3_det_xnode.sbatch (4 tasks, 2 pretended nodes) and g3_xnode_compare.sbatch")
    hs = os.path.join(A, "hspist3"); JX = 910003
    procs = []
    for i, host in ((1, "mac-pretend-a"), (2, "mac-pretend-a"), (3, "mac-pretend-b"), (4, "mac-pretend-b")):
        e = env_for(R, H, shim, SLURM_SUBMIT_DIR=hs, SLURM_ARRAY_JOB_ID=JX, SLURM_ARRAY_TASK_ID=i, SLURM_JOB_ID=JX + 100 + i); e["REHEARSAL_HOST"] = host
        lg = os.path.join(hs, "logs", f"gen3-xnode_{JX}_{i}.out")
        procs.append((i, subprocess.Popen(["/bin/bash", f"{K}/g3_det_xnode.sbatch"], cwd=hs, env=e, stdout=open(lg, "w"), stderr=subprocess.STDOUT), lg))
    oks = []
    for i, p, lg in procs:
        rc = p.wait(); t = open(lg).read()
        cells = re.findall(r"^== (\S+): exit 0; clean=1 hash=([0-9a-f]{16})", t, re.M)
        oks.append(rc == 0 and len(cells) == 2 and f"task {i} done" in t); print(f"  task {i}: exit {rc}; {cells}")
    res("g3_det_xnode.sbatch", f"array {JX}, tasks 1-4", "every task exit 0, both cells exit 0 and clean=1, done line", all(oks), f"{sum(oks)}/4 tasks")
    JC = 910004
    rc, t = job("compare", f"{K}/g3_xnode_compare.sbatch", [str(JX)], hs, env_for(R, H, shim, SLURM_SUBMIT_DIR=hs, SLURM_JOB_ID=JC), os.path.join(hs, "logs", f"gen3-xnode-compare_{JC}.out"))
    res("g3_xnode_compare.sbatch / .py", f"job {JC}", "RESULT: across nodes IDENTICAL (2 nodes); within nodes IDENTICAL; tasks used 4 of 4",
        rc == 0 and "RESULT: across nodes IDENTICAL (2 node(s)); within nodes IDENTICAL; tasks used 4 of 4" in t, re.findall(r"^RESULT: .*$", t, re.M)[:1])


def testg(R, A, K2C, H, shim, PARTS):
    hs = os.path.join(A, "hspist3")
    print("== 4. Test G on the mini list: stage 1 (chunk 1 alone), stage 2 (chunks 2-3), a chunk beyond the list, the merge (K2: from the K2 clone)")
    def chunk(jid, j, lg):
        e = env_for(R, H, shim, SLURM_SUBMIT_DIR=hs, SLURM_ARRAY_JOB_ID=jid, SLURM_ARRAY_TASK_ID=j, SLURM_JOB_ID=jid + 100 + j)
        return subprocess.Popen(["/bin/bash", f"{K}/testG_koa.sbatch"], cwd=hs, env=e, stdout=open(lg, "w"), stderr=subprocess.STDOUT)
    J1, J2 = 910005, 910006
    lg1 = os.path.join(hs, "logs", f"gen3-testG_{J1}_1.out"); t0 = time.time(); rc = chunk(J1, 1, lg1).wait(); t = open(lg1).read()
    print(f"  stage 1, chunk 1: exit {rc}, {time.time() - t0:.0f} s")
    ok1 = rc == 0 and "chunk 1 done" in t and "FAILED" not in t
    d1 = os.path.join(R, "lustre", "harddisks_gen3", "testG_koa", "chunks", "chunk_1", "testG")
    n2 = len(glob.glob(os.path.join(d1, "gen2", "*", "m_50", "wall_x_positions_*"))); n3 = len(glob.glob(os.path.join(d1, "gen3", "*", "m_50", "wall_x_positions_*")))
    rl3 = open(glob.glob(os.path.join(d1, "gen3", "*", "m_50", "run.log"))[0]).read() if n3 else ""
    c3 = len(re.findall(r"^\[EDMD3-HEALTH\] .*: clean=1 ", rl3, re.M))
    res("testG_koa.sbatch + testG_worker_koa.sh (stage 1)", f"array {J1}, task 1", "exit 0; 'chunk 1 done'; no FAILED; 13 gen2 + 12 gen3 traces at M = 50; 12 gen3 records clean=1",
        ok1 and n2 == 13 and n3 == 12 and c3 == 12, f"traces gen2 {n2}, gen3 {n3}; gen3 clean records {c3}")
    lgs = [os.path.join(hs, "logs", f"gen3-testG_{J2}_{j}.out") for j in (2, 3)]; t0 = time.time()
    ps = [chunk(J2, j, lg) for j, lg in zip((2, 3), lgs)]; rcs = [p.wait() for p in ps]; ts = [open(lg).read() for lg in lgs]
    print(f"  stage 2, chunks 2-3: exits {rcs}, {time.time() - t0:.0f} s")
    res("testG_koa.sbatch + testG_worker_koa.sh (stage 2)", f"array {J2}, tasks 2-3", "exit 0, 'chunk j done', no FAILED in both",
        rcs == [0, 0] and all(f"chunk {j} done" in x and "FAILED" not in x for j, x in zip((2, 3), ts)), f"exits {rcs}")
    lg4 = os.path.join(R, "logs", "testG_beyond.out"); rc4 = chunk(910008, 4, lg4).wait(); t4 = open(lg4).read()
    res("testG_koa.sbatch (STOP rule)", "array 910008, task 4", "exit 2: 'STOP: chunk 4 is beyond the task list'", rc4 == 2 and "STOP: chunk 4 is beyond the task list" in t4, f"exit {rc4}")
    JM = 910007; hk = os.path.join(K2C, "hspist3")
    rc, t = job("merge (K2)", f"{K}/merge_testG_k2.sbatch", [], hk, env_for(R, H, shim, SLURM_SUBMIT_DIR=hk, SLURM_JOB_ID=JM), os.path.join(hk, "logs", f"gen3-testG-merge-k2_{JM}.out"))
    red = re.findall(r"^reduce_B (\S+): exit (\d+)", t, re.M)
    res("merge_testG_k2.sbatch + merge_testG_koa.py (fix commit; + reduce_B.py)", f"job {JM}", "exit 0; 'merged 3 chunks'; one build; reduce_B exit 0 for every B cell; traces compressed; done line",
        rc == 0 and "merged 3 chunks" in t and red and all(x[1] == "0" for x in red) and "merge done" in t and "STOP" not in t,
        f"exit {rc}; {re.findall(r'^merged .*$', t, re.M)[:1]}; reduce_B {red}")
    # ---------------------------------------------------------------- 5. the fetch (Mac) and 6. the verdict (Mac), in the K2 clone
    print("== 5. the fetch (local rsync source) and 6. the registered verdict, in the K2 clone")
    hs = os.path.join(K2C, "hspist3")
    rc, t = job("fetch", f"{K}/fetch_gen3_koa.sh", ["910001", "910003", "910002"], hs, env_for(R, H, shim), os.path.join(R, "logs", "fetch.out"))
    F = os.path.join(hs, "experiments_gen3_gate_koa_261009")
    need = ("testG", "gate", "xnode", "rate_ld_k2") if {"gate", "xnode", "k2"} <= PARTS else ("testG",)
    have = {s: os.path.isdir(os.path.join(F, s)) for s in need}
    res("fetch_gen3_koa.sh (Mac)", "gate 910001, xnode 910003, K2-2 910002", f"exit 0; {', '.join(need)} fetched", rc == 0 and all(have.values()), f"{have}")
    lgA = os.path.join(R, "logs", "verdict_A.out")
    p = sh([sys.executable, f"{K}/testG_koa_verdict_261009.py"], cwd=hs, log=lgA, check=False); tA = open(lgA).read()
    inv = table_rows(tA, "| engine | M | trajectories matched")
    want = {("gen3", "50"): 14, ("gen2", "50"): 14, ("gen3", "300"): 2, ("gen2", "300"): 2, ("gen3", "1500"): 2, ("gen2", "1500"): 2}
    got = {(r[0], r[1]): int(r[2].split()[0]) for r in inv if len(r) > 2 and r[2].split()[0].isdigit()}
    res("testG_koa_verdict_261009.py (A, as committed)", "Mac python3", "exit 0; the inventory reads the mini data (matched 14/2/2 per engine of 400 expected) and ends INVALID",
        p.returncode == 0 and got == want and "TEST G (KOA): INVALID" in tA, f"matched {got}")
    lgB = os.path.join(R, "logs", "verdict_B.out")
    code = ("import sys; sys.argv = ['testG_koa_verdict_261009.py']; sys.path.insert(0, 'cluster/gen3_koa_261009'); import testG_koa_verdict_261009 as KV\n"
            "orig = KV.V.inventory\n"
            "def inv(D):\n    ok = orig(D)\n    print('[rehearsal B] inventory: ' + ('clean' if ok else 'NOT CLEAN') + ' -- forced clean, so the rest of main() runs on the mini data')\n    return True\n"
            "KV.V.inventory = inv\nKV.main()\n")
    p = sh([sys.executable, "-c", code], cwd=hs, log=lgB, check=False); tB = open(lgB).read()
    rule = table_rows(tB, "| mass | judged on (per engine)")
    res("testG_koa_verdict_261009.py (B, component check)", "Mac python3", "exit 0; the rule table (3 masses), the verdict line and the information rows print",
        p.returncode == 0 and len(rule) >= 3 and "TEST G (KOA):" in tB and "Information (no verdict): the other rows" in tB, f"rule rows {len(rule)}; {re.findall(r'^TEST G .*$', tB, re.M)[:1]}")


def table(R, t00):
    print(f"\n## Rehearsal result ({time.time() - t00:.0f} s; data {sh(['du', '-sh', R]).stdout.split()[0]})\n")
    print("| script | run | criterion | result | evidence |\n|---|---|---|---|---|")
    for r in RES: print("| " + " | ".join(str(x).replace("|", "/") for x in r) + " |")
    n_fail = sum(r[3] != "PASS" for r in RES)
    print(f"\nREHEARSAL: {'PASS' if n_fail == 0 else f'FAIL ({n_fail} of {len(RES)})'}")
    sys.exit(0 if n_fail == 0 else 1)


if __name__ == "__main__":
    main()
