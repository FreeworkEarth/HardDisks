#!/usr/bin/env python3
"""##CHRIS 2026-10-11 (plan-author decision 28, task R0; 261012 sec. 4.7.58): THE FULL MAC REHEARSAL of Test R on KOA -- every script the
KOA jobs run, from the runsheet's paste blocks as Chris types them, before anything is typed on KOA. Nothing outside --root is written;
the clone borrows the repository's objects (git clone --shared). At most one simulation process at a time (C1 runs beside it: <= 8).
THE HOMES (as on KOA today):
  home/harddisks_gen3        only what the Test R jobs read there: hspist3/00ALLINONE and hspist3/logs/BUILD_KOA_LAST.hash. The binary is a
                             MAC build of engine-gen3 a643bc7 with the koa target's name ("00ALLINONE  git a643bc7  target koa"): the
                             driver's sources from git archive (as experiments_gen3_c1_261011/build_mac_a643bc7.sh), kissfft from the
                             worktree's checkout, -O2 -ffp-contract=off (the koa flags without -march/-mtune, which arm64 does not take);
                             the hash file written as build_gen3_koa.sh writes it.
  home/harddisks_gen3_k2     engine-gen3 fd8d8a2 (KOA's ~/harddisks_gen3_k2 since K2-4 block C) with K2-1's sparse set; its logs/ holds a
                             pretended Test G merge log and a Test G chunk log (as KOA's does). Block B UPDATES IT IN PLACE, as written.
THE BLOCKS (--runsheet, default runsheet_testR_261011.txt; "TEST R BLOCK <X> (" to the first empty line), each run as one bash with
HOME = the rehearsal home and the Lustre path /mnt/lustre/koa/scratch/charing rewritten to <root>/lustre: A (squeue, srun: shims),
B (in a new shell: the update in place, the SHA lines, the manifests, the binaries, the unit tests), then the REHEARSAL SUBSTITUTIONS,
C (chunk 1 alone), D (its checks), E (the rest; its first attempt stops in a start-up check: the sha256sum shim fails), F (the checks
and the resubmission line), the printed resubmission line typed as printed, F again, the merge job's own STOP rules (a pretended
pending gen3-testR job; a chunk's done log held aside), G (the merge), H (its checks), a second G (exit 3), the fetch (local), and the
verdict on the mini data (a component check) with its one-printing guard.
THE RESIDUAL DIFFERENCES (every one printed in the table):
  - shims: the Slurm commands (squeue: empty or REHEARSAL_SQUEUE's lines; sacct, srun: not run; sbatch: runs the job script at once
    with a pretended job id and its #SBATCH output name, --array lists run one after the other, --exclude recorded; sinfo: a pretended
    line), module (no-op), gcc --version (prints 14.3.0), lscpu (the Mac's CPU), hostname (a pretended node name REHEARSAL_HOST, else the
    Mac's), sha256sum (the system's; fails when REHEARSAL_STARTUP_FAIL is set: a pretended faulty node);
  - the Mac binary (above): its SHA-256 replaces KOA's in the delivered build record (after block B), and block B's check of KOA's four
    binaries fails here (they are not on the Mac);
  - after block B, in the updated clone: the task list replaced by a MINI LIST, the first 12 registered lines in their order (r = 0
    and 1 at M = 50, 100 and 300, both engines), its SHA-256 file, the chunk size 6 (two chunks: every cell folder receives from both),
    the build record (the Mac binary's SHA-256), the two manifests rewritten by make_testR_manifests.py; the Lustre path rewritten in the
    package's job scripts; the fetch's DTN source made local;
  - block E's --array=2-216%64 becomes --array=2 (the mini list has 2 chunks); block F's counts are the mini list's;
  - the fetch runs from the updated clone, into it (on KOA's data it runs from the Mac's repository); the verdict runs from the Mac's
    repository (as it will) on that fetched folder and the mini list: its test hooks HD_TESTR_LOC, HD_TESTR_TASKS, HD_TESTR_TASKS_SHA.
usage (from the repository root): python3 hspist3/cluster/gen3_testR_261011/rehearse_testR_mac.py --root <empty dir> --pkg <engine-gen3 commit>
"""
import argparse, hashlib, os, re, shutil, subprocess, sys, time
HERE = os.path.dirname(os.path.abspath(__file__)); HS = os.path.dirname(os.path.dirname(HERE)); REPO = os.path.dirname(HS)
GEN3_WT = os.path.join(os.path.dirname(REPO), "HardDisks_gen3")
A643, K2_NOW = "a643bc7edb1b41574c0f9af0c5fea199dc797680", "fd8d8a28037678f7cb3cf7c891ca52ecc60e67a9"
SPARSE_K2 = ["/hspist3/*.py", "/hspist3/cluster/", "/hspist3/validation/"]
LUSTRE = "/mnt/lustre/koa/scratch/charing"
P = "cluster/gen3_testR_261011"
MINI_N, MINI_K = 12, 6
RES = []


def sh(cmd, cwd=None, env=None, log=None, check=True, timeout=None):
    if log:
        with open(log, "w") as out: p = subprocess.run(cmd, cwd=cwd, env=env, stdout=out, stderr=subprocess.STDOUT, timeout=timeout)
    else: p = subprocess.run(cmd, cwd=cwd, env=env, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, timeout=timeout)
    if check and p.returncode != 0:
        sys.exit(f"STOP: {' '.join(cmd)} (in {cwd}) exit {p.returncode}" + (f"; see {log}" if log else f": {(p.stdout or '')[-800:]}"))
    return p


def sha(p): return hashlib.sha256(open(p, "rb").read()).hexdigest()


def L(out):
    """the output's lines, whitespace collapsed (macOS wc -l pads its count; KOA's does not)"""
    return [re.sub(r"[ \t]+", " ", l).strip() for l in out.splitlines()]


def check_D(out):
    """block D's good lines on the mini list. ##CHRIS 2026-10-11 (run 2): the expected lines are whitespace-collapsed like the output
    (run 2 compared collapsed output lines with expected lines that kept the build line's double spaces: 0/6, a fault of this check)"""
    want = [re.sub(r"[ \t]+", " ", f"{e} M={M} traces 1 runs 1 gen3-records-clean {int(e == 'gen3')} build '00ALLINONE  git a643bc7  target koa'")
            for e in ("gen2", "gen3") for M in (50, 100, 300)]
    n = sum(w in L(out) for w in want)
    return (n == 6 and "chunk 1 done" in out and "lines 1-6; 00ALLINONE  git a643bc7  target koa" in out
            and not re.search(r"FAILED|REFUSED|^STOP", out, re.M)), n


def res(script, crit, ok, ev):
    RES.append((script, crit, "PASS" if ok else "**FAIL**", ev)); print(f"    -> {script}: {'PASS' if ok else '**FAIL**'} ({ev})", flush=True)


SBATCH_SHIM = r"""#!/bin/bash
# rehearsal shim: runs the job script AT ONCE, as Slurm would run it: SLURM_SUBMIT_DIR = this folder, a pretended job id, the output file
# of its #SBATCH --output line (%x the job name, %j the id, %A %a the array id and task). --array=<list> runs its tasks one after the other
# (the throttle %N dropped); --exclude=<nodes> is recorded (the Mac is one node); any other option is refused
ARR=""; EXC=""
while [ $# -gt 1 ]; do
  case "$1" in --array=*) ARR=${1#--array=}; ARR=${ARR%%\%*};; --exclude=*) EXC=${1#--exclude=};; *) echo "rehearsal sbatch shim: option not supported: $1" >&2; exit 97;; esac
  shift
done
S=$1; C="$REHEARSAL_ROOT/logs/sbatch_shim_counter"
id=$(( $(cat "$C" 2>/dev/null || echo 930000) + 1 )); echo "$id" > "$C"
name=$(sed -n 's/^#SBATCH --job-name=//p' "$S" | head -1); out=$(sed -n 's/^#SBATCH --output=//p' "$S" | head -1)
echo "Submitted batch job $id"
[ -z "$EXC" ] || echo "$id exclude $EXC (the Mac is one node: recorded only) 0" >> "$REHEARSAL_ROOT/logs/sbatch_shim_exits.txt"
if [ -z "$ARR" ]; then
  o=${out//%x/$name}; o=${o//%j/$id}
  SLURM_SUBMIT_DIR=$PWD SLURM_JOB_ID=$id SLURM_JOB_NAME=$name bash "$S" > "$o" 2>&1; rc=$?
  echo "$id $PWD/$S $rc" >> "$REHEARSAL_ROOT/logs/sbatch_shim_exits.txt"
else
  for t in $(echo "$ARR" | tr ',' '\n' | awk -F- '{ if (NF == 2) for (i = $1; i <= $2; i++) print i; else print $1 }'); do
    o=${out//%x/$name}; o=${o//%A/$id}; o=${o//%a/$t}
    SLURM_SUBMIT_DIR=$PWD SLURM_JOB_ID=$((id * 1000 + t)) SLURM_ARRAY_JOB_ID=$id SLURM_ARRAY_TASK_ID=$t SLURM_JOB_NAME=$name bash "$S" > "$o" 2>&1; rc=$?
    echo "${id}_$t $PWD/$S $rc" >> "$REHEARSAL_ROOT/logs/sbatch_shim_exits.txt"
  done
fi
exit 0
"""


def shims(R):
    b = os.path.join(R, "shim", "bin"); os.makedirs(b)
    files = {
        "module": "#!/bin/bash\n# rehearsal shim: Lmod is not on the Mac\nexit 0\n",
        "gcc": "#!/bin/bash\n# rehearsal shim: koa_env.sh checks 'gcc --version' for 14.3; no job compiles\n"
               "if [ \"${1:-}\" = --version ]; then echo \"gcc (GCC) 14.3.0 [rehearsal shim; here: $(/usr/bin/cc --version | head -1)]\"; exit 0; fi\n"
               "echo 'rehearsal shim: gcc is not used inside a job' >&2; exit 97\n",
        "lscpu": "#!/bin/bash\n# rehearsal shim: Linux's lscpu\necho \"Model name:          $(sysctl -n machdep.cpu.brand_string) (Mac rehearsal)\"\n",
        "hostname": "#!/bin/bash\n# rehearsal shim: a pretended node name (REHEARSAL_HOST), else the Mac's\n"
                    "if [ -n \"${REHEARSAL_HOST:-}\" ]; then echo \"$REHEARSAL_HOST\"; else exec /bin/hostname -s; fi\n",
        "squeue": "#!/bin/bash\n# rehearsal shim: Slurm's queue -- empty, or the lines of the file REHEARSAL_SQUEUE (a pretended pending job)\n"
                  "[ -n \"${REHEARSAL_SQUEUE:-}\" ] && cat \"$REHEARSAL_SQUEUE\"\nexit 0\n",
        "sacct": "#!/bin/bash\n# rehearsal shim: no accounting on the Mac\nprintf '%20s %10s %10s %8s\\n' JobID State Elapsed ExitCode\n",
        "srun": "#!/bin/bash\n# rehearsal shim: the sandbox shell is not opened (the next block runs in a new shell)\necho \"rehearsal shim: srun $* (not run)\"\nexit 0\n",
        "sinfo": "#!/bin/bash\n# rehearsal shim: a pretended summary line of the shared partition\n"
                 "echo 'PARTITION AVAIL  TIMELIMIT   NODES(A/I/O/T) NODELIST'; echo 'shared       up 3-00:00:00        70/12/8/90 cn-[...] (rehearsal shim)'\n",
        "sbatch": SBATCH_SHIM,
        "sha256sum": "#!/bin/bash\n# rehearsal shim: the system sha256sum; with REHEARSAL_STARTUP_FAIL set it fails (exit 1), as on a faulty node\n"
                     "[ -n \"${REHEARSAL_STARTUP_FAIL:-}\" ] && { echo 'rehearsal shim: sha256sum fails (a pretended faulty node)' >&2; exit 1; }\n"
                     f"exec {shutil.which('sha256sum')} \"$@\"\n",
    }
    for n, t in files.items():
        p = os.path.join(b, n); open(p, "w").write(t); os.chmod(p, 0o755)
    return b


def env_for(R, home, shim, **extra):
    e = {k: v for k, v in os.environ.items() if not k.startswith(("SLURM_", "REHEARSAL_"))}
    e.update(HOME=home, PATH=f"{shim}:{os.environ['PATH']}", REHEARSAL_ROOT=R, **{k: str(v) for k, v in extra.items()})
    return e


def build_binary(R, home):
    """the Mac build of a643bc7 with the koa target's name, and the hash file build_gen3_koa.sh writes"""
    src = os.path.join(R, "build_src"); os.makedirs(src)
    a = subprocess.run(["git", "-C", GEN3_WT, "archive", A643, "hspist3/00ALLINONE.c", "hspist3/edmd_core", "hspist3/experiment_validation.c",
                        "hspist3/experiment_validation.h"], capture_output=True, check=True).stdout
    subprocess.run(["tar", "-x", "-C", src], input=a, check=True)
    k = os.path.join(src, "hspist3", "kissfft"); os.makedirs(k)
    for f in ("kiss_fft.c", "kiss_fftr.c", "kiss_fft.h", "kiss_fftr.h", "_kiss_fft_guts.h", "kiss_fft_log.h"): shutil.copy(os.path.join(GEN3_WT, "hspist3", "kissfft", f), k)
    hs = os.path.join(home, "harddisks_gen3", "hspist3"); os.makedirs(os.path.join(hs, "logs"))
    cf = sh(["pkg-config", "--cflags", "sdl2", "SDL2_ttf"]).stdout.split(); lf = sh(["pkg-config", "--libs", "sdl2", "SDL2_ttf"]).stdout.split()
    drv = ["00ALLINONE.c", "kissfft/kiss_fft.c", "kissfft/kiss_fftr.c", "edmd_core/edmd.c", "edmd_core/edmd_accelerated.c", "edmd_core/edmd_gen3.c", "experiment_validation.c"]
    sh(["cc", "-o", os.path.join(hs, "00ALLINONE")] + drv + ["-O2", "-ffp-contract=off"] + cf + ['-DBUILD_GIT="a643bc7"', '-DBUILD_TARGET="koa"',
       '-DBUILD_CFLAGS="-O2 -ffp-contract=off [mac rehearsal]"'] + lf + ["-L/opt/homebrew/lib", "-lGLEW", "-framework", "OpenGL", "-framework", "CoreGraphics"],
       cwd=os.path.join(src, "hspist3"), log=os.path.join(R, "logs", "build_binary.txt"))
    open(os.path.join(hs, "logs", "BUILD_KOA_LAST.hash"), "w").write(f"{sha(os.path.join(hs, '00ALLINONE'))}  00ALLINONE\n")
    return sh([os.path.join(hs, "00ALLINONE"), "--version"]).stdout.splitlines()[0], sha(os.path.join(hs, "00ALLINONE"))


def k2_clone(home):
    d = os.path.join(home, "harddisks_gen3_k2")
    sh(["git", "clone", "-q", "--no-checkout", "--shared", REPO, d]); sh(["git", "-C", d, "config", "core.abbrev", "7"])
    sh(["git", "-C", d, "config", "core.sparseCheckout", "true"])
    open(os.path.join(d, ".git", "info", "sparse-checkout"), "w").write("\n".join(SPARSE_K2) + "\n")
    sh(["git", "-C", d, "-c", "advice.detachedHead=false", "checkout", "-q", K2_NOW])
    g = os.path.join(d, "hspist3", "logs"); os.makedirs(g)
    open(os.path.join(g, "gen3-testG-merge-k2_15044800.out"), "w").write("merge of ... by the scripts of clone HEAD fd8d8a2...\nchunk 1 done\n")
    open(os.path.join(g, "gen3-testG_15043031_1.out"), "w").write("env: ...\nchunk 1 on cn-03-38-07 (E5): lines 1-25; 00ALLINONE  git a643bc7  target koa\nchunk 1 done\n")
    return d


def runsheet_blocks(path):
    out, cur, shas, last = {}, None, {}, None
    for l in open(path).read().splitlines():
        m = re.match(r"^TEST R BLOCK ([A-Z]) \(", l)
        if m: cur = last = m.group(1); out[cur] = []; shas[cur] = []; continue
        if cur is not None and not l.strip(): cur = None; continue
        if cur is not None: out[cur].append(l); continue
        q = re.match(r"^\s+([0-9a-f]{64})  (\S+)$", l)
        if last and q: shas[last].append(f"{q.group(1)}  {q.group(2)}")
    return out, shas


def run_block(R, X, lines, env, tag=""):
    p, lg = os.path.join(R, "logs", f"block_{X}{tag}.sh"), os.path.join(R, "logs", f"block_{X}{tag}.out")
    open(p, "w").write("".join(l.replace(LUSTRE, os.path.join(R, "lustre")) + "\n" for l in lines))
    t0 = time.time(); q = subprocess.run(["/bin/bash", p], cwd=env["HOME"], env=env, stdout=open(lg, "w"), stderr=subprocess.STDOUT)
    print(f"  block {X}{tag}: exit {q.returncode}, {time.time() - t0:.0f} s, log logs/block_{X}{tag}.out", flush=True)
    return q.returncode, open(lg, errors="ignore").read()


def substitute(R, hs, mac_sha):
    """the rehearsal substitutions in the updated clone (listed in the table)"""
    reg = open(os.path.join(hs, P, "tasks_testR_koa_261011.txt")).read().splitlines(keepends=True)
    mini = "".join(reg[:MINI_N])
    open(os.path.join(hs, P, "tasks_testR_koa_261011.txt"), "w").write(mini)
    open(os.path.join(hs, P, "tasks_testR_koa_261011.sha256"), "w").write(hashlib.sha256(mini.encode()).hexdigest() + "\n")
    open(os.path.join(hs, P, "testR_chunk_size.txt"), "w").write(f"{MINI_K}\n")
    rec = os.path.join(hs, P, "build_record_a643bc7.txt"); s = open(rec).read()
    open(rec, "w").write(re.sub(r"^sha256 [0-9a-f]{64}$", f"sha256 {mac_sha}", s, flags=re.M))
    done = []
    for f in ("testR_koa.sbatch", "testR_merge.sbatch", "fetch_testR_koa.sh"):
        p = os.path.join(hs, P, f); t = open(p).read()
        t2 = t.replace(LUSTRE, os.path.join(R, "lustre")).replace("charing@koa-dtn.its.hawaii.edu:", "")
        if t2 != t: open(p, "w").write(t2); done.append(f)
    out = sh([sys.executable, os.path.join(P, "make_testR_manifests.py")], cwd=hs).stdout
    return mini.splitlines(), done, out


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--root", required=True); ap.add_argument("--pkg", required=True)
    ap.add_argument("--runsheet", default=os.path.join(HERE, "runsheet_testR_261011.txt"))
    ap.add_argument("--reuse-binary", default=None, help="the Mac koa-name binary of an earlier run of this driver, symlinked instead of rebuilt (the data cap); listed")
    ap.add_argument("--recheck-block-D", default=None, help="apply the block D check to a saved block_D.out and stop")
    a = ap.parse_args(); R = os.path.abspath(a.root)
    if a.recheck_block_D:
        ok_d, n_d = check_D(open(a.recheck_block_D, errors="ignore").read())
        print(f"block D re-checked on {a.recheck_block_D}: {'PASS' if ok_d else '**FAIL**'} ({n_d}/6 cell lines)"); return
    if not os.path.isdir(R) or os.listdir(R): sys.exit(f"STOP: {R} must be an existing empty folder")
    PKG = sh(["git", "-C", REPO, "rev-parse", a.pkg]).stdout.strip()
    for d in ("home", "lustre", "logs"): os.makedirs(os.path.join(R, d))
    home = os.path.join(R, "home"); shim = shims(R)
    blocks, shas = runsheet_blocks(a.runsheet)
    print(f"# Test R rehearsal: package {PKG} (engine-gen3), runsheet {a.runsheet} (blocks {''.join(sorted(blocks))}), root {R}\n", flush=True)
    m = re.search(r"git fetch --depth 1 origin ([0-9a-f]{40})", "\n".join(blocks.get("B", [])))
    res("runsheet block B", "fetches the package commit", bool(m) and m.group(1) == PKG, f"block B fetches {m.group(1) if m else None}")
    if a.reuse_binary:
        hs0 = os.path.join(home, "harddisks_gen3", "hspist3"); os.makedirs(os.path.join(hs0, "logs"))
        os.symlink(os.path.abspath(a.reuse_binary), os.path.join(hs0, "00ALLINONE"))
        mac_sha = sha(a.reuse_binary); open(os.path.join(hs0, "logs", "BUILD_KOA_LAST.hash"), "w").write(f"{mac_sha}  00ALLINONE\n")
        ver = sh([os.path.join(hs0, "00ALLINONE"), "--version"]).stdout.splitlines()[0]
        print(f"  the binary: {a.reuse_binary} (symlinked, not rebuilt; SHA-256 {mac_sha[:12]}...)", flush=True)
    else:
        ver, mac_sha = build_binary(R, home)
    res("the Mac binary (a643bc7, koa target name)", "version line '00ALLINONE  git a643bc7  target koa'", ver == "00ALLINONE  git a643bc7  target koa", ver)
    k2 = k2_clone(home); k2hs = os.path.join(k2, "hspist3")
    E = env_for(R, home, shim)
    # ---- A, B
    rc, out = run_block(R, "A", blocks["A"], E)
    res("block A", "squeue, srun (shims)", rc == 0 and "rehearsal shim: srun" in out, f"exit {rc}")
    rc, out = run_block(R, "B", blocks["B"], env_for(R, home, shim, SLURM_JOB_ID=950001))
    want_sha = shas["B"]
    got = re.findall(r"^[0-9a-f]{64}  cluster/gen3_testR_261011/\S+$", out, re.M)
    ok_b = ["0" in L(out), f"HEAD is now at {PKG[:7]}" in out, got == want_sha, "20" in L(out),
            "00ALLINONE  git a643bc7  target koa" in out, out.count("ALL PASS") == 3, "the nu_d imports: OK" in out,
            "No such file or directory" in out]
    res("block B (the update in place, the checks, the unit tests)", "as the runsheet's good lines, but KOA's four binaries (not on the Mac)",
        all(ok_b), f"checks {ok_b}; the eight SHA lines {'equal' if got == want_sha else 'DIFFER'}; KOA's binaries: "
        f"{len(re.findall(r': FAILED', out))} FAILED here (expected: the Mac build and three absent files)")
    # ---- substitutions
    mini, rew, man = substitute(R, k2hs, mac_sha)
    res("the rehearsal substitutions", "mini list (12 lines, chunk size 6), build record, manifests, Lustre path, local fetch",
        len(mini) == MINI_N and len(rew) == 3, f"rewritten {rew}; manifests: {man.count('files)')} written")
    # ---- C, D: chunk 1
    rc, out = run_block(R, "C", blocks["C"], E)
    res("block C (chunk 1 alone)", "Submitted batch job", rc == 0 and "Submitted batch job" in out, f"exit {rc}")
    rc, out = run_block(R, "D", blocks["D"], E)
    ok_d, n_d = check_D(out)
    res("block D (chunk 1's checks)", "chunk 1 on ... lines 1-6; chunk 1 done; one trace and run per cell; no FAILED/REFUSED/STOP",
        ok_d, f"{n_d}/6 cell lines")
    # ---- E: the rest, its first attempt stops at start-up
    linesE = [l.replace("--array=2-216%64", "--array=2") for l in blocks["E"]]
    rc, out = run_block(R, "E", linesE, env_for(R, home, shim, REHEARSAL_STARTUP_FAIL=1, REHEARSAL_HOST="cn-07-07-07"))
    lg = [p for p in os.listdir(os.path.join(k2hs, "logs")) if p.startswith("gen3-testR_") and p.endswith("_2.out")]
    l2 = open(os.path.join(k2hs, "logs", lg[0])).read() if len(lg) == 1 else ""
    res("block E (the rest, --array=2 here; the first attempt on a pretended faulty node)", "STOP: [cn-07-07-07] ... before 'chunk 2 on'; no folder",
        "STOP: [cn-07-07-07] the Test R inputs are not the delivered ones" in l2 and "chunk 2 on" not in l2
        and not os.path.exists(os.path.join(R, "lustre", "harddisks_gen3", "testR_koa", "chunks", "chunk_2")), l2.strip().splitlines()[-1] if l2 else "no log")
    # ---- F: the checks and the resubmission line
    linesF = [l.replace("testR_stopped_chunks.sh logs gen3-testR 216", "testR_stopped_chunks.sh logs gen3-testR 2") for l in blocks["F"]]
    rc, out = run_block(R, "F", linesF, E)
    rl = re.findall(r"^sbatch --array=\S+ --exclude=\S+ cluster/gen3_testR_261011/testR_koa\.sbatch$", out, re.M)
    res("block F (checks, the resubmission line)", "2 logs, 1 done, 0 FAILED; the line --array=2 --exclude=cn-05-08-07,cn-07-07-07",
        rl == ["sbatch --array=2 --exclude=cn-05-08-07,cn-07-07-07 cluster/gen3_testR_261011/testR_koa.sbatch"], f"{rl}")
    rc, out = run_block(R, "R", ["squeue -u charing", f"cd ~/harddisks_gen3_k2/hspist3 && {rl[0]}"] if rl else ["true"], E)
    rc, out = run_block(R, "F", linesF, E, "2")
    res("block F again (after the resubmission)", "3 logs, 2 done, 0 FAILED, nothing to resubmit",
        "3" in L(out) and "2" in L(out) and "0" in L(out) and "nothing to resubmit" in out, "")
    # ---- the merge job's own STOP rules
    q = os.path.join(R, "logs", "pretend_squeue.txt"); open(q, "w").write("930999_7 PENDING\n")
    rc, out = run_block(R, "G", blocks["G"], env_for(R, home, shim, REHEARSAL_SQUEUE=q), "_pending")
    ex = open(os.path.join(R, "logs", "sbatch_shim_exits.txt")).read().splitlines()[-1]
    mlog = sorted(p for p in os.listdir(os.path.join(k2hs, "logs")) if p.startswith("gen3-testR-merge_"))
    t = open(os.path.join(k2hs, "logs", mlog[-1])).read() if mlog else ""
    res("the merge job: a pending gen3-testR job", "exit 2, nothing moved", ex.endswith(" 2") and "still pending or running" in t
        and not os.path.exists(os.path.join(R, "lustre", "harddisks_gen3", "testR_koa", "merged")), ex)
    lg2 = sorted(p for p in os.listdir(os.path.join(k2hs, "logs")) if re.match(r"^gen3-testR_\d+_2\.out$", p))
    done2 = [p for p in lg2 if "chunk 2 done" in open(os.path.join(k2hs, "logs", p)).read()]
    hold = os.path.join(R, "logs", "held_aside_" + done2[0]); shutil.move(os.path.join(k2hs, "logs", done2[0]), hold)
    rc, out = run_block(R, "G", blocks["G"], E, "_heldaside")
    ex = open(os.path.join(R, "logs", "sbatch_shim_exits.txt")).read().splitlines()[-1]
    mlog = sorted(p for p in os.listdir(os.path.join(k2hs, "logs")) if p.startswith("gen3-testR-merge_"))
    t = open(os.path.join(k2hs, "logs", mlog[-1])).read()
    shutil.move(hold, os.path.join(k2hs, "logs", done2[0]))
    res("the merge job: chunk 2's done log held aside", "exit 2, 'NOT DONE: 2(no log with its done line)', nothing moved",
        ex.endswith(" 2") and "NOT DONE: 2(no log with its done line)" in t
        and not os.path.exists(os.path.join(R, "lustre", "harddisks_gen3", "testR_koa", "merged")), ex)
    # ---- G, H: the merge
    rc, out = run_block(R, "G", blocks["G"], E)
    ex = open(os.path.join(R, "logs", "sbatch_shim_exits.txt")).read().splitlines()[-1]
    rc, out = run_block(R, "H", blocks["H"], E)
    plan = "plan entries by treatment: APPEND 18, EQUAL 12, MOVE 18, PER CHUNK 12"
    okH = [ex.endswith(" 0"), "merge done; next: the fetch" in out, "scan: 2 chunks, 60 plan entries, 0 problem(s)" in out, plan in out,
           "merged 2 chunks: 18 files moved, 18 files appended, 12 kept per chunk (0 failed-run folders), 0 incomplete run folders, "
           "0 locks or temp files left in place, 12 build files equal; build '00ALLINONE  git a643bc7  target koa'" in out,
           out.count(": exit 0 2 trajectories: nud_runs.csv (0 fit failures), runs_health.csv (2 ##RUN sections), failed_runs.csv (0 failed runs)") == 6,
           out.count("2 traces compressed") == 6, "earlier attempt of chunk 2, stopped in the start-up checks (ran nothing)" in out,
           "4" in L(out)]
    res("blocks G and H (the merge, the reduction, the compression; its record)", "exit 0; 60 plan entries; the merged line; 6 reductions; "
        "6 compressions; chunk 2's earlier attempt listed; 4 job logs copied", all(okH), f"checks {okH}")
    rc, out = run_block(R, "G", blocks["G"], E, "_second")
    ex = open(os.path.join(R, "logs", "sbatch_shim_exits.txt")).read().splitlines()[-1]
    res("a second merge submission", "exit 3 (merged/ exists)", ex.endswith(" 3"), ex)
    # ---- the fetch (local) and the verdict (component check)
    p = sh(["bash", os.path.join(P, "fetch_testR_koa.sh")], cwd=k2hs, env=E, log=os.path.join(R, "logs", "fetch.out"), check=False)
    F = os.path.join(k2hs, "experiments_gen3_testR_koa_261011")
    fetched = [os.path.relpath(os.path.join(dp, f), F) for dp, dn, fn in os.walk(F) for f in fn]
    raw = [f for f in fetched if "wall_x_positions_" in f or "psi6_t_" in f or ".done_runs" in f or f.startswith("chunks/")]
    res("the fetch (local rsync)", "exit 0; the merge record, the tables, the job logs; no trace, psi6 or run folder",
        p.returncode == 0 and "merged/.merge_complete" in fetched and not raw and sum(f.endswith("nud_runs.csv") for f in fetched) == 6
        and sum(f.startswith("job_logs/") for f in fetched) == 4, f"{len(fetched)} files; raw {len(raw)}")
    VE = dict(E, HD_TESTR_LOC=F, HD_TESTR_TASKS=os.path.join(k2hs, P, "tasks_testR_koa_261011.txt"),
              HD_TESTR_TASKS_SHA=os.path.join(k2hs, P, "tasks_testR_koa_261011.sha256"))
    # the verdict runs on the Mac from the repository (as after KOA's fetch), on the rehearsal's fetched folder (test hooks)
    v1 = sh([sys.executable, os.path.join(HS, P, "testR_koa_verdict_261011.py")], cwd=HS, env=VE, log=os.path.join(R, "logs", "verdict_first.out"), check=False)
    v1t = open(os.path.join(R, "logs", "verdict_first.out")).read()
    v2 = sh([sys.executable, os.path.join(HS, P, "testR_koa_verdict_261011.py")], cwd=HS, env=VE, log=os.path.join(R, "logs", "verdict_second.out"), check=False)
    v2t = open(os.path.join(R, "logs", "verdict_second.out")).read()
    res("the verdict on the mini data (component check)", "the inventory clean (2 per arm, 0 excluded); the numbers; a verdict line; information skipped",
        v1.returncode == 0 and "the registered one" in v1t and "TEST R: " in v1t and "INVALID" not in v1t and "information rows need at least 20" in v1t,
        [l for l in v1t.splitlines() if l.startswith("TEST R:")][0] if "TEST R:" in v1t else "no verdict line")
    res("the verdict's one-printing guard", "a second run stops", v2.returncode != 0 and "one printing" in v2t, v2t.strip().splitlines()[-1] if v2t.strip() else "")
    # ---- the table
    print("\n| step (script) | criterion | result | evidence |\n|---|---|---|---|")
    for r in RES: print("| " + " | ".join(str(x) for x in r) + " |")
    print("\nRESIDUAL DIFFERENCES: the shims (Slurm commands, module, gcc --version, lscpu, hostname, sha256sum); the Mac binary (its SHA-256 "
          f"in the build record; KOA's four binaries absent{'; symlinked from an earlier run, not rebuilt' if a.reuse_binary else ''}); the mini list (12 lines, chunk size 6) with its SHA-256 file and the manifests "
          "rewritten; the Lustre path; the local fetch; block E's --array=2; block F's chunk count 2; the verdict's test hooks")
    nb = sum(r[2] != "PASS" for r in RES)
    print(f"\nREHEARSAL: {'PASS' if nb == 0 else 'FAIL'} ({len(RES) - nb} of {len(RES)})")


if __name__ == "__main__":
    main()
