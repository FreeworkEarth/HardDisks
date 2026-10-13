#!/usr/bin/env python3
"""##CHRIS 2026-10-12 (plan-author decision 39 item 5 f; registration 261012 sec. 4.8.10, gate G4; decision 20's rule): THE MAC REHEARSAL OF
STAGE 1, Part B: "every script a KOA job executes runs end to end on the Mac with a mini task list that holds every block type and at
least two chunks that deliver into the same (cell, mass) folder: chunks, merge, reduction, fetch, tables". It is Test R's rehearsal
(cluster/gen3_testR_261011/rehearse_testR_mac.py) adapted: the runsheet's own paste blocks are run as Chris types them. Nothing outside
--root is written; the clones borrow the repository's objects (git clone --shared); one simulation process at a time.
THE STEPS:
  1. the build (gate block B's build lines): cluster/gen3_s1_gate_261012/rehearse_s1_build_mac.py --root <root>/build --commit <pkg>, which
     clones f72d764 into home/harddisks_gen3_s1, builds it with build_s1_koa.sh and writes the build-time record with record_s1_build.sh,
     and clones the scripts clone home/harddisks_gen3_k2 at <pkg>; this home is used from here on;
  2. the steps CC takes after the gate, in the scripts clone: build_record_f72d764.txt = a copy of the build-time record; the two input
     manifests by make_s1_manifests.py;
  3. the REHEARSAL SUBSTITUTIONS (each listed): the MINI LIST (37 registered task lines, unchanged, in their registered order: chunk 1's 25
     lines, B_N0100_eta0.392699 M 50 s 0-1 and M 500 s 0-3 (the snapshot mass), D_N0100_eta0.075000 M 500 s 0-1, S_N0100_eta0.600000 M 100
     s 0-1 in both starts), its SHA-256 file and a MINI CHUNK TABLE (7 chunks in the groups B100, D, S; chunks 3 and 4 deliver into the same
     folder B_N0100_eta0.392699/m_500; two time classes in B100 and in S) before the manifests; the Lustre path in the two job scripts; the
     fetch's DTN source made local; all of it committed in the scratch clone as a local commit (the stage 1 commit's stand-in);
  4. make_s1_runsheet.py on the mini table writes the mini runsheet; its blocks run in order: A, B (without its two update lines: the
     stand-in commit is checked out already), C, D, then per group E, F, G, H, I; group S's first attempt stops in the start-up checks
     (a pretended faulty node: the sha256sum shim fails) and its printed resubmission lines are typed as printed; group B100's merge is
     first refused twice (a pretended pending job; a chunk's done log held aside: exit 2, nothing moved) and, after the merge, refused once
     more (exit 3);
  5. the verdict script (G5) on the fetched folder: --data <fetched> --interim (the tables of the groups present).
##CHRIS 2026-10-12 (plan-author decision 42 item 5 b): the runsheet's submit blocks (every E, every G, block F's resubmission line) are guarded
one-line commands; each is typed once (it submits) and a second time (it refuses: "REFUSED: ...", "THIS PASTE DID NOT SUBMIT ...", no sbatch
call); E-B100 after chunk 1; group D's block E also with a pretended pending job and with a failing squeue (REHEARSAL_SQUEUE_FAIL); group S's
guarded resubmission line once and again; B100's merge job is still shown to refuse a second submission without the guard (exit 3).
SHIMS (residual differences, each printed in the table): squeue (empty, or the lines of REHEARSAL_SQUEUE), sinfo (a pretended line), srun
(not run), sbatch (runs the job at once with a pretended job id and its #SBATCH output name; --array lists one after the other, the
throttle %N dropped; --job-name, --time, --exclude, --dependency taken or recorded; the script's arguments passed), module (nothing), gcc
--version (14.3.0), lscpu (the Mac's CPU), hostname (REHEARSAL_HOST, else the Mac's), sha256sum (fails when REHEARSAL_STARTUP_FAIL is set);
HOME; SDL's dummy drivers; the Mac binaries of step 1 (f72d764 with the koa target's name, Mac flags).
usage (from hspist3/): python3 cluster/gen3_stage1_261012/rehearse_s1_mac.py --root <empty dir> --pkg <engine-gen3 commit>
"""
import argparse, hashlib, os, re, shutil, subprocess, sys, time
HERE = os.path.dirname(os.path.abspath(__file__)); HS = os.path.dirname(os.path.dirname(HERE))
P = "cluster/gen3_stage1_261012"; LUSTRE = "/mnt/lustre/koa/scratch/charing"
BUILD = "00ALLINONE  git f72d764  target koa"
MINI = [("B_N0100_eta0.050000", 50, range(25)), ("B_N0100_eta0.392699", 50, range(2)), ("B_N0100_eta0.392699", 500, range(4)),
        ("D_N0100_eta0.075000", 500, range(2)), ("S_N0100_eta0.600000_start0", 100, range(2)), ("S_N0100_eta0.600000_turned", 100, range(2))]
# the mini chunks: (group, mini lines first-last, class): chunk 1 = the registered chunk 1; 3 and 4 into one folder
MCH = [("B100", 1, 25, 30), ("B100", 26, 27, 30), ("B100", 28, 29, 30), ("B100", 30, 31, 60), ("D", 32, 33, 30), ("S", 34, 35, 30), ("S", 36, 37, 60)]
RES = []


def sh(cmd, cwd=None, env=None, log=None, check=True):
    if log:
        with open(log, "w") as out: p = subprocess.run(cmd, cwd=cwd, env=env, stdout=out, stderr=subprocess.STDOUT)
    else: p = subprocess.run(cmd, cwd=cwd, env=env, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
    if check and p.returncode != 0:
        sys.exit(f"STOP: {' '.join(cmd)} (in {cwd}) exit {p.returncode}" + (f"; see {log}" if log else f": {(p.stdout or '')[-800:]}"))
    return p


def sha(p): return hashlib.sha256(open(p, "rb").read()).hexdigest()


def L(out): return [re.sub(r"[ \t]+", " ", l).strip() for l in out.splitlines()]


def res(step, crit, ok, ev):
    RES.append((step, crit, "PASS" if ok else "**FAIL**", ev)); print(f"    -> {step}: {'PASS' if ok else '**FAIL**'} ({ev})", flush=True)


SBATCH_SHIM = r"""#!/bin/bash
# rehearsal shim: runs the job script AT ONCE, as Slurm would run it: SLURM_SUBMIT_DIR = this folder, a pretended job id, the output file
# of its #SBATCH --output line (%x the job name, %j the id, %A %a the array id and task). --array=<list>[%N] runs its tasks one after the
# other (the throttle dropped); --job-name sets the name; --time, --exclude, --dependency are recorded (the Mac is one node); any other
# option is refused; the words after the script are its arguments
ARR=""; EXC=""; NAME=""; TIM=""; DEP=""
while [ $# -gt 0 ]; do
  case "$1" in
    --array=*) ARR=${1#--array=}; ARR=${ARR%%\%*};; --exclude=*) EXC=${1#--exclude=};; --job-name=*) NAME=${1#--job-name=};;
    --time=*) TIM=${1#--time=};; --dependency=*) DEP=${1#--dependency=};;
    --*) echo "rehearsal sbatch shim: option not supported: $1" >&2; exit 97;;
    *) break;;
  esac
  shift
done
S=$1; shift; C="$REHEARSAL_ROOT/logs/sbatch_shim_counter"
id=$(( $(cat "$C" 2>/dev/null || echo 940000) + 1 )); echo "$id" > "$C"
[ -n "$NAME" ] || NAME=$(sed -n 's/^#SBATCH --job-name=//p' "$S" | head -1)
out=$(sed -n 's/^#SBATCH --output=//p' "$S" | head -1)
echo "Submitted batch job $id"
echo "$id submit $S $* name=$NAME time=${TIM:-default} exclude=${EXC:-default} dependency=${DEP:-none} array=${ARR:-none}" >> "$REHEARSAL_ROOT/logs/sbatch_shim_exits.txt"
if [ -z "$ARR" ]; then
  o=${out//%x/$NAME}; o=${o//%j/$id}
  SLURM_SUBMIT_DIR=$PWD SLURM_JOB_ID=$id SLURM_JOB_NAME=$NAME bash "$S" "$@" > "$o" 2>&1; rc=$?
  echo "$id $S $rc" >> "$REHEARSAL_ROOT/logs/sbatch_shim_exits.txt"
else
  for t in $(echo "$ARR" | tr ',' '\n' | awk -F- '{ if (NF == 2) for (i = $1; i <= $2; i++) print i; else print $1 }'); do
    o=${out//%x/$NAME}; o=${o//%A/$id}; o=${o//%a/$t}
    SLURM_SUBMIT_DIR=$PWD SLURM_JOB_ID=$((id * 1000 + t)) SLURM_ARRAY_JOB_ID=$id SLURM_ARRAY_TASK_ID=$t SLURM_JOB_NAME=$NAME bash "$S" "$@" > "$o" 2>&1; rc=$?
    echo "${id}_$t $S $rc" >> "$REHEARSAL_ROOT/logs/sbatch_shim_exits.txt"
  done
fi
exit 0
"""


def shims(R):
    b = os.path.join(R, "shim", "bin"); os.makedirs(b)
    files = {
        "module": "#!/bin/bash\n# rehearsal shim: Lmod is not on the Mac\nexit 0\n",
        "gcc": "#!/bin/bash\n# rehearsal shim: koa_env.sh checks 'gcc --version' for 14.3; no stage 1 job compiles\n"
               "if [ \"${1:-}\" = --version ]; then echo \"gcc (GCC) 14.3.0 [rehearsal shim; here: $(/usr/bin/cc --version | head -1)]\"; exit 0; fi\n"
               "echo 'rehearsal shim: gcc is not used inside a job' >&2; exit 97\n",
        "lscpu": "#!/bin/bash\n# rehearsal shim: Linux's lscpu\necho \"Model name:          $(sysctl -n machdep.cpu.brand_string) (Mac rehearsal)\"\n",
        "hostname": "#!/bin/bash\n# rehearsal shim: a pretended node name (REHEARSAL_HOST), else the Mac's\n"
                    "if [ -n \"${REHEARSAL_HOST:-}\" ]; then echo \"$REHEARSAL_HOST\"; else exec /bin/hostname -s; fi\n",
        "squeue": "#!/bin/bash\n# rehearsal shim: Slurm's queue -- empty, or the lines of the file REHEARSAL_SQUEUE (a pretended pending job); with\n"
                  "# REHEARSAL_SQUEUE_FAIL set it fails as a slurmctld time-out would (decision 42: a failed squeue refuses)\n"
                  "[ -n \"${REHEARSAL_SQUEUE_FAIL:-}\" ] && { echo 'slurm_load_jobs error: Socket timed out on send/recv operation (rehearsal shim)' >&2; exit 1; }\n"
                  "[ -n \"${REHEARSAL_SQUEUE:-}\" ] && cat \"$REHEARSAL_SQUEUE\"\nexit 0\n",
        "sinfo": "#!/bin/bash\n# rehearsal shim: a pretended summary line of the shared partition\n"
                 "echo 'PARTITION AVAIL  TIMELIMIT   NODES(A/I/O/T) NODELIST'; echo 'shared       up 3-00:00:00        70/12/8/90 cn-[...] (rehearsal shim)'\n",
        "srun": "#!/bin/bash\n# rehearsal shim: the sandbox shell is not opened (the next block runs in a new shell)\necho \"rehearsal shim: srun $* (not run)\"\nexit 0\n",
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
    e.update(HOME=home, PATH=f"{shim}:{os.environ['PATH']}", REHEARSAL_ROOT=R, SDL_VIDEODRIVER="dummy", SDL_RENDER_DRIVER="software",
             **{k: str(v) for k, v in extra.items()})
    return e


def blocks_of(path):
    out, cur = {}, None
    for l in open(path).read().splitlines():
        m = re.match(r"^S1 BLOCK ([A-Z](?:-\w+)?) \(", l)
        if m: cur = m.group(1); out[cur] = []; continue
        if cur is not None and not l.strip(): cur = None; continue
        if cur is not None: out[cur].append(l)
    return out


def run_block(R, X, lines, env, tag=""):
    p, lg = os.path.join(R, "logs", f"block_{X}{tag}.sh"), os.path.join(R, "logs", f"block_{X}{tag}.out")
    open(p, "w").write("".join(l.replace(LUSTRE, os.path.join(R, "lustre")) + "\n" for l in lines))
    t0 = time.time(); q = subprocess.run(["/bin/bash", p], cwd=env["HOME"], env=env, stdout=open(lg, "w"), stderr=subprocess.STDOUT)
    print(f"  block {X}{tag}: exit {q.returncode}, {time.time() - t0:.0f} s, log logs/block_{X}{tag}.out", flush=True)
    return q.returncode, open(lg, errors="ignore").read()


def last_exit(R): return open(os.path.join(R, "logs", "sbatch_shim_exits.txt")).read().splitlines()[-1]


def substitute(R, hs, home):
    """the steps after the gate and the rehearsal substitutions in the scripts clone (step 2 and 3)"""
    reg = [l for l in open(os.path.join(hs, P, "tasks_stage1_koa_261012.txt")).read().splitlines(keepends=True)]
    want = []
    for cid, M, ss in MINI:
        for s in ss:
            hit = [i for i, l in enumerate(reg) if l.split()[1] == f"stage1/gen3/{cid}/m_{M}" and l.split()[3] == str(s)]
            if len(hit) != 1: sys.exit(f"STOP: the registered list has {len(hit)} lines for {cid} M {M} s {s}")
            want.append(hit[0])
    if want != sorted(want): sys.exit("STOP: the mini lines are not in the registered order")
    mini = "".join(reg[i] for i in want)
    open(os.path.join(hs, P, "tasks_stage1_koa_261012.txt"), "w").write(mini)
    open(os.path.join(hs, P, "tasks_stage1_koa_261012.sha256"), "w").write(hashlib.sha256(mini.encode()).hexdigest() + "\n")
    CT = {(c.split("\t")[5], c.split("\t")[6]): c.rstrip("\n").split("\t") for c in open(os.path.join(hs, P, "chunks_stage1_koa_261012.tsv")).readlines()[1:]}
    ML = mini.splitlines(); rows = ["chunk\tfirst\tlast\tlines\tgroup\tcell\tM\tpredicted_s_E5\tlimit_min\tclass_min"]
    for j, (g, a, b, cl) in enumerate(MCH, 1):
        f = ML[a - 1].split(); cid, M = f[1].split("/")[2], f[2]; reg_c = [c for k, c in CT.items() if k == (cid, M)][0]
        per = float(reg_c[7]) / int(reg_c[3]); pred = per * (b - a + 1)
        rows.append(f"{j}\t{a}\t{b}\t{b - a + 1}\t{g}\t{cid}\t{M}\t{pred:.1f}\t{max(30.0, 3 * pred / 60):.1f}\t{cl}")
    open(os.path.join(hs, P, "chunks_stage1_koa_261012.tsv"), "w").write("\n".join(rows) + "\n")
    shutil.copy(os.path.join(home, "harddisks_gen3_s1", "hspist3", "logs", "BUILD_RECORD_S1.txt"), os.path.join(hs, P, "build_record_f72d764.txt"))
    done = []
    for f in ("s1_koa.sbatch", "s1_merge.sbatch", "fetch_s1_koa.sh"):
        p = os.path.join(hs, P, f); t = open(p).read()
        t2 = t.replace(LUSTRE, os.path.join(R, "lustre")).replace("charing@koa-dtn.its.hawaii.edu:", "")
        if t2 != t: open(p, "w").write(t2); done.append(f)
    man = sh([sys.executable, os.path.join(P, "make_s1_manifests.py")], cwd=hs).stdout
    return mini.splitlines(), rows, done, man


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--root", required=True); ap.add_argument("--pkg", required=True); a = ap.parse_args()
    R = os.path.abspath(a.root)
    if not os.path.isdir(R) or os.listdir(R): sys.exit(f"STOP: {R} must be an existing empty folder")
    PKG = sh(["git", "-C", HS, "rev-parse", a.pkg]).stdout.strip()
    for d in ("lustre", "logs"): os.makedirs(os.path.join(R, d))
    print(f"# stage 1 rehearsal (gate G4): package {PKG} (engine-gen3), root {R}\n", flush=True)
    # ---- 1. the build
    b = sh([sys.executable, os.path.join(HS, "cluster/gen3_s1_gate_261012/rehearse_s1_build_mac.py"), "--root", os.path.join(R, "build"),
            "--commit", PKG], cwd=HS, log=os.path.join(R, "logs", "build_rehearsal.out"), check=False)
    bt = open(os.path.join(R, "logs", "build_rehearsal.out")).read()
    res("1. the build (rehearse_s1_build_mac.py: build_s1_koa.sh, record_s1_build.sh on a fresh f72d764 clone)", "BUILD REHEARSAL: PASS",
        b.returncode == 0 and "BUILD REHEARSAL: PASS" in bt, [l for l in bt.splitlines() if l.startswith("BUILD REHEARSAL")][:1])
    home = os.path.join(R, "build", "home"); k2 = os.path.join(home, "harddisks_gen3_k2"); hs = os.path.join(k2, "hspist3")
    shim = shims(R); E = env_for(R, home, shim); os.makedirs(os.path.join(hs, "logs"), exist_ok=True)
    # ---- 2, 3. after the gate; the substitutions; the stand-in commit
    mini, rows, rew, man = substitute(R, hs, home)
    sh(["git", "-C", k2, "add", "-u"]); sh(["git", "-C", k2, "add", f"hspist3/{P}/build_record_f72d764.txt", f"hspist3/{P}/s1_inputs.sha256", f"hspist3/{P}/s1_merge_inputs.sha256"])
    sh(["git", "-C", k2, "-c", "user.name=rehearsal", "-c", "user.email=rehearsal@localhost", "commit", "-q", "-m",
        "REHEARSAL stand-in of the stage 1 commit: build record, manifests, mini list and table, Lustre path, local fetch"])
    Q = sh(["git", "-C", k2, "rev-parse", "HEAD"]).stdout.strip()
    res("2-3. after the gate (the build record, the manifests) and the substitutions", "mini list 37 lines, mini table 7 chunks, 3 scripts rewritten, "
        "2 manifests; one local commit", len(mini) == 37 and len(rows) == 8 and len(rew) == 3 and man.count("files)") == 2,
        f"rewritten {rew}; stand-in commit {Q[:12]}")
    rs = os.path.join(R, "logs", "runsheet_stage1_mini.txt")
    sh([sys.executable, os.path.join(P, "make_s1_runsheet.py"), "--commit", Q, "--out", rs], cwd=hs)
    B = blocks_of(rs)
    res("4. the mini runsheet (make_s1_runsheet.py on the mini table)", "blocks A-D and E-I for B100, D, S",
        sorted(B) == sorted(["A", "B", "C", "D"] + [f"{x}-{g}" for g in ("B100", "D", "S") for x in "EFGHI"]), ", ".join(B))
    # ---- A, B
    rc, out = run_block(R, "A", B["A"], E)
    res("block A", "squeue, srun (shims)", rc == 0 and "rehearsal shim: srun" in out, f"exit {rc}")
    linesB = [l for l in B["B"] if not l.startswith(("git fetch ", "git -c advice.detachedHead=false checkout "))]
    rc, out = run_block(R, "B", linesB, env_for(R, home, shim, SLURM_JOB_ID=950001))
    want = re.findall(r"^\s+([0-9a-f]{64}  \S+)$", "\n".join(open(rs).read().split("S1 BLOCK B (")[1].split("S1 BLOCK C (")[0].splitlines()), re.M)
    got = re.findall(r"^[0-9a-f]{64}  cluster/gen3_stage1_261012/\S+$", out, re.M)
    nman = sum(1 for f in ("s1_inputs.sha256", "s1_merge_inputs.sha256") for l in open(os.path.join(hs, P, f)) if l.strip())
    rec = open(os.path.join(home, "harddisks_gen3_s1", "hspist3", "logs", "BUILD_RECORD_S1.txt")).read()
    bins = re.findall(r"^([0-9a-f]{64})  (00ALLINONE|bin/\S+)$", out, re.M)
    okB = ["0" in L(out), got == want and len(got) == 13, str(nman) in L(out), "the committed build record equals the build-time record" in out,
           len(bins) == 4 and all(f"sha256 {h} {f}" in rec for h, f in bins), BUILD in out, "ALL PASS (8 of 8 cases)" in out,
           "the merge and reducer imports: OK" in out, "No such file or directory" in out]
    res("block B (without its two update lines)", "the runsheet's good lines", all(okB), f"checks {okB}; SHA lines {len(got)} of {len(want)}")
    # ---- C, D: chunk 1 alone
    rc, out = run_block(R, "C", B["C"], E)
    res("block C (chunk 1 alone)", "Submitted batch job; chunk 1 exit 0", rc == 0 and "Submitted batch job" in out and last_exit(R).endswith(" 0"), last_exit(R))
    rc, out = run_block(R, "D", B["D"], E)
    okD = [bool(re.search(r"^chunk 1 on \S+ \(.*\): group B100, lines 1-25; " + re.escape(BUILD) + "$", out, re.M)), "chunk 1 done" in out,
           not re.search(r"FAILED|REFUSED|^STOP", out, re.M),
           re.sub(r"[ \t]+", " ", f"traces 25 psi6 25 virial 25 snapshots 0 runs 25 clean 25 build '{BUILD}'") in L(out),   # macOS wc -l pads
           bool(re.search(r"^\s*25 box=3769\.9111843077517 x 240 px exact_box=1 eta=\S+ N=100 N_L=50 N_R=50$", out, re.M))]
    res("block D (chunk 1's checks, G6 (e))", "the runsheet's good lines", all(okD), f"checks {okD}")
    # ---- the groups. ##CHRIS 2026-10-12 (plan-author decision 42 item 5 b): every submit block is ONE guarded line; each is typed once
    # (it submits) and a second time (it refuses: "REFUSED: ...", "THIS PASTE DID NOT SUBMIT ...", no sbatch call); E-B100 after chunk 1;
    # group D's block E also with a pretended pending job and with a failing squeue; group S's guarded resubmission line once and again
    EX = os.path.join(R, "logs", "sbatch_shim_exits.txt")
    nex = lambda: len(open(EX).read().splitlines())
    tasks_since = lambda n0: [l for l in open(EX).read().splitlines()[n0:] if " submit " not in l]

    def refused(X, lines, env, tag, why):
        n0 = nex(); rc, out = run_block(R, X, lines, env, tag)
        ok = (rc == 0 and "REFUSED: " in out and "THIS PASTE DID NOT SUBMIT" in out and "Submitted batch job" not in out
              and "CHECKS GOOD" not in out and nex() == n0)
        res(f"block {X} typed {why}", "REFUSED, nothing submitted", ok, (re.findall(r"^REFUSED: .*$", out, re.M) or ["no REFUSED line"])[0][:160])
        return out

    for g in ("B100", "D", "S"):
        ncg = sum(1 for c in MCH if c[0] == g); nsb = sum(l.count(" sbatch --job-name=gen3-s1-") for l in B[f"E-{g}"])
        root = os.path.join(R, "lustre", "harddisks_gen3", "stage1", g)
        if g == "D":
            q = os.path.join(R, "logs", "pretend_squeue_E.txt"); open(q, "w").write("940998_5 PENDING\n")
            refused(f"E-{g}", B[f"E-{g}"], env_for(R, home, shim, REHEARSAL_SQUEUE=q), "_pending", "with a pretended pending job of the group")
            refused(f"E-{g}", B[f"E-{g}"], env_for(R, home, shim, REHEARSAL_SQUEUE_FAIL=1), "_squeuefail", "with a failing squeue")
        envE = env_for(R, home, shim, REHEARSAL_STARTUP_FAIL=1, REHEARSAL_HOST="cn-07-07-07") if g == "S" else E
        n0 = nex(); rc, out = run_block(R, f"E-{g}", B[f"E-{g}"], envE, "_faulty" if g == "S" else "")
        ex = tasks_since(n0)
        okE = [rc == 0, "CHECKS GOOD" in out, out.count("Submitted batch job") == nsb, f"{g}: ALL {nsb} ARRAY LINE(S) SUBMITTED" in out, "REFUSED" not in out]
        if g == "S":
            logs = [p for p in os.listdir(os.path.join(hs, "logs")) if p.startswith(f"gen3-s1-{g}_")]
            st = [open(os.path.join(hs, "logs", p)).read() for p in logs]
            okE += [len(st) == ncg, all("STOP: [cn-07-07-07] the stage 1 inputs are not the delivered ones" in s and not re.search(r"^chunk \d+ on ", s, re.M) for s in st),
                    not os.path.exists(os.path.join(root, "chunks"))]
            res(f"block E-{g} typed once, on a pretended faulty node", f"CHECKS GOOD, {nsb} array lines submitted; every chunk STOP: [cn-07-07-07] before "
                "'chunk j on'; no chunk folder", all(okE), f"checks {okE}")
        else:
            okE.append(len(ex) == ncg - (1 if g == "B100" else 0) and all(l.endswith(" 0") for l in ex))
            res(f"block E-{g} typed once" + (" (after chunk 1: the guard allows exactly its folder and its done log)" if g == "B100" else ""),
                f"CHECKS GOOD, {nsb} array line(s) submitted, every task exit 0", all(okE), f"checks {okE}")
        refused(f"E-{g}", B[f"E-{g}"], E, "_second", "a second time")
        if g == "S":
            rc, out = run_block(R, f"F-{g}", B[f"F-{g}"], E, "_first")
            rl = re.findall(r"^sbatch --job-name=gen3-s1-S --time=\S+ --array=\S+ --exclude=\S+ cluster/gen3_stage1_261012/s1_koa\.sbatch$", out, re.M)
            gs = out.split("THE GUARDED RESUBMISSION LINE (one line for all of them; type it once, after the go of the plan author):\n")
            gline = gs[1].splitlines()[0] if len(gs) == 2 else ""
            okF1 = [rl == ["sbatch --job-name=gen3-s1-S --time=0:30:00 --array=6 --exclude=cn-05-08-07,cn-07-07-07 cluster/gen3_stage1_261012/s1_koa.sbatch",
                           "sbatch --job-name=gen3-s1-S --time=1:00:00 --array=7 --exclude=cn-05-08-07,cn-07-07-07 cluster/gen3_stage1_261012/s1_koa.sbatch"],
                    gline.startswith("cd ~/harddisks_gen3_k2/hspist3 && { Q=$(squeue "), all(x + " && " in gline for x in rl), '"chunk (6|7) done"' in gline]
            res(f"block F-{g} after the faulty attempt", "the two plain lines, one per class, and ONE guarded line holding both, unchanged and in order",
                all(okF1), f"checks {okF1}")
            n0 = nex(); rc, out = run_block(R, f"R-{g}", [gline], E)
            ex = tasks_since(n0)
            okR = [rc == 0, "CHECKS GOOD" in out, out.count("Submitted batch job") == 2, f"{g}: ALL 2 RESUBMISSION LINE(S) SUBMITTED" in out,
                   len(ex) == 2 and all(l.endswith(" 0") for l in ex)]
            res(f"the guarded resubmission line of F-{g} typed once", "CHECKS GOOD, 2 lines submitted, both chunks exit 0", all(okR), f"checks {okR}")
            refused(f"R-{g}", [gline], E, "_second", "a second time")
            rc, out = run_block(R, f"F-{g}", B[f"F-{g}"], E)
            nlog = 2 * ncg
        else:
            rc, out = run_block(R, f"F-{g}", B[f"F-{g}"], E)
            nlog = ncg
        okF = [str(nlog) in L(out), str(ncg) in L(out), "0" in L(out), f"chunk logs: {nlog}; chunks done: {ncg} of {ncg} (group {g};" in out,
               f"earlier attempts: {ncg if g == 'S' else 0} stopped in the start-up checks, 0 started and not finished" in out,
               f"group {g} without a done log: stopped once 0, stopped twice or more 0, started and not finished 0, no log 0" in out,
               "nothing to resubmit" in out, "THE GUARDED RESUBMISSION LINE" not in out]
        res(f"block F-{g}" + (" again, after the resubmission" if g == "S" else ""), "the counts, the chunk check, nothing to resubmit", all(okF), f"checks {okF}")
        if g == "B100":
            q = os.path.join(R, "logs", "pretend_squeue.txt"); open(q, "w").write("940999_7 PENDING\n")
            refused(f"G-{g}", B[f"G-{g}"], env_for(R, home, shim, REHEARSAL_SQUEUE=q), "_pending", "with a pretended pending job")
            lg = sorted(p for p in os.listdir(os.path.join(hs, "logs")) if re.match(rf"^gen3-s1-{g}_\d+_4\.out$", p) and "chunk 4 done" in open(os.path.join(hs, "logs", p)).read())
            hold = os.path.join(R, "logs", "held_aside_" + lg[0]); shutil.move(os.path.join(hs, "logs", lg[0]), hold)
            rc, out = run_block(R, f"G-{g}", B[f"G-{g}"], E, "_heldaside")
            ml = sorted(p for p in os.listdir(os.path.join(hs, "logs")) if p.startswith(f"gen3-s1-merge-{g}_"))
            tl = open(os.path.join(hs, "logs", ml[-1])).read()
            shutil.move(hold, os.path.join(hs, "logs", lg[0]))
            res(f"the merge of {g}: chunk 4's done log held aside", "the guard passes (CHECKS GOOD); the merge job itself: exit 2, 'NOT DONE: 4(no log with "
                "its done line)', nothing moved", "CHECKS GOOD" in out and last_exit(R).endswith(" 2") and "NOT DONE: 4(no log with its done line)" in tl
                and not os.path.exists(os.path.join(root, "merged")), last_exit(R))
        rc, out = run_block(R, f"G-{g}", B[f"G-{g}"], E)
        gx = last_exit(R); okG = ["CHECKS GOOD" in out, f"{g}: THE MERGE SUBMITTED" in out, out.count("Submitted batch job") == 1]
        rc, out = run_block(R, f"H-{g}", B[f"H-{g}"], E)
        rs_txt = open(rs).read().split(f"S1 BLOCK H-{g} (")[1].split(f"S1 BLOCK I-{g} (")[0]
        exp_plan = re.search(r'"plan entries by treatment: ([^"]+)"', rs_txt).group(1)
        exp_merged = re.search(r'"(merged \d+ chunks: [^"]+?) incomplete\n\s+(run folders[^"]+)"', rs_txt)
        exp_m = exp_merged.group(1) + " incomplete " + exp_merged.group(2)
        exp_cnt = re.search(r'^\s+"(destinations \d+; reduce_B exit 0: \d+; reduce_s1 exit 0: \d+; compressed: \d+)"$', rs_txt, re.M).group(1)
        okH = okG + [gx.endswith(" 0"), f"plan entries by treatment: {exp_plan}" in out, exp_m in out, exp_cnt in out,
                     f"merge done; next: the fetch (fetch_s1_koa.sh {g})" in out, not re.search(r"^STOP|problem\(s\): NOTHING", out, re.M),
                     "fit failures)" in out, not re.search(r"trajectories reduced *$", out, re.M)]          # the filter hides every reduce_B line
        res(f"blocks G-{g} (typed once) and H-{g} (the merge, the reduction, the compression; its record)", f"CHECKS GOOD, 1 submitted, exit 0; "
            f"'{exp_plan}'; the merged line; '{exp_cnt}'", all(okH), f"checks {okH}")
        refused(f"G-{g}", B[f"G-{g}"], E, "_second", "a second time")
        if g == "B100":
            raw = re.search(r"(sbatch --job-name=gen3-s1-merge-B100 \S+ B100)", "\n".join(B[f"G-{g}"])).group(1)
            rc, out = run_block(R, f"G-{g}", ["cd ~/harddisks_gen3_k2/hspist3", raw], E, "_raw_second")
            res(f"the merge job of {g} submitted a second time without its guard", "exit 3 (merged/ exists)", last_exit(R).endswith(" 3"), last_exit(R))
        linesI = [l.replace("cd ~/Desktop/CCS_complex_coupled_systems/Repo/HardDisks_gen3 && ", "cd ~/harddisks_gen3_k2 && ") for l in B[f"I-{g}"]]
        rc, out = run_block(R, f"I-{g}", linesI, E)
        F = os.path.join(hs, "experiments_gen3_stage1_koa_261012", g)
        fetched = [os.path.relpath(os.path.join(dp, f), F) for dp, dn, fn in os.walk(F) for f in fn]
        pi8 = "merged/stage1/gen3/B_N0100_eta0.392699/"
        raw = [f for f in fetched if not f.startswith(pi8) and ("wall_x_positions_" in f or "psi6_t_" in f or "virial_blocks_" in f or ".done_runs" in f)]
        okI = [rc == 0, "merged/.merge_complete" in fetched, not [f for f in fetched if f.startswith("chunks/")], not raw,
               sum(f.endswith("/s1_runs.csv") for f in fetched) == {"B100": 3, "D": 1, "S": 2}[g],
               sum("snapshots_" in f for f in fetched) == {"B100": 2, "D": 2, "S": 4}[g]]
        if g == "B100": okI.append(sum(f.startswith(pi8) and "wall_x_positions_" in f for f in fetched) == 6 and any(f.startswith(pi8 + "m_500/.done_runs") for f in fetched))
        res(f"block I-{g} (the fetch, local)", "the merge record and the tables; no chunk folder; no trace, psi6, virial or run folder outside the "
            "(100, pi/8) cell, which comes whole; every snapshot", all(okI), f"{len(fetched)} files; checks {okI}")
    # ---- 5. the verdict (G5) on the fetched folder
    FD = os.path.join(hs, "experiments_gen3_stage1_koa_261012")
    v = sh([sys.executable, os.path.join(HS, P, "s1_verdict_261012.py"), "--data", FD, "--interim"], cwd=HS, env=E, log=os.path.join(R, "logs", "verdict_interim.out"), check=False)
    vt = open(os.path.join(R, "logs", "verdict_interim.out")).read()
    cells = ["B_N0100_eta0.050000", "B_N0100_eta0.392699", "D_N0100_eta0.075000", "S_N0100_eta0.600000_start0", "S_N0100_eta0.600000_turned"]
    res("5. the verdict script (G5) on the fetched folder, --interim", "exit 0; INTERIM; the five mini cells in the per-cell table with their "
        "crystal orientation", v.returncode == 0 and "INTERIM" in vt and all(re.search(rf"^\| {c} \| [xy]/[xy] \|", vt, re.M) for c in cells),
        f"exit {v.returncode}; " + "; ".join(f"{c}: {o}" for c, o in re.findall(r"^\| ((?:B|D|S)_N0100_\S+) \| ([^|]+?) \|", vt, re.M)))
    # ---- the table
    print("\n| step | criterion | result | evidence |\n|---|---|---|---|")
    for r in RES: print("| " + " | ".join(str(x) for x in r) + " |")
    print("\nRESIDUAL DIFFERENCES: the shims (Slurm commands, module, gcc --version, lscpu, hostname, sha256sum); the Mac binaries (f72d764, "
          "the koa target's name, Mac flags); the mini list (37 registered lines) with its SHA-256 file and the mini chunk table (7 chunks); the "
          "build record and the manifests written in the scratch clone and committed there (the stage 1 commit's stand-in); block B without its "
          "two update lines; the Lustre path; the local fetch; group S's first attempt on a pretended faulty node")
    nb = sum(r[2] != "PASS" for r in RES)
    print(f"\nREHEARSAL G4: {'PASS' if nb == 0 else 'FAIL'} ({len(RES) - nb} of {len(RES)})")
    print(f"data under {R}: {sh(['du', '-sh', R]).stdout.split()[0]}")


if __name__ == "__main__":
    main()
