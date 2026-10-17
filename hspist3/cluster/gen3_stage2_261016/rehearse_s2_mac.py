#!/usr/bin/env python3
"""##CHRIS 2026-10-16 (plan-author decision 58 item 3; notes 261012 sec. 4.7.142): THE MAC REHEARSAL OF THE S2 PACKAGE, as P2's
(rehearse_p2_mac.py): the runsheet's pastes run as Chris types them (its own lines, read from runsheet_s2_261017.txt), with P2's shims and these
REHEARSAL SUBSTITUTIONS, each listed in the output:
  STAND-IN of the package commit: a shared scratch clone of <pkg> with (i) the chunk table cut to the four rehearsal chunks -- N 900 turned
    (0.76, the middle height k 697), N 900 pin 0 (0.78, k 685), N 1600 at 0.72 (k 866) and at 0.80 (k 1024), each alpha 5, s = 0 (decision
    58 item 3) -- (ii) the Lustre root made local and the fetch's DTN source made local (rehearse_k2_mac.rewrite_lustre), (iii)
    s2_inputs.sha256 made again (make_s2_manifests.py), committed there as a local commit; the 2100-line task list is unchanged;
  PASTE A0+A: GitHub is not used: the clone becomes a shared clone of this repository (--no-checkout, core.abbrev 7, sparse init), the fetch
    of the build commit is not needed (local objects), kissfft comes from the main checkout's hspist3/kissfft (febd4ca), the fetch of the
    package commit fetches the stand-in and checks it out; BUILD SHIMS (P2's: module, gcc, make with the Mac CFLAGS, ldd); the SRUN SHIM runs
    srun's command here with its stdin (the here-document) and SLURM_JOB_ID set, as a sandbox job would;
  PASTE B: --array=1-2100%64 becomes the four chunks, the guard's package commit the stand-in's; the SBATCH SHIM runs the array's tasks one
    after the other, and the dependent merge at once when every task of the array ended with exit 0 (--dependency=afterok), else HOLDS it
    (Slurm's DependencyNeverSatisfied); --parsable prints the job id only; the paste is typed twice (the second must refuse);
  CHECK C, CHECK D2: as written (the Lustre root local); then the merge submitted again (it must refuse: exit 3), and a merge held on a job
    with a failed task (a pretended exits record) must not run;
  PASTE E: the fetch runs from the stand-in clone (its DST is the stand-in's own hspist3/experiments_gen3_stage2_koa_261016/S2).
##CHRIS 2026-10-16 (the review before Report D58-1; notes 261012 sec. 4.7.143), the runsheet's second version: CHECK 0 in the stand-in
  clone (its origin is this repository: the package commit is in origin/main; the stand-in commit is not); paste A0+A with a pretended
  failed build in a second home (it must stop at the build: no package checkout, no record); paste B with the root absent but the logs
  present (the log guard must refuse); CHECK D2 on the refused second merge's log (it must show the count 0 and the STOP line).
Then s2_tables_mac.py on the fetched folder. Full-length runs (200 predicted periods). Writes into --root (scratch: clones, data) and --out
(a new folder: every paste's transcript, the fetched tables, the S2 tables).
usage (from hspist3/): python3 cluster/gen3_stage2_261016/rehearse_s2_mac.py --pkg <40 hex> --root <scratch folder> --out <new folder>
"""
import argparse, os, re, shlex, shutil, subprocess, sys, time
HERE = os.path.dirname(os.path.abspath(__file__)); HS = os.path.dirname(os.path.dirname(HERE)); REPO = os.path.dirname(HS)
sys.path.insert(0, os.path.join(HS, "cluster", "gen3_koa_261009")); sys.path.insert(0, HERE)
import rehearse_k2_mac as K2R
import rehearse_p2_mac as P2R
P = "cluster/gen3_stage2_261016"; BUILD = "8308fb004fdcb7e9c82197451cc04ab49ecd2f13"
PICK = (("S2_N0900_eta0.760000_k0697_pinturned", 4500), ("S2_N0900_eta0.780000_k0685_pin0", 4500),
        ("S2_N1600_eta0.720000_k0866_pinturned", 8000), ("S2_N1600_eta0.800000_k1024_pinturned", 8000))
STEPS = []
SBATCH_S2 = r"""#!/bin/bash
# rehearsal shim (S2): rehearse_k2_mac's sbatch shim (the job runs AT ONCE as Slurm would run it: SLURM_SUBMIT_DIR, a pretended job id, the
# output file of its #SBATCH --output line; --array=<list> runs its tasks one after the other; --exclude recorded) plus: --parsable prints the
# job id only; --dependency=afterok:<id> runs the job at once if every task of <id> in the exits record ended with exit 0, else HOLDS it (it is
# not run: Slurm keeps it pending, DependencyNeverSatisfied); any other option is refused
ARR=""; EXC=""; PARS=""; DEP=""
while [ $# -gt 1 ]; do
  case "$1" in --array=*) ARR=${1#--array=}; ARR=${ARR%%\%*};; --exclude=*) EXC=${1#--exclude=};; --parsable) PARS=1;;
    --dependency=afterok:*) DEP=${1#--dependency=afterok:};; *) echo "rehearsal sbatch shim: option not supported: $1" >&2; exit 97;; esac
  shift
done
S=$1; C="$REHEARSAL_ROOT/logs/sbatch_shim_counter"; X="$REHEARSAL_ROOT/logs/sbatch_shim_exits.txt"
id=$(( $(cat "$C" 2>/dev/null || echo 920000) + 1 )); echo "$id" > "$C"
name=$(sed -n 's/^#SBATCH --job-name=//p' "$S" | head -1); out=$(sed -n 's/^#SBATCH --output=//p' "$S" | head -1)
if [ -n "$PARS" ]; then echo "$id"; else echo "Submitted batch job $id"; fi
[ -z "$EXC" ] || echo "$id exclude $EXC (the Mac is one node: recorded only) 0" >> "$X"
if [ -n "$DEP" ]; then
  n=$(grep -c "^${DEP}_[0-9]* " "$X" 2>/dev/null); bad=$(grep "^${DEP}_[0-9]* " "$X" 2>/dev/null | awk '$NF != 0' | wc -l | tr -d ' ')
  if [ "${n:-0}" -eq 0 ] || [ "${bad:-0}" -ne 0 ]; then
    echo "$id held afterok:$DEP ($bad of ${n:-0} tasks with exit != 0) -" >> "$X"
    echo "[rehearsal sbatch shim] job $id HELD: afterok:$DEP is never satisfied ($bad of ${n:-0} tasks ended with exit != 0); Slurm keeps it pending" >&2
    exit 0
  fi
  echo "[rehearsal sbatch shim] job $id: afterok:$DEP satisfied (${n} tasks, all exit 0): it runs now" >&2
fi
if [ -z "$ARR" ]; then
  o=${out//%x/$name}; o=${o//%j/$id}
  SLURM_SUBMIT_DIR=$PWD SLURM_JOB_ID=$id SLURM_JOB_NAME=$name bash "$S" > "$o" 2>&1; rc=$?
  echo "$id $PWD/$S $rc" >> "$X"
else
  for t in $(echo "$ARR" | tr ',' '\n' | awk -F- '{ if (NF == 2) for (i = $1; i <= $2; i++) print i; else print $1 }'); do
    o=${out//%x/$name}; o=${o//%A/$id}; o=${o//%a/$t}
    SLURM_SUBMIT_DIR=$PWD SLURM_JOB_ID=$((id * 1000 + t)) SLURM_ARRAY_JOB_ID=$id SLURM_ARRAY_TASK_ID=$t SLURM_JOB_NAME=$name bash "$S" > "$o" 2>&1; rc=$?
    echo "${id}_$t $PWD/$S $rc" >> "$X"
  done
fi
exit 0
"""
SRUN_S2 = r"""#!/bin/bash
# rehearsal shim (S2): the sandbox job is this Mac -- srun's options are dropped, its command runs here with srun's stdin (the paste's
# here-document) and SLURM_JOB_ID set, as Slurm sets it in the job
while [ $# -gt 0 ]; do case "$1" in -p|-t|-c|-A|-N|-n) shift 2;; --*=*|--pty) shift;; -*) shift;; *) break;; esac; done
echo "[rehearsal srun shim] the sandbox job runs here: $*" >&2
export SLURM_JOB_ID=${REHEARSAL_SRUN_JOB_ID:-930001}
exec "$@"
"""


def step(lab, ok, det): STEPS.append((lab, bool(ok), det)); print(f"[{'PASS' if ok else 'FAIL'}] {lab}: {det}", flush=True)


def blocks(path):
    """the runsheet's paste blocks: header 'S2 PASTE|CHECK X (...)' -> its command lines (up to the first empty line)"""
    out, cur = {}, None
    for l in open(path).read().split("\n"):
        m = re.match(r"^S2 (?:PASTE|CHECK) (0|A0\+A|B|C|D2|E) \(", l)
        if m: cur = m.group(1); out[cur] = []; continue
        if cur is not None:
            if not l.strip(): cur = None if out[cur] else cur; continue
            if l.startswith("  "): cur = None; continue
            out[cur].append(l)
    return out


def units(lines):
    """the paste's commands: one line each, except a here-document (from its '<<' line to its delimiter), which is one command"""
    out, i = [], 0
    while i < len(lines):
        m = re.search(r"<<'(\w+)'\s*$", lines[i])
        if m:
            j = lines.index(m.group(1), i + 1); out.append(lines[i:j + 1]); i = j + 1
        else: out.append([lines[i]]); i += 1
    return out


def run_paste(lines, env, cwd, log, title):
    """the paste in one bash; every command shown as '$ <its lines>' before it runs (a terminal's echo of the paste); stderr into stdout"""
    script = ""
    for u in units(lines):
        script += "".join(f"printf '%s\\n' {shlex.quote('$ ' + l)}\n" for l in u) + "\n".join(u) + "\n"
    t0 = time.time(); p = subprocess.run(["/bin/bash"], input=script, env=env, cwd=cwd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
    open(log, "w").write(f"### {title}\n### exit {p.returncode}, {time.time() - t0:.0f} s\n" + p.stdout)
    print(f"  {title}: exit {p.returncode}, {time.time() - t0:.0f} s, log {os.path.basename(log)}", flush=True)
    return p.returncode, p.stdout


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--pkg", required=True); ap.add_argument("--root", required=True); ap.add_argument("--out", required=True)
    a = ap.parse_args(); R = os.path.abspath(a.root); O = os.path.abspath(a.out)
    for d in (R, O):
        if os.path.exists(d): sys.exit(f"REFUSED: {d} exists (a new folder only)")
    os.makedirs(os.path.join(R, "logs")); os.makedirs(O); home = os.path.join(R, "home"); os.makedirs(home)
    LUS = os.path.join(R, "lustre"); B = blocks(os.path.join(HS, P, "runsheet_s2_261017.txt"))
    step("the runsheet's pastes and checks", set(B) == {"0", "A0+A", "B", "C", "D2", "E"}, ", ".join(f"{k}: {len(v)} lines" for k, v in B.items()))
    NOSUB = "S2: NO SUBMISSION CONFIRMED BY THIS PASTE"
    TL = open(os.path.join(HS, P, "tasks_s2_koa_261017.txt")).read().split("\n")
    MINI = tuple(next(j for j, l in enumerate(TL, 1) if l.split()[1:4] == [f"stage2/gen3/{c}/m_{M}", str(M), "0"]) for c, M in PICK)
    # the stand-in of the package commit
    S = os.path.join(R, "standin"); sh_ = K2R.clone(S, a.pkg, K2R.SPARSE_KOA); P2R.kissfft(S)
    ct = os.path.join(S, "hspist3", P, "chunks_s2_koa_261017.tsv"); L = open(ct).read().split("\n")
    keep = [L[0]] + [l for l in L[1:] if l and int(l.split("\t")[0]) in MINI]; open(ct, "w").write("\n".join(keep) + "\n")
    done = K2R.rewrite_lustre(S, R, local_fetch=True)
    mk = subprocess.run([sys.executable, "-B", os.path.join(S, "hspist3", P, "make_s2_manifests.py")], cwd=os.path.join(S, "hspist3"), capture_output=True, text=True)
    for c in (["config", "user.email", "rehearsal@localhost"], ["config", "user.name", "S2 rehearsal"], ["add", "-A", "hspist3/" + P],
              ["commit", "-q", "-m", "S2 rehearsal stand-in: the mini chunk table, local Lustre and fetch, the input manifest again"], ["branch", "standin"]):
        K2R.sh(["git", "-C", S] + c)
    SC = K2R.sh(["git", "-C", S, "rev-parse", "HEAD"]).stdout.strip()
    step("the stand-in commit", sh_ == a.pkg and len(keep) == 5 and mk.returncode == 0,
         f"{SC[:12]} on {a.pkg[:12]}: chunk table {len(keep) - 1} chunks {MINI}; Lustre/fetch made local in {len(done)} files (the S2 ones: "
         f"{', '.join(d for d in done if 's2_' in d)}); manifests: {mk.stdout.strip().splitlines()[1] if mk.returncode == 0 else mk.stderr[-200:]}")
    shim = K2R.shims(R); open(os.path.join(shim, "sbatch"), "w").write(SBATCH_S2); open(os.path.join(shim, "srun"), "w").write(SRUN_S2)
    for n in ("sbatch", "srun"): os.chmod(os.path.join(shim, n), 0o755)
    bsh = P2R.build_shims(R); env = K2R.env_for(R, home, shim, REHEARSAL_ROOT=R); envA = dict(env, PATH=f"{bsh}:{env['PATH']}")
    # check 0 (Mac): in the stand-in clone, whose origin is this repository
    c0 = [l.replace("cd ~/Desktop/CCS_complex_coupled_systems/Repo/HardDisks", f"cd {S}") for l in B["0"]]
    rc, out = run_paste(c0, env, home, os.path.join(O, "check_0.txt"), "check 0 (the package commit in origin/main)")
    rc, out2 = run_paste([l.replace(a.pkg, SC).replace(a.pkg[:8], SC[:8]) for l in c0], env, home, os.path.join(O, "check_0_negative.txt"),
                         "check 0 for the stand-in commit (not in origin/main)")
    # every check below reads PRINTED lines only (anchored, re.M): run_paste echoes each command as '$ <line>', and the commands hold the texts they echo
    step("check 0: the package commit is found in origin/main; a commit that is not is reported (the stand-in's origin/main is current: "
         "the fetch updating a stale one is not tested here)", re.search(rf"^S2: {a.pkg[:8]} IS ON GITHUB \(in origin/main\)$", out, re.M) is not None
         and not re.search(r"^S2: \w+ IS NOT ON GITHUB|^S2: git fetch failed", out, re.M) and re.search(rf"^S2: {SC[:8]} IS NOT ON GITHUB YET", out2, re.M) is not None
         and not re.search(r"^S2: \w+ IS ON GITHUB|^S2: git fetch failed", out2, re.M), "; ".join(l for l in (out + out2).splitlines() if l.startswith("S2:")))
    # paste A0+A
    subs = [("git clone --depth 1 --branch main --filter=blob:none --sparse https://github.com/FreeworkEarth/HardDisks.git harddisks_gen3_s2",
             f"git clone -q --no-checkout --shared {REPO} harddisks_gen3_s2 && git -C harddisks_gen3_s2 config core.abbrev 7 && git -C harddisks_gen3_s2 sparse-checkout init --no-cone"),
            (f"git fetch --depth 1 origin {BUILD}", "true"), ("git clone https://github.com/mborgerding/kissfft.git hspist3/kissfft", f"git clone -q {os.path.join(HS, 'kissfft')} hspist3/kissfft"),
            (f"git fetch --depth 1 origin {a.pkg}", f"git fetch -q {S} standin"), (f"checkout {a.pkg}", f"checkout {SC}"), ("/mnt/lustre/koa/scratch/charing", LUS)]
    linesA, nsub = [], 0
    for l in B["A0+A"]:
        for x, y in subs:
            if x in l: l = l.replace(x, y); nsub += 1
        linesA.append(l)
    rc, out = run_paste(linesA, envA, home, os.path.join(O, "paste_A0A.txt"), "paste A0+A (the sandbox job through the srun shim)")
    tail = out.split("\nthe merge and reducer imports: OK\n")[-1]
    ok = (re.search(r"^BUILD OK$", out, re.M) and re.search(r"^BUILD RECORD OK: ", out, re.M) and re.search(r"^the merge and reducer imports: OK$", out, re.M)
          and re.search(rf"^version +00ALLINONE  git {BUILD[:7]}  target koa$", out, re.M) and re.search(r"^build_inputs (\d+) of \1 equal", out, re.M)
          and re.search(r"^ls: .*S2: No such file or directory$", tail, re.M) and re.search(r"^S2 PASTE A0\+A: END OF THE SANDBOX JOB$", out, re.M)
          and not re.search(r"^STOP", out, re.M))   # set -v echoes the commands, whose text holds 'STOP: ...'
    noecho = not re.search(r"^module purge$", out, re.M)   # koa_env.sh is sourced with set -v off: its lines are not echoed
    step("paste A0+A: clone, build at the package's parent, package checkout, record, checks", rc == 0 and bool(ok) and noecho, f"{nsub} substitutions; exit {rc}; "
         f"koa_env.sh's lines echoed: {'no' if noecho else 'YES'}; " + "; ".join(x for x in re.findall(r"^(commit_build \S+|version .*|build_inputs .*|\d+$)", out, re.M)[:6]))
    # paste A0+A with a pretended failed build, in a second home: it must stop at the build
    home2 = os.path.join(R, "home_failed_build"); os.makedirs(home2)
    fb = "bash -c \"echo 'STOP: make koa failed (the rehearsal pretends a failed build)'; exit 4\""
    linesF = [l.replace("bash hspist3/cluster/build_koa.sh", fb) for l in linesA]
    rc, out = run_paste(linesF, dict(envA, HOME=home2), home2, os.path.join(O, "paste_A0A_failed_build.txt"), "paste A0+A with a pretended failed build (a second home)")
    step("paste A0+A with a failed build: STOP at the build, no package checkout, no record", re.search(r"^STOP: the build failed$", out, re.M) is not None
         and not re.search(rf"^HEAD is now at {SC[:7]}", out, re.M) and not re.search(r"^BUILD RECORD OK", out, re.M)
         and not re.search(r"^S2 PASTE A0\+A: END OF THE SANDBOX JOB$", out, re.M)
         and not os.path.exists(os.path.join(home2, "harddisks_gen3_s2", "hspist3", "logs", "BUILD_RECORD_S2.txt")),
         f"exit {rc}; " + "; ".join(l for l in out.splitlines() if l.startswith("STOP")))
    H2 = os.path.join(home, "harddisks_gen3_s2", "hspist3")
    # paste B, twice
    linesB = [l.replace("--array=1-2100%64", "--array=" + ",".join(map(str, MINI))).replace("/mnt/lustre/koa/scratch/charing", LUS).replace(a.pkg, SC).replace(a.pkg[:8], SC[:8]) for l in B["B"]]
    rc, out = run_paste(linesB, env, home, os.path.join(O, "paste_B_first.txt"), "paste B (the array: 4 tasks one after the other; then the dependent merge)")
    X = open(os.path.join(R, "logs", "sbatch_shim_exits.txt")).read() if os.path.exists(os.path.join(R, "logs", "sbatch_shim_exits.txt")) else ""
    tex = [l.split()[-1] for l in X.splitlines() if re.match(r"^\d+_\d+ ", l)]
    step("paste B: the guarded submit, the array and the dependent merge", re.search(r"^CHECKS GOOD$", out, re.M) and re.search(r"^S2: THE ARRAY IS SUBMITTED: job \d+$", out, re.M)
         and re.search(r"^Submitted batch job \d+$", out, re.M) and re.search(r"^S2: THE MERGE IS SUBMITTED, held by Slurm", out, re.M)
         and tex == ["0"] * 4 and re.search(r"^\[rehearsal sbatch shim\] job \d+: afterok:\d+ satisfied", out, re.M) and not re.search(r"^REFUSED", out, re.M),
         "; ".join(l for l in out.splitlines() if l.startswith(("CHECKS", "Submitted", "S2:", "[rehearsal sbatch shim]"))) + f"; task exits {tex}")
    rc, out = run_paste(linesB, env, home, os.path.join(O, "paste_B_second.txt"), "paste B pasted again")
    step("paste B a second time: refused (by the root guard)", re.search(r"^REFUSED: .*/stage2/S2 exists: the array ran before$", out, re.M) and len(re.findall(r"^REFUSED", out, re.M)) == 1
         and re.search(rf"^{NOSUB}", out, re.M) and not re.search(r"^Submitted|^CHECKS GOOD", out, re.M),
         "; ".join(l for l in out.splitlines() if l.startswith(("REFUSED", "S2:"))))
    linesB3 = [l.replace("--array=1-2100%64", "--array=" + ",".join(map(str, MINI))).replace("/mnt/lustre/koa/scratch/charing", LUS + "_absent").replace(a.pkg, SC).replace(a.pkg[:8], SC[:8]) for l in B["B"]]
    rc, out = run_paste(linesB3, env, home, os.path.join(O, "paste_B_root_absent.txt"), "paste B with the root absent but the logs present")
    step("paste B with the root absent and the logs present: refused by the log guard", re.search(r"^REFUSED: logs/gen3-s2\*\.out exist", out, re.M)
         and len(re.findall(r"^REFUSED", out, re.M)) == 1 and re.search(rf"^{NOSUB}", out, re.M) and not re.search(r"^Submitted|^CHECKS GOOD", out, re.M),
         "; ".join(l for l in out.splitlines() if l.startswith(("REFUSED", "S2:"))))
    # check C, check D2
    rc, out = run_paste([l.replace("/mnt/lustre/koa/scratch/charing", LUS) for l in B["C"]], env, home, os.path.join(O, "check_C.txt"), "check C (read only)")
    step("check C: 4 logs, 4 done, 0 FAILED/REFUSED/STOP, the chunk check", re.search(r"^chunk logs: 4; chunks done: 4 of 4 ", out, re.M) is not None and re.findall(r"^\s*(\d+)\s*$", out, re.M)[:3] == ["4", "4", "0"],
         " | ".join(l.strip() for l in out.splitlines() if l.strip() and not l.startswith("$")))
    rc, out = run_paste(B["D2"], env, home, os.path.join(O, "check_D2.txt"), "check D2 (the merge's log; read only)")
    step("check D2: the merge and p2_reduce", re.search(r"^S2 merge and reduction done", out, re.M) and re.search(r"^p2_reduce: 4 folders \(0 without the stage 1 tables\), 4 trajectories, 0 with an error", out, re.M)
         and re.search(r"^merged 4 chunks: ", out, re.M) and re.search(r"^4$", out, re.M) and not re.search(r"^STOP", out, re.M), "; ".join(l[:130] for l in out.splitlines() if l.startswith(("merged ", "p2_reduce:", "S2 merge", "4"))))
    rc, out = run_paste([f"cd {H2} && sbatch {P}/s2_merge.sbatch", f"tail -2 \"$(ls -t {H2}/logs/gen3-s2-merge_*.out | head -1)\""], env, home,
                        os.path.join(O, "merge_again.txt"), "the merge submitted again (no guard, no dependency)")
    step("the merge a second time: refused (exit 3)", re.search(r"^STOP: .*/merged exists -- a merge ran before; not merging twice$", out, re.M) is not None, "; ".join(l[:160] for l in out.splitlines() if "STOP" in l))
    rc, out = run_paste(B["D2"], env, home, os.path.join(O, "check_D2_on_the_refused_merge.txt"), "check D2 on the newest merge log (the refused one)")
    step("check D2 on a stopped merge's log: the count 0 and its STOP lines shown", re.search(r"^0$", out, re.M) is not None and re.search(r"^STOP", out, re.M) is not None,
         "; ".join(l[:160] for l in out.splitlines() if re.match(r"^(0|STOP)", l)))
    nlog = len([f for f in os.listdir(os.path.join(H2, "logs")) if f.startswith("gen3-s2-merge_")])
    open(os.path.join(R, "logs", "sbatch_shim_exits.txt"), "a").write(f"990001_1 pretended-task 0\n990001_2 pretended-task 2\n")
    rc, out = run_paste([f"cd {H2} && sbatch --dependency=afterok:990001 {P}/s2_merge.sbatch"], env, home, os.path.join(O, "merge_held.txt"),
                        "a merge held on a job with a failed task (pretended exits record: task 2 exit 2)")
    nlog2 = len([f for f in os.listdir(os.path.join(H2, "logs")) if f.startswith("gen3-s2-merge_")])
    step("the dependency: a merge on a job with a failed task is held, not run", re.search(r"^\[rehearsal sbatch shim\] job \d+ HELD: afterok:990001", out, re.M) and nlog2 == nlog, f"merge logs before {nlog}, after {nlog2}; "
         + "; ".join(l[:170] for l in out.splitlines() if "HELD" in l))
    # paste E
    rc, out = run_paste([l.replace("cd ~/Desktop/CCS_complex_coupled_systems/Repo/HardDisks", f"cd {S}") for l in B["E"]], env, home, os.path.join(O, "paste_E.txt"), "paste E (the fetch, local)")
    F = os.path.join(S, "hspist3", "experiments_gen3_stage2_koa_261016", "S2")
    step("paste E: the fetch", re.search(r"^fetched into .*p2_reduce: p2_reduce: 4 folders", out, re.M) is not None, out.strip().splitlines()[-1][:220])
    snaps = [f for dp, dn, fn in os.walk(F) for f in fn if f.startswith("snapshots_")]
    traces = [f for dp, dn, fn in os.walk(F) for f in fn if f.startswith("wall_x_positions_")]
    dumps = [f for dp, dn, fn in os.walk(F) for f in fn if f == "initial_state.txt"]
    step("the fetch's filter: 4 gzipped traces, no snapshot, the 4 initial states", not snaps and len(traces) == 4 and all(t.endswith(".csv.gz") for t in traces) and len(dumps) == 4,
         f"snapshots {len(snaps)}; traces {len(traces)} ({', '.join(sorted({t.rsplit('.', 2)[-2] + '.' + t.rsplit('.', 1)[-1] for t in traces}))}); initial states {len(dumps)}")
    shutil.copytree(os.path.join(F, "merged"), os.path.join(O, "fetched_merged"), ignore=shutil.ignore_patterns(".done_runs", "wall_x_positions_*"))
    t = subprocess.run([sys.executable, "-B", os.path.join(HS, P, "s2_tables_mac.py"), "--data", F, "--out", os.path.join(O, "s2_tables")], capture_output=True, text=True)
    open(os.path.join(O, "s2_tables_stdout.txt"), "w").write(t.stdout + t.stderr)
    step("s2_tables_mac.py on the fetched folder", t.returncode == 0, f"exit {t.returncode}; {len(t.stdout.splitlines())} lines" + ("" if t.returncode == 0 else f"; {t.stderr[-300:]}"))
    print("\n== the rehearsal\n\n| step | result | detail |\n|---|---|---|")
    rows = [f"| {lab} | {'PASS' if ok else 'FAIL'} | {det} |" for lab, ok, det in STEPS]; print("\n".join(rows))
    tail = f"\nS2 REHEARSAL: {'PASS' if all(s[1] for s in STEPS) else 'FAIL'} ({sum(s[1] for s in STEPS)} of {len(STEPS)})"
    print(tail); open(os.path.join(O, "REHEARSAL.md"), "w").write("| step | result | detail |\n|---|---|---|\n" + "\n".join(rows) + "\n" + tail + "\n")
    print(f"data under {R}: {K2R.sh(['du', '-sh', R]).stdout.split()[0]}; under {O}: {K2R.sh(['du', '-sh', O]).stdout.split()[0]}")


if __name__ == "__main__":
    main()
