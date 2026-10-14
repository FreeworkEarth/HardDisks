#!/usr/bin/env python3
"""##CHRIS 2026-10-13 (plan-author decision 47 item 3): RC BLOCK C's guarded line, typed as the contingency runsheet prints it, in bash with
shims (nothing is submitted): a fake home with ~/harddisks_gen3_rc (.git/HEAD at the named commit) and ~/harddisks_gen3_k2's merge log; D's
Lustre root replaced by a scratch folder (the only substitution); squeue (a pretended queue filtered by -n), sbatch (records, prints a job).
usage (from hspist3/): python3 experiments_gen3_s1_261012/complete_rehearsal/test_rc_block_c.py <contingency runsheet> <scratch dir>
"""
import os, subprocess, sys
rs, tmp = sys.argv[1], os.path.abspath(sys.argv[2])
T = open(rs).read(); C = T.split("RC BLOCK C (")[1].split("RC BLOCK D (")[0]
line = [l for l in C.splitlines() if l.startswith("cd ~/harddisks_gen3_rc/hspist3 && ")][0]
full = T.split("git fetch --depth 1 origin ")[1].split()[0]
LUS = "/mnt/lustre/koa/scratch/charing/harddisks_gen3/stage1/D"; lroot = os.path.join(tmp, "lustre_D")
line = line.replace(LUS, lroot)
H = os.path.join(tmp, "home"); os.makedirs(os.path.join(H, "harddisks_gen3_rc", ".git")); os.makedirs(os.path.join(H, "harddisks_gen3_rc", "hspist3", "logs"))
os.makedirs(os.path.join(H, "harddisks_gen3_k2", "hspist3", "logs"))
b = os.path.join(tmp, "bin"); os.makedirs(b)
for n, t in {"squeue": "#!/bin/bash\nN=''; while [ $# -gt 0 ]; do case \"$1\" in -n) N=$2; shift;; esac; shift; done\n"
                       "[ -n \"${QUEUE:-}\" ] && awk -v n=\",$N,\" 'index(n, \",\" $3 \",\")' \"$QUEUE\"\nexit 0\n",
             "sbatch": "#!/bin/bash\necho \"$*\" >> \"$SBATCH_LOG\"; echo 'Submitted batch job 9600001'\n"}.items():
    open(os.path.join(b, n), "w").write(t); os.chmod(os.path.join(b, n), 0o755)
def setup(head=full, mc=True, rc=False, mlog=True):
    open(os.path.join(H, "harddisks_gen3_rc", ".git", "HEAD"), "w").write(head + "\n")
    os.makedirs(os.path.join(lroot, "merged"), exist_ok=True)
    for f, want in ((".merge_complete", mc), (".reduction_completed", rc)):
        p = os.path.join(lroot, "merged", f)
        if want and not os.path.exists(p): open(p, "w").write("x\n")
        if not want and os.path.exists(p): os.rename(p, p + f".aside_{len(os.listdir(os.path.join(lroot, 'merged')))}")
    p = os.path.join(H, "harddisks_gen3_k2", "hspist3", "logs", "gen3-s1-merge-D_15055588.out")
    if mlog and not os.path.exists(p): open(p, "w").write(f"merge of {LUS} by the scripts of ~/harddisks_gen3_k2\n")
    if not mlog and os.path.exists(p): os.rename(p, p + ".aside")
cases = [("all checks hold", {}, [], True, None),
         ("the merge job still in squeue", {}, ["9100 RUNNING gen3-s1-merge-D"], False, "gen3-s1-merge-D"),
         ("a completion already queued (a second paste)", {}, ["9101 PENDING gen3-s1-complete-D"], False, "gen3-s1-complete-D"),
         ("the rc clone at another commit", dict(head="0" * 40), [], False, "the rc clone is not at"),
         ("no .merge_complete (the merge stopped during its moves)", dict(mc=False), [], False, ".merge_complete is missing"),
         (".reduction_completed exists", dict(rc=True), [], False, "completed before"),
         ("the merge log missing", dict(mlog=False), [], False, "the merge log")]
print(f"# RC block C's guard, from {os.path.basename(rs)} (commit {full[:7]}; D's root replaced by a scratch folder)\n")
print("| case | expected | submitted | printed | result |\n|---|---|---|---|---|")
allok = True
for i, (lab, st, q, submit, why) in enumerate(cases):
    setup(**st); qf = os.path.join(tmp, f"q{i}.txt"); open(qf, "w").write("".join(x + "\n" for x in q)); sl = os.path.join(tmp, f"sbatch{i}.txt")
    e = dict(os.environ, HOME=H, PATH=b + ":" + os.environ["PATH"], QUEUE=qf, SBATCH_LOG=sl)
    out = subprocess.run(["bash", "-c", line], env=e, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True).stdout
    sub = open(sl).read().splitlines() if os.path.exists(sl) else []
    ok = (len(sub) == 1 and "CHECKS GOOD" in out and "D: THE COMPLETION SUBMITTED" in out and sub[0].startswith("--time=24:00:00 --job-name=gen3-s1-complete-D "
          "cluster/gen3_s1_complete_261013/s1_complete.sbatch D " + H + "/harddisks_gen3_k2/hspist3/logs/gen3-s1-merge-D_15055588.out")) if submit else \
         (not sub and "REFUSED: " in out and why in out and "D: THIS PASTE DID NOT SUBMIT THE COMPLETION" in out and "CHECKS GOOD" not in out)
    allok &= ok; setup()
    shown = "; ".join(l for l in out.splitlines() if l.startswith(("CHECKS GOOD", "REFUSED", "D: ")))
    print(f"| {lab} | {'submits' if submit else 'refuses: ' + why} | {len(sub)} | {shown[:200]} | {'PASS' if ok else '**FAIL**'} |")
print(f"\nRC BLOCK C GUARD TEST: {'PASS' if allok else 'FAIL'} ({len(cases)} cases)")
