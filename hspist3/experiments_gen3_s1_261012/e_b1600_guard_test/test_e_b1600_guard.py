#!/usr/bin/env python3
"""##CHRIS 2026-10-13 (plan-author decision 47 addendum 1, amendment A4): THE GUARD OF BLOCK E-B1600, typed as the runsheet prints it, in bash
with shims (nothing is submitted, nothing runs): a fake home with ~/harddisks_gen3_k2 (its .git/HEAD at 7bdfb67e, an empty logs/), and on PATH
squeue (the lines of a pretended queue whose third field is in -n), sinfo (prints SINFO_F for -o %F, fails with SINFO_FAIL), sbatch (records its
arguments, prints "Submitted batch job <n>"). The Lustre root of B1600 does not exist on the Mac (its check passes, as on KOA before E-B1600).
usage (from hspist3/): python3 experiments_gen3_s1_261012/e_b1600_guard_test/test_e_b1600_guard.py <runsheet> <scratch dir>
"""
import os, re, subprocess, sys

rs, tmp = sys.argv[1], os.path.abspath(sys.argv[2])
L = open(rs).read().split("S1 BLOCK E-B1600 (")[1].split("\nS1 BLOCK ")[0].splitlines()
line = [l for l in L if l.startswith("cd ~/harddisks_gen3_k2/hspist3 && ")][0]
narr = line.count(" sbatch --job-name=gen3-s1-B1600 ")
os.makedirs(os.path.join(tmp, "home", "harddisks_gen3_k2", ".git")); os.makedirs(os.path.join(tmp, "home", "harddisks_gen3_k2", "hspist3", "logs"))
open(os.path.join(tmp, "home", "harddisks_gen3_k2", ".git", "HEAD"), "w").write("7bdfb67ec99500bd901e2a138a1af74e5be18ebb\n")
b = os.path.join(tmp, "bin"); os.makedirs(b)
shims = {
    "squeue": "#!/bin/bash\nN=''; while [ $# -gt 0 ]; do case \"$1\" in -n) N=$2; shift;; esac; shift; done\n"
              "[ -n \"${QUEUE:-}\" ] && awk -v n=\",$N,\" 'index(n, \",\" $3 \",\")' \"$QUEUE\"\nexit 0\n",
    "sinfo": "#!/bin/bash\n[ -n \"${SINFO_FAIL:-}\" ] && { echo 'sinfo: error: Unable to contact slurm controller (test shim)' >&2; exit 1; }\n"
             "echo \"${SINFO_F:-70/12/8/90}\"\n",
    "sbatch": "#!/bin/bash\nC=\"$SBATCH_LOG\"; n=$(( $(wc -l < \"$C\" 2>/dev/null || echo 0) + 1 )); echo \"$*\" >> \"$C\"; echo \"Submitted batch job 9500$n\"\n",
}
for n, t in shims.items():
    open(os.path.join(b, n), "w").write(t); os.chmod(os.path.join(b, n), 0o755)
qdir = os.path.join(tmp, "queues"); os.makedirs(qdir)
cases = [("idle 12 of 90, a B900 array and merge job running", dict(SINFO_F="70/12/8/90"), ["9001_7 RUNNING gen3-s1-B900", "9002 RUNNING gen3-s1-merge-B900"], True, None),
         ("idle 9 of 90 (exactly 10 %)", dict(SINFO_F="73/9/8/90"), [], True, None),
         ("idle 8 of 90 (below 10 %)", dict(SINFO_F="74/8/8/90"), [], False, "the 128-core rule"),
         ("a gen3-s1-S array job queued", dict(SINFO_F="70/12/8/90"), ["9003_2 PENDING gen3-s1-S"], False, "gen3-s1-S job"),
         ("only S's merge job running (not its arrays)", dict(SINFO_F="70/12/8/90"), ["9004 RUNNING gen3-s1-merge-S"], True, None),
         ("sinfo fails", dict(SINFO_FAIL="1"), [], False, "the 128-core rule"),
         ("sinfo prints no A/I/O/T line", dict(SINFO_F="n/a"), [], False, "the 128-core rule"),
         ("a second paste: B1600's arrays queued", dict(SINFO_F="70/12/8/90"), ["9005_1 PENDING gen3-s1-B1600"], False, "gen3-s1-B1600 job")]
print(f"# the guard of E-B1600 (amendment A4), from {os.path.relpath(rs)}: {narr} array line(s)\n")
print("| case | expected | submitted | printed (CHECKS GOOD / REFUSED / sinfo) | result |\n|---|---|---|---|---|")
allok = True
for i, (lab, env, q, submit, why) in enumerate(cases):
    qf = os.path.join(qdir, f"q{i}.txt"); open(qf, "w").write("".join(x + "\n" for x in q))
    slog = os.path.join(tmp, f"sbatch_{i}.txt")
    e = dict(os.environ, HOME=os.path.join(tmp, "home"), PATH=b + ":" + os.environ["PATH"], QUEUE=qf, SBATCH_LOG=slog, **env)
    p = subprocess.run(["bash", "-c", line], env=e, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
    nsub = len(open(slog).read().splitlines()) if os.path.exists(slog) else 0
    out = p.stdout
    shown = "; ".join(l for l in out.splitlines() if l.startswith(("CHECKS GOOD", "REFUSED", "sinfo -p", "B1600:")))
    ok = (nsub == narr and "CHECKS GOOD" in out and f"B1600: ALL {narr} ARRAY LINE(S) SUBMITTED" in out) if submit else \
         (nsub == 0 and "CHECKS GOOD" not in out and "REFUSED: " in out and why in out and "THIS PASTE DID NOT SUBMIT" in out)
    allok &= ok
    print(f"| {lab} | {'submits' if submit else 'refuses: ' + why} | {nsub} | {shown[:230]} | {'PASS' if ok else '**FAIL**'} |")
print(f"\nE-B1600 GUARD TEST: {'PASS' if allok else 'FAIL'} ({len(cases)} cases)")
