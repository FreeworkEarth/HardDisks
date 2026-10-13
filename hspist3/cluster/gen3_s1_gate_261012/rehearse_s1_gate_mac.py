#!/usr/bin/env python3
"""##CHRIS 2026-10-12 (plan-author decision 39 item 5 c; decision 20's rule): THE MAC REHEARSAL OF THE STAGE 1 GATE PACKAGE, end to end, before
Chris types it: the build record, the gate job, the cross-node array and its compare, the fetch and the gate verdict, in scratch stand-ins
of the four KOA folders. Nothing outside --root is written. Its stand-ins also serve the stage 1 rehearsal (rehearse_s1_mac.py, gate G4).
  home/harddisks_gen3      a symlink to --gen3-home: the a643bc7 Mac builds of rehearse_k2_mac.py (koa, koa-ld, the harnesses, their
                           hash files) -- the gate binary; read only.
  home/harddisks_resched2  7b08827 (the E0 reference), a shared clone with KOA's step-1 sparse set and the pinned kissfft, built here with
                           rehearse_k2_mac.build() (the Makefile's koa target with the Mac flags, so the version line reads "target koa").
  home/harddisks_gen3_s1   f72d764 (the campaign build), the same clone and builds (koa, koa-ld, the harnesses, the hash files) -- the
                           NEW clone of decision 39 item 5 a; record_s1_build.sh then runs as in the runsheet's block B.
  home/harddisks_gen3_k2   --commit (the gate package's commit) with the K2 clone's sparse set: runs the jobs, the fetch and the verdict.
SHIMS (each a listed residual difference): the four rewritten files of the scripts clone -- k2_binaries_a643bc7.sha256 (the Mac a643bc7
builds' SHA-256), koa_harness_ref_a643bc7.sha256 (the SHA-256 of the a643bc7 stand-in's own harness outputs, run here once: the KOA files
name GCC 14.3 and glibc's libm), s1_gate_inputs.sha256 regenerated for those two, and the Lustre root in the three job scripts and the
fetch (<root>/lustre); E0's case list in s1_gate_koa.sbatch cut to ctrl_min,ctrl_leg (and the verdict told so: --e0-cases) (afix and afix_leg write six event logs of 17 MB
per build: about 200 MB, which this message's data cap does not hold; the runner is the a643bc7 gate's, unchanged, and ran all four on
KOA in job 15040861); sbatch emulated here (run at once, the #SBATCH --output name, the Slurm variables; --array runs its tasks one after
another; --dependency and --exclude recorded); squeue prints nothing; module, gcc --version (14.3.0), lscpu, and hostname (the array's
odd tasks "cn-r1", even "cn-r2", so the compare sees two nodes; otherwise "cn-r0"); HOME. On arm64 the long double is the double.
usage (from hspist3/): python3 cluster/gen3_s1_gate_261012/rehearse_s1_gate_mac.py --root R --commit C --gen3-home H
"""
import argparse, glob, os, re, shutil, subprocess, sys
HERE = os.path.dirname(os.path.abspath(__file__)); HS = os.path.dirname(os.path.dirname(HERE))
sys.path.insert(0, os.path.join(HS, "cluster", "gen3_koa_261009"))
import rehearse_k2_mac as K2R
F72 = "f72d764a8a7e837644c26bd8abc175d9e2cac268"; P = "cluster/gen3_s1_gate_261012"
STEPS = []


def step(lab, ok, det): STEPS.append((lab, bool(ok), det)); print(f"[{'PASS' if ok else 'FAIL'}] {lab}: {det}", flush=True)


def rewrite(path, old, new, count=1):
    t = open(path).read()
    if t.count(old) != count: sys.exit(f"STOP: {path} holds '{old}' {t.count(old)} times, not {count}")
    open(path, "w").write(t.replace(old, new))


def shims(R):
    d = os.path.join(R, "shims"); os.makedirs(d, exist_ok=True)
    cpu = K2R.sh(["sysctl", "-n", "machdep.cpu.brand_string"]).stdout.strip()
    for n, body in (("module", "exit 0\n"), ("gcc", 'if [ "${1:-}" = --version ]; then echo "gcc (GCC) 14.3.0"; else exec cc "$@"; fi\n'),
                    ("lscpu", f'echo "Model name:            {cpu}"\n'), ("squeue", "exit 0\n"),
                    ("hostname", 'if [ -n "${SLURM_ARRAY_TASK_ID:-}" ]; then [ $((SLURM_ARRAY_TASK_ID % 2)) -eq 1 ] && echo cn-r1 || echo cn-r2; else echo cn-r0; fi\n')):
        p = os.path.join(d, n); open(p, "w").write("#!/bin/bash\n" + body); os.chmod(p, 0o755)
    return d


def sbatch(R, env, cwd, script, args=(), array=None, extra=()):
    """emulated sbatch: runs the script at once (each array task in turn) with the Slurm variables and its #SBATCH --output file"""
    t = open(os.path.join(cwd, script)).read()
    name = re.search(r"^#SBATCH --job-name=(\S+)", t, re.M).group(1); outp = re.search(r"^#SBATCH --output=(\S+)", t, re.M).group(1)
    cf = os.path.join(R, "logs", "sbatch_counter"); os.makedirs(os.path.dirname(cf), exist_ok=True)
    jid = int(open(cf).read()) + 1 if os.path.exists(cf) else 930001; open(cf, "w").write(str(jid))
    rcs = []
    for task in (array or [None]):
        e = dict(env, SLURM_SUBMIT_DIR=cwd, SLURM_JOB_ID=str(jid), SLURM_JOB_NAME=name)
        o = outp.replace("%x", name).replace("%j", str(jid)).replace("%A", str(jid))
        if task is not None: e.update(SLURM_ARRAY_JOB_ID=str(jid), SLURM_ARRAY_TASK_ID=str(task)); o = o.replace("%a", str(task))
        else: e.pop("SLURM_ARRAY_TASK_ID", None); e.pop("SLURM_ARRAY_JOB_ID", None)
        with open(os.path.join(cwd, o), "w") as fh: rcs.append(subprocess.run(["bash", script] + list(args), cwd=cwd, env=e, stdout=fh, stderr=subprocess.STDOUT).returncode)
    with open(os.path.join(R, "logs", "sbatch_shim.txt"), "a") as fh: fh.write(f"{jid} {name} {script} {' '.join(args)} array={array} extra={' '.join(extra)} exits={rcs}\n")
    return jid, rcs


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--root", required=True); ap.add_argument("--commit", required=True)
    ap.add_argument("--gen3-home", required=True); a = ap.parse_args()
    R = os.path.abspath(a.root); G = os.path.abspath(a.gen3_home)
    if os.path.exists(R) and os.listdir(R): sys.exit(f"STOP: {R} exists and is not empty")
    home = os.path.join(R, "home"); os.makedirs(os.path.join(R, "logs")); os.makedirs(home)
    # ---------------------------------------------------------------- the stand-ins
    os.symlink(G, os.path.join(home, "harddisks_gen3")); GH = os.path.join(G, "hspist3")
    vg = K2R.sh([os.path.join(GH, "00ALLINONE"), "--version"]).stdout.splitlines()[0]
    step("the gate binary stand-in (a643bc7 Mac builds)", vg == "00ALLINONE  git a643bc7  target koa", vg)
    REF = os.path.join(home, "harddisks_resched2"); K2R.clone(REF, K2R.REF_C, K2R.SPARSE_KOA); K2R.kissfft(REF)
    vr = K2R.build(os.path.join(REF, "hspist3"), R, ld=False, harness=False)
    step("the 7b08827 reference (clone and Mac koa build)", vr[0] == "00ALLINONE  git 7b08827  target koa", vr[0])
    S1 = os.path.join(home, "harddisks_gen3_s1"); head = K2R.clone(S1, F72, K2R.SPARSE_KOA); K2R.kissfft(S1)
    vs = K2R.build(os.path.join(S1, "hspist3"), R)
    step("the campaign stand-in, the NEW clone at f72d764 (clone and Mac builds: koa, koa-ld, the harnesses)", head == F72 and vs == ["00ALLINONE  git f72d764  target koa", "00ALLINONE  git f72d764  target koa-ld"], f"{head[:12]}; {vs}")
    C = os.path.join(home, "harddisks_gen3_k2"); ch = K2R.clone(C, a.commit, K2R.SPARSE_K2); H2 = os.path.join(C, "hspist3"); os.makedirs(os.path.join(H2, "logs"), exist_ok=True)
    # ---------------------------------------------------------------- shims
    with open(os.path.join(H2, "cluster/gen3_koa_261009/k2_binaries_a643bc7.sha256"), "w") as fh:
        for f in ("00ALLINONE", "bin/00ALLINONE_ld", "bin/gen3_m1", "bin/gen3_m2"): fh.write(f"{K2R.sha(os.path.join(GH, f))}  {f}\n")
    rh = os.path.join(R, "ref_harness"); os.makedirs(rh)
    for n, cmd in (("m1_audit_output.txt", ["bin/gen3_m1", "audit"]), ("m2_audit_quick_output.txt", ["bin/gen3_m2", "audit", "--quick"]), ("m2_audit_output.txt", ["bin/gen3_m2", "audit"])):
        with open(os.path.join(rh, n), "w") as fh: subprocess.run([os.path.join(GH, cmd[0])] + cmd[1:], cwd=rh, stdout=fh, stderr=subprocess.STDOUT)
    with open(os.path.join(H2, P, "koa_harness_ref_a643bc7.sha256"), "w") as fh:
        fh.write("# REHEARSAL SHIM: the a643bc7 Mac stand-in's own harness outputs (the KOA file names GCC 14.3 and glibc)\n")
        for n in ("m1_audit_output.txt", "m2_audit_quick_output.txt", "m2_audit_output.txt"): fh.write(f"{K2R.sha(os.path.join(rh, n))}  {n}\n")
    for f in ("s1_gate_koa.sbatch", "s1_xnode.sbatch", "s1_xnode_compare.sbatch"):
        t = open(os.path.join(H2, P, f)).read(); n = t.count(K2R.LUSTRE)
        open(os.path.join(H2, P, f), "w").write(t.replace(K2R.LUSTRE, os.path.join(R, "lustre")))
    rewrite(os.path.join(H2, P, "s1_gate_koa.sbatch"), "--cases ctrl_min,ctrl_leg,afix,afix_leg", "--cases ctrl_min,ctrl_leg")
    rewrite(os.path.join(H2, P, "fetch_s1_gate_koa.sh"), 'SRC="charing@koa-dtn.its.hawaii.edu:/mnt/lustre/koa/scratch/charing/harddisks_gen3"',
            f'SRC="{os.path.join(R, "lustre", "harddisks_gen3")}"')
    K2R.sh([sys.executable, os.path.join(P, "make_s1_gate_manifest.py")], cwd=H2)
    changed = [l.strip() for l in K2R.sh(["git", "-C", C, "status", "--porcelain", "--untracked-files=no"]).stdout.splitlines() if l.strip()]
    step(f"the scripts clone at {ch[:7]} and its shims", len(changed) == 7, "; ".join(changed))
    env = dict(os.environ, HOME=home, PATH=shims(R) + ":" + os.environ["PATH"], REHEARSAL_ROOT=R)
    # ---------------------------------------------------------------- block B: the build record
    r = subprocess.run(["bash", os.path.join(H2, P, "record_s1_build.sh")], env=env, capture_output=True, text=True)
    step("block B: record_s1_build.sh", r.returncode == 0 and "BUILD RECORD OK" in r.stdout, r.stdout.strip().splitlines()[-1] if r.stdout.strip() else r.stderr[-200:])
    print(r.stdout)
    r2 = subprocess.run(["bash", os.path.join(H2, P, "record_s1_build.sh")], env=env, capture_output=True, text=True)
    step("STOP check: the record a second time", r2.returncode == 2 and "never rebuilt" in r2.stdout, r2.stdout.strip()[-150:])
    # ---------------------------------------------------------------- block C: the jobs
    jg, rg = sbatch(R, env, H2, os.path.join(P, "s1_gate_koa.sbatch"))
    gl = open(glob.glob(os.path.join(H2, "logs", f"gen3-s1-gate_{jg}.out"))[0]).read()
    step("block C: the gate job", rg == [0] and "s1 gate done" in gl and "STOP" not in gl, f"exit {rg}")
    jx, rx = sbatch(R, env, H2, os.path.join(P, "s1_xnode.sbatch"), array=list(range(1, 9)))
    step("block C: the cross-node array (8 tasks)", rx == [0] * 8, f"exits {rx}")
    jc, rc = sbatch(R, env, H2, os.path.join(P, "s1_xnode_compare.sbatch"), extra=("--dependency=singleton",))
    cl = open(glob.glob(os.path.join(H2, "logs", f"gen3-s1-xnode-compare_{jc}.out"))[0]).read()
    step("block C: the compare job", rc == [0] and "cross-node compare done" in cl, f"exit {rc}")
    jg2, rg2 = sbatch(R, env, H2, os.path.join(P, "s1_gate_koa.sbatch"))
    g2l = open(glob.glob(os.path.join(H2, "logs", f"gen3-s1-gate_{jg2}.out"))[0]).read()
    step("STOP check: the gate job a second time", rg2 == [2] and "never writes into an existing root" in g2l, g2l.strip().splitlines()[-1][:150])
    print("== the gate job's log\n" + gl + "\n== the compare job's log\n" + cl)
    # ---------------------------------------------------------------- block E: the fetch and the verdict
    f = subprocess.run(["bash", os.path.join(H2, P, "fetch_s1_gate_koa.sh")], cwd=C, env=env, capture_output=True, text=True)
    step("block E: the fetch (source made local)", f.returncode == 0 and "build record: f72d764" in f.stdout, f.stdout.strip().splitlines()[-1][:200])
    v = subprocess.run([sys.executable, os.path.join(P, "s1_gate_verdict_261012.py"), "--e0-cases", "ctrl_min,ctrl_leg"], cwd=H2, env=env, capture_output=True, text=True)
    print("== the gate verdict on the fetched folder\n" + v.stdout + v.stderr)
    step("the gate verdict: G6 (c) and G6 (d)", v.returncode == 0 and "G6 (c): PASS" in v.stdout and "G6 (d): PASS" in v.stdout,
         "; ".join(l for l in v.stdout.splitlines() if l.startswith("G6 (")))
    print("== the rehearsal\n\n| step | result | detail |\n|---|---|---|")
    for lab, ok, det in STEPS: print(f"| {lab} | {'PASS' if ok else 'FAIL'} | {det} |")
    print(f"\nGATE REHEARSAL: {'PASS' if all(s[1] for s in STEPS) else 'FAIL'} ({sum(s[1] for s in STEPS)} of {len(STEPS)})")
    print(f"data under {R}: {K2R.sh(['du', '-sh', R]).stdout.split()[0]}")


if __name__ == "__main__":
    main()
